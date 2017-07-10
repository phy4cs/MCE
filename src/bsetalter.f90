MODULE bsetalter
  
  use globvars
  use Ham
  use alarrays
  use outputs

!***********************************************************************************!
!*
!*         Basis Set alteration module Module
!*           
!*   Contains subroutines for:
!*
!*      1) Relocating the wavefunction on a static grid
!*      2) Leaking the basis set, used for large q values in the Hennon-Heiles model
!*      3) Cloning subroutine, which increases the basis set, used mainly for the
!*              spin boson model
!*      4) Retriving cloning information from previous calculations when using the
!*              AIMC-MCE propagation system
!*      
!***********************************************************************************!

contains


!-----------------------------------------------------------------------------------

  subroutine reloc_basis(bs, z0, nbf, x, time, gridsp, mup, muq)

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    complex(kind=8), dimension(:,:), intent(in) :: z0   ! This is the initial grid.
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: time, gridsp
    integer, intent(inout) :: nbf, x
    
    type(basisfn) :: bf
    complex(kind=8), dimension (:,:), allocatable :: ovrlp_mat, zt, znew    
    complex(kind=8), dimension (:), allocatable :: bigA,cnew,dnew,C_k
    real(kind=8), dimension (:), allocatable :: DC, qstrt, pstrt
    real(kind=8) :: realz, imagz
    integer, dimension (:,:), allocatable :: coordinates
    integer, dimension (:), allocatable :: qsize,psize,columnsize,columnrep,remove
    integer :: j, k, p, q, l, m, n, z, y, ierr

    if (errorflag .ne. 0) return

    if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
      write (0,"(a)") "Basis Set relocation called, but initial basis set was not a grid!"
      errorflag=1
      return
    end if

    ! Create a smaller grid for the first step, max size 1600 bfs.     
    ! This section has been disabled, and should be used with care as it will cause
    ! errors if used with the coulomb potential and not using the even grid+1 layout
    ! due to the absence of D values for the basis set
    if ((x==1).and.(basis=="GRID").and.(.false.)) then

      allocate (qstrt(ndim), stat=ierr)
      if (ierr==0) allocate (pstrt(ndim), stat=ierr)
      if (ierr==0) allocate (qsize(ndim), stat=ierr)
      if (ierr==0) allocate (psize(ndim), stat=ierr)
      if (ierr/=0) then
        write (0,"(a)") "Error in allocation starting p and q values for grid"
        errorflag=1
        return
      end if
      
      ! Set the grid dimensions for the smaller initial grid
      if (ndim==3) then
        qsize(1) = 4
        qsize(2) = 2
        qsize(3) = 2 
        psize(1) = 4
        psize(2) = 2
        psize(3) = 2 
      else if (ndim==1) then
        qsize(1) = 4
        psize(1) = 4
      else
        write (0,"(a)") "Using a grid of neither 3 nor 1 dimensions!"
        errorflag=1
        return
      end if

      nbf = product(qsize(1:ndim))*product(psize(1:ndim))+1
      if ((mod(qsizez,2)==0).and.(sys.ne."CP")) nbf = nbf + 1   ! Adds central point at end   
      allocate (coordinates(nbf,ndim*2), stat=ierr)
      if (ierr==0) allocate (columnsize(ndim*2), stat=ierr)
      if (ierr==0) allocate (columnrep(ndim*2), stat=ierr)
      if (ierr/=0) then
        write (0,"(a)") "Error in deallocation of column and coordinate arrays in reloc"
        errorflag=1
        return
      end if
      coordinates(:,:) = 0
      
      ! Once the number of basis functions has been changed, the basis set is
      ! de- and re-allocated to the correct new size
      call deallocbs(bs)
      call allocbs(bs,nbf)
      
      ! A set of coordinates are set up to convert the grid dimensions into a single
      ! array of values to be used to assign the correct positions in phase space to
      ! each member of the new smaller initial basis set, with columns assigned in 
      ! order of dimension and in (q,p) pairs, ie q_x, p_x, q_y, p_y, q_z, p_z.
      ! The columnreps array determines how many times the coordiantes for a
      ! particular dimension are repeated over the entire coordinates array in
      ! order to list all coordinates systematically
      n=0
      do m=1,ndim*2
        if (mod(m,2)==1) then
          n=n+1
          columnsize(m)=qsize(n)
        else
          columnsize(m)=psize(n)
        end if
        columnrep(m) = 1
        if (m/=1) then
          do l=1,m-1
            columnrep(m) = columnrep(m)*columnsize(l)
          end do
        end if
      end do
      
      ! Write to the coordinates array all the coordinates in order.
      ! These are the grid point numbers, not the values themselves, with each
      ! block in each column listing from 1 to the size of that particular dimension
      do m=1,ndim*2
        z=nbf/(columnrep(m)*columnsize(m))
        y=nbf/columnrep(m)     
        do n=1,nbf
          if (n==1) then
            coordinates(n,m) = 1
          else if (mod(n-1,z)==0) then 
            if (mod(n-1,y)==0) then
              coordinates(n,m) = 1
            else
              coordinates(n,m) = coordinates(n-1,m) + 1
            end if
          else
            coordinates(n,m) = coordinates(n-1,m)
          end if
        end do
      end do
      
      ! Calculate the bottom left corner of the grid in actual position
      ! and momentum coordinates
      do m=1,ndim
         qstrt(m) = muq(m)-((dble(qsize(m))-1.0d0)*sigq*gridsp/2.0d0)
         pstrt(m) = mup(m)-((dble(psize(m))-1.0d0)*sigp*gridsp/2.0d0)
      end do
      
      ! Uses the coordinates array to calculate the positions in phase space of each
      ! basis function, and writes that position along with the initial values for 
      ! the other quantities to the basis set array (D is calculated later)
      do k=1,nbf
        bf=bs(k)
        do m=1,ndim
          realz = qstrt(m)+(gridsp*((coordinates(k,(2*m)-1)-1)*sigq))
          imagz = pstrt(m)+(gridsp*((coordinates(k,(2*m)  )-1)*sigp))
          bf%z(m)=cmplx(realz, imagz,kind=8)
          bf%s_pes(1) = 0.0d0
          bf%d_pes(1) = (1.0d0,0.0d0)
          bf%a_pes(1) = (1.0d0,0.0d0)
          bf%D_big = (0.0d0,0.0d0)
        end do
        bs(k)=bf
      end do 
      
      ! Adds the central point of the grid to the basis set if needed.
      ! It should be noted that this will not work for the coulomb potential
      ! as the centre of the grid is over the singularity, hence why this
      ! possibility is not allowed by the condional
      if ((mod(qsizez,2)==0).and.(sys.ne."CP")) then
        do m=1,ndim
          bf%z(m)=cmplx(muq(m),mup(m),kind=8)
          bf%s_pes(1) = 0.0d0
          bf%d_pes(1) = (1.0d0,0.0d0)
          bf%a_pes(1) = (1.0d0,0.0d0)
          bf%D_big = (1.0d0,0.0d0)
        end do
        bs(nbf)=bf
      end if
      
      deallocate (qstrt, stat=ierr)
      if (ierr==0) deallocate (pstrt, stat=ierr)
      if (ierr==0) deallocate (qsize, stat=ierr)
      if (ierr==0) deallocate (psize, stat=ierr)
      if (ierr==0) deallocate (coordinates, stat=ierr)
      if (ierr==0) deallocate (columnsize, stat=ierr)
      if (ierr==0) deallocate (columnrep, stat=ierr)
      if (ierr/=0) then
        write (0,"(a,a)") "Error in deallocation of starting p and q values ", &
                               "and column and coordinate arrays in reloc"
        errorflag=1
        return
      end if
      
    end if

    !!!!!!!! Beginning of the reprojection section!!!!!!!

    allocate (zt(nbf,ndim), stat=ierr)
    if (ierr==0) allocate (bigA(nbf), stat=ierr)
    if (ierr==0) allocate (cnew(in_nbf), stat=ierr)
    if (ierr==0) allocate (ovrlp_mat(in_nbf,nbf), stat=ierr)
    if (ierr/=0) then
      write (0,"(a)") "Error in allocating temporary z values for reprojection subroutine"
      errorflag=1
      return
    end if 

    ! Build the arrays containing all needed basis set information
    do j=1,nbf
      zt(j,1:ndim)=bs(j)%z(1:ndim)
      bigA(j) = bs(j)%D_big * bs(j)%d_pes(1) * exp(i*bs(j)%s_pes(1))
    end do

    call deallocbs(bs)

    ! Ovrlp_mat is the overlap with the initial wavepacket
    do j=1,nbf
      do k=1,in_nbf
        cnew(k) = (0.0d0, 0.0d0)
        ovrlp_mat(k,j) = ovrlpij(z0(k,:), zt(j,:)) 
      end do
    end do

    cnew = matmul(ovrlp_mat,bigA)

    deallocate(ovrlp_mat, stat=ierr)
    if (ierr==0) deallocate(bigA, stat=ierr)
    if (ierr==0) deallocate(zt, stat=ierr)
    if (ierr==0) allocate(remove(in_nbf), stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in allocation of remove array or deallocation of first", &
               " overlap matrix in reprojection subroutine"
      errorflag = 1
      return
    end if
 
    remove=0   ! Array of flags to say whether or not a grid point needs removing

    if ((nbfadapt.eq."YES")) then 
      do k=1,in_nbf
        if (abs(cnew(k)).lt.bfeps) then
          remove(k) = 1
        else
          remove(k) = 0
        end if
      end do
      nbf = in_nbf - sum(remove(1:in_nbf))
      write (6,"(i0,a,i0,a,i0)") sum(remove(1:in_nbf)), " of ", in_nbf, &
                " bfs were removed. nbf now = ", nbf
    end if

    allocate (C_k(nbf), stat = ierr)  
    if (ierr==0) allocate (ovrlp_mat(nbf,nbf), stat=ierr)
    if (ierr/=0) then
      write (0,"(a)") "Error in allocation of overlap matrix or C_k in reprojection sub"
      errorflag=1
      return
    end if

    p=1
    do k=1,in_nbf
      if (remove(k)==0) then
        C_k(p) = cnew(k)
        q=1
        do j=1,in_nbf
          if (remove(j)==0) then
            ovrlp_mat(q,p) = ovrlpij(z0(j,:), z0(k,:))
            q=q+1
          end if
        end do           
        p=p+1
      end if
    end do

    allocate (znew(nbf,ndim))

    p=1
    do k=1,in_nbf
      if (remove(k)==0) then
        znew(p,:) = z0(k,:)          
        p=p+1
      end if
    end do

    deallocate (cnew, stat = ierr)
    if (ierr==0) allocate (cnew(nbf), stat=ierr)
    if (ierr==0) allocate (dnew(nbf), stat=ierr)
    if (ierr==0) deallocate (znew, stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in de- and re-allocation of cnew or in allocation of dnew", &
                   " in reprojection subroutine"
      errorflag=1
      return
    end if    

    cnew = C_k

    if (matfun.eq.'zgesv') then
      call lineq(ovrlp_mat, cnew, dnew)
    else if (matfun.eq.'zheev') then
      call matinv2(ovrlp_mat, cnew, dnew)
    else
      write(0,"(a)") "Error! Matrix function not recognised! Value is ", matfun
      errorflag = 1
      return
    end if

    deallocate(ovrlp_mat, stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error in allocation of Linear Algebra check array or deallocation",&
                  " of input overlap"
      errorflag=1
      return
    end if

    call allocbs(bs, nbf)
    allocate (DC(in_nbf), stat=ierr)
    if (ierr==0) deallocate (cnew, stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error allocating DC array or deallocating Btemp, cnew or chkmat ",& 
                    "in relocation subroutine"
      errorflag = 1
      return
    end if

    do j=1,in_nbf
      DC(j) = -1.0d0
    end do  

    j=1
    do k=1,in_nbf
      if (remove(k)==0) then
        bs(j)%D_big = dnew(j)
        bs(j)%s_pes(1) = 0.0d0
        bs(j)%z(1:ndim) = z0(k,1:ndim)
        bs(j)%a_pes(1) = (1.0d0,0.0d0)
        bs(j)%d_pes(1) = (1.0d0,0.0d0)
        DC(k) = abs(dble(dconjg(dnew(j))*C_k(j)))
        j=j+1
      end if
    end do

    if (j-1/=nbf) then
      write (0,"(a)") "Error! Mismatch in (new) basis set size!"
      errorflag = 1
      return
    end if

    if (debug==1) call graphwavefn(DC, x-1, z0)

    deallocate(remove, stat=ierr)
    if (ierr==0) deallocate(dnew, stat=ierr)
    if (ierr==0) deallocate(C_k, stat=ierr)
    if (ierr==0) deallocate(DC, stat=ierr)
    if (ierr/=0) then
      write (0,"(a)"),"Error deallocating remove, dnew, C_k or DC arrays in reloc"
      errorflag = 1
      return
    end if

    if (nbfadapt=="YES") then
      open(unit=4532,file="nbf.dat",status="old",access="append",iostat=ierr)
      if (ierr/=0) then
        write (0,"(a)"),"Error opening nbf file in reloc"
        errorflag = 1
        return
      end if
      write(4532,'(e12.5,i5)') time, nbf
      close(4532)
!      if (x==1) then
!        open(unit=1,file="trajat1.dat")
!        do k=1,nbf
!          write(1,*) k, (bs(k)%z(m), m=1,ndim)
!        end do
!        close(1)
!      end if 
    end if

    return

  end subroutine reloc_basis

!--------------------------------------------------------------------------------------------------

  subroutine leaking(bs, nbf, x)

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    integer, intent(inout) :: nbf, x  
    complex(kind=8), dimension (:,:), allocatable :: d_pes, zqp, a_pes
    complex(kind=8), dimension(:), allocatable :: D_k
    real(kind=8), dimension(:,:), allocatable :: s_pes
    integer, dimension (:), allocatable :: remove
    integer :: k, j, m, r, nbfnew, ierr

    if (errorflag/=0) return

    if (size(bs).ne.nbf) then
      write (0,"(a)") "Error with nbf. Does not match the size of the basis set"
      errorflag = 1
      return
    end if  

    allocate (remove(nbf), stat=ierr)
    if (ierr==0) allocate (zqp(nbf,ndim), stat=ierr)
    if (ierr==0) allocate (d_pes(nbf,npes), stat=ierr)
    if (ierr==0) allocate (a_pes(nbf,npes), stat=ierr)
    if (ierr==0) allocate (s_pes(nbf,npes), stat=ierr)
    if (ierr==0) allocate (D_k(nbf), stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error allocating the basis set variables in wavefunction leaking",&
                  " subroutine"
      errorflag = 1
      return
    end if

    remove = 0

    do k=1,nbf
      do m=1,ndim
        zqp(k,m) = bs(k)%z(m)
      end do
      do r=1,npes
        d_pes(k,r) = bs(k)%d_pes(r)
        a_pes(k,r) = bs(k)%a_pes(r)
        s_pes(k,r) = bs(k)%s_pes(r)
      end do
      D_k(k) = bs(k)%D_big
      if (sqrt(sum(dble(zqp(k,1:ndim)))**2.0d0).gt.40.0d0) then
        remove(k) = 1
      end if
    end do

    nbfnew = nbf - sum(remove(1:nbf))
    if (sum(remove(1:nbf)).ne.0) then
      write (0,"(i0,a,i0,a,i0,a,i0)") sum(remove(1:nbf)), " of ", nbf, &
          " bfs were removed at step ", x, ". nbf now = ", nbfnew
    end if

    call deallocbs(bs)
    call allocbs(bs,nbfnew)

    j=1
    do k=1,nbf
      if (remove(k)==0) then
        bs(j)%D_big = D_k(k)
        bs(j)%s_pes(1:npes) = s_pes(k,1:npes)
        bs(j)%z(1:ndim) = zqp(k,1:ndim)
        bs(j)%a_pes(1:npes) = a_pes(k,1:npes)
        bs(j)%d_pes(1:npes) = d_pes(k,1:npes)
        j=j+1
      end if
    end do  

    nbf = nbfnew  
        
    deallocate (remove, stat=ierr)
    if (ierr==0) deallocate (zqp, stat=ierr)
    if (ierr==0) deallocate (d_pes, stat=ierr)
    if (ierr==0) deallocate (a_pes, stat=ierr)
    if (ierr==0) deallocate (s_pes, stat=ierr)
    if (ierr==0) deallocate (D_k, stat=ierr)
    if (ierr/=0) then
      write (0,"(a,a)") "Error deallocating the basis set variables in wavefunction ",& 
                 "leaking subroutine"
      errorflag = 1
      return
    end if

    return

  end subroutine leaking

!***********************************************************************************!

  subroutine cloning(bs,nbf,x,time,clone, clonenum, reps)
  
    !NOTE: Cloning is only set up for 2 PESs. This needs to be generalised!

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    type(basisfn), dimension(:), allocatable :: bsnew
    real(kind=8), intent(in) :: time
    integer, dimension(:), allocatable, intent(inout) :: clone, clonenum
    integer, intent (inout) :: nbf
    integer, intent (in) :: x, reps
    real(kind=8) :: brforce, normar
    integer, dimension(:), allocatable :: clonehere, clonecopy, clonecopy2 
    integer :: k, m, j, n, nbfnew, ierr, r, clonetype
    character(LEN=3)::rep

    if (errorflag==1) return

    if ((cloneflg=="YES").or.((cloneflg=="BLIND+").and.(x.ne.0))) then
      clonetype = 1 ! 1=conditional cloning, 2=blind cloning
    else if ((cloneflg=="BLIND").or.((cloneflg=="BLIND+").and.(x.eq.0))) then
      clonetype = 2
    else
      write (0,"(2a)") "Cloneflg is invalid. Should be 'YES' or 'BLIND', but was ", cloneflg
      errorflag = 1
      return
    end if

    allocate (clonehere(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy2(nbf), stat=ierr)
    if (ierr/=0) then
      write (0,"(a)") "Error allocating the clonehere array"
      errorflag = 1
      return
    end if
    
    if (npes.ne.2) then
      write(6,*) "Error. Cloning currently only valid for npes=2"
      errorflag = 1
      return
    end if

    do k=1,nbf
      if (clone(k).lt.x-clonefreq) clone(k)=0
      clonehere(k) = 0
    end do

    if (clonetype==1) then
      do k=1,nbf
        normar = 0.0d0
        do r=1,npes
          normar = normar + dconjg(bs(k)%a_pes(r))*bs(k)%a_pes(r)
        end do
        brforce = ((abs(bs(k)%a_pes(1)*bs(k)%a_pes(2))**2.0)/(normar**2.0))
        if ((brforce.gt.0.249).and.(clone(k)==0).and.(clonenum(k).lt.clonemax)) then
          clone(k) = x
          clonehere(k) = 1
        end if 
        clonecopy(k) = clone(k)
        clonecopy2(k) = clonenum(k)
      end do
    else if (clonetype==2) then
      if (mod(x,clonefreq)==0) then
        do k=1,nbf
          if (clonenum(k).lt.clonemax) then
            clone(k) = x
            clonehere(k) = 1
          end if 
        end do 
      end if
      do k=1,nbf
        clonecopy(k) = clone(k)
        clonecopy2(k) = clonenum(k)
      end do
    end if          

    nbfnew = nbf + sum(clonehere(:))
    
    deallocate (clone, stat=ierr)
    if (ierr==0) deallocate (clonenum, stat=ierr)
    if (ierr==0) allocate (clone(nbfnew), stat=ierr)
    if (ierr==0) allocate (clonenum(nbfnew), stat=ierr)
    if (ierr/=0) then
      write (0,"(a)") "Error in de- and re-allocation of clone arrays"
      errorflag = 1
      return
    end if
    do k=1,nbf
      clone(k) = clonecopy(k)
      clonenum(k) = clonecopy2(k)
    end do
    deallocate(clonecopy, stat=ierr)
    if (ierr==0) deallocate(clonecopy2, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating the cloning copy arrays"
      errorflag = 1
      return
    end if 

    if (nbfnew/=nbf) then

      call allocbs(bsnew, nbfnew)

      j=1
      
      write(rep,"(i3.3)") reps
      open(unit=47756,file="Clonetrack-"//trim(rep)//".out",status="old",access="append",iostat=ierr)
      
      do k=1,nbf
        do m=1,ndim
          bsnew(k)%z(m) = bs(k)%z(m)
        end do
        do r=1,npes
          bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
        end do
        if (clonehere(k) == 1) then
          clone(k) = x
          clone(nbf+j) = x
          clonenum(k) = clonenum(k) + 1
          clonenum(nbf+j) = clonenum(k)
          if (method=="MCEv1") then
            bsnew(k)%D_big = bs(k)%D_big
            bsnew(k)%d_pes(in_pes) = bs(k)%d_pes(in_pes)
          else
            bsnew(k)%D_big = bs(k)%D_big * abs(bs(k)%a_pes(in_pes))
            bsnew(k)%d_pes(in_pes) = bs(k)%d_pes(in_pes)/abs(bs(k)%a_pes(in_pes))
          end if
          do r=1,npes
            if (r.ne.in_pes) then
              bsnew(k)%d_pes(r) = (0.0d0,0.0d0)
              bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
            end if
          end do
          if (method=="MCEv1") then
            bsnew(nbf+j)%D_big = bs(k)%D_big
          else
            bsnew(nbf+j)%D_big = bs(k)%D_big * &
                                  sqrt(1.-(dconjg(bs(k)%a_pes(in_pes))*bs(k)%a_pes(in_pes)))
          end if
          bsnew(nbf+j)%d_pes(in_pes) = (0.0d0,0.0d0)         
          bsnew(nbf+j)%s_pes(in_pes) = bs(k)%s_pes(in_pes)   
          do r=1,npes
            if (r.ne.in_pes) then
              if (x.eq.0) then
                bsnew(nbf+j)%d_pes(r) = (1.0d0,0.0d0)
              else
                if (method=="MCEv1") then
                  bsnew(nbf+j)%d_pes(r) = bs(k)%d_pes(r)
                else
                  bsnew(nbf+j)%d_pes(r) = bs(k)%d_pes(r)/&
                                  sqrt(1.-(dconjg(bs(k)%a_pes(1))*bs(k)%a_pes(1)))
                end if                  
              end if
              bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
            end if
          end do
          do m=1,ndim
            bsnew(nbf+j)%z(m) = bs(k)%z(m)
          end do
          write(47756,"(3i5,2es25.17e3)") x, k, nbf+j, abs(bs(k)%a_pes(1)), abs(bs(k)%a_pes(2))
          j = j+1
        else
          bsnew(k)%D_big = bs(k)%D_big
          do r=1,npes
            bsnew(k)%d_pes(r) = bs(k)%d_pes(r)
          end do
        end if
      end do
      
      close (47756)
   
      call deallocbs(bs)
      call allocbs(bs, nbfnew)
   
      do k=1,nbfnew
        bs(k)%D_big = bsnew(k)%D_big
        do r=1,npes
          bs(k)%a_pes(r) = bsnew(k)%d_pes(r) * cdexp(i*bsnew(k)%s_pes(r))
          bs(k)%d_pes(r) = bsnew(k)%d_pes(r)
          bs(k)%s_pes(r) = bsnew(k)%s_pes(r)
        end do
        do m=1,ndim
          bs(k)%z(m) = bsnew(k)%z(m)
        end do
      end do
      
      call deallocbs(bsnew)          
   
      n = nbfnew-nbf
   
      if (n.ne.0) then
        do k=1,n
          clone(nbf+k) = x
        end do
      end if
   
      write (6,"(i0,a,i0,a,i0)") sum(clonehere(:)), " bfs cloned in step ", x, &
            ". nbf now = ", nbfnew
   
      nbf = nbfnew

    end if
    
    deallocate(clonehere, stat=ierr)
    if (ierr/=0) then
      write(0,"(a)") "Error deallocating the clone-here array"
      errorflag = 1
      return
    end if 

  end subroutine cloning   
  
!------------------------------------------------------------------------------------

  subroutine retrieveclon (dummybs, bset, reps, x, time, nbf, map_bfs)   ! Currently only works for 2pes.
                                                     ! and for swarm trains - no aimce for swarms
    implicit none
  
    type(basisfn), dimension (:), intent(inout), allocatable :: bset, dummybs 
    real(kind=8), intent(in) :: time
    integer, dimension(:,:), intent(in) :: map_bfs
    integer, intent(inout) :: nbf
    integer, intent (in) :: reps, x

    real(kind=8), dimension(ndim) ::dummyarr  
    real(kind=8) :: amp1, amp2
    integer, dimension(:), allocatable :: carriages
    integer :: k, l, m, r, ierr, step, inbf, finbf, stepfwd, p, q, locnbf
    character(LEN=3) :: rep 
    
    if (errorflag.ne.0) return
    
    write(rep,"(i3.3)") reps
    p=0
    
    locnbf = size(bset)
    dummyarr = 0.0d0
    
    if (locnbf.eq.size(dummybs)) then
      do k=1,locnbf
        do r=1,npes
          if ((abs(bset(k)%d_pes(r)) - abs(dummybs(k)%d_pes(r))).gt.1.0d-6) then
            write (0,*) "Difference in expected and calculated a values too large in step ", x, " bf ", k
            write (0,*) (bset(k)%a_pes(r)), (dummybs(k)%a_pes(r))
          end if
!          bset(k)%a_pes(r) = dummybs(k)%a_pes(r)
!          bset(k)%d_pes(r) = dummybs(k)%d_pes(r)
!          bset(k)%s_pes(r) = dummybs(k)%s_pes(r)
        end do
!        do m=1,ndim
!          bset(k)%z(m) = dummybs(k)%z(m)
!        end do
      end do
    else
      do k=1,locnbf
        dummybs(k)%D_big = bset(k)%D_big
      end do
      open (unit=354+reps, file="Clonetrack-"//rep//".out", status="old",iostat=ierr)
      if (ierr/=0) then
        write(0,"(3a,i0)") "Error opening Clonetrack-", rep, ".out file. ierr was ", ierr
        errorflag = 1
        return
      end if
      allocate (carriages(def_stp), stat = ierr)
      if (ierr/=0) then
        write(0,"(2(a,i0))") "Error allocating carriages array for step ", x, ". ierr was ", ierr
        errorflag = 1
        return
      end if
            
      stepfwd = x + ((def_stp-1)/2)*trainsp
      do q=1,def_stp
        carriages(q) = stepfwd - (q-1)*trainsp
      end do
      
      do while (ierr==0)
        read (354+reps,"(3i5,2es25.17e3)", iostat=ierr) step, inbf, finbf, amp1, amp2
        do q=1,def_stp
          if ((step==carriages(q)).and.(ierr==0)) then
            p=p+1
            k=map_bfs(q,inbf)
            l=map_bfs(q,finbf)
            dummybs(k)%D_big = bset(k)%D_big * amp1
            dummybs(l)%D_big = bset(k)%D_big * amp2
            write (6,"(3(a,i4))") "Cloned basis set ", k, " to basis set ", l, " in step ", x
            if (dble(dummybs(l)%a_pes(2))-(dble(bset(k)%a_pes(2))/amp2).gt.1.0d-6) then
              write (0,"(4(a,es25.17e3),a)") "read a_pes values     : (", dble(dummybs(k)%a_pes(1)),",", &
                         dimag(dummybs(k)%a_pes(1)),") , (", dble(dummybs(k)%a_pes(2)),",",dimag(dummybs(k)%a_pes(2)),")"
              write (0,"(4(a,es25.17e3),a)") "                        (", dble(dummybs(l)%a_pes(1)),",", &
                               dimag(dummybs(l)%a_pes(1)),") , (", dble(dummybs(l)%a_pes(2)),",",&
                               dimag(dummybs(l)%a_pes(2)),")"
              write (0,"(4(a,es25.17e3),a)") "expected a_pes values : (", dble(bset(k)%a_pes(1))/amp1,",", &
                         dimag(bset(k)%a_pes(1))/amp1,") , (", 0.0d0,",",0.0d0,")"
              write (0,"(4(a,es25.17e3),a)") "                        (", 0.0d0,",", &
                               0.0d0,") , (", dble(bset(k)%a_pes(2))/amp2,",",&
                               dimag(bset(k)%a_pes(2))/amp2,")"
              call outbs(dummybs,reps,dummyarr, dummyarr, time,x,0)
              errorflag = 1
              return
            end if
          end if
        end do
      end do
      close (354+reps)
      deallocate (carriages, stat = ierr)
      if (ierr/=0) then
        write(0,"(2(a,i0))") "Error deallocating carriages array for step ", x, ". ierr was ", ierr
        errorflag = 1
        return
      end if 
      
      if (p==0) then
        write (0,"(3a)") "Expected cloning, but no cloning was found in the Cloningtrack-", rep,".out file"
        write (0,"(a,i0)") "This should not have happened. Stepfwd was ", stepfwd
        write (0,"(2(a,i0))") "Size of bset was ", locnbf, " and size of dummybs was ", size(dummybs)
        do k=1,size(carriages)
          write (0,*) carriages(k)
        end do
        errorflag = 1
        return
      else if (size(dummybs)-locnbf.ne.p) then
        write(0,"(2a,i0)")"The difference between the basis set sizes was not equal to the number of ",&
                    "cloning events found for step ", x+stepfwd
        write(0,"(2(a,i0))") "Difference was ",size(dummybs)-locnbf, " and number of cloning events was ", p 
        errorflag = 1
        return
      end if
      
      call deallocbs(bset)
      call allocbs(bset,size(dummybs))
      do k=1,size(dummybs)
        bset(k)%D_big = dummybs(k)%D_big
        do r=1,npes
          bset(k)%a_pes(r) = dummybs(k)%a_pes(r)
          bset(k)%d_pes(r) = dummybs(k)%d_pes(r)
          bset(k)%s_pes(r) = dummybs(k)%s_pes(r)
        end do
        do m=1,ndim
          bset(k)%z(m) = dummybs(k)%z(m)
        end do
      end do
      
    end if 
    
    nbf = size(dummybs)
    call deallocbs(dummybs)       
          
    return
    
  end subroutine retrieveclon

!***********************************************************************************!
end module bsetalter
