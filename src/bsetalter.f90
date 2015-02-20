MODULE bsetalter
  
  use globvars
  use Ham
  use alarrays
  use outputs

!*************************************************************************************************!
!*
!*         Basis Set alteration module Module
!*           
!*   Contains subroutines for:
!*
!*      1) Relocating the wavefunction on a static grid
!*      2) Leaking the basis set, used for large q values in the Hennon-Heiles model
!*      3) Cloning subroutine, which increases the basis set, used mainly for the spin boson model
!*      
!*************************************************************************************************!

contains


!--------------------------------------------------------------------------------------------------

  subroutine reloc_basis(bs, z0, nbf, x, time, gridsp, mup, muq)

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    complex(kind=8), dimension(:,:), intent(in) :: z0
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: time, gridsp
    integer, intent(inout) :: nbf, x
    type(basisfn) :: bf
    complex(kind=8), dimension (:), allocatable :: bigA, cnew, dnew, dnew2, Btemp, C_k
    complex(kind=8), dimension (:,:), allocatable :: ovrlp_mat, chkmat, zt, znew
    complex(kind=8) :: nrm
    real(kind=8), dimension (:), allocatable :: DC, qstrt, pstrt
    real(kind=8) :: absB, realz, imagz
    integer, dimension (:,:), allocatable :: coordinates
    integer, dimension (:), allocatable :: qsize, psize, columnsize, columnrep, remove
    integer :: j, k, p, q, l, m, n, z, y, ierr, nbfnew, nbfgrd
    character (LEN=15) :: filenm
    character (LEN=5)  :: step

    if (errorflag .ne. 0) return

    if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
       print *, "Basis Set relocation called, but initial basis set was not a regular grid!"
       errorflag=1
       return
    end if

    !!!!!!z0 is the initial grid

    !!!!!!! Create a smaller grid for the first step, max size 1600 bfs (ie 10x10x2x2x2x2).

    if ((x==1).and.(basis=="GRID")) then

       allocate (qstrt(ndim), stat=ierr)
       if (ierr==0) allocate (pstrt(ndim), stat=ierr)
       if (ierr==0) allocate (qsize(ndim), stat=ierr)
       if (ierr==0) allocate (psize(ndim), stat=ierr)
       if (ierr/=0) then
          print *, "Error in allocation starting p and q values for grid"
          errorflag=1
          return
       end if
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
          print *, "Using a grid of neither 3 nor 1 dimensions! How did that happen?"
          errorflag=1
          return
      end if
       nbf = product(qsize(1:ndim))*product(psize(1:ndim))
       allocate (coordinates(nbf,ndim*2), stat=ierr)
       if (ierr==0) allocate (columnsize(ndim*2), stat=ierr)
       if (ierr==0) allocate (columnrep(ndim*2), stat=ierr)
       if (ierr/=0) then
          print *, "Error in deallocation of column and coordinate arrays in reloc"
          errorflag=1
          return
       end if
       call deallocbs(bs)
       call allocbs(bs,nbf)
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
       do m=1,ndim*2
          columnrep(m) = 1
          if (m/=1) then
             do l=1,m-1
                columnrep(m) = columnrep(m)*columnsize(l)
             end do
          end if
       end do
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
       do m=1,ndim
          qstrt(m) = muq(m)-((dble(qsize(m))-1.0d0)*sigq*gridsp/2.0d0)
          pstrt(m) = mup(m)-((dble(psize(m))-1.0d0)*sigp*gridsp/2.0d0)
       end do
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
       if (mod(qsizez,2)==0) then
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
          print *, "Error in deallocation of starting p and q values and column and coordinate arrays in reloc"
          errorflag=1
          return
       end if
    end if

    allocate (zt(nbf,ndim), stat=ierr)
    if (ierr==0) allocate (bigA(nbf), stat=ierr)
    if (ierr==0) allocate (cnew(in_nbf), stat=ierr)
    if (ierr==0) allocate (ovrlp_mat(in_nbf,nbf), stat=ierr)
    if (ierr/=0) then
       print *, "Error in allocating temporary z values for reprojection subroutine"
       errorflag=1
       return
    end if 

    do j=1,nbf
       zt(j,1:ndim)=bs(j)%z(1:ndim)
       bigA(j) = bs(j)%D_big * bs(j)%d_pes(1) * exp(i*bs(j)%s_pes(1))
    end do

    call deallocbs(bs)

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
       print *, "Error in allocation of remove array or deallocation of first overlap matrix in reprojection subroutine"
       errorflag = 1
       return
    end if

    remove=0

    if ((nbfadapt.eq."YES")) then 
       do k=1,in_nbf
          if (abs(cnew(k)).lt.bfeps) then
             remove(k) = 1
          else
             remove(k) = 0
          end if
       end do
       nbf = in_nbf - sum(remove(1:in_nbf))
       print "(i0,a,i0,a,i0)", sum(remove(1:in_nbf)), " of ", in_nbf, " bfs were removed. nbf now = ", nbf
    end if

    allocate (C_k(nbf), stat = ierr)  
    if (ierr==0) allocate (ovrlp_mat(nbf,nbf), stat=ierr)
    if (ierr/=0) then
       print *, "Error in allocation of overlap matrix or C_k in reprojection subroutine"
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
       print *, "Error in de- and re-allocation of cnew or allocation of dnew in reprojection subroutine"
       errorflag=1
       return
    end if    

    cnew = C_k

    if (matfun.eq.'zgesv') then
       call lineq(ovrlp_mat, cnew, dnew)
    else if (matfun.eq.'zheev') then
       call matinv2(ovrlp_mat, cnew, dnew)
    else
       print *, "Error! Matrix function not recognised! Value is ", matfun
       errorflag = 1
       return
    end if

    deallocate(ovrlp_mat, stat=ierr)
    if (ierr/=0) then
       print *, "Error in allocation of Linear Algebra check array or deallocation of input overlap"
       errorflag=1
       return
    end if

    call allocbs(bs, nbf)
    allocate (DC(in_nbf), stat=ierr)
    if (ierr==0) deallocate (cnew, stat=ierr)
    if (ierr/=0) then
       print *, "Error allocating DC array or deallocating Btemp, cnew or chkmat in reloc"
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

    if (x==1) then
       do k=1,in_nbf
          if (abs(C_k(k)).lt.bfeps) then
             DC(k) = -1.0d0
          end if
       end do
    end if
      
    if (j-1/=nbf) then
       print *, "Error! Mismatch in (new) basis set size!"
       errorflag = 1
       return
    end if

    if (debug==1) call graphwavefn(DC, x-1, z0)

    deallocate(remove, stat=ierr)
    if (ierr==0) deallocate(dnew, stat=ierr)
    if (ierr==0) deallocate(C_k, stat=ierr)
    if (ierr==0) deallocate(DC, stat=ierr)
    if (ierr/=0) then
       print *,"Error deallocating remove, dnew, C_k or DC arrays in reloc"
       errorflag = 1
       return
    end if

    if (nbfadapt=="YES") then
       open(unit=4532,file="nbf.dat",status="old",access="append",iostat=ierr)
       if (ierr/=0) then
          print *,"Error opening nbf file in reloc"
          errorflag = 1
          return
       end if
       write(4532,'(e12.5,i5)') time, nbf
       close(4532)
       if (x==1) then
          open(unit=1,file="trajat1.dat")
          do k=1,nbf
             write(1,*) k, (bs(k)%z(m), m=1,ndim)
          end do
          close(1)
       end if 
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
       print *, "Error with nbf. Does not match the size of the basis set"
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
       print *, "Error allocating the basis set variables in wavefunction leaking subroutine"
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
       print "(i0,a,i0,a,i0,a,i0)", sum(remove(1:nbf)), " of ", nbf, " bfs were removed at step ", x, ". nbf now = ", nbfnew
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
       print *, "Error deallocating the basis set variables in wavefunction leaking subroutine"
       errorflag = 1
       return
    end if

    return

  end subroutine leaking

!*************************************************************************************************!

  subroutine cloning(bs,nbf,x,time,clone)

    implicit none

    type(basisfn), dimension(:), allocatable, intent(inout) :: bs
    type(basisfn), dimension(:), allocatable :: bsnew
    real(kind=8), intent(in) :: time
    integer, dimension(:), allocatable, intent(inout) :: clone
    integer, intent (inout) :: nbf
    integer, intent (in) :: x
    real(kind=8) :: brforce
    integer, dimension(:), allocatable :: clonehere, clonecopy 
    integer :: k, m, j, n, nbfnew, ierr, r, clonetype

    if (errorflag==1) return

    clonetype = 1 ! 1=conditional cloning, 2=blind cloning

    allocate (clonehere(nbf), stat=ierr)
    if (ierr==0) allocate(clonecopy(nbf), stat=ierr)
    if (ierr/=0) then
       print *, "Error allocating the clonehere array"
       errorflag = 1
       return
    end if

    do k=1,nbf
       if (clone(k).lt.x-50) clone(k)=0
       clonehere(k) = 0
    end do

    if (clonetype==1) then
       do k=1,nbf
          brforce = ((abs((bs(k)%a_pes(1))*abs(bs(k)%a_pes(2))))**2.0)
          if ((brforce.gt.0.249).and.(clone(k)==0)) then
             clone(k) = x
             clonehere(k) = 1
          end if 
          clonecopy(k) = clone(k)
       end do
    else if (clonetype==2) then
       if (mod(x,150)==0) then
          do k=1,nbf
             if (clone(k)==0) then
                clone(k) = x
                clonehere(k) = 1
             end if 
             clonecopy(k) = clone(k)
          end do 
       end if
    end if          

    nbfnew = nbf + sum(clonehere(:))

    if (nbfnew/=nbf) then

       call allocbs(bsnew, nbfnew)

       j=1

       do k=1,nbf
          do m=1,ndim
             bsnew(k)%z(m) = bs(k)%z(m)
          end do
          do r=1,npes
             bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
          end do
          if (clonehere(k) == 1) then
             bsnew(k)%D_big = bs(k)%D_big * abs(bs(k)%a_pes(1))
             bsnew(k)%d_pes(1) = bs(k)%d_pes(1)/abs(bs(k)%a_pes(1))
             do r=2,npes
                bsnew(k)%d_pes(r) = (0.0d0,0.0d0)
                bsnew(k)%s_pes(r) = bs(k)%s_pes(r)
             end do
             bsnew(nbf+j)%D_big = bs(k)%D_big * sqrt(1-(dconjg(bs(k)%a_pes(1)*bs(k)%a_pes(1))))
             bsnew(nbf+j)%d_pes(1) = (0.0d0,0.0d0)
             bsnew(nbf+j)%s_pes(1) = bs(k)%s_pes(1)
             do r=2,npes
                bsnew(nbf+j)%d_pes(r) = bs(k)%d_pes(r)/sqrt(1-(dconjg(bs(k)%a_pes(1)*bs(k)%a_pes(1))))
                bsnew(nbf+j)%s_pes(r) = bs(k)%s_pes(r)
             end do
             do m=1,ndim
                bsnew(nbf+j)%z(m) = bs(k)%z(m)
             end do
             j = j+1
          else
             bsnew(k)%D_big = bs(k)%D_big
             do r=1,npes
                bsnew(k)%d_pes(r) = bs(k)%d_pes(r)
             end do
          end if
       end do
   
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
   
       deallocate (clone, stat=ierr)
       if (ierr==0) allocate (clone(nbfnew), stat=ierr)
       if (ierr/=0) then
          print *, "Error in de- and re-allocation of clone array"
          errorflag = 1
          return
       end if
   
       n = nbfnew-nbf
   
       do k=1,nbf
          clone(k) = clonecopy(k)
       end do
       if (n.ne.0) then
          do k=1,n
             clone(nbf+k) = x
          end do
       end if
   
       print "(1x,i0,a,i0,a,i0)", sum(clonehere(:)), " bfs cloned in step ", x, ". nbf now = ", nbfnew
   
       nbf = nbfnew

    end if

    deallocate(clonehere, stat=ierr)
    if (ierr==0) deallocate(clonecopy, stat=ierr)
    if (ierr/=0) then
       print *, "Error deallocating the cloning arrays"
       errorflag = 1
       return
    end if        

  end subroutine cloning     

!*************************************************************************************************!
end module bsetalter
