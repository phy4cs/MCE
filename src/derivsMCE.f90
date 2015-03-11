MODULE derivsMCE

	use globvars
	use Ham
	use alarrays
	use outputs
	use redirect

!***********************************************************************************!
!*
!*         MCE Derivatives Module
!*           
!*   Contains subroutines for:
!*
!*      1) Calling the time derivative subroutines for each variable and returning
!*        the results to the timestep subroutine
!*      2) Calculating the time derivative of the z values
!*      3) Calculating the derivative of the single configuration MCEv1 d prefactor
!*      4) Calculating the derivative of the single configuration MCEv2 d prefactor
!*      5) Calculating the time derivative of the classical action
!*      6) Calculating the derivative of the multi-configuration MCEv1 D prefactor
!*      7) Calculating the derivative of the multi-configuration MCEv2D prefactor
!*      8) Calculating the first component of the Hamiltonian used in 7)
!*      9) Calculating the second component of the Hamiltonian used in 7)
!*     10) Calculating the third component of the Hamiltonian used in 7) and in 3)
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

	subroutine deriv(bsin, bsout, x, time, genflg)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(basisfn), dimension (:), intent (inout) :: bsout
		real(kind=8), intent(in) :: time
		complex(kind=8), dimension(:,:), allocatable :: dz, dd
		real(kind=8), dimension(:,:), allocatable :: ds
		complex(kind=8), dimension(:), allocatable ::dD_big
		integer, intent(in) :: x, genflg
		integer :: k, m, r, ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate (dz(size(bsin),ndim), stat=ierr)
		if (ierr==0) allocate (dd(size(bsin),npes),stat=ierr)
		if (ierr==0) allocate (ds(size(bsin),npes),stat=ierr)
		if (ierr==0) allocate (dD_big(size(bsin)),stat=ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in derivatives array allocation in deriv"
			errorflag=1
			return
		end if

		dz=zdot(bsin,time)
		ds=sdot(bsin,dz,time)

		if (method=="MCEv1") then
			dd=ddotv1(bsin,time)
			dD_big=bigDdotv1(size(bsin))
		else if ((method=="MCEv2").or.(trim(method)=="CCS")) then
			dd=ddotv2(bsin,time)
			if (((basis=="TRAIN").or.(basis=="SWTRN")).and.(genflg==1)) then
				dD_big=bigDdotv1(size(bsin))
			else
				dD_big=bigDdotv2 (bsin,x,time)
			end if
		else
			write(0,"(a)"), "Error! Method unrecognised!"
			write(0,"(a)"), "How did you even get this far?"
			errorflag = 1
			return
		end if  

		if (errorflag .ne. 0) return     

		do k = 1,size(bsin)
			do m = 1,ndim
				bsout(k)%z(m) = dz(k,m)
			end do
			do r = 1,npes
				bsout(k)%d_pes(r) = dd(k,r)
				bsout(k)%s_pes(r) = ds(k,r)
			end do
			bsout(k)%D_big = dD_big(k)
		end do

		deallocate (dz,dd,ds,dD_big, stat=ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in derivatives array deallocation in deriv"
			errorflag=1
			return
		end if

		return

	end subroutine deriv    

!------------------------------------------------------------------------------------

	function zdot(bsin,t)

		implicit none

		complex(kind=8), dimension (:) ,allocatable:: zdottemp, z
		complex(kind=8), dimension (:,:,:),allocatable :: dhdz
		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension(size(bsin),ndim) :: zdot
		complex(kind=8), dimension(:),allocatable :: a, ac
		complex(kind=8) :: asum
		real(kind=8), intent (in) :: t
		integer :: k, m, r, s, ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate(zdottemp(ndim), stat = ierr)
		if (ierr==0) allocate (z(ndim), stat = ierr)
		if (ierr==0) allocate (dhdz(npes,npes,ndim), stat = ierr)
		if (ierr==0) allocate (a(npes), stat = ierr)
		if (ierr==0) allocate (ac(npes), stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in array allocation in zdot"
			errorflag=1
			return
		end if

		do k=1,size(bsin)
			zdottemp = (0.0d0, 0.0d0)
			z = bsin(k)%z
			do r=1,npes
				a(r) = bsin(k)%a_pes(r)
				ac(r) = dconjg(a(r))
			end do
			call dh_dz(dhdz, z, t)
			asum = (0.0d0, 0.0d0)
			do r=1,npes
				asum = asum + (ac(r)*a(r))
				do s=1,npes
					do m=1,ndim
						zdottemp(m) = zdottemp(m) + (dhdz(r,s,m)*ac(r)*a(s))
					end do
				end do
			end do
			do m=1,ndim
				zdot(k,m)=(-1.0d0*i*zdottemp(m))/asum
			end do
		end do

		deallocate(zdottemp,z,dhdz,a,ac, stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in array deallocation in zdot"
			errorflag=1
			return
		end if

		return

	end function zdot

!------------------------------------------------------------------------------------

	function ddotv1(bsin,time)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), allocatable :: H
		complex(kind=8), dimension (:,:),allocatable :: ovrlp, ovrlpin, d2H_temp
		complex(kind=8), dimension (:,:),allocatable :: zczd, ddot_temppes
		complex(kind=8), dimension (:,:,:),allocatable :: d2H
		complex(kind=8), dimension (:),allocatable :: atemp, ddot_temp, ddot_out, ddot_in
		complex(kind=8), dimension(size(bsin),npes) :: ddotv1
		real(kind=8), intent (in) :: time
		real(kind=8) :: absB
		integer :: j,k,r,s,nbf,ierr

		if (errorflag .ne. 0) return

		ierr = 0

		nbf = size(bsin)

		allocate(ovrlp(nbf,nbf),stat=ierr)
		if (ierr==0) allocate (ovrlpin(nbf,nbf),stat=ierr)
		if (ierr==0) allocate (zczd(nbf,nbf),stat=ierr)
		if (ierr==0) allocate (d2H_temp(nbf,nbf),stat=ierr)
		if (ierr==0) allocate (d2H(nbf,nbf,npes),stat=ierr)
		if (ierr==0) allocate (atemp(nbf),stat=ierr)
		if (ierr==0) allocate (ddot_temp(nbf),stat=ierr)
		if (ierr==0) allocate (ddot_out(nbf),stat=ierr)
		if (ierr==0) allocate (ddot_in(nbf),stat=ierr)
		if (ierr==0) allocate (ddot_temppes(nbf, npes),stat=ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of matrices or arrays in ddotv1"
			errorflag=1
			return
		end if

		ovrlp = ovrlpmat(bsin)

		call allocham(H,nbf)
		call Hord(bsin,H,time)

		zczd = zczdot(bsin, time)

		do r=1,npes
			do j=1,nbf
				do k=1,nbf
					if (j==k) then
						d2H(k,j,r) = (0.0d0,0.0d0)
					else
						d2H(k,j,r) = H(k,j)%Hjk(r,r) - H(j,j)%Hjk(r,r) - zczd(k,j)
					end if
				end do
			end do
		end do

		do r=1,npes
			do j=1,nbf
				do k=1,nbf
					if (k.ne.j) then
						d2H_temp(k,j) = d2H(k,j,r) * ovrlp(k,j)
					else
						d2H_temp(k,j) = (0.0d0,0.0d0)
					end if
				end do
				atemp(j) = bsin(j)%A_pes(r)
			end do
			ddot_temp = matmul(d2H_temp,atemp)
			do k=1,nbf
				ddot_temppes(k,r) = ddot_temp(k)
			end do
		end do
	 

		do r=1,npes
			ddot_temp = (0.0d0,0.0d0)
			do s=1,npes
				if (s.ne.r) then
					do k=1,nbf
						do j=1,nbf
							ddot_temp(k) = ddot_temp(k)+H(k,j)%Hjk(r,s)*ovrlp(k,j)*&
															bsin(j)%a_pes(s)
						end do
					end do
				end if
			end do
			do k=1,nbf
				ddot_temppes(k,r) = ddot_temppes(k,r) + ddot_temp(k)
			end do
		end do

		call deallocham(H)
 
		do r=1,npes   

			do k=1,nbf
				ddot_in(k) = ddot_temppes(k,r)
				do j=1,nbf
					ovrlpin(j,k) = ovrlp(j,k)    !!!! ovrlp is changed by the matinv subroutine
				end do                         !!!! so must be redefined for each pes
			end do

			if (matfun.eq.'zgesv') then
				call lineq (ovrlpin, ddot_in, ddot_out)
			else if (matfun.eq.'zheev') then
				call matinv2(ovrlpin, ddot_in, ddot_out)
			else
				write(0,"(a)"), "Error! Matrix function not recognised! Value is ", matfun
				errorflag = 1
				return
			end if 

			do k=1,nbf
				ddotv1(k,r)=-1.0d0*i*ddot_out(k)*cdexp(-1.0d0*i*(bsin(k)%s_pes(r)))
			end do
		
		end do

		deallocate(ovrlp,ovrlpin,zczd,d2H_temp,d2H,atemp,stat=ierr)
		if(ierr==0) deallocate(ddot_temp,ddot_out,ddot_in,ddot_temppes,stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in deallocation of matrices or arrays in ddotv1"
			errorflag=1
			return
		end if
 
		return
 
	end function ddotv1  

!------------------------------------------------------------------------------------

	function ddotv2(bsin,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		real(kind=8), intent (in) :: t
		complex(kind=8), dimension(:,:),allocatable :: Hkk
		complex(kind=8), dimension(:),allocatable :: dk
		complex(kind=8), dimension(size(bsin),npes) :: ddotv2
		complex(kind=8), dimension(:),allocatable :: z
		real(kind=8), dimension(:),allocatable :: Sk
		integer :: k, r, s, ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate(Hkk(npes,npes),dk(npes),z(ndim),Sk(npes),stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of matrices or arrays in ddotv2"
			errorflag=1
			return
		end if

		do r=1,size(ddotv2,2)
			do k=1,size(ddotv2,1)
				ddotv2(k,r) = (0.0d0, 0.0d0)
			end do
		end do

		do k=1,size(bsin)
			z = bsin(k)%z  
			call Hij(Hkk,z,z,t)
			dk = bsin(k)%d_pes
			Sk = bsin(k)%S_pes
			do r=1,npes
				do s=1,npes
					if (r.ne.s) then
						ddotv2(k,r) = ddotv2(k,r) + Hkk(r,s) * ovrlpij(z,z) * dk(s) &
						                                     * cdexp(i*(Sk(s)-Sk(r)))
					end if
				end do
			end do
		end do

		do r=1,size(ddotv2,2)
			do k=1,size(ddotv2,1)
				ddotv2(k,r) = -1.0d0*i*ddotv2(k,r)
			end do
		end do

		deallocate(Hkk,dk,z,Sk,stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in deallocation of matrices or arrays in ddotv2"
			errorflag=1
			return
		end if

		return

	end function ddotv2

!------------------------------------------------------------------------------------

	function sdot(bsin,dz,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension (:,:), intent(in) :: dz
		real(kind=8), intent (in) :: t
		real(kind=8), dimension (size(bsin),npes) :: sdot
		complex(kind=8), dimension(:), allocatable :: zk, zkc, zkdot, zkdotc
		complex(kind=8), dimension(:,:), allocatable :: Hkk
		integer :: k, r, m, ierr
		complex(kind=8) :: zsum

		if (errorflag .ne. 0) return
		
		ierr = 0

		allocate(zk(ndim),zkc(ndim),zkdot(ndim),zkdotc(ndim),Hkk(npes,npes),stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of matrices or arrays in sdot"
			errorflag=1
			return
		end if

		do k = 1,size(bsin)
			do m = 1,ndim
				zk(m) = bsin(k)%z(m)
				zkdot(m) = dz(k,m)
				zkc(m) = dconjg(zk(m))
				zkdotc(m) = dconjg(zkdot(m))
			end do
			call Hij(Hkk,zk,zk,t)
			zsum = (0.0d0, 0.0d0)
			do m = 1,ndim
				zsum = zsum + i*0.5d0*(zkdot(m)*zkc(m)-zkdotc(m)*zk(m))
			end do
			do r = 1,npes
				sdot(k,r) = zsum - dble(Hkk(r,r))
			end do
		end do

		deallocate(zk,zkc,zkdot,zkdotc,Hkk,stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of matrices or arrays in sdot"
			errorflag=1
			return
		end if

		return

	end function sdot

!------------------------------------------------------------------------------------

	function bigDdotv1(nbf)
 
		implicit none

		integer, intent (in) :: nbf
		complex(kind=8), dimension(nbf)::bigDdotv1
		integer :: j

		if (errorflag .ne. 0) return

		do j=1,nbf
			bigDdotv1 = (0.0d0,0.0d0)
		end do

		return

	end function bigDdotv1    

!------------------------------------------------------------------------------------

	function bigDdotv2(bsin,x, time)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin 
		type(hamiltonian), dimension (:,:), allocatable :: H
		complex(kind=8), dimension (:,:), allocatable :: ovrlp, ovrlpphi, ovrlpin
		complex(kind=8), dimension (:,:), allocatable :: h_av_jj, h_av_jk, zconjzdot
		complex(kind=8), dimension (:,:), allocatable :: d2H, d2H2, d2Hdiff
		complex(kind=8), dimension(:), allocatable::tempD, chk, Ddot, d2HD, Dtemp
		complex(kind=8), dimension(size(bsin)) :: bigDdotv2
		complex(kind=8) :: ovrlpdif
		real(kind=8), intent(in) :: time
		integer, intent(in) :: x
		integer :: j, k, nbf, ierr
		real(kind=8)::absB

		if (errorflag .ne. 0) return

		ierr = 0

		nbf = size(bsin)

		allocate (ovrlp(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (ovrlpphi(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (h_av_jj(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (h_av_jk(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (zconjzdot(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (d2H(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (d2H2(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (d2Hdiff(nbf,nbf),stat=ierr)
		if (ierr == 0) allocate (tempD(nbf),stat=ierr)
		if (ierr == 0) allocate (chk(nbf),stat=ierr)
		if (ierr == 0) allocate (Ddot(nbf),stat=ierr)
		if (ierr == 0) allocate (d2HD(nbf),stat=ierr)
		if (ierr == 0) allocate (Dtemp(nbf),stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of matrices or arrays in bigDdotv2"
			errorflag=1
			return
		end if    

		ovrlp = ovrlpmat(bsin)
		ovrlpphi = ovrlpphimat(bsin)

		if ((x.eq.1).and.(time.eq.0.0d0)) then
			do j=1,nbf
				do k=1,nbf
					ovrlpdif = ovrlpphi(j,k) - ovrlp(j,k) 
					if ((ovrlpdif.ne.0.0d0).and.(basis.ne."TRAIN").and.(basis.ne."SWTRN")) then
						write(0,"(a)"), "Error! Initial phi-overlap has disimilarilies to z-overlap"
						write(0,'(a,a,i0,a,i0,a)'), "These matricies should be identical ",&
																	"but differences found at coordinate ", j,",",k,"."
						write(0,'(a,4(e15.8,a))'), "Expected (",dimag(i*ovrlp(j,k)),","&
						                     ,dimag(ovrlp(j,k)),") but got (",&
																  dimag(i*ovrlpphi(j,k)),",",dimag(ovrlpphi(j,k)),")"
						errorflag = 1
					end if
				end do
			end do
			if (errorflag == 1) return
		end if

		call allocham(H, nbf)
		call Hord(bsin, H, time)

		h_av_jj = Hjk_avrg(H,ovrlp,bsin)
		h_av_jk = phiHphi(H,ovrlp,bsin)
		zconjzdot = zczdot(bsin,time)
		if (errorflag .ne. 0) return    
		do j=1,nbf
			do k=1,nbf
				zconjzdot(j,k) = zconjzdot(j,k)*ovrlpphi(j,k)
			end do
		end do

		do k=1,nbf
			Dtemp(k) = bsin(k)%D_big
			do j=1,nbf
				d2H(j,k) = h_av_jk(j,k) - h_av_jj(j,k) - zconjzdot(j,k)
				if (d2H(j,k)/=d2H(j,k)) then
					write(0,"(1x,a,i4,a,i4,a)"), "d2H(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		tempD = matmul(d2H, Dtemp)

		do k=1,nbf
			tempD(k)=-1.0d0*i*tempD(k)
		end do


		do j=1,nbf
			Ddot(j) = (0.0d0, 0.0d0)
			if (tempD(j)/=tempD(j)) then
				write(0,"(1x,a,i4,a)"), "tempD(", j, ") is NaN. Terminating"
				errorflag = 1
			return
			end if
		end do

		d2hD = tempD

		if (matfun.eq.'zgesv') then
			call lineq (ovrlpphi, tempD, Ddot)
		else if (matfun.eq.'zheev') then
			call matinv2(ovrlpphi, tempD, Ddot)
		else
			write(0,"(a)"), "Error! Matrix function not recognised! Value is ", matfun
			errorflag = 1
			return
		end if

		do k=1,nbf
			bigDdotv2(k)=Ddot(k)
		end do

	 call deallocham(H)

		deallocate (ovrlp,ovrlpphi,h_av_jj,h_av_jk,zconjzdot,d2H,d2H2,d2Hdiff,stat=ierr)
		if (ierr == 0 ) deallocate (tempD,chk,Ddot,d2HD,Dtemp,stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in deallocation of matrices or arrays in bigDdotv2"
			errorflag=1
			return
		end if   

		return

	end function bigDdotv2

!------------------------------------------------------------------------------------

	function phiHphi(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: phiHphi
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(phiHphi,2)
			do j=1,size(phiHphi,1)
				phiHphi(j,k) = (0.0d0, 0.0d0)
				do r=1,npes
					do s=1,npes
						phiHphi(j,k) = phiHphi(j,k) + (dconjg(bsin(j)%a_pes(r)) * ovrlp(j,k) * &
																				H(j,k)%Hjk(r,s) * bsin(k)%a_pes(s))
					end do
				end do
				if (phiHphi(j,k)/=phiHphi(j,k)) then
					write(0,"(1x,a,i4,a,i4,a)"), "phiHphi(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function phiHphi

!------------------------------------------------------------------------------------

	function Hjk_avrg(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: Hjk_avrg
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(Hjk_avrg,2)
			do j=1,size(Hjk_avrg,1)
				Hjk_avrg(j,k) = (0.0d0, 0.0d0)
				do r=1,npes
					do s=1,npes
						Hjk_avrg(j,k) = Hjk_avrg(j,k) + (dconjg(bsin(j)%a_pes(r)) * ovrlp(k,k)*&
																							H(k,k)%Hjk(r,s) * bsin(k)%a_pes(s))
					end do
				end do
				Hjk_avrg(j,k) = Hjk_avrg(j,k) * ovrlp(j,k)
				if (Hjk_avrg(j,k)/=Hjk_avrg(j,k)) then
					write(0,"(1x,a,i4,a,i4,a)"), "Hjk_avrg(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function Hjk_avrg

!------------------------------------------------------------------------------------

	function zczdot(bsin,t)

		type(basisfn), dimension (:), intent (in) :: bsin 
		complex(kind=8), dimension (size(bsin),size(bsin)) :: zczdot
		complex(kind=8), dimension (:,:), allocatable :: dz
		integer :: j, k, m, ierr
		real(kind=8), intent (in) :: t

		if (errorflag .ne. 0) return

		ierr = 0
		allocate(dz(size(bsin),ndim), stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in allocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		dz = zdot(bsin, t)

		do k=1,size(zczdot,2)
			do j=1,size(zczdot,1)
				zczdot(j,k) = (0.0d0, 0.0d0)
				do m=1,ndim
					zczdot(j,k) = zczdot(j,k) + ((dconjg(bsin(j)%z(m)-bsin(k)%z(m)))*dz(k,m))
				end do
				if (zczdot(j,k)/=zczdot(j,k)) then
					write(0,"(1x,a,i4,a,i4,a)"), "zczdot(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
				zczdot(j,k) = i*zczdot(j,k)
			end do
		end do

		deallocate(dz, stat = ierr)
		if (ierr/=0) then
			write(0,"(a)"), "Error in deallocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		return

	end function zczdot

!***********************************************************************************!

end module derivsMCE
