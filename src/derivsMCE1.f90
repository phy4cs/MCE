MODULE derivsMCE1

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
!*      1) Calculating the time derivative of the z values
!*      2) Calculating the derivative of the single configuration MCEv1 d prefactor
!*      3) Calculating the time derivative of the classical action
!*      4) Calculating the derivative of the multi-configuration MCEv1 D prefactor
!*      5) Calculating the third component of the Hamiltonian used in 2)
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

	function zdot_MCE1(bsin,t)

		implicit none

		complex(kind=8), dimension (:) ,allocatable:: zdottemp, z
		complex(kind=8), dimension (:,:,:),allocatable :: dhdz
		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension(size(bsin),ndim) :: zdot_MCE1
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
			print *, "Error in array allocation in zdot"
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
				zdot_MCE1(k,m)=(-1.0d0*i*zdottemp(m))/asum
			end do
		end do

		deallocate(zdottemp,z,dhdz,a,ac, stat = ierr)
		if (ierr/=0) then
			print *, "Error in array deallocation in zdot"
			errorflag=1
			return
		end if

		return

	end function zdot_MCE1

!------------------------------------------------------------------------------------

	function ddot_MCE1(bsin,time)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), allocatable :: H
		complex(kind=8), dimension (:,:),allocatable :: ovrlp, ovrlpin, d2H_temp
		complex(kind=8), dimension (:,:),allocatable :: zczd, ddot_temppes
		complex(kind=8), dimension (:,:,:),allocatable :: d2H
		complex(kind=8), dimension (:),allocatable :: atemp, ddot_temp, ddot_out, ddot_in
		complex(kind=8), dimension(size(bsin),npes) :: ddot_MCE1
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
			print *, "Error in allocation of matrices or arrays in ddotv1"
			errorflag=1
			return
		end if

		ovrlp = ovrlpmat(bsin)

		call allocham(H,nbf)
		call Hord(bsin,H,time)

		zczd = zczdot_MCE1(bsin, time)

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
				print *, "Error! Matrix function not recognised! Value is ", matfun
				errorflag = 1
				return
			end if 

			do k=1,nbf
				ddot_MCE1(k,r)=-1.0d0*i*ddot_out(k)*cdexp(-1.0d0*i*(bsin(k)%s_pes(r)))
			end do
		
		end do

		deallocate(ovrlp,ovrlpin,zczd,d2H_temp,d2H,atemp,stat=ierr)
		if(ierr==0) deallocate(ddot_temp,ddot_out,ddot_in,ddot_temppes,stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of matrices or arrays in ddotv1"
			errorflag=1
			return
		end if
 
		return
 
	end function ddot_MCE1  

!------------------------------------------------------------------------------------

	function sdot_MCE1(bsin,dz,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension (:,:), intent(in) :: dz
		real(kind=8), intent (in) :: t
		real(kind=8), dimension (size(bsin),npes) :: sdot_MCE1
		complex(kind=8), dimension(:), allocatable :: zk, zkc, zkdot, zkdotc
		complex(kind=8), dimension(:,:), allocatable :: Hkk
		integer :: k, r, m, ierr
		complex(kind=8) :: zsum

		if (errorflag .ne. 0) return
		
		ierr = 0

		allocate(zk(ndim),zkc(ndim),zkdot(ndim),zkdotc(ndim),Hkk(npes,npes),stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of matrices or arrays in sdot"
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
				sdot_MCE1(k,r) = zsum - dble(Hkk(r,r))
			end do
		end do

		deallocate(zk,zkc,zkdot,zkdotc,Hkk,stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of matrices or arrays in sdot"
			errorflag=1
			return
		end if

		return

	end function sdot_MCE1

!------------------------------------------------------------------------------------

	function bigDdot_MCE1(nbf)
 
		implicit none

		integer, intent (in) :: nbf
		complex(kind=8), dimension(nbf)::bigDdot_MCE1
		integer :: j

		if (errorflag .ne. 0) return

		do j=1,nbf
			bigDdot_MCE1 = (0.0d0,0.0d0)
		end do

		return

	end function bigDdot_MCE1    

!------------------------------------------------------------------------------------

	function zczdot_MCE1(bsin,t)

		type(basisfn), dimension (:), intent (in) :: bsin 
		complex(kind=8), dimension (size(bsin),size(bsin)) :: zczdot_MCE1
		complex(kind=8), dimension (:,:), allocatable :: dz
		integer :: j, k, m, ierr
		real(kind=8), intent (in) :: t

		if (errorflag .ne. 0) return

		ierr = 0
		allocate(dz(size(bsin),ndim), stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		dz = zdot_MCE1(bsin, t)

		do k=1,size(zczdot_MCE1,2)
			do j=1,size(zczdot_MCE1,1)
				zczdot_MCE1(j,k) = (0.0d0, 0.0d0)
				do m=1,ndim
					zczdot_MCE1(j,k) = zczdot_MCE1(j,k) + ((dconjg(bsin(j)%z(m)-bsin(k)%z(m)))*dz(k,m))
				end do
				if (zczdot_MCE1(j,k)/=zczdot_MCE1(j,k)) then
					print "(1x,a,i4,a,i4,a)", "zczdot(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
				zczdot_MCE1(j,k) = i*zczdot_MCE1(j,k)
			end do
		end do

		deallocate(dz, stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		return

	end function zczdot_MCE1

!***********************************************************************************!

end module derivsMCE1