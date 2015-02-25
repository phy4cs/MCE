MODULE derivs

	use globvars
	use derivsMCE1
	use derivsMCE2
	use derivsCCS

!*************************************************************************************************!
!*
!*         Deriv Redirection Module
!*           
!*   Contains subroutines for:
!*
!*   1) Redirecting the subroutine which calls the time derivative functions
!*      
!*************************************************************************************************!

contains

!--------------------------------------------------------------------------------------------------

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
			print *, "Error in derivatives array allocation in deriv"
			errorflag=1
			return
		end if
		
		select case (method)
			case ("MCEv1")
				dz=zdot_MCE1(bsin,time)
				ds=sdot_MCE1(bsin,dz,time)
				dd=ddot_MCE1(bsin,time)
				dD_big=bigDdot_MCE1(size(bsin))
			case ("MCEv2")
				dz=zdot_MCE2(bsin,time)
				ds=sdot_MCE2(bsin,dz,time)
				dd=ddot_MCE2(bsin,time)
				dD_big=bigDdot_MCE2 (bsin,x,time,genflg)
			case ("CCS")
				dz=zdot_CCS(bsin,time)
				ds=sdot_CCS(bsin,dz,time)
				dd=ddot_CCS(bsin)
				dD_big=bigDdot_CCS (bsin,x,time,genflg)
			case default
				print *, "Error! The propagation method was not recognised!"
				print *, "If you are seeing this something is terribly wrong"
				errorflag=1
		end select			

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
			print *, "Error in derivatives array deallocation in deriv"
			errorflag=1
			return
		end if

		return

	end subroutine deriv
	  
!*************************************************************************************************!

end module derivs
