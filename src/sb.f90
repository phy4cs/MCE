MODULE sb

	use globvars

!*************************************************************************************************!
!                                                                                                 !
!      This Module contains the subroutines and functions for the Spin Boson model.               !
!                                                                                                 !
!      It can be used as a template for other models, as all other subroutines and functions      !
!      in the program are very general, allowing different Hamiltonians or numbers of PESs        !
!      to be used.                                                                                !
!                                                                                                 !
!      The subroutines / functions which must be included in any future                           !
!      adaptations of this module are:                                                            !
!                                                                                                 !
!         1) ***subroutine readparams*** This reads the system specific parameters                !
!         2) ***subroutine genzinit*** This creates the initial CS                                !
!         3) ***function Hij(z1,z2)*** This builds a npes x npes matrix which represent a single  !
!                  basis function combination for the Hamiltonian                                 !
!         4) ***function dh_dz(z)***   This calculates the derivative of the Ehrenfest            !
!                  Hamiltonian for use in calculating the time derivative for z                   !
!                                                                                                 !
!       All of the subroutines and functions which are dependent upon the above should be         !
!       included in this module also. The subroutines should be referenced in the redirect        !
!       switchboard module and identified by a two character code, ie spin boson = sb             !
!                                                                                                 !
!*************************************************************************************************!

!*************************************************************************************************!
!*
!*         Spin Boson Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the Spin Boson Parameters
!*      2) Calculating the array for the frequency distribution wm
!*      3) Calculating the array for the amplitudes Cm 
!*      4) Calculating the array for the width distribution sig 
!*      5) Calculating the real and imaginary parts of the m dimensional initial wavefunction zinit
!*      6) Calculating a single 2 x 2 matrix which is equivalent to a single element of the 
!*         Hamiltonian matrix
!*      7) Calculating the Bath Hamiltonian for a pair of basis functions
!*      8) Calculating the Coupling Hamiltonian for a pair of basis functions
!*      9) Combining the elements of the differential of the hamiltonan for each combination of 2
!*         PESs with respect to z
!*     10) Calculating the differential of the hamiltonian for the PES combination (1,1)
!*     11) Calculating the differential of the hamiltonian for the PES combination (1,2)
!*     12) Calculating the differential of the hamiltonian for the PES combination (2,1)
!*     13) Calculating the differential of the hamiltonian for the PES combination (2,2)
!*      
!*************************************************************************************************!

contains

!--------------------------------------------------------------------------------------------------

	subroutine readparams_sb
		implicit none
		character(LEN=100)::LINE
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		open(unit=128, file='inham.dat', status='old', iostat=ierr)

		if (ierr.ne.0) then
			print *, 'Error in opening inham.dat file'
			errorflag = 1
			return
		end if

		read(128,*,iostat=ierr)LINE

		do while (ierr==0)

			if(LINE=='SBDelta') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,delta_sb
				if (ierr.ne.0) then
					print *, "Error reading Delta value"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='SBEps') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,eps_sb
				if (ierr.ne.0) then
					print *, "Error reading Epsilon value"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='SBw') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,wc_sb
				if (ierr.ne.0) then
					print *, "Error reading wc value"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='SBkondo') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,kondo_sb
				if (ierr.ne.0) then
					print *, "Error reading kondo parameter value"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='SBwmax') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,wmax_sb
				if (ierr.ne.0) then
					print *, "Error reading wmax value"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='SBBeta') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,beta_sb
				if (ierr.ne.0) then
					print *, "Error reading beta value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='SBupnorm') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,uplimnorm
				if (ierr.ne.0) then
					print *, "Error reading upper limit of the norm"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='SBdownnorm') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,lowlimnorm
				if (ierr.ne.0) then
					print *, "Error reading lower limit of the norm"
					errorflag = 1
					return
				end if
				n = n+1
			end if

			read(128,*,iostat=ierr) LINE

		end do

		if (wmax_sb .ne. 5.0d0 * wc_sb) then
			wmax_sb = wc_sb * 5.0d0
		end if

		close (128)

		if (n.ne.8) then
			print *, "Not all required variables read in readparams_sb subroutine"
			errorflag = 1
			return
		end if

		return

	end subroutine readparams_sb

!--------------------------------------------------------------------------------------------------

	subroutine genwm_sb(wm_sb)   !   Level 2 Subroutine

		implicit none
		integer::m, ierr
		real(kind=8), dimension(:), allocatable, intent(inout) :: wm_sb

		if (errorflag .ne. 0) return
		ierr = 0

		if (ndim.le.0) then
			print *, "Number of dimensions not correctly read or stored"
			errorflag = 1
			return
		else if (.not.allocated(wm_sb)) then
			print *, "genwm subroutine called but wm not allocated"
			errorflag=1
			return
		else
			do m=1,size(wm_sb)
			  wm_sb(m)=-1.0d0*wc_sb*dlog(1.0d0-((dble(m)*(1.0d0-dexp(-1.0d0*wmax_sb/wc_sb)))/dble(ndim)))
			end do
		end if

		return

	end subroutine genwm_sb

!--------------------------------------------------------------------------------------------------

	subroutine genCm_sb(Cm_sb, wm_sb)   !   Level 2 Subroutine

		implicit none
		integer::m, ierr
		real(kind=8), dimension(:), intent(in) :: wm_sb 
		real(kind=8), dimension(:), allocatable, intent(inout) :: Cm_sb  

		if (errorflag .ne. 0) return
		ierr = 0

		if (ndim.le.0) then
			print *, 'Number of dimensions not correctly read or stored'
			errorflag = 1
			return
		else if (.not.allocated(Cm_sb)) then
			print *, "genCm subroutine called but Cm not allocated"
			errorflag=1
			return
		else
			do m=1,size(Cm_sb)
				Cm_sb(m)=wm_sb(m)*sqrt(kondo_sb*wc_sb*(1.0d0-dexp(-1.0d0*wmax_sb/wc_sb))/dble(ndim))
			end do
		end if

		return

	end subroutine genCm_sb

!--------------------------------------------------------------------------------------------------

	subroutine gensig_sb(sig_sb, wm_sb)   !   Level 2 Subroutine

		implicit none
		real(kind=8), dimension (:), allocatable, intent (inout) :: sig_sb
		real(kind=8), dimension (:), intent (in) :: wm_sb
		integer::m, ierr

		if (errorflag .ne. 0) return
		ierr = 0

		if (ndim.le.0) then
			print *, 'Number of dimensions not correctly read or stored'
			errorflag = 1
			return
		else if (.not.allocated(sig_sb)) then
			print *, 'Allocation error in sig'
			errorflag = 1
			return
		else
			do m=1,size(sig_sb)
				sig_sb(m) = 1.0d0/(sqrt(dexp(beta_sb*wm_sb(m))-1.0d0))
			end do
		end if

		return

	end subroutine gensig_sb

!--------------------------------------------------------------------------------------------------

	subroutine genzinit_sb(mup, muq)   !   Level 1 Subroutine

		implicit none
		integer::m, n, recalc, ierr
		real(kind=8) :: Ezin
		real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq
		real(kind=8), dimension(:), allocatable :: sig_sb, wm_sb, Cm_sb
		complex(kind=8), dimension (:), allocatable:: zin
		complex(kind=8),dimension(:,:), allocatable::H

		if (errorflag .ne. 0) return

		n = 0
		recalc = 1
		ierr = 0

		allocate(zin(ndim), H(npes,npes), stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of zin and H in genzinit"
			errorflag=1
			return
		end if  

		if (ndim.le.0) then
			print *, 'Number of dimensions not correctly read or stored'
			errorflag = 1
			return
		end if

		if ((.not.allocated(mup)).or.(.not.allocated(muq)).or.(size(mup).ne.size(muq)) &
					 .or.(size(mup).ne.ndim)) then
			print *, 'Allocation error in mup or muq'
			errorflag = 1
			return
		end if

		allocate (sig_sb(ndim), stat = ierr)
		if (ierr==0) allocate (wm_sb(ndim), stat=ierr)
		if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
		if (ierr/=0) then
			print *, "Error in sig or wm allocation in genzinit"
			errorflag=1
			return
		end if

		call genwm_sb(wm_sb)
		call gensig_sb(sig_sb, wm_sb)    

		do while (recalc == 1)
			do m=1,ndim
				mup(m)=ZBQLNOR(mu,sig_sb(m)*sigp)
				muq(m)=ZBQLNOR(mu,sig_sb(m)*sigq)
				zin(m)=cmplx(muq(m),mup(m),kind=8)
			end do
			if (ECheck.eq."YES") then
				call Hij_sb(H,zin,zin)
				Ezin = dble(H(in_pes,in_pes))
				if ((Ezin.gt.Ebfmax).or.(Ezin.lt.Ebfmin)) then
					if (n.lt.Ntries) then
						n = n+1
						print *,"Initial Basis did not meet energy requirements. Recalculating..."
						print *,"Ebfmax = ", Ebfmax, "Ebfmin = ", Ebfmin, "Ezin = ", Ezin
						recalc = 1
					else
						print *,"Initial basis recalculated ", n, "times but still outside acceptable region."
						Print *,"Terminating calculation"
						errorflag = 1
						recalc = 0
					end if
				else
					recalc = 0           
				end if
			else
				recalc = 0
			end if
		end do

		deallocate(zin, H, sig_sb, wm_sb, stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of zin, H, wm or sig in genzinit"
			errorflag=1
			return
		end if 

		return

	end subroutine genzinit_sb

!--------------------------------------------------------------------------------------------------

	subroutine Hij_sb(H,z1,z2)

		implicit none
		complex(kind=8), dimension (:), intent(in)::z1,z2
		complex(kind=8), dimension(:,:), intent (inout)::H
		complex(kind=8):: Hbath, Hcoup
		real(kind=8), dimension(:), allocatable :: wm_sb, Cm_sb
		real(kind=8) :: chk
		integer :: ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate (wm_sb(ndim), stat=ierr)
		if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
		if (ierr/=0) then
			print *, "Error in wm or Cm allocation in Hij"
			errorflag=1
			return
		end if

		call genwm_sb(wm_sb)
		call genCm_sb(Cm_sb, wm_sb)

		Hbath = Hb_sb(z1,z2,wm_sb)
		Hcoup = Hc_sb(z1,z2,wm_sb, Cm_sb)

		chk = sum((abs(z1(1:ndim)-z2(1:ndim)))**2)

		if (chk.le.20000) then
			H(1,1) = Hbath+Hcoup+eps_sb
			H(1,2) = cmplx(delta_sb,0.0d0,kind=8)
			H(2,1) = cmplx(delta_sb,0.0d0,kind=8)
			H(2,2) = Hbath-Hcoup-eps_sb
		else
			H = (0.0d0, 0.0d0)
		end if

		deallocate(wm_sb, Cm_sb, stat=ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of wm or Cm in Hij"
			errorflag=1
			return
		end if 

		return   

	end subroutine Hij_sb

!--------------------------------------------------------------------------------------------------

	function Hb_sb(z1,z2,wm_sb)

		implicit none
		complex(kind=8), dimension (:), intent(in)::z1,z2
		real(kind=8), dimension(:), intent(in) :: wm_sb
		complex(kind=8)::Hb_sb
		integer :: m

		if (errorflag .ne. 0) return

		Hb_sb = (0.0d0, 0.0d0)

		do m=1,ndim
			Hb_sb = Hb_sb + (wm_sb(m)*(dconjg(z1(m))*z2(m)+0.5d0))
		end do

		return

	end function Hb_sb
		
!--------------------------------------------------------------------------------------------------

	function Hc_sb(z1,z2,wm_sb, Cm_sb)

		implicit none
		complex(kind=8), dimension (:), intent(in)::z1,z2
		real(kind=8), dimension(:), intent(in) :: wm_sb, Cm_sb
		complex(kind=8)::Hc_sb
		integer :: m

		if (errorflag .ne. 0) return

		Hc_sb = (0.0d0, 0.0d0)

		do m=1,ndim
			Hc_sb = Hc_sb + ((Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))*(dconjg(z1(m))+z2(m)))
		end do

		return

	end function Hc_sb

!--------------------------------------------------------------------------------------------------

	function dh_dz_sb(z)

		implicit none
		complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_sb
		complex(kind=8),dimension(:),intent(in)::z
		real(kind=8), dimension(:), allocatable :: wm_sb, Cm_sb
		integer :: ierr

		if (errorflag .ne. 0) return

		if (npes.ne.2) then
			print *, "No. PES is not equal to 2, but Spin Boson 2 PES derivative called"
			errorflag = 1
			return
		end if

		allocate (wm_sb(ndim), stat=ierr)
		if (ierr==0) allocate (Cm_sb(ndim), stat=ierr)
		if (ierr/=0) then
			print *, "Error in wm or Cm allocation in dh_dz"
			errorflag=1
			return
		end if 

		call genwm_sb(wm_sb)
		call genCm_sb(Cm_sb, wm_sb)

		dh_dz_sb(1,1,:) = dhdz_sb_11(z,wm_sb,Cm_sb)
		dh_dz_sb(1,2,:) = dhdz_sb_12()
		dh_dz_sb(2,1,:) = dhdz_sb_21()
		dh_dz_sb(2,2,:) = dhdz_sb_22(z,wm_sb,Cm_sb)

		deallocate (wm_sb, Cm_sb, stat=ierr)
		if (ierr/=0) then
			print *, "Error in wm or Cm allocation in dh_dz"
			errorflag=1
			return
		end if 

		return

	end function dh_dz_sb

!--------------------------------------------------------------------------------------------------

	function dhdz_sb_11(z,wm_sb,Cm_sb)

		implicit none
		complex(kind=8),dimension(ndim)::dhdz_sb_11
		complex(kind=8),dimension(:),intent(in)::z
		real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
		integer::m

		if (errorflag .ne. 0) return

		do m=1,ndim
			dhdz_sb_11(m) = wm_sb(m)*z(m)+(Cm_sb(m)/sqrt(2.0d0*wm_sb(m)))
		end do

		return

	end function dhdz_sb_11

!--------------------------------------------------------------------------------------------------

	function dhdz_sb_12()

		implicit none
		complex(kind=8),dimension(ndim)::dhdz_sb_12
		integer::m

		if (errorflag .ne. 0) return

		do m=1,ndim
			dhdz_sb_12 = (0.0d0,0.0d0)
		end do

		return

	end function dhdz_sb_12

!--------------------------------------------------------------------------------------------------

	function dhdz_sb_21()

		implicit none
		complex(kind=8),dimension(ndim)::dhdz_sb_21
		integer::m

		if (errorflag .ne. 0) return

		do m=1,ndim
			dhdz_sb_21(m) = (0.0d0,0.0d0)
		end do

		return

	end function dhdz_sb_21

!--------------------------------------------------------------------------------------------------

	function dhdz_sb_22(z,wm_sb,Cm_sb)

		implicit none
		complex(kind=8),dimension(ndim)::dhdz_sb_22
		complex(kind=8),dimension(:),intent(in)::z
		real(kind=8),dimension(:),intent(in)::wm_sb, Cm_sb
		integer::m

		if (errorflag .ne. 0) return

		do m=1,ndim
			dhdz_sb_22(m) = wm_sb(m)*z(m)-Cm_sb(m)/sqrt(2.0d0*wm_sb(m))
		end do

		return

	end function dhdz_sb_22

!*************************************************************************************************!

end module sb

