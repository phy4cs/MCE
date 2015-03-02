MODULE Chks

	use globvars
	use Ham

!***********************************************************************************!
!*
!*         Checking Module
!*           
!*   Contains subroutines for:
!*
!*      1) Checking the initial norm and population sum falls in acceptable range, 
!*         and if desired changes the compression parameter and restarts the basis
!*         set generation.
!*      2) Checking the position components of the coherent state are not too 
!*         widely spaced
!*      3) Checking that the conserved quantites do not deviate too much from their
!*         initial values
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

	subroutine initnormchk(bs,recalcs,restart, alcmprss, gridsp, absnorm, popsum)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bs
		integer, intent (inout) :: recalcs, restart
		real(kind=8), intent(inout) :: alcmprss, gridsp, absnorm, popsum
		complex(kind=8)::normtemp
		integer :: r, istat

		if (errorflag .ne. 0) return

		if (((basis.eq."GRID").and.(mod(in_nbf,2)==1)).or.(basis.eq."GRSWM").or.(basis.eq."TRAIN")) then
			uplimnorm = 1.000001d0
		end if
		popsum = 0.0d0

		normtemp = norm(bs)

		absnorm = abs(normtemp)

		do r=1,npes
			popsum = popsum + pop(bs, r)
		end do
 
		if(size(bs).ne.1) then
			 
			if (abs(popsum-absnorm).gt.1.0d-10) then
				print *, "Error! Difference between norm and population sum is too high"
				print *, ""
				print *, "ABS(Norm)  ", absnorm
				print *, "Popsum     ", popsum
				print *, "Difference ", abs(popsum-absnorm)
				print *, ""
				restart = 1
			end if

			if (cmprss.eq."N") then      
				if ((absnorm.gt.uplimnorm).or.(absnorm.lt.lowlimnorm)) then
					print "(a,a,e13.5e3)", "Warning. Initial Norm outside established ",&
											 "parameters, with a value of ", absnorm
					print *, ""
					if ((basis.eq."SWARM").or.(basis.eq."SWTRN")) restart = 1
				end if
			else
				if (absnorm.lt.lowlimnorm) then
					print '(a,es16.8e3)', " Initial Norm too low with a value of ", absnorm
					if (basis.eq."GRID") then
						print *, "Reducing grid spacing to ", (gridsp * 0.95d0)/sqrt(2.)
						print *, ""
						gridsp = gridsp * 0.95d0
					else if ((basis.eq."SWARM").or.(basis.eq."SWTRN")) then
						print *, "Increasing compression parameter to", 1/(alcmprss * 0.95d0)
						print *, ""
						alcmprss = alcmprss * 0.95d0
					end if 
					restart = 1
				else if (absnorm.gt.uplimnorm) then
					print '(a,es16.8e3)', " Initial Norm too high with a value of ", absnorm
					if (basis.eq."GRID") then
						print *, "Increasing grid spacing to ", (gridsp * 1.05d0)/sqrt(2.)
						print *, ""
						gridsp = gridsp * 1.05d0
					else if ((basis.eq."SWARM").or.(basis.eq."SWTRN")) then
						print *, "Reducing compression parameter to", 1/(alcmprss * 1.05d0)
						print *, ""
						alcmprss = alcmprss * 1.05d0
					end if 
					restart = 1
				end if
			end if
			
			if ((restart.eq.1).and.(recalcs.lt.Ntries)) then
				if (cmprss.eq."N") then
					recalcs = recalcs + 1
				end if
				print *, "Recalculating..."
				print *, ""
				return
			else
				return
			end if
	 
		else
	 
			return
	 
		end if

	end subroutine initnormchk

!------------------------------------------------------------------------------------

	subroutine enchk(bf,t,n,redo,k)

		implicit none
		type(basisfn), intent(inout) :: bf
		real(kind=8), intent(in) :: t
		integer, intent(inout) :: n, redo
		integer, intent(in) :: k
		complex(kind=8),dimension(:,:), allocatable::H
		real(kind=8)::Echk
		integer :: ierr

		if (errorflag/=0) return

		if (ECheck.eq."YES") then
			allocate(H(npes,npes), stat = ierr)
			if (ierr/=0) then
				print *, "Error in H allocation in genbasis"
				errorflag=1
				return
			end if
			call Hij(H, bf%z, bf%z, t)
			Echk = dble(H(in_pes,in_pes))
			if ((Echk.gt.Ebfmax).or.(Echk.lt.Ebfmin)) then
				if (n.lt.Ntries) then
					n = n+1
					print *,"Basis ", k, " did not meet energy requirements. ",&
								 "Recalculating..."
					redo=1
				else
					print *,"Basis ", k, " recalculated ", n, "times but still outside ",&
									"acceptable region."
					print *,"Terminating calculation"
					errorflag = 1
					redo=0
				end if
			else
				redo=0
			end if
		else
			redo=0
		end if

		if ((redo/=1).and.(redo/=0)) then
			print *, "Error! Somehow, the redo flag is not 1 or 0"
			errorflag = 1
			return
		end if

		return

	end subroutine enchk

!------------------------------------------------------------------------------------

	subroutine trajchk(bs)   !   Level 1 Subroutine

		implicit none

		type(basisfn),  dimension(:), intent(in) :: bs
		complex(kind=8), dimension(ndim) :: z
		complex(kind=8) :: trajq
		integer :: j, m, flag

		if (errorflag .ne. 0) return

		do j = 1,size(bs)
			flag = 0
			z = bs(j)%z
			do m=1,ndim
				trajq = dble(z(m))*dsqrt(2.0d0/gam)
				if (abs(trajq).gt.20000) then
					print '(a,i0,a,i0,a,e16.8)', "Trajectory ", j, " in dof ", m, &
												" is equal to ", dble(z(m))*dsqrt(2.0d0)
					flag = flag + 1
				end if
			end do
		end do

		if (flag.gt.0) then
			errorflag = 1
			print "(1x,i0,a,a)", flag, " trajectories have position components outside ",&
											"acceptable range."
			return
		end if

		return

	end subroutine trajchk

!------------------------------------------------------------------------------------

	subroutine conservchk(initehr, initnorm, absehr, absnorm, reps)   !   Level 1 Subroutine

		implicit none

		real(kind=8), intent (in) :: initehr, initnorm, absehr, absnorm
		character(LEN=15) :: filenm
		integer, intent(in) :: reps
		character(LEN=3):: rep

		if (errorflag .ne. 0) return

		write(rep,"(i3.3)") reps

		filenm = "normpop-"//trim(rep)//".out"

		if (abs(initnorm-absnorm).ge.1.0d-3) then
			write (*,*) ""
			write (*,*) "*************************Simulation Failed*************************"
			write (*,*) "*************Norm deviated too much from initial value*************"
			errorflag = 1
			return
		end if

		if ((abs(1.0d0-(absehr/initehr)).ge.1.0d-2).and.(method=="MCEv2")) then
			write (*,*) ""
			write (*,*) "*************************Simulation Failed*************************"
			write (*,*) "*******Ehrenfest Energy deviated too much from initial value*******"
			errorflag = 1
			return
		end if

		return

	end subroutine conservchk    

!*************************************************************************************************!  

end module chks


