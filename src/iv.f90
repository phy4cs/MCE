MODULE iv

	use globvars

!***********************************************************************************!
!*
!*		Inverted Gaussian Module
!*					 
!*		Contains subroutines for:
!*
!*		1) Reading the Morse Potential Parameters
!*		2) Calculating the real and imaginary parts of the m dimensional initial 
!*				wavefunction zinit
!*		3) Calculating a single element of the Hamiltonian matrix
!*		4) Calculating the derivative of the hamiltonian
!*			
!***********************************************************************************!

contains

!------------------------------------------------------------------------------------

	subroutine readparams_iv
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
			if(LINE=='IVm') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,mass_iv
				if (ierr.ne.0) then
					print *, "Error reading mass value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='IVw') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,freq_iv
				if (ierr.ne.0) then
					print *, "Error reading frequency value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='IVlambda') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,lambda_iv
				if (ierr.ne.0) then
					print *, "Error reading lambda value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='IVIntensity') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,inten_iv
				if (ierr.ne.0) then
					print *, "Error reading laser intensity value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='IVupnorm') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,uplimnorm
				if (ierr.ne.0) then
					print *, "Error reading upper limit of the norm"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='IVdownnorm') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,lowlimnorm
				if (ierr.ne.0) then
					print *, "Error reading lower limit of the norm"
					errorflag = 1
					return
				end if
				n = n+1
			end if
			read (128,*,iostat=ierr)LINE
		end do

		close (128)

		if (n.ne.6) then
			print *, "Not all required variables read in readparams_iv subroutine"
			errorflag = 1
			return
		end if

		return

	end subroutine readparams_iv

!------------------------------------------------------------------------------------

	subroutine genzinit_iv(mup, muq)	 !	 Level 1 Subroutine

		implicit none
		real(kind=8), dimension(:), allocatable, intent(inout) :: mup, muq

		if (errorflag .ne. 0) return

		muq(1:ndim) = 0.0d0*sigp
		mup(1:ndim) = 0.0d0*sigq

		return

	end subroutine genzinit_iv

!------------------------------------------------------------------------------------

	subroutine Hij_iv(H,z1,z2,t)

		implicit none
		complex(kind=8), dimension (:), intent(in)::z1,z2
		complex(kind=8), dimension(:,:), intent (inout)::H
		real(kind=8), intent (in) :: t
		integer :: m, ierr
		complex(kind=8), dimension (:), allocatable :: Htemp, z1c, rho
		real (kind=8) :: rt2, eta

		if (errorflag .ne. 0) return

		if (npes.ne.1) then
			print *, "Error! There is more than 1 pes for the Inverse Gaussian"
			errorflag = 1
			return
		end if

		allocate (Htemp(ndim), stat=ierr)
		if (ierr/=0) then
			print *, "Error allocating Htemp in Hij_iv"
			errorflag=1
		end if

		eta = lambda_iv*gam/(gam+lambda_iv)

		allocate(z1c(size(z1)))
		allocate(rho(size(z1)))
		z1c(1:ndim)=dconjg(z1(1:ndim))

		do m=1,ndim
			rho(m) = (z1c(m)+z2(m))/sqrt(2.0*gam)
		end do

		do m=1,ndim
			Htemp(m) = (0.0d0,0.0d0)
			Htemp(m) = Htemp(m) - ((hbar**2)*gam/(4.0d0*mass_iv))&
			             *(z1c(m)**2.0d0+z2(m)**2.0d0-2.0d0*z1c(m)*z2(m)-1.0d0)
			Htemp(m) = Htemp(m) - sqrt(eta/lambda_iv)*cdexp(-1.0d0*eta*rho(m)**2.0d0)
			if (m==1) Htemp(m) = Htemp(m) + inten_iv*rho(m)*dcos(freq_iv*t)
		end do

		H(1,1) = sum(Htemp(1:ndim))
		
		if (H(1,1)/=H(1,1)) then
			print *, "Error! Hamiltonian element NaN"
			errorflag = 1
			return
		end if

		deallocate(Htemp, z1c)

		return	 

	end subroutine Hij_iv

!------------------------------------------------------------------------------------

	function dh_dz_iv(z,t)

		implicit none
		complex(kind=8),dimension(npes,npes,ndim) :: dh_dz_iv
		complex(kind=8),dimension(:),intent(in)::z
		complex(kind=8), dimension(:), allocatable :: zc, rho
		real(kind=8), intent (in) :: t
		complex(kind=8) :: dhdztmp
		real(kind=8) :: rt2, eta, eta2
		integer :: m

		if (errorflag .ne. 0) return

		allocate (zc(size(z)))
		allocate (rho(size(z)))
		zc(1:ndim)=dconjg(z(1:ndim))

		eta = lambda_iv*gam/(gam+lambda_iv)
		do m=1,ndim
			rho(m) = (zc(m)+z(m))/sqrt(2.0*gam)
		end do

		do m=1,ndim		
			dhdztmp = (0.0d0,0.0d0) 
			dhdztmp = dhdztmp - ((hbar**2)*gam/(2.0d0*mass_iv))*(zc(m)-z(m))
			dhdztmp = dhdztmp + sqrt(2.0*(eta**3.0)*gam/lambda_iv)*rho(m)*&
			                    exp(-1.0d0*eta*rho(m)**2.0)
			if (m==1) dhdztmp = dhdztmp + inten_iv*(1.0/sqrt(2.0d0*gam))*dcos(freq_iv*t)
			dh_dz_iv (1,1,m) = dhdztmp
		end do

		return

	end function dh_dz_iv

!------------------------------------------------------------------------------------

	function dipole_iv(bs, x)	 !	 Level 1 Function

		implicit none
		type(basisfn),dimension(:),intent(in)::bs
		integer, intent(in) :: x
		complex(kind=8) :: dipole_iv, zsum, ovrlp
		complex(kind=8), dimension(:), allocatable :: D, Dc, zk, zj, rho
		real(kind=8), dimension (:), allocatable :: s
		real(kind=8) :: eta, eta2, rt2
		integer::k,j,m,ierr

		if (errorflag .ne. 0) return

		ierr = 0

		dipole_iv = (0.0d0,0.0d0)
		eta = gam*lambda_iv/(gam+lambda_iv)

		if ((ndim.ne.1).and.(ndim.ne.3)) then
			print *, "Error! ndim should be 1 or 3 but is not."
			errorflag=1
			return
		end if

		if ((npes.ne.1).or.(in_pes.ne.1)) then
			print *, "Error! npes and in_pes should be 1 but are not."
			errorflag=1
			return
		end if		

		allocate (D(size(bs)), stat = ierr)
		if (ierr==0) allocate (Dc(size(bs)), stat=ierr)
		if (ierr==0) allocate (zk(ndim), stat=ierr)
		if (ierr==0) allocate (zj(ndim), stat=ierr)
		if (ierr==0) allocate (rho(ndim), stat=ierr)
		if (ierr==0) allocate (s(size(bs)), stat=ierr)
		if (ierr/=0) then
			print *, "Error allocating basis set variable arrays for dipole calculation"
			errorflag = 1
			return
		end if

		do k=1,size(bs)
			D(k)=bs(k)%D_big
			Dc(k)=dconjg(D(k))
			s(k)=bs(k)%s_pes(1)
		end do		 

		do k=1,size(bs)
			do j=1,size(bs)
				do m=1,ndim
					zk(m) = bs(k)%z(m)
					zj(m) = bs(j)%z(m)
					rho(m) = (dconjg(zj(m))+zk(m))/sqrt(2.0*gam)
				end do
				ovrlp = product(cdexp((dconjg(zj(1:ndim))*zk(1:ndim))&
				      -(0.5d0*dconjg(zj(1:ndim))*zj(1:ndim))&
							-(0.5d0*dconjg(zk(1:ndim))*zk(1:ndim))))
				do m=1,ndim
					dipole_iv = dipole_iv + (-1.0d0*Dc(j)*D(k)*cdexp(i*(s(k)-s(j)))* ovrlp*&
						                        2.0d0*sqrt((eta**3.0d0)/lambda_iv)*rho(m)&
						                                    *cdexp(-1.0d0*eta*(rho(m)**2.0d0)))
				end do
			end do
		end do

		deallocate (D, stat = ierr)
		if (ierr==0) deallocate (Dc, stat=ierr)
		if (ierr==0) deallocate (zk, stat=ierr)
		if (ierr==0) deallocate (zj, stat=ierr)
		if (ierr==0) deallocate (rho, stat=ierr)
		if (ierr==0) deallocate (s, stat=ierr)
		if (ierr/=0) then
			print *, "Error deallocating basis set variable arrays ",&
			           "for dipole acceleration calculation"	 
			errorflag = 1	 
			return
		end if

		return

	end function dipole_iv

!***********************************************************************************!

end module iv