MODULE derivsCCS

	use globvars
	use Ham
	use alarrays
	use outputs
	use redirect

!***********************************************************************************!
!*
!*         CCS Derivatives Module
!*           
!*   Contains subroutines for:
!*
!*      1) Calculating the time derivative of the z values
!*      2) Calculating the derivative of the single configuration MCEv2 d prefactor
!*      3) Calculating the time derivative of the classical action
!*      4) Calculating the derivative of the multi-configuration MCEv2D prefactor
!*      5) Calculating the first component of the Hamiltonian used in 4)
!*      6) Calculating the second component of the Hamiltonian used in 4)
!*      7) Calculating the third component of the Hamiltonian used in 4)
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

	function zdot_CCS(bsin,t)

		implicit none

		complex(kind=8), dimension (:) ,allocatable:: zdottemp, z
		complex(kind=8), dimension (:,:,:),allocatable :: dhdz
		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension(size(bsin),ndim) :: zdot_CCS
		real(kind=8), intent (in) :: t
		integer :: k, m, ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate(zdottemp(ndim), stat = ierr)
		if (ierr==0) allocate (z(ndim), stat = ierr)
		if (ierr==0) allocate (dhdz(npes,npes,ndim), stat = ierr)
		if (ierr/=0) then
			print *, "Error in array allocation in zdot"
			errorflag=1
			return
		end if

		do k=1,size(bsin)
			zdottemp = (0.0d0, 0.0d0)
			z = bsin(k)%z
			call dh_dz(dhdz, z, t)
			do m=1,ndim
				zdottemp(m) = zdottemp(m) + dhdz(1,1,m)
			end do
			do m=1,ndim
				zdot_CCS(k,m)=(-1.0d0*i*zdottemp(m))
			end do
		end do

		deallocate(zdottemp,z,dhdz, stat = ierr)
		if (ierr/=0) then
			print *, "Error in array deallocation in zdot"
			errorflag=1
			return
		end if

		return

	end function zdot_CCS

!------------------------------------------------------------------------------------

	function ddot_CCS(bsin)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension(size(bsin),npes) :: ddot_CCS
		integer :: k, r
		
		if (errorflag .ne. 0) return

		do k=1,size(ddot_CCS,1)
			ddot_CCS(k,1) = (0.0d0, 0.0d0)
		end do

		return

	end function ddot_CCS

!------------------------------------------------------------------------------------

	function sdot_CCS(bsin,dz,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension (:,:), intent(in) :: dz
		real(kind=8), intent (in) :: t
		real(kind=8), dimension (size(bsin),npes) :: sdot_CCS
		complex(kind=8), dimension(:), allocatable :: zk, zkc, zkdot, zkdotc
		complex(kind=8), dimension(:,:), allocatable :: Hkk
		integer :: k, m, ierr
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
			sdot_CCS(k,1) = zsum - dble(Hkk(1,1))
		end do

		deallocate(zk,zkc,zkdot,zkdotc,Hkk,stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of matrices or arrays in sdot"
			errorflag=1
			return
		end if

		return

	end function sdot_CCS

!------------------------------------------------------------------------------------

	function bigDdot_CCS(bsin,x,time,genflg)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin 
		type(hamiltonian), dimension (:,:), allocatable :: H
		complex(kind=8), dimension (:,:), allocatable :: ovrlp, ovrlpphi, ovrlpin
		complex(kind=8), dimension (:,:), allocatable :: h_av_jj, h_av_jk, zconjzdot
		complex(kind=8), dimension (:,:), allocatable :: d2H, d2H2, d2Hdiff
		complex(kind=8), dimension(:), allocatable::tempD, chk, Ddot, d2HD, Dtemp
		complex(kind=8), dimension(size(bsin)) :: bigDdot_CCS
		complex(kind=8) :: ovrlpdif
		real(kind=8), intent(in) :: time
		integer, intent(in) :: x, genflg
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
			print *, "Error in allocation of matrices or arrays in bigDdotv2"
			errorflag=1
			return
		end if    

		ovrlp = ovrlpmat(bsin)
		ovrlpphi = ovrlpphimat(bsin)

		if ((x.eq.1).and.(time.eq.0.0d0)) then
			do j=1,nbf
				do k=1,nbf
					ovrlpdif = ovrlpphi(j,k) - ovrlp(j,k) 
					if (ovrlpdif.ne.0.0d0) then
						print '(a)', "Error! Initial phi-overlap has disimilarilies to z-overlap"
						print '(a,a,i0,a,i0,a)', "These matricies should be identical ",&
																	"but differences found at coordinate ", j,",",k,"."
						print '(a,4(e15.8,a))', "Expected (",dimag(i*ovrlp(j,k)),","&
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

		h_av_jj = Hjk_avrg_CCS(H,ovrlp,bsin)
		h_av_jk = phiHphi_CCS(H,ovrlp,bsin)
		zconjzdot = zczdot_CCS(bsin,time)
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
					print " (1x,a,i4,a,i4,a)", "d2H(", j, ",", k, ") is NaN. Terminating"
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
				print " (1x,a,i4,a)", "tempD(", j, ") is NaN. Terminating"
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
			print *, "Error! Matrix function not recognised! Value is ", matfun
			errorflag = 1
			return
		end if

		do k=1,nbf
			bigDdot_CCS(k)=Ddot(k)
		end do

		call deallocham(H)

		deallocate (ovrlp,ovrlpphi,h_av_jj,h_av_jk,zconjzdot,d2H,d2H2,d2Hdiff,stat=ierr)
		if (ierr == 0 ) deallocate (tempD,chk,Ddot,d2HD,Dtemp,stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of matrices or arrays in bigDdotv2"
			errorflag=1
			return
		end if   

		return

	end function bigDdot_CCS

!------------------------------------------------------------------------------------

	function phiHphi_CCS(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: phiHphi_CCS
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(phiHphi_CCS,2)
			do j=1,size(phiHphi_CCS,1)
				phiHphi_CCS(j,k) = (dconjg(bsin(j)%a_pes(1)) * ovrlp(j,k) * &
														H(j,k)%Hjk(1,1) * bsin(k)%a_pes(1))
				if (phiHphi_CCS(j,k)/=phiHphi_CCS(j,k)) then
					print " (1x,a,i4,a,i4,a)", "phiHphi(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function phiHphi_CCS

!------------------------------------------------------------------------------------

	function Hjk_avrg_CCS(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: Hjk_avrg_CCS
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(Hjk_avrg_CCS,2)
			do j=1,size(Hjk_avrg_CCS,1)
				Hjk_avrg_CCS(j,k) = (dconjg(bsin(j)%a_pes(1)) * H(k,k)%Hjk(r,s) * bsin(k)%a_pes(1))
				Hjk_avrg_CCS(j,k) = Hjk_avrg_CCS(j,k) * ovrlp(j,k)
				if (Hjk_avrg_CCS(j,k)/=Hjk_avrg_CCS(j,k)) then
					print " (1x,a,i4,a,i4,a)", "Hjk_avrg(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function Hjk_avrg_CCS

!------------------------------------------------------------------------------------

	function zczdot_CCS(bsin,t)

		type(basisfn), dimension (:), intent (in) :: bsin 
		real(kind=8), intent (in) :: t
		complex(kind=8), dimension (size(bsin),size(bsin)) :: zczdot_CCS
		complex(kind=8), dimension (:,:), allocatable :: dz
		integer :: j, k, m, ierr

		if (errorflag .ne. 0) return

		ierr = 0
		allocate(dz(size(bsin),ndim), stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		dz = zdot_CCS(bsin, t)

		do k=1,size(zczdot_CCS,2)
			do j=1,size(zczdot_CCS,1)
				zczdot_CCS(j,k) = (0.0d0, 0.0d0)
				do m=1,ndim
					zczdot_CCS(j,k) = zczdot_CCS(j,k) + &
															((dconjg(bsin(j)%z(m)-bsin(k)%z(m)))*dz(k,m))
				end do
				if (zczdot_CCS(j,k)/=zczdot_CCS(j,k)) then
					print "(1x,a,i4,a,i4,a)", "zczdot(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
				zczdot_CCS(j,k) = i*zczdot_CCS(j,k)
			end do
		end do

		deallocate(dz, stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		return

	end function zczdot_CCS

!***********************************************************************************!

end module derivsCCS
