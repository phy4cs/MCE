MODULE derivsMCE2

	use globvars
	use Ham
	use alarrays
	use outputs
	use redirect

!***********************************************************************************!
!*
!*         MCE2 Derivatives Module
!*           
!*   Contains subroutines for:
!*
!*      1) Calculating the time derivative of the z values
!*      2) Calculating the derivative of the single configuration MCEv2 d prefactor
!*      3) Calculating the time derivative of the classical action
!*      4) Calculating the derivative of the multi-configuration MCEv2 D prefactor
!*      5) Calculating the first component of the Hamiltonian used in 4)
!*      6) Calculating the second component of the Hamiltonian used in 4)
!*      7) Calculating the third component of the Hamiltonian used in 4)
!*      
!***********************************************************************************!

contains

!***********************************************************************************!

	function zdot_MCE2(bsin,t)

		implicit none

		complex(kind=8), dimension (:) ,allocatable:: zdottemp, z
		complex(kind=8), dimension (:,:,:),allocatable :: dhdz
		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension(size(bsin),ndim) :: zdot_MCE2
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
				zdot_MCE2(k,m)=(-1.0d0*i*zdottemp(m))/asum
			end do
		end do

		deallocate(zdottemp,z,dhdz,a,ac, stat = ierr)
		if (ierr/=0) then
			print *, "Error in array deallocation in zdot"
			errorflag=1
			return
		end if

		return

	end function zdot_MCE2

!------------------------------------------------------------------------------------

	function ddot_MCE2(bsin,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		real(kind=8), intent (in) :: t
		complex(kind=8), dimension(:,:),allocatable :: Hkk
		complex(kind=8), dimension(:),allocatable :: dk
		complex(kind=8), dimension(size(bsin),npes) :: ddot_MCE2
		complex(kind=8), dimension(:),allocatable :: z
		real(kind=8), dimension(:),allocatable :: Sk
		integer :: k, r, s, ierr

		if (errorflag .ne. 0) return

		ierr = 0

		allocate(Hkk(npes,npes),dk(npes),z(ndim),Sk(npes),stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of matrices or arrays in ddotv2"
			errorflag=1
			return
		end if

		do r=1,size(ddot_MCE2,2)
			do k=1,size(ddot_MCE2,1)
				ddot_MCE2(k,r) = (0.0d0, 0.0d0)
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
						ddot_MCE2(k,r) = ddot_MCE2(k,r) + Hkk(r,s) * ovrlpij(z,z) * dk(s) &
						                                     * cdexp(i*(Sk(s)-Sk(r)))
					end if
				end do
			end do
		end do

		do r=1,size(ddot_MCE2,2)
			do k=1,size(ddot_MCE2,1)
				ddot_MCE2(k,r) = -1.0d0*i*ddot_MCE2(k,r)
			end do
		end do

		deallocate(Hkk,dk,z,Sk,stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of matrices or arrays in ddotv2"
			errorflag=1
			return
		end if

		return

	end function ddot_MCE2

!------------------------------------------------------------------------------------

	function sdot_MCE2(bsin,dz,t)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin
		complex(kind=8), dimension (:,:), intent(in) :: dz
		real(kind=8), intent (in) :: t
		real(kind=8), dimension (size(bsin),npes) :: sdot_MCE2
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
				sdot_MCE2(k,r) = zsum - dble(Hkk(r,r))
			end do
		end do

		deallocate(zk,zkc,zkdot,zkdotc,Hkk,stat = ierr)
		if (ierr/=0) then
			print *, "Error in allocation of matrices or arrays in sdot"
			errorflag=1
			return
		end if

		return

	end function sdot_MCE2

!------------------------------------------------------------------------------------

	function bigDdot_MCE2(bsin,x,time,genflg)

		implicit none

		type(basisfn), dimension (:), intent (in) :: bsin 
		type(hamiltonian), dimension (:,:), allocatable :: H
		complex(kind=8), dimension (:,:), allocatable :: ovrlp, ovrlpphi, ovrlpin
		complex(kind=8), dimension (:,:), allocatable :: h_av_jj, h_av_jk, zconjzdot
		complex(kind=8), dimension (:,:), allocatable :: d2H, d2H2, d2Hdiff
		complex(kind=8), dimension(:), allocatable::tempD, chk, Ddot, d2HD, Dtemp
		complex(kind=8), dimension(size(bsin)) :: bigDdot_MCE2
		complex(kind=8) :: ovrlpdif
		real(kind=8), intent(in) :: time
		integer, intent(in) :: x, genflg
		integer :: j, k, nbf, ierr
		real(kind=8)::absB

		if (errorflag .ne. 0) return

		ierr = 0
		
		if (((basis=="TRAIN").or.(basis=="SWTRN")).and.(genflg==1)) then
			do j=1,nbf
				bigDdot_MCE2 = (0.0d0,0.0d0)
			end do
			return
		end if

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
					if ((ovrlpdif.ne.0.0d0).and.(basis.ne."TRAIN").and.(basis.ne."SWTRN")) then
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

		h_av_jj = Hjk_avrg_MCE2(H,ovrlp,bsin)
		h_av_jk = phiHphi_MCE2(H,ovrlp,bsin)
		zconjzdot = zczdot_MCE2(bsin,time)
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
			bigDdot_MCE2(k)=Ddot(k)
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

	end function bigDdot_MCE2

!------------------------------------------------------------------------------------

	function phiHphi_MCE2(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: phiHphi_MCE2
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(phiHphi_MCE2,2)
			do j=1,size(phiHphi_MCE2,1)
				phiHphi_MCE2(j,k) = (0.0d0, 0.0d0)
				do r=1,npes
					do s=1,npes
						phiHphi_MCE2(j,k) = phiHphi_MCE2(j,k) + (dconjg(bsin(j)%a_pes(r)) * &
																ovrlp(j,k) * H(j,k)%Hjk(r,s) * bsin(k)%a_pes(s))
					end do
				end do
				if (phiHphi_MCE2(j,k)/=phiHphi_MCE2(j,k)) then
					print " (1x,a,i4,a,i4,a)", "phiHphi(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function phiHphi_MCE2

!------------------------------------------------------------------------------------

	function Hjk_avrg_MCE2(H,ovrlp,bsin)

		implicit none
		type(basisfn), dimension (:), intent (in) :: bsin
		type(hamiltonian), dimension (:,:), intent(in) :: H
		complex(kind=8), dimension (:,:), intent (in) :: ovrlp
		complex(kind=8), dimension (size(bsin),size(bsin)) :: Hjk_avrg_MCE2
		integer :: j, k , r, s

		if (errorflag .ne. 0) return

		do k=1,size(Hjk_avrg_MCE2,2)
			do j=1,size(Hjk_avrg_MCE2,1)
				Hjk_avrg_MCE2(j,k) = (0.0d0, 0.0d0)
				do r=1,npes
					do s=1,npes
						Hjk_avrg_MCE2(j,k) = Hjk_avrg_MCE2(j,k) + (dconjg(bsin(j)%a_pes(r)) * &
																		ovrlp(k,k) * H(k,k)%Hjk(r,s) * bsin(k)%a_pes(s))
					end do
				end do
				Hjk_avrg_MCE2(j,k) = Hjk_avrg_MCE2(j,k) * ovrlp(j,k)
				if (Hjk_avrg_MCE2(j,k)/=Hjk_avrg_MCE2(j,k)) then
					print " (1x,a,i4,a,i4,a)", "Hjk_avrg(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
			end do
		end do

		return

	end function Hjk_avrg_MCE2

!------------------------------------------------------------------------------------

	function zczdot_MCE2(bsin,t)

		type(basisfn), dimension (:), intent (in) :: bsin 
		complex(kind=8), dimension (size(bsin),size(bsin)) :: zczdot_MCE2
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

		dz = zdot_MCE2(bsin, t)

		do k=1,size(zczdot_MCE2,2)
			do j=1,size(zczdot_MCE2,1)
				zczdot_MCE2(j,k) = (0.0d0, 0.0d0)
				do m=1,ndim
					zczdot_MCE2(j,k) = zczdot_MCE2(j,k) + &
															((dconjg(bsin(j)%z(m)-bsin(k)%z(m)))*dz(k,m))
				end do
				if (zczdot_MCE2(j,k)/=zczdot_MCE2(j,k)) then
					print "(1x,a,i4,a,i4,a)", "zczdot(", j, ",", k, ") is NaN. Terminating"
					errorflag = 1
					return
				end if
				zczdot_MCE2(j,k) = i*zczdot_MCE2(j,k)
			end do
		end do

		deallocate(dz, stat = ierr)
		if (ierr/=0) then
			print *, "Error in deallocation of dz arrays in zczdot"
			errorflag=1
			return
		end if   

		return

	end function zczdot_MCE2

!***********************************************************************************!

end module derivsMCE2