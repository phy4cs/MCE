MODULE readpars

	use globvars
	use alarrays
	use Ham
	use redirect

!*************************************************************************************************!
!*
!*         Input File Reading Module
!*           
!*   Contains subroutines for:
!*
!*      1) Reading the run conditions (debug,gen,prop,cmprss,method,reptot,conjflg)
!*      2) Reading the system type (currently only Spin Boson supported)
!*      3) Reading the energy cutoff parameters (ECheck,Ntries,Ebfmin,Ebfmax)
!*      4) Reading the basis set parameters (ndim,in_nbf,matfun,npes,in_pes,grid)
!*      5) Reading initial wavefunction parameters (initialcmprss,gam,mu,hbar and calculation of sigp & sigq)
!*      6) Reading a pre-calculated basis set (inc. bs params and all bs values)
!*      7) Reading time propagation parameters (dtmin,dtmax,dtinit,timeend,timestrt,step)
!*      
!*************************************************************************************************!

contains

!*************************************************************************************************!
!          Shared Reading Subroutines
!*************************************************************************************************!

	subroutine readrunconds   !   Level 1 Subroutine

		implicit none
		character(LEN=100)::LINE, LINE2
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		open(unit=140, file='input.dat', status='old', iostat=ierr)

		if (ierr.ne.0) then
			print *, 'Error in opening runconds.dat file'
			errorflag = 1
			return
		end if

		read(140,*,iostat=ierr)LINE

		do while (ierr==0)

			if(LINE=='debug') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,debug
				if (ierr.ne.0) then
					print *, "Error reading debug status"
					errorflag = 1
					return
				end if
				if ((debug.ne.0).and.(debug.ne.1)) then
					print *, "Error in debug state. Debug value should be only 0 (for off) or 1 (for on)"
					print "(a,i3)", "Debug value read was ", debug
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=='gen') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,LINE2
				if (ierr.ne.0) then
					print *, "Error reading basis set generation status"
					errorflag = 1
					return
				end if
				if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
					gen = "Y"
				else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
					gen = "N"
				else
					print *, "Error. gen value must be YES/NO. Read ", trim(LINE2)
				end if
				n=n+1
			else if (LINE=='prop') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,LINE2
				if (ierr.ne.0) then
					print *, "Error reading basis set propagation status"
					errorflag = 1
					return
				end if
				if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
					prop = "Y"
				else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
					prop = "N"
				else
					print *, "Error. prop value must be YES/NO. Read ", trim(LINE2)
				end if
				n=n+1              
			else if (LINE=='cmprss') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,LINE2
				if (ierr.ne.0) then
					print *, "Error reading compression parameter change status"
					errorflag = 1
					return
				end if
				if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
					cmprss = "Y"
				else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
					cmprss = "N"
				else
					print *, "Error. cmprss value must be YES/NO. Read ", trim(LINE2)
				end if
				n=n+1 
			else if (LINE=='method') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,LINE2
				if (ierr.ne.0) then
					print *, "Error reading basis set propagation method"
					errorflag = 1
					return
				end if
				if ((LINE2(1:5).eq.'mcev1').or.(LINE2(1:5).eq.'MCEv1')) then
					method = "MCEv1"
				else if ((LINE2(1:5).eq.'mcev2').or.(LINE2(1:5).eq.'MCEv2')) then
					method = "MCEv2"
				else if ((LINE2(1:3).eq.'ccs').or.(LINE2(1:3).eq.'CCS')) then
					method = "CCS"
				else
					print *, "Error. Method must be MCEv1,MCEv2 or CCS. Read ", trim(LINE2)
				end if
				n=n+1
			else if (LINE=='Repeats') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,reptot
				if (ierr.ne.0) then
					print *, "Error reading number of repeats"
					errorflag = 1
					return
				end if
				n=n+1                    
			else if (LINE=='Conjugate_Repeats') then
				backspace(140)
				read(140,*,iostat=ierr)LINE,LINE2
				if (ierr.ne.0) then
					print *, "Error reading conjugate repeats flag"
					errorflag = 1
					return
				end if
				if ((LINE2(1:1).eq.'y').or.(LINE2(1:1).eq.'Y')) then
					conjflg = 1
				else if ((LINE2(1:1).eq.'n').or.(LINE2(1:1).eq.'N')) then
					conjflg = 0
				else
					print *, "Error. Conjugate repeats flag must be Yes or No. Read ", trim(LINE2)
				end if
				n=n+1                     
			end if

			read(140,*,iostat=ierr) LINE

		end do

		close(140)

		if ((gen.eq."N").and.(prop.eq."N")) then
			print *, "Error! Run conditions are for no basis set generation or propagation. So what now genius?"
			errorflag=1
			return
		end if

		if ((conjflg==1).and.(mod(reptot,2).ne.0)) then
			print *,"Warning! An odd number of repeats selected but conjugate repetition chosen!"
			print *,"Incrementing repeat total to even number"
			reptot = reptot + 1
		end if

		if ((conjflg==1).and.(gen=="N")) then
			print *,"Error! Conjugate repetition is not compatible for simulations with pre-calculated basis set"
			errorflag=1
			return
		end if

		if (n.ne.7) then
			print *, "Not all required variables read in readrunconds subroutine"
			errorflag = 1
			return
		end if

		return

	end subroutine readrunconds

!--------------------------------------------------------------------------------------------------

	subroutine readsys   !   Level 1 Subroutine

		implicit none
		character(LEN=100)::LINE
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		open(unit=127, file='input.dat', status='old', iostat=ierr)

		if (ierr.ne.0) then
			print *, 'Error in opening input.dat file'
			errorflag = 1
			return
		end if

		read(127,*,iostat=ierr)LINE

		do while (ierr==0)

			if(LINE=='System:') then
				backspace(127)
				read(127,*,iostat=ierr)LINE,sys
				if (ierr.ne.0) then
					print *, "Error reading System Name"
					errorflag = 1
					return
				end if
				n = n+1
			end if
			read(127,*,iostat=ierr) LINE

		end do

		close(127)

		select case (sys)
			case ("SB")
				if (npes.lt.2) then
					print *, "Spin Boson model must have at least 2 pes'"
					errorflag = 1
					return
				end if
				if ((method.ne."MCEv1").and.(method.ne."MCEv2")) then
					print *, "Spin Boson model can only be simulated by MCEv1 or MCEv2"
					errorflag = 1
					return
				end if
				if (basis.eq."GRID") then
					print *, "This method must not use a static grid."
					errorflag = 1
					return
				end if 
			case ("HP")
				if (npes.ne.1) then
					print *, "Harmonic Potential only valid for 1 PES"
					errorflag = 1
					return
				end if
				if (trim(method).ne."CCS") then
					print *, "Harmonic Potential can only be simulated by CCS"
					errorflag = 1
					return  
				end if
			case ("FP")
				if (npes.ne.1) then
					print *, "Free Particle only valid for 1 PES"
					errorflag = 1
					return
				end if
				if (trim(method).ne."CCS") then
					print *, "Free Particle can only be simulated by CCS"
					errorflag = 1
					return  
				end if
			case ("MP")
				if (npes.ne.1) then
					print *, "Morse Potential only valid for 1 PES"
					errorflag = 1
					return
				end if
				if (trim(method).ne."CCS") then
					print *, "Morse Potential can only be simulated by CCS"
					errorflag = 1
					return  
				end if 
			case ("IV")
				if (npes.ne.1) then
					print *, "Inverted Gaussian only valid for 1 PES"
					errorflag = 1
					return
				end if
				if ((ndim.ne.1).and.(ndim.ne.3)) then
					print  *, "Inverted Gaussian is only valid for 1 or 3 dimensional"
					errorflag = 1
					return 
				end if
				if (trim(method).ne."CCS") then
					print *, "Inverted Gaussian can only be simulated by CCS"
					errorflag = 1
					return  
				end if 
				if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
					print *, "This method must use a static grid."
					errorflag = 1
					return
				end if
			case ("CP")
				if (npes.ne.1) then
					print *, "Coulomb Potential only valid for 1 PES"
					errorflag = 1
					return
				end if
				if (mod(qsizez,2)==1) then
					print *, "Odd parity for the grid parameters will result in a point on the singularity"
					print *, "This would make reprojection impossible, as well as seriously affecting results"
					print *, "Change the grid parameters and restart."
					errorflag=1
					return
				end if
				if (in_nbf.eq.(qsizez*psizez+1)) then
					print *, "The in_nbf value has been reset so that there is a point on the singularity."
					print *, "This will make reprojection impossible. Resetting to even parity."
					in_nbf = qsizez*psizez
				end if
				if (ndim.ne.3) then
					print  *, "Coulomb Potential is only valid in 3 dimensions"
					errorflag = 1
					return 
				end if
				if (trim(method).ne."CCS") then
					print *, "Coulomb Potential can only be simulated by CCS"
					errorflag = 1
					return  
				end if 
				if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
					print *, "This method must use a static grid."
					errorflag = 1
					return
				end if
			case ("HH")
				if (npes.ne.1) then
					print *, "Inverted Gaussian only valid for 1 PES"
					errorflag = 1
					return
				end if
				if (trim(method).ne."CCS") then
					print *, "Henon-Heiles Potential can only be simulated by CCS"
					errorflag = 1
					return  
				end if 
				if ((ndim.ne.2).and.(ndim.ne.6).and.(ndim.ne.10)) then
					print  *, "Henon-Heiles potential is only valid currently for the 2,6,or 10 dimensional systems"
					errorflag = 1
					return 
				end if
				if (basis.eq."GRID") then
					print *, "This method must not use a static grid."
					errorflag = 1
					return
				end if 
			case default
				print *, "System is not recognised. Value is ", sys
				errorflag = 1
				return
		end select             

		if (n.ne.1) then
			print *, "Not all required variables read in readsys subroutine"
			errorflag = 1
			return
		end if

		call readparams

		return

	end subroutine readsys

!--------------------------------------------------------------------------------------------------

	subroutine readecut   !   Level 1 Subroutine

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

			if(LINE=='ECheck') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,ECheck
				if (ierr.ne.0) then
					print *, "Error reading ECheck value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='Ntries') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,Ntries
				if (ierr.ne.0) then
					print *, "Error reading Ntries value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='Ebfmin') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,Ebfmin
				if (ierr.ne.0) then
					print *, "Error reading Ebfmin value"
					errorflag = 1
					return
				end if
				n = n+1
			else if(LINE=='Ebfmax') then
				backspace(128)
				read(128,*,iostat=ierr)LINE,Ebfmax
				if (ierr.ne.0) then
					print *, "Error reading Ebfmax value"
					errorflag = 1
					return
				end if
				n = n+1
			end if

			read(128,*,iostat=ierr) LINE

		end do

		close (128)

		if ((ECheck.ne.'NO').and.(ECheck.ne.'YES').and.(ECheck.ne.'No').and.(ECheck.ne.'Yes') &
			  .and.(ECheck.ne.'yes').and.(ECheck.ne.'no').and.(ECheck.ne.'Y').and.(ECheck.ne.'N') &
			  .and.(ECheck.ne.'y').and.(ECheck.ne.'n')) then
			print *, "Invalid value for ECheck. Must be YES/NO/Yes/No/yes/no/Y/N/y/n. Value is ", ECheck
			errorflag = 1
			return
		else if ((ECheck.eq.'NO').or.(ECheck.eq.'No').or.(ECheck.eq.'no').or.(ECheck.eq.'N') &
					 .or.(ECheck.eq.'n')) then
			ECheck = 'NO'
		else if ((ECheck.eq.'YES').or.(ECheck.eq.'Yes').or.(ECheck.eq.'yes').or.(ECheck.eq.'Y') &
					 .or.(ECheck.eq.'y')) then
			ECheck = 'YES'
		end if

		if (Ebfmin.ge.Ebfmax) then
			print *, "Invalid values for Ebfmin and/or Ebfmax. Max must be higher than min"
			errorflag = 1
			return
		end if

		if (n.ne.4) then
			print *, "Not all required variables read in readecut subroutine"
			errorflag = 1
			return
		end if

		return

	end subroutine readecut


!*************************************************************************************************!
!          Reading Subroutines for Basis Set Generation
!*************************************************************************************************!

	subroutine readbsparams   !   Level 1 Subroutine

		IMPLICIT NONE
		character(LEN=100)::LINE
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		OPEN(UNIT=127, FILE='input.dat',STATUS='OLD', iostat=ierr)

		if (ierr .ne. 0) then
			print *, 'error in opening input.dat file'
			errorflag = 1
			return
		end if

		read(127,*,iostat=ierr)LINE

		do while (ierr==0)

			if(LINE== "ndim") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,ndim 
				if(ierr.ne.0) then
					Print *,  "Error reading ndim"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="in_nbf") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,in_nbf
				if(ierr.ne.0) then
					Print *,  "Error reading in_nbf"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="matfun") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,matfun
				if(ierr.ne.0) then
					print *,  "Error reading matrix function"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="npes") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,npes
				if(ierr.ne.0) then
					Print *,  "Error reading npes"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="in_PES") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,in_pes
				if(ierr.ne.0) then
					print *,  "Error reading in_PES"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="basis") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,basis
				if(ierr.ne.0) then
					Print *,  "Error reading basis option"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="gridsp") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,initsp
				if(ierr.ne.0) then
					print *,  "Error reading grid spacing"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="psizex") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,psizex
				if(ierr.ne.0) then
					print *,  "Error reading grid size in p coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="qsizex") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,qsizex
				if(ierr.ne.0) then
					print *,  "Error reading grid size in q coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="psizey") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,psizey
				if(ierr.ne.0) then
					print *,  "Error reading grid size in p coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="qsizey") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,qsizey
				if(ierr.ne.0) then
					print *,  "Error reading grid size in q coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="psizez") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,psizez
				if(ierr.ne.0) then
					print *,  "Error reading grid size in p coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="qsizez") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,qsizez
				if(ierr.ne.0) then
					print *,  "Error reading grid size in q coordinate"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="trainsp") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,trainsp
				if(ierr.ne.0) then
					print *,  "Error reading train spacing"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="def_stp") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,def_stp
				if(ierr.ne.0) then
					print *,  "Error reading default number of basis functions per train"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="nbfadapt") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,nbfadapt
				if(ierr.ne.0) then
					print *,  "Error reading bf adapt flag"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="nbfepsilon") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,bfeps
				if(ierr.ne.0) then
					print *,  "Error reading bf adapt cutoff parameter"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="Cloning") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,cloneflg
				if(ierr.ne.0) then
					print *,  "Error reading cloning flag"
					errorflag = 1
					return
				end if
				n = n+1
			end if

			read(127,*,iostat=ierr) LINE

		end do

		close(127) 
	  
		if ((in_pes.gt.npes).or.(in_pes.le.0)) then
			print *, "Initial PES does not exist"
			ierr=-1
			errorflag = 1
			return
		end if

		if (in_nbf.le.0) then
			print *, "Number of Basis Functions <= 0"
			ierr=-1
			errorflag = 1
			return
		end if

		if (npes.le.0) then
			print *, "Number of PES <= 0"
			ierr=-1
			errorflag = 1
			return
		end if

		if ((npes.eq.1).and.(trim(method)/="CCS")) then
			print *, "Only one PES selected, but propagation method is not CCS!"
			ierr=-1
			errorflag=1
			return
		end if

		if (ndim.le.0) then
			print *, "Number of Degrees of Freedom <= 0"
			ierr=-1
			errorflag = 1
			return
		end if

		if ((basis.ne.'TRAIN').and.(basis.ne.'train').and.(basis.ne.'SWARM').and.(basis.ne.'swarm').and.(basis.ne.'SWTRN')&
			.and.(basis.ne.'swtrn').and.(basis.ne.'GRID').and.(basis.ne.'grid').and.(basis.ne.'GRSWM').and.(basis.ne.'grswm')) then
			print *, "Invalid value for basis. Must be TRAIN/SWARM/GRID/SWTRN/GRSWM and all upper/lower case Value is ", basis
			errorflag = 1
			return
		else if ((basis.eq.'TRAIN').or.(basis.eq.'train'))then
			basis = 'TRAIN'
		else if ((basis.eq.'SWARM').or.(basis.eq.'swarm'))then
			basis = 'SWARM'
		else if ((basis.eq.'GRID').or.(basis.eq.'grid'))then
			basis = 'GRID'
		else if ((basis.eq.'SWTRN').or.(basis.eq.'swtrn'))then
			basis = 'SWTRN'
		else if ((basis.eq.'GRSWM').or.(basis.eq.'grswm'))then
			basis = 'GRSWM'
		end if

		if ((matfun.ne.'zgesv').and.(matfun.ne.'ZGESV').and.(matfun.ne.'zheev').and.(matfun.ne.'ZHEEV')) then
			print *, "Invalid value for matrix function. Must be ZGESV/zgesv or ZHEEV/zheev. Value is ", matfun
			errorflag = 1
			return
		else if ((matfun.eq.'zgesv').or.(matfun.eq.'ZGESV')) then
			matfun = 'zgesv'
		else
			matfun = 'zheev'
		end if

		if (basis.eq."GRID") then
			if (initsp .lt. 0.8d0) then
				print *, "Error! Grid points are too close together"
				errorflag=1
				return
			end if
			if ((ndim.ne.1).and.(ndim.ne.3)) then
				print *, "ndim is neither 1 nor 3. This is currently invalid."
				errorflag=1
				return
			else if (ndim==1) then
				if ((qsizex.gt.qsizey).and.(qsizex.gt.qsizez).and.(psizex.gt.psizey).and.(psizex.gt.psizez)) then
					qsizey=0
					qsizez=0
					psizey=0
					psizez=0
				else if ((qsizey.gt.qsizex).and.(qsizey.gt.qsizez).and.(psizey.gt.psizex).and.(psizey.gt.psizez)) then
					qsizex=0
					qsizez=0
					psizex=0
					psizez=0
				else if ((qsizez.gt.qsizex).and.(qsizez.gt.qsizey).and.(psizez.gt.psizex).and.(psizez.gt.psizey)) then
					qsizex=0
					qsizey=0
					psizex=0
					psizey=0
				else 
					print *, "Error! The largest grid dimensions do not match. Check the values."
					errorflag = 1
					return
				end if
			end if
			if ((mod(qsizez,2).ne.mod(psizez,2)).and.(mod(qsizey,2).ne.mod(psizey,2)).and.(mod(qsizex,2).ne.mod(psizex,2))) then
				print *, "Error! Grid is not symmetrically spaced around initial CS."
				print *, "There should be equal distance between initial CS and the four closest grid points"
				print *, "This equates to  all psize and qsize values being either all even or all odd."
				errorflag=1
				return
			end if
			if ((max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)).ne.in_nbf) then
				print *, "Warning! Grid size does not match in_nbf. in_nbf should be product of all psize and qsize"
				if (mod(qsizez,2)==1) then
					in_nbf = (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1))
					print "(a,i0,a)", "Altering in_nbf to ", &
							 (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)), "..."
				else
					in_nbf = (max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)) + 1
					print "(a,i0,a)", "Altering in_nbf to ", &
									(max(qsizex,1)*max(psizex,1)*max(qsizey,1)*max(psizey,1)*max(qsizez,1)*max(psizez,1)) + 1, "..."
				end if
			end if
		end if

		if (basis.eq."GRSWM") then
			if (initsp .lt. 0.8d0) then
				print *, "Error! Grid points are too close together"
				errorflag=1
				return
			end if
			if (ndim.ne.3) then
				print *, "ndim is not 3. This is currently invalid."
				errorflag=1
				return
			end if
			if (mod(qsizez,2).ne.mod(psizez,2)) then
				print *, "Error! Grid is not symmetrically spaced around initial CS."
				print *, "There should be equal distance between initial CS and the four closest grid points in the grid dimension (z)"
				print *, "This equates to both psizez and qsizez values being even or odd."
				errorflag=1
				return
			end if
			if ((qsizez*psizez).ne.in_nbf) then
				print *, "Warning! Grid size does not match in_nbf. in_nbf should be product of psizez and qsizez"
				if (mod(qsizez,2)==1) then
					in_nbf = qsizez*psizez
					print "(a,i0,a)", "Altering in_nbf to ", qsizez*psizez, "..."
				else
					in_nbf = qsizez*psizez + 1
					print "(a,i0,a)", "Altering in_nbf to ", qsizez*psizez+1, "..."
				end if
			end if
		end if

		if ((basis.eq."TRAIN").or.(basis.eq."SWTRN")) then
			if (step.eq."A") then
				print *, "Trains are only valid for static stepsizes."
				errorflag = 1
				return
			end if
			if (trainsp.le.0) then
				print *, "Error! Spacing between basis set train elements is <=0"
				errorflag = 1
				return
			end if
			if ((method.ne."MCEv2").and.(method.ne."CCS")) then
				print *, "Error! Trains can only work with MCEv2 or CCS. MCEv1 cannot calculate the amplitudes correctly"
				errorflag = 1
				return
			end if        
		end if

		if ((nbfadapt.ne.'NO').and.(nbfadapt.ne.'YES').and.(nbfadapt.ne.'No').and.(nbfadapt.ne.'Yes') &
			  .and.(nbfadapt.ne.'yes').and.(nbfadapt.ne.'no').and.(nbfadapt.ne.'Y').and.(nbfadapt.ne.'N') &
			  .and.(nbfadapt.ne.'y').and.(nbfadapt.ne.'n')) then
			print *, "Invalid value for nbfadapt. Must be YES/Yes/yes/Y/y or NO/No/no/N/n. Value is ", nbfadapt
			errorflag = 1
			return
		else if ((nbfadapt.eq.'NO').or.(nbfadapt.eq.'No').or.(nbfadapt.eq.'no').or.(nbfadapt.eq.'N') &
					 .or.(nbfadapt.eq.'n')) then
			nbfadapt = 'NO'
		else if ((nbfadapt.eq.'YES').or.(nbfadapt.eq.'Yes').or.(nbfadapt.eq.'yes').or.(nbfadapt.eq.'Y') &
					 .or.(nbfadapt.eq.'y')) then
			nbfadapt = 'YES'
		end if

		if ((cloneflg.ne.'NO').and.(cloneflg.ne.'YES').and.(cloneflg.ne.'No').and.(cloneflg.ne.'Yes') &
			  .and.(cloneflg.ne.'yes').and.(cloneflg.ne.'no').and.(cloneflg.ne.'Y').and.(cloneflg.ne.'N') &
			  .and.(cloneflg.ne.'y').and.(cloneflg.ne.'n')) then
			print *, "Invalid value for nbfadapt. Must be YES/Yes/yes/Y/y or NO/No/no/N/n. Value is ", nbfadapt
			errorflag = 1
			return
		else if ((cloneflg.eq.'NO').or.(cloneflg.eq.'No').or.(cloneflg.eq.'no').or.(cloneflg.eq.'N') &
					 .or.(cloneflg.eq.'n')) then
			cloneflg = 'NO'
		else if ((cloneflg.eq.'YES').or.(cloneflg.eq.'Yes').or.(cloneflg.eq.'yes').or.(cloneflg.eq.'Y') &
					 .or.(cloneflg.eq.'y')) then
			cloneflg = 'YES'
		end if

		if (nbfadapt.eq."YES") then
			if ((basis.ne."GRID").and.(basis.ne."GRSWM")) then
				print *, "Adaptive basis set chosen but gridding disabled. This is currently an invalid combination."
				print *, "basis should be GRID or GRSWM. Value was ", basis 
				errorflag=1
				return
			end if
			if (bfeps.lt.0.0d0) then
				print *, "Basis set adaptive cutoff parameter less than zero. This is not valid"
				errorflag=1
				return
			end if
		end if     

		initsp=initsp*dsqrt(2.0d0)      

		if (n.ne.18) then
			print *, "Not all required variables read in readbsparams subroutine. n=", n
			errorflag = 1
			return
		end if

		return

	end subroutine readbsparams

!--------------------------------------------------------------------------------------------------

	subroutine readzparams   !   Level 1 Subroutine

		IMPLICIT NONE
		character(LEN=100)::LINE
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		OPEN(UNIT=127, FILE='input.dat',STATUS='OLD', iostat=ierr)

		if (ierr .ne. 0) then
			print *, 'error in opening input.dat file'
			errorflag = 1
			return
		end if

		read(127,*,iostat=ierr)LINE

		do while (ierr==0)

			if (LINE=="ALCMP") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,initalcmprss
				if(ierr.ne.0) then
					print *,  "Error reading compression parameter"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="gamma") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,gam
				if(ierr.ne.0) then
					print *,  "Error reading gamma factor"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="mu") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,mu
				if(ierr.ne.0) then
					print *,  "Error reading mu (centre of initial random gaussian)"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="hbar") then
				backspace(127)
				read(127,*,iostat=ierr)LINE,hbar
				if(ierr.ne.0) then
					print *,  "Error reading hbar value. This value is optional, defaulting to 1"
					hbar = 1.0d0
				else if (hbar.ne.1.0d0) then 
					print *, "hbar changed from default to ", hbar
				end if
			end if

			read(127,*,iostat=ierr) LINE

		end do

		close(127)

		sigp = sqrt(1.0d0/(2.0d0*gam))
		sigq = sqrt(gam/2.0d0)

		if (n.ne.3) then
			print *, "Not all required variables read in readzparams subroutine."
			errorflag = 1
			return
		end if

		return

	end subroutine readzparams

!*************************************************************************************************!
!          Reading Subroutines for Basis Set Propagation
!*************************************************************************************************!

	subroutine readbasis(bs, mup, muq, rep, t)   !   Level 1 Subroutine

		implicit none
		type(basisfn), dimension (:), allocatable, intent(inout) :: bs
		integer::ierr, n, j, k, m, r, cflg
		real(kind=8), dimension(:), allocatable, intent(out) :: mup, muq 
		real(kind=8), intent(inout) :: t
		integer, intent(in) :: rep  
		character(LEN=100)::LINE
		character(LEN=13)::filename
		real(kind=8)::rl, im
		complex(kind = 8) :: dsum1

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0
		cflg = 0

		print "(a)","Starting read subroutine"

		if (rep.lt.10) then
			write (filename,"(a,i1,a)") "Outbs-00", rep, ".out"
		else if (rep.lt.100) then
			write(filename,"(a,i2,a)") "Outbs-0", rep, ".out"
		else
			write(filename,"(a,i3,a)") "Outbs-", rep, ".out"
		end if

		print *, "Opening file ", trim(filename)

		open(unit=200, file=filename, status="old", iostat=ierr)

		if (ierr .ne. 0) then
			print *, 'error in opening Outbs.out file'
			errorflag = 1
			return
		end if

		read(200,*,iostat=ierr)LINE

		do while ((LINE.ne."zinit").and.(ierr==0))
			if (LINE=="ndof") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,ndim
				if(ierr.ne.0) then
					print *,  "Error reading ndim"
					errorflag = 1
					return
				end if
				print *, "ndim   = ", ndim
				n = n+1
			else if (LINE=="nconf") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,npes
				if(ierr.ne.0) then
					print *,  "Error reading npes"
					errorflag = 1
					return
				end if
				print *, "npes   = ", npes
				n = n+1
			else if (LINE=="nbasisfns") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,in_nbf
				if(ierr.ne.0) then
					print *,  "Error reading in_nbf"
					errorflag = 1
					return
				end if
				print *, "in_nbf    = ", in_nbf
				n = n+1
			else if (LINE=="initial_PES") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,in_pes
				if(ierr.ne.0) then
					print *,  "Error reading in_PES"
					errorflag = 1
					return
				end if
				print *, "in_pes = ", in_pes
				n = n+1
			else if (LINE=="matfun") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,matfun
				if(ierr.ne.0) then
					print *,  "Error reading matrix function"
					errorflag = 1
					return
				end if
				print *, "matfun = ", matfun
				n = n+1
			else if (LINE=="time") then
				backspace(200)
				read(200,*,iostat=ierr)LINE,t
				if(ierr.ne.0) then
					print *,  "Error reading time"
					errorflag = 1
					return
				end if
				print *, "time = ", t
				n = n+1
			end if
			read(200,*,iostat=ierr)LINE
		end do

		if (n.ne.6) then
			print "(a,i2,a)", "Error in reading parameters. Only ", n, " of 6 parameters read."
			errorflag = 1
			return
		end if 

		allocate (mup(ndim), stat=ierr)
		if (ierr == 0) allocate (muq(ndim), stat=ierr)
		if (ierr/=0) then
			print *, "Error in allocation of mup and muq"
			errorflag=1
		end if

		if (LINE.ne."zinit") then
			print *, "Error! Expected zinit, but read ", trim(LINE)
		end if
		backspace(200)

		do m=1,ndim
			read(200,*,iostat=ierr)LINE, j, muq(j), mup(j)
			if(ierr.ne.0) then
				print *,  "Error reading zinit value ", m
				errorflag = 1
				return
			end if
			if(m.ne.j) then
				print *,  "Error! Count mismatch in zinit. Expected ", m, "but read ", j
				errorflag = 1
				return
			end if
		end do
		 
		read(200,*,iostat=ierr)LINE

		if (LINE.ne."basis") then
			print *, "Error! Expected basis, but read ", trim(LINE)
		end if

		backspace(200)

		if (size(bs).ne.in_nbf) then
			print *, "Basis set size has changed. Reallocating basis set."
			call deallocbs(bs)
			call allocbs(bs, in_nbf)
		end if

		do j=1,in_nbf
			read(200,*,iostat=ierr)LINE,k
			if(k.ne.j) then
				print "(a,i2,a,i2)", "Error. Expected basis function ", j, " but got ", k
			end if
			read (200,*,iostat=ierr)LINE
			if (LINE.ne."D") then
				print *, "Error! Expected D but read ", trim(LINE)
			end if
			backspace(200)
			read(200,*,iostat=ierr)LINE,rl, im
			if (LINE.eq."D") then
				bs(j)%D_big=cmplx(rl,im,kind=8)
			else
				print *, "Something has gone very wrong here"
				errorflag = 1
				return
			end if
			do r=1,npes
				read(200,*,iostat=ierr)LINE
				if (LINE.ne."a") then
					print *, "Error! Expected a but read ", trim(LINE)
				end if
				backspace(200)
				read(200,*,iostat=ierr)LINE,k,rl, im
				if (k.ne.r) then
					print "(a,i2,a,i2)", "Error. Expected a from pes ", r, "but got ", k
				end if
				bs(j)%a_pes(r)=cmplx(rl,im,kind=8)
			end do
			do r=1,npes
				read(200,*,iostat=ierr)LINE
				if (LINE.ne."d") then
					print *, "Error! Expected d but read ", trim(LINE)
				end if
				backspace(200)
				read(200,*,iostat=ierr)LINE,k,rl, im
				if (k.ne.r) then
					print "(a,i2,a,i2)", "Error. Expected d from pes ", r, "but got ", k
				end if
				bs(j)%d_pes(r)=cmplx(rl,im,kind=8)
			end do
			do r=1,npes
				read(200,*,iostat=ierr)LINE
				if (LINE.ne."s") then
					print *, "Error! Expected s but read ", trim(LINE)
				end if
				backspace(200)
				read(200,*,iostat=ierr)LINE,k,rl
				if (k.ne.r) then
					print "(a,i2,a,i2)", "Error. Expected s from pes ", r, "but got ", k
				end if
				bs(j)%s_pes(r)=rl
			end do
			do m=1,ndim
				read(200,*,iostat=ierr)LINE
				if (LINE.ne."z") then
					print *, "Error! Expected z, but read ", trim(LINE)
				end if
				backspace(200)
				read(200,*,iostat=ierr)LINE,k,rl, im
				if (k.ne.m) then
					print "(a,i2,a,i2)", "Error. Expected z dimension ", m, "but got ", k
				end if
				bs(j)%z(m)=cmplx(rl,im,kind=8)
			end do
		end do

		if (t==0.0d0) then
			if (method.eq."MCEv1") then
				do j=1,in_nbf
					do r=1,npes
						bs(j)%d_pes(r) = bs(j)%d_pes(r) * bs(j)%D_big
						bs(j)%a_pes(r) = bs(j)%d_pes(r) * cdexp(i*bs(j)%s_pes(r))
					end do
					bs(j)%D_big = (1.0d0,0.0d0)
				end do
			else
				do j=1,in_nbf
					dsum1 = (0.0d0,0.0d0)
					do r=1,npes
						dsum1 = dsum1 + bs(j)%d_pes(r)
						if (r.eq.in_pes) then
							bs(j)%d_pes(r) = (1.0d0,0.0d0)
						else
							bs(j)%d_pes(r) = (0.0d0,0.0d0)
						end if
						bs(j)%a_pes(r) = bs(j)%d_pes(r) * cdexp(i*bs(j)%s_pes(r))
					end do
					bs(j)%D_big = dsum1 * bs(j)%D_big
				end do       
			end if
		else
			do j=1,in_nbf
				if ((dble(bs(j)%D_big) /= 1.0d0).and.(method.eq."MCEv1")) then
					print *, "The D_big amplitudes are not compatible with MCEv1 propagation for t/=0"
					errorflag=1
					return
				else if ((dble(bs(j)%d_pes(1)) /= 1.0d0).and.(method.eq."CCS")) then
					print *, "The d_pes amplitudes are not compatible with CCS propagation for t/=0"
					errorflag=1
					return
				end if
			end do
		end if   
			 
	end subroutine readbasis

!--------------------------------------------------------------------------------------------------

	subroutine readtimepar   !   Level 1 Subroutine

		implicit none
		character(LEN=100)::LINE, LINE2
		integer::ierr, n

		if (errorflag .ne. 0) return

		ierr = 0
		n = 0

		OPEN(UNIT=135, FILE='prop.dat',STATUS='OLD', iostat=ierr)

		if (ierr .ne. 0) then
			print *, 'error in opening prop.dat file'
			errorflag = 1
			return
		end if

		read(135,*,iostat=ierr)LINE

		do while (ierr==0)

			if(LINE== "dtmin") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,dtmin 
				if(ierr.ne.0) then
					Print *,  "Error reading minimum dt"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="dtmax") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,dtmax
				if(ierr.ne.0) then
					Print *,  "Error reading maximum dt"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="dtinit") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,dtinit
				if(ierr.ne.0) then
					print *,  "Error reading initial dt"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="time_end") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,timeend
				if(ierr.ne.0) then
					print *,  "Error reading end time for propagation"
					errorflag = 1
					return
				end if
				n = n+1    
			else if (LINE=="time_start") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,timestrt
				if(ierr.ne.0) then
					print *,  "Error reading starting time for propagation"
					errorflag = 1
					return
				end if
				n = n+1
			else if (LINE=="step") then
				backspace(135)
				read(135,*,iostat=ierr)LINE,LINE2
				if(ierr.ne.0) then
					print *,  "Error reading step type for propagation"
					errorflag = 1
					return
				end if
				if ((LINE2(1:1).eq.'S').or.(LINE2(1:1).eq.'s')) then
					step = "S"
				else if ((LINE2(1:1).eq.'A').or.(LINE2(1:1).eq.'a')) then
					step = "A"
				else
					print *, "Error reading step type. Expected 'Static' or 'Adaptive', but got ", trim(LINE2)
					errorflag = 1
				end if
				n = n+1
			end if    

			read(135,*,iostat=ierr)LINE

		end do

		if (n.ne.6) then
			print *, "Not all required variables read in readtimepar subroutine."
			errorflag = 1
			return
		end if

		return

	end subroutine readtimepar
		
!*************************************************************************************************!

END MODULE readpars
