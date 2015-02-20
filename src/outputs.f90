MODULE outputs

  use globvars

!*************************************************************************************************!
!*
!*         Output Module
!*           
!*   Contains subroutines for:
!*
!*      1) Generating and outputting a histogram using a complex variable $$
!*      2) Generating and outputting a histogram using a real variable
!*      3) Generating and outputting a two-dimensional histogram using a complex variable $$
!*      4) Outputting the basis set values and parameters for later use
!*      5) Creating a file with headers for outputting the simulation data for a single thread 
!*      6) Outputting the simulation data at a single timestep for a single thread 
!*      7) Outputting the simulation data for all timesteps over all threads
!*      8) Creating a file with headers and plotting scripts for outputting the wavefunction values for a single thread
!*      9) Outputting the wavefunction values at a single timestep for a single thread
!*     10) Creating a file with headers and plotting scripts for outputting the trajectories for a single thread
!*     11) Outputting the trajectory values at a single timestep for a single thread
!*     12) Interpolation controller subroutine for averaging all adaptive timestep data for a single multi-thread run
!*     13) Neville Algorithm subroutine which controls the polynomial interpolation for a single thread
!*     14) Polynomial Interpolation subroutine which calculates the value of a single datapoint **
!*     15) Data seeking algorithm which finds the array index of a particular time value in an ordered list **
!*     
!*      $$ - this subroutine is not currently used and was created for debugging purposes
!*      ** - this subroutine was copied verbatim from Press et al Numerical Recipes.
!*      
!*************************************************************************************************!

contains

!*************************************************************************************************!
!        Output Subroutines
!*************************************************************************************************!

  subroutine histogram(A, n, filename, cutup, cutdown)   !   Disused

    implicit none

    integer:: n_box_p, n_box_q, boxn, u, ierr
    integer, dimension(:), allocatable::p_dist, q_dist
    real:: boxn_rl
    real(kind=8),intent(inout)::cutup, cutdown
    complex(kind=8),dimension(:),intent(in)::A
    integer,intent(in)::n
    character(LEN=10),intent(in)::filename

    if (errorflag .ne. 0) return

    ierr = 0
    boxn = 200

    boxn_rl = real(boxn)
    allocate(p_dist(boxn), stat = ierr)
    allocate(q_dist(boxn), stat = ierr)
    if (ierr/=0) then
       print *, "Error in p or q distribution allocation for histogram"
       errorflag=1
       return
    end if


    if ((cutup == 0.0).and.(cutdown == 0.0)) then
       if (maxval(dble(A)).gt.maxval(dimag(A))) then
          cutup = maxval(dble(A))
       else
          cutup = maxval(dimag(A))
       end if

       if (minval(dble(A)).lt.minval(dimag(A))) then
          cutdown = minval(dble(A))
       else
          cutdown = minval(dimag(A))
       end if
    end if

    open(unit=135,file=filename,status='unknown',iostat=ierr)

    do u=1,boxn
       p_dist(u) = 0
       q_dist(u) = 0
    end do

    do u=1,n
       n_box_p = 0
       n_box_q = 0       
       if ((dimag(A(u)).ge.(cutdown-((cutup-cutdown)/(2*boxn_rl)))) &
           .and.(dimag(A(u)).le.(cutup+((cutup-cutdown)/(2*boxn_rl))))) then
          n_box_p=int((((dimag(A(u))-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
          if ((n_box_p.gt.boxn).or.(n_box_p.lt.0)) then
             print *,"Error! Invalid Box Calculated. n_box_p = ", n_box_p
          else if (n_box_p.ne.0) then
             p_dist(n_box_p)=p_dist(n_box_p) + 1
          end if
       end if
       if ((dble(A(u)).ge.(cutdown-((cutup-cutdown)/(2*boxn_rl)))) &
           .and.(dble(A(u)).le.(cutup+((cutup-cutdown)/(2*boxn_rl))))) then
          n_box_q=int((((dble(A(u))-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
          if ((n_box_q.gt.boxn).or.(n_box_q.lt.0)) then
             print *,"Error! Invalid Box Calculated. n_box_q = ", n_box_q
          else if (n_box_q.ne.0) then
             q_dist(n_box_q)=q_dist(n_box_q) + 1
          end if
       end if
    end do

    write (135,*) 'Bin Real_Part Imag_Part'

    do u=1,boxn
       write (135,"(3(1x,e13.5e3))"), ((cutup-cutdown)*real(u)/boxn_rl)+cutdown, &
                         real(q_dist(u)), real(p_dist(u))
    end do

    close(135)

    return

  end subroutine histogram

!--------------------------------------------------------------------------------------------------

  subroutine histogram2(A, n, filename, cutup, cutdown) ! Level 1 Subroutine

    implicit none

    integer:: n_box_p, boxn, u, ierr
    integer, dimension(:), allocatable::p_dist
    real:: boxn_rl
    real(kind=8),intent(inout)::cutup, cutdown
    real(kind=8),dimension(:),intent(in)::A
    integer,intent(in)::n
    character(LEN=12),intent(in)::filename

    if (errorflag .ne. 0) return

    ierr = 0
    boxn = 200

    boxn_rl = real(boxn)
    allocate(p_dist(boxn), stat = ierr)
    if (ierr/=0) then
       print *, "Error in p distribution allocation for histogram"
       errorflag=1
       return
    end if

    if ((cutup == 0.0).and.(cutdown == 0.0)) then
       cutup = maxval(A)
       cutdown = minval(A)
    end if

    open(unit=135,file=filename,status='unknown',iostat=ierr)

    do u=1,boxn
       p_dist(u) = 0
    end do

    do u=1,n
       n_box_p = 0      
       if ((A(u).ge.(cutdown-((cutup-cutdown)/(2*boxn_rl)))) &
           .and.(A(u).le.(cutup+((cutup-cutdown)/(2*boxn_rl))))) then
          n_box_p=int((((A(u)-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
          if ((n_box_p.gt.boxn).or.(n_box_p.lt.0)) then
             print *,"Error! Invalid Box Calculated. n_box_p = ", n_box_p
          else if (n_box_p.ne.0) then
             p_dist(n_box_p)=p_dist(n_box_p) + 1
          end if
       end if
    end do

    write (135,*) 'Bin Freq'

    do u=1,boxn
       write (135,"(2(1x,e13.5e3))"), ((cutup-cutdown)*real(u)/boxn_rl)+cutdown, &
                         real(p_dist(u))
    end do

    close(135)

    return

  end subroutine histogram2

!--------------------------------------------------------------------------------------------------

  subroutine hist2d(A, n, filename, cutup, cutdown)   !   Disused
    implicit none

    integer:: n_box_p, n_box_q, boxn, u, v, ierr
    integer, dimension(:,:), allocatable::dist
    real:: boxn_rl
    real(kind=8),intent(inout)::cutup, cutdown
    real(kind=8)::lowcut, highcut, step
    real(kind=8),dimension(:),allocatable::boxvals
    complex(kind=8),dimension(:),intent(in)::A
    integer,intent(in)::n
    character(LEN=10),intent(in)::filename
    character(LEN=18)::myfmt, myfmt2

    if (errorflag .ne. 0) return

    ierr = 0
    boxn = 100

    if (size(A).ne.n) then
       print *,"Error! Size of array for histogram mismatched!"
       errorflag = 1
       return
    end if

    boxn_rl = real(boxn)
    allocate(dist(boxn,boxn), stat = ierr)
    allocate(boxvals(boxn), stat = ierr)
    if (ierr/=0) then
       print *, "Error in distribution or boxval allocation for 2D histogram"
       errorflag=1
       return
    end if

    if ((cutup == 0.0).and.(cutdown == 0.0)) then
       if (maxval(dble(A)).gt.maxval(dimag(A))) then
          cutup = maxval(dble(A))
       else
          cutup = maxval(dimag(A))
       end if

       if (minval(dble(A)).lt.minval(dimag(A))) then
          cutdown = minval(dble(A))
       else
          cutdown = minval(dimag(A))
       end if
    end if

    if (abs(cutup).gt.abs(cutdown)) then
       cutdown = cutup * (-1.0d0)
    else
       cutup = (-1.0d0) * cutdown
    end if

    do u=1,boxn
       boxvals(u) = ((cutup-cutdown)*real(u)/boxn_rl)+cutdown
    end do

    open(unit=135,file=filename,status='unknown',iostat=ierr)

    dist = 0

    step = ((cutup-cutdown)/(2*boxn_rl))
    lowcut = cutdown-step
    highcut = cutup+step

    do u=1,n
       n_box_p = 0
       n_box_q = 0
       if ((dimag(A(u)).ge.lowcut) .and.(dimag(A(u)).le.highcut)) then
          n_box_p=int((((dimag(A(u))-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
       end if
       if ((dble(A(u)).ge.lowcut) .and.(dble(A(u)).le.highcut)) then
          n_box_q=int((((dble(A(u))-cutdown)*(boxn_rl-1.0))/(cutup-cutdown))+1.5)
       end if
       if ((n_box_p.gt.boxn).or.(n_box_q.gt.boxn).or.(n_box_p.lt.0).or.(n_box_q.lt.0)) then
          print *,"Error! Invalid Box Calculated. n_box_p = ", n_box_p, " n_box_q = ", n_box_q
       else if ((n_box_p.ne.0).and.(n_box_q.ne.0)) then
          dist(n_box_p,n_box_q) = dist(n_box_p,n_box_q) + 1
       end if
    end do

    write (135,"(a,i3,1x,a,1x,e13.5e3,1x,a,1x,e13.5e3,1x,a)") "# ", boxn, "equal boxes from", cutdown, "to", cutup,&
                       ". Real in y direction"

    write (myfmt,'(a,i0,a)') '(', 3, '(1x,e13.5e3))'
    write (myfmt2, '(a, i0, a)'), '(13x', boxn, '(1x,e13.5e3))'

!    write (135,myfmt2) boxvals(1:boxn)    

    do u=1,boxn
       do v=1,boxn
          write (135,myfmt), boxvals(u), boxvals(v), real(dist(u,v))
       end do
       write(135,*)""
    end do

    close(135)

    return

  end subroutine hist2d

!--------------------------------------------------------------------------------------------------

  subroutine outbs(bs, reps, mup, muq, t)   !   Level 1 Subroutine
    implicit none
    type(basisfn), dimension (:), intent(in) :: bs
    real(kind=8), dimension(:), intent(in) :: mup, muq
    real(kind=8), intent(in) :: t
    integer::m, j, r, ierr
    character(LEN=13)::filename
    integer, intent(in) :: reps
    character(LEN=3):: rep

    ierr = 0

    write(rep,"(i3.3)") reps

    if (errorflag.eq.0) then
       filename = "Outbs-"//trim(rep)//".out" 
    else
       filename = "Errbs-"//trim(rep)//".out" 
    end if

    open(unit=135, file=filename, status='unknown', iostat=ierr)

    write(135,"(a,1x,i4)"     ), 'ndof'       , ndim
    write(135,"(a,1x,i4)"     ), 'nconf'      , npes
    write(135,"(a,1x,i4)"     ), 'nbasisfns'  , size(bs)
    write(135,"(a,1x,i4)"     ), 'initial_PES', in_pes
    write(135,"(a,1x,a)"      ), 'matfun'     , matfun
    write(135,"(a,1x,es25.17e3)"), 'time'       , t
    write(135,*),""
    do m=1,ndim
       write(135,"(a,1x,i3,2(1x,es25.17e3))"), 'zinit', m, muq(m), mup(m)
    end do

    do j=1,size(bs)
       write(135,*),""
       write(135,"(a,1x,i5)"),'basis', j
       write(135,"(a,2(1x,es25.17e3))"), "D", dble(bs(j)%D_big), dimag(bs(j)%D_big)
       do r=1,npes
          write(135,"(a,i4,2(1x,es25.17e3))"), "a",r, dble(bs(j)%a_pes(r)), dimag(bs(j)%a_pes(r)) 
       end do  
       do r=1,npes
          write(135,"(a,i4,2(1x,es25.17e3))"), "d",r, dble(bs(j)%d_pes(r)), dimag(bs(j)%d_pes(r)) 
       end do  
       do r=1,npes
          write(135,"(a,i4,1x,es25.17e3)"), "s",r, (bs(j)%s_pes(r))
       end do  
       do m=1,ndim
          write(135,"(a,i4,2(1x,es25.17e3))"), "z",m, dble(bs(j)%z(m)), dimag(bs(j)%z(m))
       end do
    end do

    close(135)

    if (((basis=="GRID").or.(basis=="GRSWM")).and.(ndim==3)) then
       open(unit=136, file="outgrd.out", status="unknown", iostat=ierr)
       do j=1,size(bs)
          if (mod(j,2)==1) write(136,"(3(1x,es25.17e3))"), dble(bs(j)%z(1)), dble(bs(j)%z(2)), dble(bs(j)%z(3))
       end do
       close(136)
    end if

    return

  end subroutine outbs

!--------------------------------------------------------------------------------------------------

  subroutine graphwavefn(DC, x, z0)

    implicit none
    
    complex(kind=8), dimension(:,:), intent(in) :: z0
    real(kind=8), dimension(:), intent(in) :: DC
    integer, intent (in) :: x
    character(LEN=5):: stepstr
    character(LEN=16) :: filenm   
    integer :: k, q, p

    if (errorflag .ne. 0) return

    write(stepstr,"(i5.5)") x

    filenm = "Wavefn-"//trim(stepstr)//".out" 

    open(unit=215,file=filenm)

!    k=0
!    do p=1,psizez
!       do q=1,qsizez
!          k=q+((p-1)*qsizez)
!          write (215,"(3g15.8)") dble(z0(k,1)), dimag(z0(k,1)), DC(k)   
!       end do
!       write (215,*),""
!    end do

    do k=1,in_nbf
       write (215,"(3g15.8)") dble(z0(k,1)), dimag(z0(k,1)), DC(k)
    end do      

    close (215)

    open (unit=216,file="wavefn")

    write(216,"(a)"), "# Set terminal and output"
!    write(216,"(a)") , "set terminal postscript eps enhanced color font 'Helvetica,10'"
    write(216,"(a)") , "set terminal png"    
    write(216,"(a,a,a)") , "set output 'wavefn",stepstr,".png'"!.eps'"

    write(216,"(a)") , "# Set various features of the plot"
    write(216,"(a)") , "set pm3d"
    write(216,"(a)") , "unset surface  # don't need surfaces"
    write(216,"(a)") , "set view map"
    write(216,"(a)") , "set key outside"
    write(216,"(a)") , "set cntrparam cubicspline  # smooth out the lines"
    write(216,"(a)") , "set cntrparam levels 50    # sets the num of contour lines"
    write(216,"(a)") , "set pm3d interpolate 20,20 # interpolate the color"
 
    write(216,"(a)") , "# Set a nice color palette"
    write(216,"(a,a)"), 'set palette model RGB defined ( 0"white", 0.01"black", 1"blue", 2"cyan"',&
                   ',3"green",4"yellow",5"red",8"purple" )'
 
    write(216,"(a)") , "# Axes"
    write(216,"(a)") , "set xlabel 'Q(t)'"
    write(216,"(a)") , "set ylabel 'P(t)'"
!    if (mod(psize,2)==0) then
!       write(216,"(a)") , "set cbrange [-0.005:0.5]"
!    else 
       write(216,"(a)") , "set cbrange [-0.005:1.0]"
!    end if
    write(216,"(a)") , "set format x '%.1f'"
    write(216,"(a)") , "set format y '%.1f'"
    write(216,"(a)") , "set format z '%.2f'"
 
    write(216,"(a)") , "# Now plot"
    write(216,"(a,a,a)"), "splot '", filenm, "' using 1:2:3 notitle with lines lt 1"

    close (216)

    call system("gnuplot wavefn")
    call system("rm wavefn")
 
  end subroutine graphwavefn

!--------------------------------------------------------------------------------------------------

  subroutine outD(bs,x,t)

    implicit none
    type(basisfn),dimension(:),intent(in)::bs
    real(kind=8),intent(in)::t
    integer,intent(in)::x
    real(kind=8),dimension(:),allocatable::absD
    integer :: ierr,j
    character(LEN=18)::myfmt

    if (errorflag .ne. 0) return

    ierr = 0

    write(myfmt,'(a,i0,a)') '(', 1+in_nbf, '(1x,e16.8e3))'

    if (x==2) then
       open (unit=1443,file="D.out",status="new",iostat=ierr)
    else
       open (unit=1443,file="D.out",access="append",status="old",iostat=ierr)
    end if

    do j=1,in_nbf
       absD(j)=sqrt(bs(j)%D_big*dconjg(bs(j)%D_big))
    end do

    write(1443,myfmt) t,(absD(j),j=1,in_nbf)

    close (1443)

  end subroutine outD
       

!--------------------------------------------------------------------------------------------------

  subroutine outnormpopadapheads(reps)   !   Level 1 Subroutine

    implicit none
    integer :: ierr, fileun
    character(LEN=15)::filenm, filenm2
    integer, intent(in) :: reps
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "normpop-"//trim(rep)//".out"  

    fileun = 610+reps

    open(unit=fileun,file=filenm,status='unknown',iostat=ierr)

    if (ierr.ne.0) then
       print *, "Error in initial opening of ", trim(filenm), " file"
       errorflag = 1
       return
    end if

    if (npes==2) then 
       write(fileun,"(a,a)"), "Time Norm Re(ACF(t)) Im(ACF(t)) |ACF(t)| Re(Extra) Im(Extra) |Extra| ",&
                               "Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
       write(fileun,*), ""
       write(fileun,*), ""
    else
       write(fileun,"(a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pops 1...n"
       write(fileun,*), ""
       write(fileun,*), ""
    end if

    close(fileun)

    filenm2="plotacf-"//rep//".gnu"

    open (unit=175,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(175,"(a)") 'set terminal png'
    write(175,"(a,a,a)") 'set output "ACF-',trim(rep),'.png"'
    write(175,"(a)") 'set title "Graph of Autocorrelation Function"'
    write(175,"(a)") 'set xlabel "Time (au)"'
    write(175,"(a)") 'set ylabel "ACF"'
    write(175,"(a,a,a)") 'plot "', filenm, '" u 1:3 t "Re(ACF)" w l, "" u 1:4 t "Im(ACF)" w l, "" u 1:5 t "|ACF|" w l'
    close(175)

    filenm2="plotnrm-"//rep//".gnu"

    open (unit=176,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(176,"(a)") 'set terminal png'
    write(176,"(a,a,a)") 'set output "Norm-',trim(rep),'.png"'
    write(176,"(a)") 'set title "Graph of Norm"'
    write(176,"(a)") 'set xlabel "Time (au)"'
    write(176,"(a)") 'set ylabel "Norm"'
    write(176,"(a,a,a)") 'plot "', filenm, '" u 1:2 t "Norm" w l'
    close(176)

    filenm2="plotext-"//rep//".gnu"

    open (unit=177,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(177,"(a)") 'set terminal png'
    if (sys=="IV") then
       write(177,"(a,a,a)") 'set output "Dipole-',trim(rep),'.png"'
       write(177,"(a)") 'set title "Graph of Dipole Acceleration"'
       write(177,"(a)") 'set ylabel "Dipole Acceleration"'
    else
       write(177,"(a,a,a)") 'set output "Extra-',trim(rep),'.png"'
       write(177,"(a)") 'set title "Graph of Extra Calculated Quantity"'
       write(177,"(a)") 'set ylabel "Extra"'
    end if
    write(177,"(a)") 'set xlabel "Time (au)"'
    write(177,"(a,a,a)") 'plot "', filenm, '" u 1:6 t "Real" w l, "" u 1:7 t "Imaginary" w l, "" u 1:8 t "Absolute" w l'
    close(177)

    filenm2="plotdif-"//rep//".gnu"

    open (unit=178,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(178,"(a)") 'set terminal png'
    write(178,"(a,a,a)") 'set output "PopDiff-',trim(rep),'.png"'
    write(178,"(a)") 'set title "Graph of Population Difference"'
    write(178,"(a)") 'set ylabel "Population Difference"'
    write(178,"(a)") 'set xlabel "Time (au)"'
    write(178,"(a,a,a)") 'plot "', filenm, '" u 1:13 t "Population Difference" w l, "" u 1:2 t "Norm" w l'
    close(178)

  end subroutine outnormpopadapheads

!--------------------------------------------------------------------------------------------------

  subroutine outnormpopadap(norm, acft, extra, ehr, popt, x, reps, time)   !   Level 1 Subroutine   

    implicit none
    real(kind=8), intent (in) :: norm, ehr, time
    complex(kind=8), intent (in) :: acft, extra
    real(kind=8), dimension (:), intent (in) :: popt
    integer, intent (in) :: x
    integer :: ierr, fileun
    character(LEN=18) :: filenm,myfmt
    integer, intent(in) :: reps
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "normpop-"//trim(rep)//".out"     

    fileun = 610+reps  

    open(unit=fileun,file=trim(filenm),status='old',access='append',iostat=ierr)

    if (ierr.ne.0) then
       print "(a,a,a,1x,i5)", "Error opening ", trim(filenm), " file in step", x
       errorflag = 1
       return
    end if

    write(myfmt,'(a,i0,a)') '(', 9+npes, '(1x,e16.8e3))'

    if (npes==2) then
       write(fileun,"(13(1x,e16.8e3))"), time, norm, dble(acft), dimag(acft), abs(acft), &
                             dble(extra), dimag(extra), abs(extra), ehr, &
                             popt(1), popt(2), popt(1)+popt(2), popt(1)-popt(2)
    else
        write(fileun,myfmt), time, norm, dble(acft), dimag(acft), abs(acft), &
                             dble(extra), dimag(extra), abs(extra), ehr, popt(:)
    end if
    

    close (fileun)
  
  end subroutine outnormpopadap

!--------------------------------------------------------------------------------------------------

  subroutine outdimacfheads(reps)   !   Level 1 Subroutine

    implicit none
    integer :: ierr, fileun, m
    character(LEN=21)::filenm, filenm2
    integer, intent(in) :: reps
    character(LEN=4):: rep, dof
    character(LEN=50000) :: lngchar1, lngchar2

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "ndimacf-"//trim(rep)//".out"  

    fileun = 1340+reps

    open(unit=fileun,file=filenm,status='unknown',iostat=ierr)

    if (ierr.ne.0) then
       print *, "Error in initial opening of ", trim(filenm), " file"
       errorflag = 1
       return
    end if

    lngchar1 = "Time Re(ACF(total))"
 
    do m=1,ndim+1
       write(dof,"(i4.4)") m
       lngchar2 = trim(lngchar1)
       lngchar1 = trim(lngchar2)//" Re(ACF("//trim(dof)//")) "
    end do

    lngchar2 = trim(lngchar1)
    lngchar1 = trim(lngchar2)//" Im(ACF(total)) "    

    do m=1,ndim+1
       write(dof,"(i4.4)") m
       lngchar2 = trim(lngchar1)
       lngchar1 = trim(lngchar2)//" Im(ACF("//trim(dof)//")) "
    end do

    lngchar2 = trim(lngchar1)
    lngchar1 = trim(lngchar2)//" |ACF(total)| "

    do m=1,ndim+1
       write(dof,"(i4.4)") m
       lngchar2 = trim(lngchar1)
       lngchar1 = trim(lngchar2)//" |ACF("//trim(dof)//")| "
    end do  

    write(fileun,"(a)") trim(lngchar1)
    write(fileun,*), ""
    write(fileun,*), ""

    close(fileun)

    filenm2="plotrendimacf-"//trim(rep)//".gnu"

    open (unit=175,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(175,"(a)") 'set terminal png'
    write(175,"(a,a,a)") 'set output "RendimACF-',trim(rep),'.png"'
    write(175,"(a)") 'set title "Graph of n 1-dimensional Autocorrelation Functions"'
    write(175,"(a)") 'set xlabel "Time (au)"'
    write(175,"(a)") 'set ylabel "ACF"'
    lngchar1 = 'plot "'//trim(filenm)//'" u 1:2 w l'
    do m=1,ndim
       write(dof,"(i0)") m+2   
       lngchar2 = lngchar1
       lngchar1 = trim(lngchar2)//', "" u 1:'//trim(dof)//' w l'
    end do
    write(175,"(a)") trim(lngchar1)
    close(175)

    filenm2="plotimndimacf-"//trim(rep)//".gnu"

    open (unit=176,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(176,"(a)") 'set terminal png'
    write(176,"(a,a,a)") 'set output "ImndimACF-',trim(rep),'.png"'
    write(176,"(a)") 'set title "Graph of n 1-dimensional Autocorrelation Functions"'
    write(176,"(a)") 'set xlabel "Time (au)"'
    write(176,"(a)") 'set ylabel "ACF"'
    write(dof,"(i0)") 3+ndim
    lngchar1 = 'plot "'//trim(filenm)//'" u 1:'//trim(dof)//' w l'
    do m=1,ndim
       write(dof,"(i0)") m+3+ndim   
       lngchar2 = lngchar1
       lngchar1 = trim(lngchar2)//', "" u 1:'//trim(dof)//' w l'
    end do
    write(176,"(a)") trim(lngchar1)
    close(176)

    filenm2="plotabndimacf-"//trim(rep)//".gnu"

    open (unit=177,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(177,"(a)") 'set terminal png'
    write(177,"(a,a,a)") 'set output "AbndimACF-',trim(rep),'.png"'
    write(177,"(a)") 'set title "Graph of n 1-dimensional Autocorrelation Functions"'
    write(177,"(a)") 'set xlabel "Time (au)"'
    write(177,"(a)") 'set ylabel "ACF"'
    write(dof,"(i0)") 4+2*ndim
    lngchar1 = 'plot "'//trim(filenm)//'" u 1:'//trim(dof)//' w l'
    do m=1,ndim
       write(dof,"(i0)") m+4+2*ndim  
       lngchar2 = lngchar1
       lngchar1 = trim(lngchar2)//', "" u 1:'//trim(dof)//' w l'
    end do
    write(177,"(a)") trim(lngchar1)
    close(177)

  end subroutine outdimacfheads

!--------------------------------------------------------------------------------------------------

  subroutine outdimacf(time,ndimacf,reps)   !   Level 1 Subroutine   

    implicit none
    real(kind=8), dimension (:), intent (in) :: ndimacf
    real(kind=8), intent(in) :: time
    integer :: ierr, fileun
    character(LEN=18) :: filenm,myfmt
    integer, intent(in) :: reps
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "ndimacf-"//trim(rep)//".out"     

    fileun = 1340+reps  

    open(unit=fileun,file=filenm,status='old',access='append',iostat=ierr)

    if (ierr.ne.0) then
       print "(a,a,a)", "Error opening ", trim(filenm), " file"
       errorflag = 1
       return
    end if

    write(myfmt,'(a,i0,a)') '(', 1+size(ndimacf), '(1x,E16.8E3))'

    write(fileun,myfmt), time, ndimacf(:)    

    close (fileun)
  
  end subroutine outdimacf
!--------------------------------------------------------------------------------------------------

  subroutine outnormpopstat(norm, acft, extra, ehr, pops)   !   Level 1 Subroutine

    implicit none
    real(kind=8),dimension(:), intent (in) :: norm, ehr
    complex(kind=8), dimension(:), intent (in) :: acft, extra
    real(kind=8),dimension(:,:), intent(in) :: pops
    real(kind=8) :: time
    integer :: ierr, t
    character(LEN=11) :: filenm
    character(LEN=17) :: myfmt

    if (errorflag .ne. 0) return 

    ierr = 0   

    filenm = "normpop.out"

    open(unit=150,file=trim(filenm),status='unknown',iostat=ierr)

    if (ierr.ne.0) then
       print "(a,1x,i5)", "Error opening ", trim(filenm), " file"
       errorflag = 1
       return
    end if

    if (npes==2) then 
       write(150,"(a,a)"), "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| ", &
                            "Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
       write(150,*), ""
       write(150,*), ""
    else
       write(150,"(a)") "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra) |Extra| Sum(HEhr) Pops 1...n"
       write(150,*), ""
       write(150,*), ""
       write(myfmt,'(a,i0,a)') '(', 9+npes, '(1x,e16.8e3))'
    end if
   
    if (npes==2) then
       do t=1,size(norm)
          time = dtinit * (t-1)
          write(150,"(13(1x,e16.8e3))"), time, norm(t), dble(acft(t)), dimag(acft(t)), abs(acft(t)), &
                   dble(extra(t)), dimag(extra(t)), abs(extra(t)), ehr(t), &
                   pops(t,1), pops(t,2), pops(t,1)+pops(t,2), pops(t,1)-pops(t,2)
       end do
    else
       do t=1,size(norm)
          time = dtinit * (t-1)
          write(150,myfmt), time, norm(t), dble(acft(t)), dimag(acft(t)), abs(acft(t)), &
                   dble(extra(t)), dimag(extra(t)), abs(extra(t)), ehr(t), pops(t,:)
       end do
    end if

    close (150)
  
  end subroutine outnormpopstat

!--------------------------------------------------------------------------------------------------

  subroutine outvarsheads(reps,nbf)

    implicit none
    integer :: ierr, k
    character(LEN=16) :: filenm
    character(LEN=21) :: filenm2
    integer, intent(in) :: reps, nbf
    character(LEN=3):: rep
    character(LEN=5):: klowstr, khighstr, kstr
    character(LEN=50000):: lngchar, lngchar2

    if (errorflag .ne. 0) return

    if (nbfadapt.eq."YES") then
       print *, "Warning! Outvarsheads subroutine called improperly."
       print *, "This output file is incompatible with a non-constant basis set size"
       print *, "File will not be created."
       return
    end if 
       
    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "PropVars-"//rep//".out"

    open(unit=151,file=filenm,status='unknown',iostat=ierr)

    if (ierr.ne.0) then
       print *, "Error in initial opening of PropVars.dat file"
       errorflag = 1
       return
    end if

    write(151,"(a)"), "Time Re(D(1:nbf)) Im(D(1:nbf)) Re(d1(1:nbf)) Im(d1(1:nbf)) Re(d2(1:nbf)) Im(d2(1:nbf)) ",&
                    "S1(1:nbf) S2(1:nbf)"
    write(151,*), ""
    write(151,*), ""

    close(151)

    do k=1,nbf

       write(kstr,"(i0)") k
       write(klowstr,"(i0)") k+1
       write(khighstr,"(i0)") k+1+nbf
    
       if (k==1) then     
          lngchar = 'plot "'//trim(filenm)//'" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       else
          lngchar2 = lngchar
          lngchar = trim(lngchar2)//', "" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       end if
  
    end do

    filenm2="plotamps-"//rep//".gnu"

    open (unit=175,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(175,"(a)") 'set terminal png'
    write(175,"(a,a,a)") 'set output "Amplitudes-',trim(rep),'.png"'
    write(175,"(a)") 'set nokey'
    write(175,"(a)") 'unset key'
    write(175,"(a)") 'set title "Graph of D amplitudes"'
    write(175,"(a)") 'set xlabel "Re(D)"'
    write(175,"(a)") 'set ylabel "Im(D)"'
    write(175,"(a)") trim(lngchar)
    close (175)

    do k=1,nbf

       write(kstr,"(i0)") k
       write(klowstr,"(i0)") k+1+(2*nbf)
       write(khighstr,"(i0)") k+1+(3*nbf)
    
       if (k==1) then     
          lngchar = 'plot "'//trim(filenm)//'" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       else
          lngchar2 = lngchar
          lngchar = trim(lngchar2)//', "" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       end if
  
    end do

    do k=1,nbf

       write(kstr,"(i0)") k
       write(klowstr,"(i0)") k+1+(4*nbf)
       write(khighstr,"(i0)") k+1+(5*nbf)
    
       lngchar2 = lngchar
       lngchar = trim(lngchar2)//', "" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
  
    end do

    filenm2="plotd-"//rep//".gnu"

    open (unit=176,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(176,"(a)") 'set terminal png'
    write(176,"(a,a,a)") 'set output "SCAmplitudes-',trim(rep),'.png"'
    write(176,"(a)") 'set title "Graph of single configuration d amplitudes"'
    write(176,"(a)") 'set nokey'
    write(176,"(a)") 'unset key'
    write(176,"(a)") 'set xlabel "Re(d)"'
    write(176,"(a)") 'set ylabel "Im(d)"'
    write(176,"(a)") trim(lngchar)
    close (176)

    do k=1,nbf

       write(kstr,"(i0)") k
       write(khighstr,"(i0)") k+1+(6*nbf)
    
       if (k==1) then     
          lngchar = 'plot "'//trim(filenm)//'" u 1:'//trim(khighstr)//' w l'
       else
          lngchar2 = lngchar
          lngchar = trim(lngchar2)//', "" u 1:'//trim(khighstr)//' w l'
       end if
  
    end do

    do k=1,nbf

       write(kstr,"(i0)") k
       write(khighstr,"(i0)") k+1+(7*nbf)
    
       lngchar2 = lngchar
       lngchar = trim(lngchar2)//', "" u 1:'//trim(khighstr)//' w l'
  
    end do


    filenm2="plotact-"//rep//".gnu"

    open (unit=177,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(177,"(a)") 'set terminal png'
    write(177,"(a,a,a)") 'set output "Actions-',trim(rep),'.png"'
    write(177,"(a)") 'set title "Graph of Action with time"'
    write(177,"(a)") 'set nokey'
    write(177,"(a)") 'unset key'
    write(177,"(a)") 'set xlabel "Time"'
    write(177,"(a)") 'set ylabel "Action"'
    write(177,"(a)") trim(lngchar)
    close (177)



  end subroutine outvarsheads

!--------------------------------------------------------------------------------------------------

  subroutine outvars(bs, x, reps, time)

    implicit none
    type(basisfn), dimension (:), intent (in) :: bs
    integer, intent(in) :: reps
    real(kind=8), intent(in) :: time
    integer, intent (in) :: x
    real(kind=8), dimension(:), allocatable :: s1ar, s2ar, imdbigar, rldbigar, rld1ar, imd1ar, rld2ar, imd2ar
    integer::ierr, k
    character(LEN=16) :: filenm
    character(LEN=18) :: myfmt
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    if (nbfadapt.eq."YES") return

    write(rep,"(i3.3)") reps

    filenm = "PropVars-"//trim(rep)//".out" 

    allocate(s1ar(size(bs)), stat = ierr)
    if (ierr==0) allocate(s2ar(size(bs)), stat = ierr)
    if (ierr==0) allocate(imdbigar(size(bs)), stat = ierr)
    if (ierr==0) allocate(rldbigar(size(bs)), stat = ierr)
    if (ierr==0) allocate(rld1ar(size(bs)), stat = ierr)
    if (ierr==0) allocate(imd1ar(size(bs)), stat = ierr)
    if (ierr==0) allocate(rld2ar(size(bs)), stat = ierr)
    if (ierr==0) allocate(imd2ar(size(bs)), stat = ierr)
    if (ierr/=0) then
       print *, "Error allocating arrays for output of basis set variables in step ", x
       errorflag = 1
       return
    end if

    do k=1,size(bs)
       rldbigar(k) = dble(bs(k)%D_big)
       imdbigar(k) = dimag(bs(k)%D_big) 
       rld1ar(k) = dble(bs(k)%d_pes(1))
       imd1ar(k) = dimag(bs(k)%d_pes(1))
       s1ar(k) = bs(k)%s_pes(1)
       if (npes.gt.1) then
          rld2ar(k) = dble(bs(k)%d_pes(2))
          imd2ar(k) = dimag(bs(k)%d_pes(2))
          s2ar(k) = bs(k)%s_pes(2)
       else if (npes==1) then
          rld2ar(k) = 0.0d0
          imd2ar(k) = 0.0d0
          s2ar(k) = 0.0d0          
       end if       
    end do

    open(unit=151,file=filenm,status='old',access='append',iostat=ierr)

    if (ierr.ne.0) then
       print "(a,a,a,1x,i5)", "Error opening ", filenm, " file for step", x
       errorflag = 1
       return
    end if

    write(myfmt,"(a,i0,a)") "(", 8*size(bs) + 1, "(1x,e13.5e3))"

    write(151,myfmt), time, rldbigar(1:size(bs)), imdbigar(1:size(bs)), rld1ar(1:size(bs)), imd1ar(1:size(bs)), &
            rld2ar(1:size(bs)), imd2ar(1:size(bs)), s1ar(1:size(bs)), s2ar(1:size(bs))   

    close(151)

    deallocate(s1ar, stat = ierr)
    if (ierr==0) deallocate(s2ar, stat = ierr)
    if (ierr==0) deallocate(imdbigar, stat = ierr)
    if (ierr==0) deallocate(rldbigar, stat = ierr)
    if (ierr==0) deallocate(rld1ar, stat = ierr)
    if (ierr==0) deallocate(imd1ar, stat = ierr)
    if (ierr==0) deallocate(rld2ar, stat = ierr)
    if (ierr==0) deallocate(imd2ar, stat = ierr)
    if (ierr/=0) then
       print *, "Error deallocating arrays for output of basis set variables in step ", x
       errorflag = 1
       return
    end if

  end subroutine outvars

!--------------------------------------------------------------------------------------------------

  subroutine outtrajheads(reps, nbf) 
 
    implicit none
    integer :: ierr, k
    character(LEN=16) :: filenm
    character(LEN=21) :: filenm2
    integer, intent(in) :: reps, nbf
    character(LEN=3):: rep
    character(LEN=5):: klowstr, khighstr, kstr
    character(LEN=50000):: lngchar, lngchar2

    if (errorflag .ne. 0) return

    if (nbfadapt.eq."YES") then
       print *, "Warning! Outtrajheads subroutine called improperly."
       print *, "This output file is incompatible with a non-constant basis set size"
       print *, "File will not be created."
       return
    end if 

    ierr = 0

    write(rep,"(i3.3)") reps

    filenm = "Traj-"//trim(rep)//".out"

    open(unit=151,file=filenm,status='unknown',iostat=ierr)

    if (ierr.ne.0) then
       print *, "Error in initial opening of Traj.dat file"
       errorflag = 1
       return
    end if

    write(151,*), "Time q(k=1..nbf) p(k=1..nbf)"
    write(151,*), ""
    write(151,*), ""

    close(151)

    do k=1,nbf

       write(kstr,"(i0)") k
       write(klowstr,"(i0)") k+1
       write(khighstr,"(i0)") k+1+nbf
    
       if (k==1) then     
          lngchar = 'plot "'//trim(filenm)//'" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       else
          lngchar2 = lngchar
          lngchar = trim(lngchar2)//', "" u '//trim(klowstr)//':'//trim(khighstr)//' w l'
       end if
  
    end do

    filenm2="plot-"//trim(filenm)

    open (unit=175,file=filenm2,status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
       print *, 'Error in opening ', trim(filenm2), ' output file'
       errorflag=1
       return
    end if

    write(175,"(a)") 'set terminal png'
    write(175,"(a,a,a)") 'set output "trajectories-',trim(rep),'.png"'
    write(175,"(a)") 'set title "Graph of convergence with different numbers of repetitions"'
    write(175,"(a)") 'set nokey'
    write(175,"(a)") 'unset key'
    write(175,"(a)") 'set xlabel "q(t)"'
    write(175,"(a)") 'set ylabel "p(t)"'
    write(175,"(a)") trim(lngchar)
    close (175)

  end subroutine outtrajheads

!--------------------------------------------------------------------------------------------------

  subroutine outtraj(bs, x, reps, time) 

    implicit none
    type(basisfn), dimension (:), intent (in) :: bs
    real(kind=8), intent(in) :: time
    integer, intent (in) :: x
    real(kind=8), dimension (:), allocatable :: p, q
    integer::ierr, k, fields
    character(LEN=17) :: filenm, myfmt
    integer, intent(in) :: reps
    character(LEN=3):: rep

    if (errorflag .ne. 0) return

    if (nbfadapt.eq."YES") return    

    allocate (p(size(bs)), q(size(bs)))

    ierr = 0

    write(rep,"(i3.3)") reps


    filenm = "Traj-"//trim(rep)//".out"  

    do k=1,size(bs)
       q(k) = sum(dble(bs(k)%z(1:ndim)))
       p(k) = sum(dimag(bs(k)%z(1:ndim)))
    end do

    open(unit=151,file=filenm,status='old',access='append',iostat=ierr)

    if (ierr.ne.0) then
       print "(a,a,a,1x,i5)", "Error opening ", filenm, " file for step", x
       errorflag = 1
       return
    end if

    fields=size(p)+size(q)+1
    write(myfmt,'(a,i0,a)') '(', fields, '(1x,e13.5e3))'

    write(151,myfmt), time, q(1:size(q)), p(1:size(p))  

    close(151)

  end subroutine outtraj

!--------------------------------------------------------------------------------------------------

  subroutine interpolate(cols,errorflag)

    implicit none
    integer, intent (inout) :: errorflag
    integer, intent (in) :: cols
    integer :: fileunin, fileunout, i, j, k, timesteps,ierr, lines, tot
    integer, dimension (:), allocatable :: valid
    real(kind=8), dimension (:), allocatable :: fileline
    real(kind=8), dimension (:,:), allocatable :: input, output, output2
    character (LEN=100) :: filename1, filename2, LINE, LINE2, LINE3
    character (LEN=3) :: repstr, folstr
    character (LEN=14) :: myfmt
    logical :: file_exists
  
    if (errorflag==1) return
  
    allocate (valid(reptot))
    valid = 0
    lines = 0
    if (mod(timeend-timestrt,dtinit).ne.0) timeend=timeend - (mod(timeend-timestrt,dtinit)) + dtinit
    timesteps = (timeend-timestrt)/dtinit + 1
    allocate (output(cols,timesteps))
    output = 0.0d0
    tot = 0
  
    allocate(fileline(cols))
  
    do i=1,reptot
       fileunin = 1150 + i
       fileunout = 7150 + i
  
       lines=0
       
       write(repstr,"(i3.3)") i  
  
       filename1 = "normpop-"//trim(repstr)//".out"
       filename2 = "normpop-"//trim(repstr)//"_interp.out"  
  
       inquire(file=trim(filename1), exist=file_exists)
       if (file_exists) then
          valid(i) = 1
       else
          valid(i) = 0
       end if
  
       if (valid(i)==0) cycle
       tot=tot+1
  
       OPEN(UNIT=fileunin, FILE=trim(filename1),STATUS='OLD', iostat=ierr)
       if (ierr .ne. 0) then
          print *, "Error in opening ", trim(filename1), " file"
          errorflag=1
          return
       end if
  
       OPEN(UNIT=fileunout, FILE=trim(filename2),STATUS='UNKNOWN', iostat=ierr)
       if (ierr .ne. 0) then
          print *, "Error in opening ", trim(filename2), " file"
          errorflag=1
          return
       end if
  
       if (i==1) then
          read(fileunin,"(a141)",iostat=ierr) LINE2
          if (ierr .ne. 0) then
             print *, "Error in reading from ", trim(filename1), " file to get headers"
             errorflag=1  
             return  
          end if
          rewind (fileunin)
       end if
  
       do while (LINE/="0.00000000E+00")
          read(fileunin,*,iostat=ierr) LINE
          if (ierr .ne. 0) then
             print *, "Error in reading from ", trim(filename1), " file to skip to data"
             errorflag=1
             return
          end if
       end do
       backspace (fileunin)
     
       do while (ierr==0)
          LINE3=LINE
          lines = lines + 1
          read(fileunin,*,iostat=ierr) LINE
          if ((LINE==LINE3).and.(LINE/="0.00000000E+00")) lines = lines-1      ! Stops any line being re-read
       end do
  
       allocate(input(cols,lines), stat=ierr)
       if (ierr/=0) then
          print *, "Allocation error in input array"
          errorflag=1
          return
       end if
       
       rewind (fileunin)
       do while (LINE/="0.00000000E+00")
          read(fileunin,*,iostat=ierr) LINE
          if (ierr .ne. 0) then
             print *, "Error in reading from ", trim(filename1), " file to skip to data again"
             errorflag=1
             return
          end if
       end do
       backspace (fileunin)
     
       ierr=0
  
       do j=1,lines
          read(fileunin,*,iostat=ierr) fileline(1:cols)
          if (ierr/=0) then
             print *, "Error writing to input array"
             errorflag=1
             return
          end if
          do k=1,cols
             input(k,j) = fileline(k)
          end do
       end do
  
       call nevillealg(input,output,cols,lines,errorflag)   ! Interpolate subroutine for single file
       if (errorflag==1) return
  
       write(myfmt,'(a,i0,a)') '(', cols, '(1x,e16.8e3))'
  
       do k=1,timesteps
          write(fileunout,myfmt) output(:,k)
       end do     
  
       rewind (fileunout)
  
       close (fileunin)
  
       deallocate(input)
       if (ierr/=0) then
          print *, "Dellocation error in input array"
          errorflag=1
          return
       end if
  
    end do  
  
    deallocate(output)
    if (ierr/=0) then
       print *, "Dellocation error in output array"
       errorflag=1
       return
    end if
  
    OPEN(UNIT=501, FILE="normpop.out",STATUS='new', iostat=ierr) 
  
    write (501,*), LINE2
    write (501,*), ""
    write (501,*), ""
  
    allocate (output2(cols,timesteps))
  
    output2 = 0.0d0
    fileline = 0.0d0  
  
    do i=1,reptot
  
      if (valid(i)==0) cycle
  
      fileunout=7150+i
  
      do k=1,timesteps
         read (fileunout,*,iostat=ierr) fileline
         do j=1,cols
            output2(j,k)=output2(j,k)+fileline(j)
         end do
      end do
  
      close(fileunout)
  
    end do
  
    do k=1,timesteps
       do j=1,cols
          output2(j,k) = output2(j,k)/tot
       end do
    end do 
  
    write(myfmt,'(a,i0,a)') '(', cols, '(1x,e16.8e3))'
  
    do k=1,timesteps
       write(501,myfmt) output2(:,k)
    end do

    close (501)
  
    deallocate(valid,output2,fileline)
  
    stop

  end subroutine interpolate

!--------------------------------------------------------------------------------------------------

  subroutine nevillealg (input,output,cols,lines,errorflag)

    implicit none

    real(kind=8), dimension (:,:), intent (in) :: input
    real(kind=8), dimension (:,:), intent (inout) :: output
    integer, intent (inout) :: errorflag
    integer, intent (in) :: cols, lines
    real(kind=8) :: time, timetmp, dataout, errout
    real(kind=8), dimension (:), allocatable :: timein, datain, timeall
    integer :: i, k, n, z, jdown, jup, jold

    if (errorflag==1) return

    time = timestrt

    allocate(timeall(lines))

    do i=1,lines
       timeall(i) = input(1,i)
    end do

    n=4

    if (lines.le.n) then
       print *, "File Error! There is not enough data in the files to interpolate!"
       errorflag=1
       return
    end if

    if (input(1,1) == timestrt) then
       output(:,1) = input (:,1)
       time = timestrt
    else
       print *, "Something is wrong. Input array does not start at start time"
       errorflag = 1
       return
    end if
 
    if (input(1,2) == timestrt + dtinit) then
       output(:,2) = input (:,2)
       time = time + dtinit
    else
       print *, "I expected the first timestep to be the same as the initial timestep"
       errorflag = 1
       return
    end if

    jold = 2
    z=2
    timetmp=time  

    do while (time.lt.timeend)

       z=z+1

       if (timeend-dtinit.lt.time) then
          timetmp = timeend
       else
          timetmp = timetmp + dtinit
       end if

       jdown = jold

       call hunt(timeall,lines,timetmp,jdown,errorflag)
       if (errorflag==1) return

       jup = jdown + 1

       if (timetmp==timeend) then
          n=3
       end if

       jold = jdown

       if (mod(n,2)==0) then
          jdown = jdown - (n/2) + 1
          jup = jup + (n/2) - 1
       else
          if (((abs(input(1,jup)-timetmp)).lt.(abs(input(1,jdown)-timetmp))).and.(timetmp/=timeend)) then
             jup = jup + (n/2)
             jdown = jdown - (n/2) + 1
          else
             jup = jup + (n/2) - 1
             jdown = jdown - (n/2)
          end if
       end if

       allocate (timein(n), datain(n))

       do i=1,n
          timein(i) = input(1,jdown+i-1)
       end do
  
       do k=2,cols
          do i=1,n
             datain(i) = input(k,jdown+i-1)
          end do

          call polint (timein, datain, n, timetmp, dataout, errout, errorflag)
          if (errorflag==1) return

          output(k,z) = dataout

       end do
     
       time = timetmp

       output(1,z) = time

       deallocate(timein, datain)

    end do  

  end subroutine nevillealg

!--------------------------------------------------------------------------------------------------

  subroutine polint (xa,ya,n,x,y,dy,errorflag)

    implicit none

    integer, intent (in) :: n
    integer, intent (inout) :: errorflag
    integer :: NMAX, i, m, na
    real(kind=8), intent (inout) :: dy, x, y
    real(kind=8), dimension (:), intent (inout) :: xa, ya
    real(kind=8), dimension (:), allocatable :: c, d
    real(kind=8) :: den, dif, dift, ho, hp, w

    if (errorflag==1) return

    na=1
    NMAX=10

    allocate (c(NMAX), d(NMAX))

    dif = abs(x-xa(1))
    do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
          na=i
          dif=dift
       end if
       c(i)=ya(i)
       d(i)=ya(i)
    end do
    y=ya(na)
    na=na-1
    do m=1,n-1
       do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if (den.eq.0) then  
             print *, "Error in polint when t=", x, "for xa values", xa(i)," and", xa(i+m), "when i and m are ",i, "and ", m
             errorflag=1
             return
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
       end do
       if (2*na.lt.n-m) then
          dy=c(na+1)
       else
          dy=d(na)
          na=na-1
       end if
       y=y+dy
    end do
    return

  end subroutine polint

!-------------------------------------------------------------------------------------------------

  subroutine hunt(timeall,lines,timetmp,jdown,errorflag)

    implicit none

    integer, intent (inout) :: jdown, errorflag
    integer, intent (in) :: lines
    real(kind=8), intent (in) :: timetmp
    real(kind=8), dimension(:), intent (in) :: timeall
    integer :: inc, jup, jm
    logical :: ascnd

    if (errorflag==1) return

    ascnd=timeall(lines).gt.timeall(1)

    if ((timetmp.ge.timeall(jdown)).or.(jdown.gt.lines)) then
       jdown = 0
       jup=lines+1
       goto 3
    end if
    inc=1
    if ((timetmp.ge.timeall(jdown)).eqv.ascnd) then
1      jup=jdown+inc
       if (jup.gt.lines) then
          jup=lines+1
       else if ((timetmp.ge.timeall(jup)).eqv.ascnd) then
          jdown=jup
          inc=inc+inc
          goto 1
       end if
    else
       jup = jdown
2      jdown = jup - inc
       if (jdown.lt.1) then
          jdown = 0
       else if ((timetmp.lt.timeall(jdown)).eqv.ascnd) then
          jup=jdown
          inc=inc+inc
          goto 2
       end if
    end if
3   if ((jup - jdown).eq.1) then
       if (timetmp.eq.timeall(lines)) jdown=lines-1
       if (timetmp.eq.timeall(1)) jdown=1
       return
    end if
    jm=(jup+jdown)/2
    if ((timetmp.ge.timeall(jm)).eqv.ascnd) then
       jdown=jm
    else
       jup=jm
    end if
    goto 3

  end subroutine hunt

!*************************************************************************************************!

END MODULE outputs
