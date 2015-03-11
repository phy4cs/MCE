Program avrgnorm

  implicit none
  integer::ierr, i=1, k=1, l=1, n, m, tot=0,folders, folreps, totreps
  real(kind=8),dimension(:),allocatable::pop1, pop2, popsum, popdiff, nrm, rlacf
  real(kind=8),dimension(:),allocatable::imacf, abacf, ehr, time, rlex, imex, abex
  real(kind=8)::pop1av, pop2av, popsumav, popdiffav, nrmav, rlacfav, imacfav
  real(kind=8)::abacfav, rlexav, imexav, abexav, ehrav, timeav
  character(LEN=100)::LINE, filename
  character(LEN=17)::myfmt
  character(LEN=6)::repstr, lstr, timestp
  character(LEN=5000)::lngchar, lngchar2
  integer, dimension(:),allocatable::valid, lines
  logical :: file_exists
           
  call getarg(0,filename)
  call getarg(1,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read first argument for ", trim(filename)
    write (0,"(a)") "This should be the number of seed folders used."
    stop
  end if
  read (LINE,*)folders
  call getarg(2,LINE)
  if (ierr==-1) then
    write (0,"(a,a)") "Error! Could not read second argument for ", trim(filename)
    write (0,"(a)") "This should be the number of repeats used."
    stop
  end if
  read (LINE,*) totreps

  folreps = totreps/folders

  allocate(valid(folders))
  allocate(lines(folders))
  lines = 0

  do i=1,folders
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"    
    inquire(file=trim(filename), exist=file_exists)
    if (file_exists) then
      valid(i) = 1
    else
      valid(i) = 0
    end if
  end do

  do i=1,folders       
    if (valid(i)==0) then
      cycle
    end if
    tot = tot+1       
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"
    OPEN(UNIT=128, FILE=trim(filename),STATUS='OLD', iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(a,a,a)") 'Error in opening ', trim(filename), ' file'
      stop
    end if
    read(128,*,iostat=ierr) LINE
    do while (ierr==0)
      lines(i) = lines(i) + 1
      read(128,*,iostat=ierr) LINE
    end do
    close (128)
    lines(i) = lines(i) - 1
  end do

  do i=1,folders
    if ((lines(i).ne.lines(1)).and.(lines(i).ne.0)) then
      write (0,"(a,i0)") "Error - Line number mismatch in repeat ", i
      write (0,"(a,i0,a,i0)") "Expected ", lines(1), " lines, but got ", lines(i)
      stop
    end if
  end do

  write (6,"(i0,a,i0,a)") tot, " files of ",folders, " were valid."
       
  allocate(time(tot))
  allocate(nrm(tot))
  allocate(rlacf(tot))
  allocate(imacf(tot))
  allocate(abacf(tot))
  allocate(rlex(tot))
  allocate(imex(tot))
  allocate(abex(tot))
  allocate(ehr(tot))
  allocate(pop1(tot))
  allocate(pop2(tot))
  allocate(popsum(tot))
  allocate(popdiff(tot))
  n = 1130
  m = 7150

  do i=1,folders
       
    if (valid(i)==0) then
      cycle
    end if     
    write (repstr,"(i0)") i
    filename = "normpop_"//trim(repstr)//".out"
    n = n+1
    OPEN(UNIT=n, FILE=trim(filename),STATUS='OLD', iostat=ierr)

    if (ierr .ne. 0) then
      write (0,"(a,a,a)") 'Error in opening ', trim(filename), ' file'
      stop
    end if

    do
      read (n,*,iostat=ierr)LINE
      if (ierr.ne.0) then
        write (0,"(a,a,a)") "Read Error in normpop_", trim(repstr), ".out"
        stop
      end if
      if (LINE=="0.00000000E+000") then
        backspace (n)
        exit
      else
        cycle
      end if
    end do

  end do

  do l=1,tot

    m=7150+l

    write(repstr,"(i0)") folreps*l

    filename= "normpop_cumul_"//trim(repstr)//".out"

    open (unit=m, file=trim(filename), status='unknown', iostat=ierr)

    if (ierr .ne. 0) then
      write (0,"(a,a,a)") 'Error in opening ', trim(filename), ' output file'
      stop
    end if

    write (m,*) "Time Norm Re(ACF(t)) Im(ACF(t) |ACF(t)| Re(Extra) Im(Extra)", &
                           " |Extra| Sum(HEhr) Pop1 Pop2 Pop1+Pop2 Pop2-Pop1"
    write (m,*) ""
    write (m,*) ""

  end do

  do k = 1,lines(1)

    do i=1,tot
      n=1130+i
      read(n,*,iostat=ierr)time(i),nrm(i),rlacf(i),imacf(i),abacf(i),rlex(i),&
                        imex(i),abex(i),ehr(i),pop1(i),pop2(i),popsum(i),popdiff(i)
      if (time(i).ne.time(1)) then
        write (0,"(a,e15.8,a,i0)") "File synchronisation error when reading at ", & 
                                    "time ", time(1) ," in unit ", n
        write (0,"(a,i0,a,i0,a,i0)") "i = ", i, " n = ", n
        write (0,"(a,e15.8,a,e15.8)") "Expected t=", time(1), " but got t=", time(i)
        stop
      end if
    end do

    do i=1,tot

      m=7150+i
       
      timeav = sum(time(1:i))/i
      nrmav = sum(nrm(1:i))/i
      rlacfav = sum(rlacf(1:i))/i
      imacfav = sum(imacf(1:i))/i
      abacfav = sum(abacf(1:i))/i
      rlexav = sum(rlex(1:i))/i
      imexav = sum(imex(1:i))/i
      abexav = sum(abex(1:i))/i
      ehrav = sum(ehr(1:i))/i
      pop1av = sum(pop1(1:i))/i
      pop2av = sum(pop2(1:i))/i
      popsumav = sum(popsum(1:i))/i
      popdiffav = sum(popdiff(1:i))/i

      write (m,"(13(1x,e15.8))") timeav,nrmav,rlacfav,imacfav,abacfav,rlexav,&
                        imexav,abexav,ehrav,pop1av,pop2av,popsumav,popdiffav

    end do

  end do

  do l=1,tot

    write(repstr,"(i0)") folreps*l
    
    if (l==1) then     
      lngchar = 'plot "normpop_cumul_'//trim(repstr)//'.out" u 1:13 t "' & 
                            //trim(repstr)//' Reps" w l'
    else
      lngchar2 = lngchar
      lngchar = trim(lngchar2)//', "normpop_cumul_'//trim(repstr)// & 
                       '.out" u 1:13 t "'//trim(repstr)//' Reps" w l'
    end if

  end do

  lngchar2 = lngchar
  lngchar = trim(lngchar2)//', "normpop_cumul_'//trim(repstr)// & 
                                      '.out" u 1:2 t "Total Av Norm" w l'

  open (unit=175,file="plotpopdiff.gpl",status="unknown",iostat=ierr)
  if (ierr .ne. 0) then
    write (0,"(a)") 'Error in opening plotpopdiff.gpl output file'
    stop
  end if

  write(175,"(a)") 'set terminal png'
  write(175,"(a)") 'set output "popsculm.png"'
  write(175,"(a)") 'set title "Graph of convergence over multiple repetitions"'
  write(175,"(a)") 'set xlabel "Time"'
  write(175,"(a)") 'set ylabel "Population difference"'
  write(175,"(a)") trim(lngchar)
  close (175)

  if (tot .gt. 1) then

    open(unit=180,file="popdiffresiduals.out",status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(a)") 'Error in opening popdiffresiduals.out output file'
      stop
    end if

    do i=1,tot
      n=7150+i
      rewind(n)
      do
        read(n,*,iostat=ierr) LINE
        if (LINE=="0.00000000E+00") then
          backspace(n)
          exit
        else
          cycle
        end if
      end do          
    end do

    do k = 1,lines(1)

      do i=1,tot
        n=7150+i
        read(n,*,iostat=ierr)time(i),nrm(i),rlacf(i),imacf(i),abacf(i),rlex(i),&
                      imex(i),abex(i),ehr(i),pop1(i),pop2(i),popsum(i),popdiff(i)
        if (time(i).ne.time(1)) then
          write (0,"(a,a,e15.8,a,i0)"),"File synchronisation error when reading at",&
                               "time ", time(1) ," in unit ", n
          write (0,"(a,i0,a,i0,a,i0)") "i = ", i, " n = ", n
          stop
        end if
      end do

      do i=1,tot
        popdiff(i)=popdiff(tot)-popdiff(i)
      end do

      write (myfmt,'(a,i0,a)') '(', tot+1, '(1x,e12.5))'
       
      write (180,trim(myfmt)) time(1), popdiff

    end do

    close (180)

    do l=1,tot

      write(repstr,"(i0)") folreps*l
      write(lstr,"(i0)") l+1
    
      if (l==1) then     
        lngchar = 'plot "popdiffresiduals.out" u 1:2 t "'//trim(repstr)//&
                     ' Reps" w l'
      else
        lngchar2 = lngchar
        lngchar = trim(lngchar2)//', "" u 1:'//trim(lstr)//' t "'//trim(repstr)//&
                     ' Reps" w l'
      end if

    end do

    open (unit=176,file="plotpopres.gpl",status="unknown",iostat=ierr)
    if (ierr .ne. 0) then
      write (0,"(a)") 'Error in opening plotpopres.gpl output file'
      stop
    end if

    write(176,"(a)") 'set terminal png'
    write(176,"(a)") 'set output "popsculmres.png"'
    write(176,"(a)") 'set title "Graph of residuals over multiple repetitions"'
    write(176,"(a)") 'set xlabel "Time"'
    write(176,"(a)") 'set ylabel "Residual Population Difference"'
    write(176,"(a)") trim(lngchar)
    close (176)
  end if

  stop

end program avrgnorm
