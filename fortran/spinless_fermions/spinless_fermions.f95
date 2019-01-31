program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites = 6
    integer :: num_particles = 3
    integer :: level

    integer :: smax
    integer :: ns

    integer :: nlines
    integer, allocatable :: sizedata(:)
    integer, allocatable :: ordinatedata(:)

    type (state) :: s1
    type (state) :: s2
    type (state) :: s3
    type (state) :: s4

    type (state) :: t_s1
    type (state) :: t_s2
    type (state) :: tstate

    integer :: i, j
    integer :: junk
    integer :: ts, tns, ta, tb, tpos
    integer :: ts1, tns1, ta1, tb1, tpos1
    integer :: ts2, tns2, ta2, tb2, tpos2
    integer :: currenthead = 1

    double complex, allocatable :: g11(:,:)
    double complex, allocatable :: g22(:,:)
    double complex, allocatable :: u12(:,:)
    double complex, allocatable :: u21(:,:)
    double complex, allocatable :: testmatrix(:,:)
    double complex, allocatable :: testmatrix1(:,:)
    double complex, allocatable :: testmatrix2(:,:)
    double complex, allocatable :: fg11(:,:)
    double complex, allocatable :: fg12(:,:)
    double complex, allocatable :: fg21(:,:)
    double complex, allocatable :: fg22(:,:)

    double complex :: bineet_test1
    double complex :: bineet_test2
    real(8) :: bineet_real_test1
    real(8) :: modsqz

    character(len = 100) :: fname

    real(8) :: w
    real(8) :: dos

    real(8), parameter :: pi = 3.1415926535

    integer, parameter :: checklevel = 3


    smax = num_sites - 1

    ns = num_particles

    ! write (*,*) inv(array)

    call setstate(s1, smax, ns, 0, 0)
    call setstate(s2, smax, ns, 1, 0)
    call setstate(s3, smax, ns-1, 0, 1)
    call setstate(s4, smax, ns-1, 1, 1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    fname = 'splitstates.dat'
    open(10, file = trim(fname))

    call splitstates(10, s1)
    call splitstates(10, s2)
    call splitstates(10, s3)
    call splitstates(10, s4)

    close(10)

    do level = smax, 1, -1

        ! fname = 'splitstatesatlevel.dat'
        write(fname, '(a, i2.2, a)') "splitstatesatlevel", level, ".dat"
        open(10, file = trim(fname))

        call splitstatesatlevel(10, s1, level)
        call splitstatesatlevel(10, s2, level)
        call splitstatesatlevel(10, s3, level)
        call splitstatesatlevel(10, s4, level)

        close(10)

        write(fname, '(a, i2.2, a)') "splitstatesatlevelraw", level, ".dat"
        open(10, file = trim(fname))

        call splitstatesatlevelraw(10, s1, level)
        call splitstatesatlevelraw(10, s2, level)
        call splitstatesatlevelraw(10, s3, level)
        call splitstatesatlevelraw(10, s4, level)

        close(10)

        nlines = 0

        open(11, file = fname)
        do
            read(11, *, end = 110)
            nlines = nlines + 1
        end do
        110 close(11)

        open(11, file = fname)
        write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
        open(10, file = fname)

        ! print *, nlines
        
        ! allocate(rawdata(nlines, 5))
        allocate(sizedata(nlines))

        
        ! rawdata = transpose(rawdata)
        
        do i = 1, nlines
            read(11, *) ts, tns, ta, tb, sizedata(i)
            if ( i .eq. 1) then
                write (10,*) ts, tns, ta, tb, 1
            else
                write (10,*) ts, tns, ta, tb, sum(sizedata(1:i-1)) + 1
            end if
        end do
        close(11)
        close(10)
        deallocate(sizedata)

    end do

    w = 0.0
    do while (w .eq. 0.0) !FOR TESTING; comment out the two lines below

    ! w = -5
    ! do while (w .le. 5)

    level = 3
    nlines = 0
    currenthead = 1
    write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
    open(11, file = fname)
    do
        read(11, *, end = 111)
        nlines = nlines + 1
    end do
    111 close(11)

    open(11, file = fname)
    i = 1
    do while (i .le. nlines)
        read(11,*) ts, tns, ta, tb, tpos
        write(fname, '(a, i2.2, a, i2.2, a)') 'g_', ts, '_',  tpos, ".dat"
        call setstate(tstate, ts, tns, ta, tb)
        open(10, file = fname)
        call matprinttofile(10, g(diagonalblock(tstate), w) )
        close(10)
        i = i + 1
    end do
    close(11)

    level = 3
    do while (level .le. smax)
        write(*,*) level
        nlines = 0
        currenthead = 1
        write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
        open(11, file = fname)
        do
            read(11, *, end = 112)
            nlines = nlines + 1
        end do
        112 close(11)
        
        do while(currenthead .lt. nlines .and. level .ge. 2)
            write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
            open(11, file = fname)
            i = 1
            do while (i .lt. currenthead)
                read(11,*) junk, junk, junk, junk, junk
                i = i + 1
            end do
            
            read(11,*) ts1, tns1, ta1, tb1, tpos1
            call setstate(t_s1, ts1, tns1, ta1, tb1)
            read(11,*) ts2, tns2, ta2, tb2, tpos2
            call setstate(t_s2, ts2, tns2, ta2, tb2)
            close(11)

            if (getletter(t_s1) .eq. 1 .and. getletter(t_s2) .eq. 3) then
                ! write(*,*) getlstate(t_s1), tpos1, getlstate(t_s2), tpos2

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1, '_',  tpos1, ".dat"
                call setfiletomatrix(fname, g11)
                ! write(*,*) fname
                ! call matprint(g11)

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2, '_',  tpos2, ".dat"
                call setfiletomatrix(fname, g22)
                ! write(*,*) fname
                ! call matprint(g22)

                allocate( u12( size(g11, 1), size(g22, 2) ) )
                allocate( u21( size(g22, 1), size(g11, 2) ) )
                
                allocate( fg11( size(g11, 1), size(g11, 2)) )
                allocate( fg22( size(g22, 1), size(g22, 2)) )
                allocate( fg12( size(u12, 1), size(u12, 2)) )
                allocate( fg21( size(u21, 1), size(u21, 2)) )
                
                u12 = connecting_u(t_s1, t_s2)
                u21 = transpose(u12)

                ! call matprint(u12)
                ! call matprint(u21)

                fg11 = inv( inv(g11) - tmm( u12, g22, u21 ) )
                fg22 = inv( inv(g22) - tmm( u21, g11, u12 ) )
                fg12 = tmm( fg11, u12, g22 )
                fg21 = tmm( fg22, u21, g11 )

                ! call matprint(fg11)

                if (level .eq. checklevel) then
                    write(*,*) 'g11'
                    call matprint(g11)
                    write(*,*) 'g11 inverse'
                    call matprint(inv(g11))
                    write(*,*) 'g22'
                    call matprint(g22)
                    write(*,*) 'g22 inverse'
                    call matprint(inv(g22))
                    write(*,*) 'u12'
                    call matprint(u12)
                    write(*,*) 'u21'
                    call matprint(u21)
                    write(*,*) 'g22*u21'
                    call matprint(matprod(g22,u21))
                    write(*,*) 'u12 * g22 * u21'
                    call matprint(tmm( u12, g22, u21 ))
                    write(*,*) 'inv(g11) - u12 * g22 * u21'
                    call matprint(inv(g11) - tmm( u12, g22, u21 ))
                    write(*,*) 'full connected fg11'
                    call matprint(fg11)
                    write(*,*) '-------------------'
                end if

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1+1, '_',  tpos1, ".dat"
                open(20, file = fname)
                call matprinttofile(20, bmat(fg11, fg12, fg21, fg22))
                close(20)

                deallocate(g11)
                deallocate(g22)
                deallocate(u12)
                deallocate(u21)
                deallocate(fg11)
                deallocate(fg22)
                deallocate(fg12)
                deallocate(fg21)

                currenthead = currenthead + 2
            else if (getletter(t_s1) .eq. 2 .and. getletter(t_s2) .eq. 4) then
                ! write(*,*) getlstate(t_s1),tpos1, getlstate(t_s2), tpos2

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1, '_',  tpos1, ".dat"
                call setfiletomatrix(fname, g11)
                ! write(*,*) fname
                ! call matprint(g11)

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2, '_',  tpos2, ".dat"
                call setfiletomatrix(fname, g22)
                ! write(*,*) fname
                ! call matprint(g22)

                allocate( u12( size(g11, 1), size(g22, 2) ) )
                allocate( u21( size(g22, 1), size(g11, 2) ) )
                
                allocate( fg11( size(g11, 1), size(g11, 2)) )
                allocate( fg22( size(g22, 1), size(g22, 2)) )
                allocate( fg12( size(u12, 1), size(u12, 2)) )
                allocate( fg21( size(u21, 1), size(u21, 2)) )
                
                u12 = connecting_u(t_s1, t_s2)
                u21 = transpose(u12)

                ! call matprint(u12)
                ! call matprint(u21)

                fg11 = inv( inv(g11) - tmm( u12, g22, u21 ) )
                fg22 = inv( inv(g22) - tmm( u21, g11, u12 ) )
                fg12 = tmm( fg11, u12, g22 )
                fg21 = tmm( fg22, u21, g11 )

                if (level .eq. checklevel) then
                    write(*,*) 'g11'
                    call matprint(g11)
                    write(*,*) 'g11 inverse'
                    call matprint(inv(g11))
                    write(*,*) 'g22'
                    call matprint(g22)
                    write(*,*) 'g22 inverse'
                    call matprint(inv(g22))
                    write(*,*) 'u12'
                    call matprint(u12)
                    write(*,*) 'u21'
                    call matprint(u21)
                    write(*,*) 'g22*u21'
                    call matprint(matprod(g22,u21))
                    write(*,*) 'u12 * g22 * u21'
                    call matprint(tmm( u12, g22, u21 ))
                    write(*,*) 'inv(g11) - u12 * g22 * u21'
                    call matprint(inv(g11) - tmm( u12, g22, u21 ))
                    write(*,*) 'full connected fg11'
                    call matprint(fg11)
                    write(*,*) '-------------------'
                end if

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1+1, '_',  tpos1, ".dat"
                open(20, file = fname)
                call matprinttofile(20, bmat(fg11, fg12, fg21, fg22))
                close(20)

                deallocate(g11)
                deallocate(g22)
                deallocate(u12)
                deallocate(u21)
                deallocate(fg11)
                deallocate(fg22)
                deallocate(fg12)
                deallocate(fg21)

                currenthead = currenthead + 2
            else
                ! write(*,*) getlstate(t_s1), tpos1

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1, '_',  tpos1, ".dat"
                call setfiletomatrix(fname, g11)

                if (level .eq. checklevel) then
                    write(*,*) 'g11'
                    call matprint(g11)
                    write(*,*) 'g11 inverse'
                    call matprint(inv(g11))
                    write(*,*) '-------------------'
                end if
                
                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1+1, '_',  tpos1, ".dat"
                open(15, file = fname)
                allocate(fg11(size(g11, 1), size(g11, 2)))
                fg11 = g11
                call matprinttofile(15, fg11)
                close(15)

                deallocate(g11)
                deallocate(fg11)


                currenthead = currenthead + 1
            end if

            if (currenthead .eq. nlines) then

                write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
                open(11, file = fname)
                i = 1
                do while (i .ne. currenthead)
                    read(11,*) junk, junk, junk, junk, junk
                    i = i + 1
                end do
            
                read(11,*) ts2, tns2, ta2, tb2, tpos2
                call setstate(t_s2, ts2, tns2, ta2, tb2)
                close(11)
                
                ! write(*,*) getlstate(t_s2), tpos2
                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2, '_',  tpos2, ".dat"
                call setfiletomatrix(fname, g22)
                
                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2+1, '_',  tpos2, ".dat"
                open(15, file = fname)
                allocate(fg22(size(g22, 1), size(g22, 2)))
                fg22 = g22
                call matprinttofile(15, fg22)

                if (level .eq. checklevel) then
                    write(*,*) 'g22'
                    call matprint(g22)
                    write(*,*) 'g22 inverse'
                    call matprint(inv(g22))
                    write(*,*) '-------------------'
                end if

                close(15)

                deallocate(g22)
                deallocate(fg22)
            end if

        end do
        level = level + 1
    end do

    call setstate(s1, smax+1, ns, 0, 0)
    call setstate(s2, smax+1, ns, 1, 0)
    
    ! write(*,*) getlstate(s1), getlstate(s2)
    write(*,*) smax+1
    tpos1 = 1
    tpos2 = 1 + getstatesize(s1)
    write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', smax+1, '_',  tpos1, ".dat"
                call setfiletomatrix(fname, g11)
    ! write(*,*) fname
    write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', smax+1, '_',  tpos2, ".dat"
                call setfiletomatrix(fname, g22)
    ! write (*,*) fname

    allocate( u12( size(g11, 1), size(g22, 2) ) )
    allocate( u21( size(g22, 1), size(g11, 2) ) )
    
    allocate( fg11( size(g11, 1), size(g11, 2)) )
    allocate( fg22( size(g22, 1), size(g22, 2)) )
    allocate( fg12( size(u12, 1), size(u12, 2)) )
    allocate( fg21( size(u21, 1), size(u21, 2)) )
    
    u12 = connecting_u(s1, s2)
    u21 = transpose(u12)

    ! write(*,*) 'The 2 input matrices at smax+1 level'
    ! call diagonal_imaginary_part(g11)
    ! write(*,*) '*****************'
    ! call diagonal_imaginary_part(g22)

    ! call matprint(u12)
    ! call matprint(u21)

    fg11 = inv( inv(g11) - tmm( u12, g22, u21 ) )
    fg22 = inv( inv(g22) - tmm( u21, g11, u12 ) )
    
    ! write(*,*) 'After final calculation'
    ! call diagonal_imaginary_part(fg11)
    ! write(*,*) '*****************'
    ! call diagonal_imaginary_part(fg22)

    dos = - imag(trace(fg11) + trace(fg22))/ &
    (getstatesize(s1) + getstatesize(s2)) / (pi)

    open(30, file = 'dosdata.dat')
    close(30)
    open(30, file = 'dosdata.dat', status = 'old', position = 'append')
    write(30,*) w, dos
    close(30)


    deallocate(g11)
    deallocate(g22)
    deallocate(u12)
    deallocate(u21)
    deallocate(fg11)
    deallocate(fg22)
    deallocate(fg12)
    deallocate(fg21)


    write(*,*) w, dos
    w = w + 10.0/1000.0

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TESTING AREA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! allocate(testmatrix1(3, 3))
    ! allocate(testmatrix2(3, 3))

    ! ! testmatrix1 = 0.0
    ! ! testmatrix2 = 0.0

    ! do i = 1, 3
    !     do j = 1, 3
    !         testmatrix1(i,j) = complex(i, j)
    !     end do
    ! end do

    ! ! call matprint(testmatrix1)

    ! do i = 1, 3
    !     do j = 1, 3
    !         testmatrix2(i,j) = complex(j, i*i)
    !     end do
    ! end do

    ! ! call matprint(testmatrix2)


    ! call matprint( matmul(testmatrix1, testmatrix2) )

    ! write(*,*) '4,1'
    ! call setfiletomatrix('g_04_01.dat', testmatrix)
    ! call diagonal_only_part(testmatrix)
    ! deallocate(testmatrix)

    ! write(*,*) '4,2'
    ! call setfiletomatrix('g_04_02.dat', testmatrix)
    ! call diagonal_only_part(testmatrix)
    ! write(*,*) '------------------'
    ! call matprint(testmatrix)
    ! deallocate(testmatrix)

    ! write(*,*) '4,3'
    ! call setfiletomatrix('g_04_05.dat', testmatrix)
    ! call diagonal_only_part(testmatrix)
    ! deallocate(testmatrix)

    ! write(*,*) '4,4'
    ! call setfiletomatrix('g_04_08.dat', testmatrix)
    ! call diagonal_only_part(testmatrix)
    ! deallocate(testmatrix)
    
    !write(*,*) real((3.0, 4.0) * conjg((3.0, 4.0) ))
    
    ! allocate (testmatrix1 (2,2))
    ! testmatrix1(1,1) = complex( 0.0, -0.09)
    ! testmatrix1(2,2) = complex( 0.0, -0.09)
    ! testmatrix1(1,2) = complex( -1.0, 0.0)
    ! testmatrix1(2,1) = complex( -1.0, 0.0)

    ! allocate (testmatrix2 (2,1))
    ! testmatrix2(1,1) = complex(1.0,0.0)
    ! testmatrix2(2,1) = complex(0.0,0.0)

    ! call matprint(testmatrix1)
    ! write(*,*)
    ! call matprint(testmatrix2)

    ! write(*,*) ' test prod '

    ! call matprint(matprod(testmatrix1, testmatrix2)) 
    
    ! bineet_test1 = complex(0.0,-0.09)
    ! bineet_test2 = complex(1.0,0.0)

    ! ! write(*,*) bineet_test1*bineet_test2

    ! write(*,*) bineet_test1
    ! write(*,*) cmplxfilterbelow(bineet_test1)
    
    ! modsqz = real(bineet_test1 * conjg(bineet_test1))
    ! write(*,*) modsqz

    ! if ( modsqz .le. 1.0E-10) then
    !     write(*,*) 'below'
    !     ! cmplxfilterbelow = complex(0.0, 0.0)
    ! else
    !     write(*,*) 'above'
    ! end if

end program spinless_fermions