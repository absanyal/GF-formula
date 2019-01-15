program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites = 8
    integer :: num_particles = 4
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

    complex, allocatable :: g11(:,:)
    complex, allocatable :: g22(:,:)
    complex, allocatable :: tau12(:,:)
    complex, allocatable :: tau21(:,:)
    complex, allocatable :: blockh(:,:)
    complex, allocatable :: blockg(:,:)
    complex, allocatable :: testmatrix(:,:)
    ! complex, allocatable :: tau(:, :)
    complex, allocatable :: fg11(:,:)
    complex, allocatable :: fg12(:,:)
    complex, allocatable :: fg21(:,:)
    complex, allocatable :: fg22(:,:)

    character(len = 100) :: fname

    real :: w
    ! integer :: level

    w = 0

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

    w = 0

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
        ! write(*,*) level
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
                write(*,*) getlstate(t_s1), tpos1, getlstate(t_s2), tpos2

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

                deallocate(g11)
                deallocate(g22)

                currenthead = currenthead + 2
            else if (getletter(t_s1) .eq. 2 .and. getletter(t_s2) .eq. 4) then
                write(*,*) getlstate(t_s1),tpos1, getlstate(t_s2), tpos2

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

                deallocate(g11)
                deallocate(g22)

                currenthead = currenthead + 2
            else
                write(*,*) getlstate(t_s1), tpos1

                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts1, '_',  tpos1, ".dat"
                call setfiletomatrix(fname, g11)
                
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
            
                ! read(11,*) ts1, tns1, ta1, tb1, tpos1
                ! call setstate(t_s1, ts1, tns1, ta1, tb1)
                read(11,*) ts2, tns2, ta2, tb2, tpos2
                call setstate(t_s2, ts2, tns2, ta2, tb2)
                close(11)
                
                write(*,*) getlstate(t_s2), tpos2
                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2, '_',  tpos2, ".dat"
                call setfiletomatrix(fname, g22)
                
                write(fname, '(a, i2.2, a, i2.2, a)') &
                'g_', ts2+1, '_',  tpos2, ".dat"
                open(15, file = fname)
                allocate(fg22(size(g22, 1), size(g22, 2)))
                fg22 = g22
                call matprinttofile(15, fg22)
                close(15)

                deallocate(g22)
                deallocate(fg22)
            end if

        end do
        level = level + 1
    end do

end program spinless_fermions