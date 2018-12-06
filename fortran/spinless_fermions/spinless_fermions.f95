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

    complex, allocatable :: g1(:,:)
    complex, allocatable :: g2(:,:)
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

    

    level = 3
    nlines = 0
    write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
    open(11, file = fname)
    do
        read(11, *, end = 111)
        nlines = nlines + 1
    end do
    111 close(11)

    ! print *, nlines
    
    do while(currenthead .lt. nlines)
    write(fname, '(a, i2.2, a)') "stateordinatesatlevel", level, ".dat"
    open(11, file = fname)
    i = 1
    do while (i .lt. currenthead)
        ! print *, 'reading junk'
        read(11,*) junk, junk, junk, junk, junk
        i = i + 1
    end do
    
    read(11,*) ts1, tns1, ta1, tb1, tpos1
    call setstate(t_s1, ts1, tns1, ta1, tb1)
    ! write(*,*) getlstate(t_s1)
    read(11,*) ts2, tns2, ta2, tb2, tpos2
    call setstate(t_s2, ts2, tns2, ta2, tb2)

    ! write(*,*) getlstate(t_s1), getlstate(t_s2)

    if (getletter(t_s1) .eq. 1 .and. getletter(t_s2) .eq. 3) then
        write(*,*) getlstate(t_s1), getlstate(t_s2)
        currenthead = currenthead + 2
    else if (getletter(t_s1) .eq. 2 .and. getletter(t_s2) .eq. 4) then
        write(*,*) getlstate(t_s1), getlstate(t_s2)
        currenthead = currenthead + 2
    else
        write(*,*) getlstate(t_s1)
        currenthead = currenthead + 1
    end if
    
    ! currenthead = currenthead + 1
    ! write(*,*) currenthead

    if (currenthead .eq. nlines) then
        ! read(11,*) ts2, tns2, ta2, tb2, tpos2
        ! call setstate(t_s2, ts2, tns2, ta2, tb2)
        write(*,*) getlstate(t_s2)
    end if

    close(11)
    end do

end program spinless_fermions