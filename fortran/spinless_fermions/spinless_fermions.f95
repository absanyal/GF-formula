program spinless_fermions
    
    use statemanip
    use gf
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

    integer :: i
    integer :: junk
    integer :: ts, tns, ta, tb, tsize

    ! complex, allocatable :: h1(:,:)
    ! complex, allocatable :: h2(:,:)
    ! complex, allocatable :: htau12(:,:)
    ! complex, allocatable :: htau21(:,:)
    complex, allocatable :: blockh(:,:)

    character(len = 100) :: fname

    complex, allocatable :: testmatrix(:,:)

    ! integer, DIMENSION(3, 3) :: array=reshape( (/ 1, 0, 0, &
    !                                           0, 2, 0, &
    !                                           3, 0, 3 /), &
    !                                        shape(array), order=(/2,1/) )
    complex, allocatable :: tau(:, :)


    
    ! allocate(eye3(3, 4))
    ! eye3 = int_zeros(3, 4)
    ! call matprint(eye3)

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

    ! fname = 'splitstatesraw.dat'
    ! open(10, file = trim(fname))

    ! call splitstatesraw(10, s1)
    ! call splitstatesraw(10, s2)
    ! call splitstatesraw(10, s3)
    ! call splitstatesraw(10, s4)

    ! close(10)

    ! fname = 'splitstatesatlevelraw.dat'
    write(fname, '(a, i2.2, a)') "splitstatesatlevelraw", level, ".dat"
    open(10, file = trim(fname))

    call splitstatesatlevelraw(10, s1, level)
    call splitstatesatlevelraw(10, s2, level)
    call splitstatesatlevelraw(10, s3, level)
    call splitstatesatlevelraw(10, s4, level)

    close(10)

    ! fname = 'findmels.dat'
    ! open(10, file = trim(fname))

    ! call findmels(10, s1, s1, 1)
    ! call findmels(10, s1, s2, 0)
    ! call findmels(10, s2, s2, 1)
    ! call findmels(10, s2, s3, 0)
    ! call findmels(10, s3, s3, 1)
    ! call findmels(10, s3, s4, 0)
    ! call findmels(10, s4, s4, 1)

    ! close(10)

    ! fname = 'findmelsatlevel.dat'
    ! open(10, file = trim(fname))

    ! call findmelsatlevel(10, s1, s1, level)
    ! call findmelsatlevel(10, s1, s2, level)
    ! call findmelsatlevel(10, s2, s2, level)
    ! call findmelsatlevel(10, s2, s3, level)
    ! call findmelsatlevel(10, s3, s3, level)
    ! call findmelsatlevel(10, s3, s4, level)
    ! call findmelsatlevel(10, s4, s4, level)

    ! close(10)

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
        read(11, *) junk, junk, junk, junk, sizedata(i)
        if ( i .eq. 1) then
            write (10,*) sizedata(i), 1
        else
            write (10,*) sizedata(i), sum(sizedata(1:i-1)) + 1
        end if
    end do
    close(11)
    close(10)
    deallocate(sizedata)

    end do

    ! call setstate(s1, 3, 2, 0, 1)
    ! call setstate(s2, 3, 2, 1, 1)
    ! allocate(tau(getstatesize(s1), getstatesize(s2)))
    ! tau = transpose(connecting_tau(s1, s2))
    ! call matprint(tau)

    nlines = 0
    fname = 'splitstatesatlevelraw03.dat'
    open(11, file = fname)
    do
        read(11, *, end = 111)
        nlines = nlines + 1
    end do
    111 close(11)
    ! allocate(sizedata(nlines))

    open(11, file = fname)
    read (11,*) ts, tns, ta, tb, tsize
    call setstate(tstate, ts, tns, ta, tb)
    rewind(11)
    if (getletter(tstate) .eq. 1) then
        i = 1
        do while (i .le. nlines)
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s1, ts, tns, ta, tb)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            i = i + 1
            write(*, *) getlstate(t_s1), getlstate(t_s2)
            allocate(blockh(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            blockh = connectedblock(t_s1, t_s2)
            call matprint(blockh)
            deallocate(blockh)
        end do
    else
        read (11,*) ts, tns, ta, tb, tsize
        call setstate(t_s1, ts, tns, ta, tb)
        write (*,*) getlstate(t_s1), getstatesize(t_s1)
        allocate(blockh(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        blockh = diagonalblock(t_s1)
        call matprint(blockh)
        deallocate(blockh)
        i = 2
        do while (i .le. nlines - 1)
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s1, ts, tns, ta, tb)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            i = i + 1
            write(*, *) getlstate(t_s1), getstatesize(t_s1), &
                        getlstate(t_s2), getstatesize(t_s2)
            allocate(blockh(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            ! write(*,*) size(blockh, 1), size(blockh, 2)
            blockh = connectedblock(t_s1, t_s2)
            call matprint(blockh)
            deallocate(blockh)
        end do
        read (11,*) ts, tns, ta, tb, tsize
        call setstate(t_s1, ts, tns, ta, tb)
        write (*,*) getlstate(t_s1), getstatesize(t_s1)
        allocate(blockh(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        blockh = diagonalblock(t_s1)
        call matprint(blockh)
        deallocate(blockh)
    end if
    close(11)

    ! write(*,*) "*********TESTING CODE *****************"

    ! call setstate(t_s1, 3, 2, 0, 1)
    ! call setstate(t_s2, 3, 2, 1, 1)
    ! write(*,*) getlstate(t_s1), getlstate(t_s2)
    ! allocate(testmatrix(getstatesize(t_s1), getstatesize(t_s2)))
    ! testmatrix = connecting_tau(t_s2, t_s1)
    ! ! testmatrix = bmat(  twobytwo(), 2*twobytwo(), &
    ! !                     twobytwo(), 3*twobytwo())
    ! call matprint(testmatrix)




end program spinless_fermions