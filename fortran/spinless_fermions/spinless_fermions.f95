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
    integer :: ts, tns, ta, tb, tsize

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
    level = 3
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
        j = 1
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
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            blockh = connectedblock(t_s1, t_s2)
            blockg = g(blockh, w)
            ! call matprint(blockg)
            write(fname, '(a, i2.2,a,i2.2  a)') "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            ! write (*,*) "wrote", fname
            deallocate(blockh)
            deallocate(blockg)
        end do
    else
        j = 1
        read (11,*) ts, tns, ta, tb, tsize
        call setstate(t_s1, ts, tns, ta, tb)
        write (*,*) getlstate(t_s1), getstatesize(t_s1)
        allocate(blockh(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        allocate(blockg(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        blockh = diagonalblock(t_s1)
        blockg = g(blockh, w)
        ! call matprint(blockg)
        write(fname, '(a, i2.2,a,i2.2  a)') "g_level_", level,"_", j, ".dat"
        open(10, file = fname)
        call matprinttofile(10, blockg)
        j = j + 1
        close(10)
        ! write (*,*) "wrote", fname
        deallocate(blockh)
        deallocate(blockg)
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
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            blockh = connectedblock(t_s1, t_s2)
            blockg = g(blockh, w)
            ! call matprint(blockg)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            ! write (*,*) "wrote", fname
            deallocate(blockh)
            deallocate(blockg)
        end do
        read (11,*) ts, tns, ta, tb, tsize
        call setstate(t_s1, ts, tns, ta, tb)
        write (*,*) getlstate(t_s1), getstatesize(t_s1)
        allocate(blockh(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        allocate(blockg(getstatesize(t_s1), &
                            getstatesize(t_s1)))
        blockh = diagonalblock(t_s1)
        blockg = g(blockh, w)
        ! call matprint(blockg)
        write(fname, '(a, i2.2,a,i2.2  a)') "g_level_", level,"_", j, ".dat"
        open(10, file = fname)
        call matprinttofile(10, blockg)
        j = j + 1
        close(10)
        ! write (*,*) "wrote", fname
        deallocate(blockh)
        deallocate(blockg)
    end if
    close(11)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    level = 4

    do while (level .le. smax)

    nlines = 0
    write(fname, '(a, i2.2, a)') "splitstatesatlevelraw", level, ".dat"
    open(11, file = trim(fname))
    do
        read(11, *, end = 112)
        nlines = nlines + 1
    end do
    112 close(11)
    ! allocate(sizedata(nlines))

    open(11, file = fname)
    read (11,*) ts, tns, ta, tb, tsize
    call setstate(tstate, ts, tns, ta, tb)
    rewind(11)
    j = 1
    if (getletter(tstate) .eq. 1) then
        i = 1
        do while (i .le. nlines)
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s1, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g1(getstatesize(t_s1), getstatesize(t_s1)))
            read(10, *) g1
            close(10)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g2(getstatesize(t_s2), getstatesize(t_s2)))
            read(10, *) g2
            close(10)
            i = i + 1
            write(*, *) getlstate(t_s1), getstatesize(t_s1), &
                         getlstate(t_s2), getstatesize(t_s2)
            ! call matprint(g1)
            ! write(*,*) "---------------------------------------"
            ! call matprint(g2)
            ! write(*,*) "//////////////////////////////////////"
            allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
            allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
            allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
            allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
            allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
            allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            tau12 = - connecting_tau(t_s1, t_s2)
            tau21 = transpose(tau12)
            fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
            fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
            fg12 = - matmul(g1, matmul(tau12, g2))
            fg21 = - matmul(g2, matmul(tau21, g1))
            blockg = bmat(fg11, fg12, fg21, fg22)
            write(fname, '(a, i2.2,a,i2.2  a)') &
                            "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            deallocate(g1)
            deallocate(g2)
            deallocate(tau12)
            deallocate(tau21)
            deallocate(fg11)
            deallocate(fg12)
            deallocate(fg21)
            deallocate(fg22)
            deallocate(blockg)
        end do
    else
        i = 1
        read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s1, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g1(getstatesize(t_s1), getstatesize(t_s1)))
            read(10, *) g1
            close(10)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g2(getstatesize(t_s2), getstatesize(t_s2)))
            read(10, *) g2
            close(10)
            i = i + 1
            write(*, *) getlstate(t_s1), getstatesize(t_s1), &
                         getlstate(t_s2), getstatesize(t_s2)
            ! call matprint(g1)
            ! write(*,*) "---------------------------------------"
            ! call matprint(g2)
            ! write(*,*) "//////////////////////////////////////"
            ! allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
            ! allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
            ! allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
            ! allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
            ! allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
            ! allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            ! tau12 = - connecting_tau(t_s1, t_s2)
            ! tau21 = transpose(tau12)
            ! fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
            ! fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
            ! fg12 = - matmul(g1, matmul(tau12, g2))
            ! fg21 = - matmul(g2, matmul(tau21, g1))
            blockg = g1
            write(fname, '(a, i2.2,a,i2.2  a)') &
                            "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            deallocate(g1)
            deallocate(g2)
            ! deallocate(tau12)
            ! deallocate(tau21)
            ! deallocate(fg11)
            ! deallocate(fg12)
            ! deallocate(fg21)
            ! deallocate(fg22)
            deallocate(blockg)
        ! i = 2
        do while (i .le. nlines - 1)
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s1, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g1(getstatesize(t_s1), getstatesize(t_s1)))
            read(10, *) g1
            close(10)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g2(getstatesize(t_s2), getstatesize(t_s2)))
            read(10, *) g2
            close(10)
            i = i + 1
            write(*, *) getlstate(t_s1), getstatesize(t_s1), &
                         getlstate(t_s2), getstatesize(t_s2)
            ! call matprint(g1)
            ! write(*,*) "---------------------------------------"
            ! call matprint(g2)
            ! write(*,*) "//////////////////////////////////////"
            allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
            allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
            allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
            allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
            allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
            allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            tau12 = - connecting_tau(t_s1, t_s2)
            tau21 = transpose(tau12)
            fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
            fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
            fg12 = - matmul(g1, matmul(tau12, g2))
            fg21 = - matmul(g2, matmul(tau21, g1))
            blockg = bmat(fg11, fg12, fg21, fg22)
            write(fname, '(a, i2.2,a,i2.2  a)') &
                            "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            deallocate(g1)
            deallocate(g2)
            deallocate(tau12)
            deallocate(tau21)
            deallocate(fg11)
            deallocate(fg12)
            deallocate(fg21)
            deallocate(fg22)
            deallocate(blockg)
        end do
        call setstate(t_s1, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g1(getstatesize(t_s1), getstatesize(t_s1)))
            read(10, *) g1
            close(10)
            i = i + 1
            read (11,*) ts, tns, ta, tb, tsize
            call setstate(t_s2, ts, tns, ta, tb)
            write(fname, '(a, i2.2,a,i2.2  a)') &
            "g_level_", level-1,"_", i, ".dat"
            open(10, file = fname)
            allocate(g2(getstatesize(t_s2), getstatesize(t_s2)))
            read(10, *) g2
            close(10)
            i = i + 1
            write(*, *) getlstate(t_s1), getlstate(t_s2)
            ! call matprint(g1)
            ! write(*,*) "---------------------------------------"
            ! call matprint(g2)
            ! write(*,*) "//////////////////////////////////////"
            ! allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
            ! allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
            ! allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
            ! allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
            ! allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
            ! allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
            allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
                            getstatesize(t_s1) + getstatesize(t_s2)))
            ! tau12 = - connecting_tau(t_s1, t_s2)
            ! tau21 = transpose(tau12)
            ! fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
            ! fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
            ! fg12 = - matmul(g1, matmul(tau12, g2))
            ! fg21 = - matmul(g2, matmul(tau21, g1))
            blockg = g1
            write(fname, '(a, i2.2,a,i2.2  a)') &
                            "g_level_", level,"_", j, ".dat"
            open(10, file = fname)
            call matprinttofile(10, blockg)
            j = j + 1
            close(10)
            deallocate(g1)
            deallocate(g2)
            ! deallocate(tau12)
            ! deallocate(tau21)
            ! deallocate(fg11)
            ! deallocate(fg12)
            ! deallocate(fg21)
            ! deallocate(fg22)
            deallocate(blockg)
    end if
    close(11)

    level = level + 1

    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! level = smax + 1
    ! i = 1
    ! j = 1

    ! call setstate(t_s1, smax, ns, 0, 0)
    ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! "g_level_", level-1,"_", i, ".dat"
    ! open(10, file = fname)
    ! allocate(g1(getstatesize(s1) + getstatesize(s2), &
    !             getstatesize(s1) + getstatesize(s2)))
    ! write(*,*) size(g1, 1), size(g1, 2)
    ! read(10, *) g1
    ! close(10)
    ! i = i + 1
    ! read (11,*) ts, tns, ta, tb, tsize
    ! call setstate(t_s2, smax, ns, 1, 0)
    ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! "g_level_", level-1,"_", i, ".dat"
    ! open(10, file = fname)
    ! allocate(g2(getstatesize(s3) + getstatesize(s4), &
    !             getstatesize(s3) + getstatesize(s4)))
    ! read(10, *) g2
    ! close(10)
    ! i = i + 1
    ! write(*, *) getlstate(t_s1), getstatesize(t_s1), &
    !                 getlstate(t_s2), getstatesize(t_s2)
    ! ! call matprint(g1)
    ! ! write(*,*) "---------------------------------------"
    ! ! call matprint(g2)
    ! ! write(*,*) "//////////////////////////////////////"
    ! allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
    ! allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
    ! allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
    ! allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
    ! allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
    ! allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
    ! allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
    !                 getstatesize(t_s1) + getstatesize(t_s2)))
    ! tau12 = - connecting_tau(t_s1, t_s2)
    ! tau21 = transpose(tau12)
    ! fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
    ! fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
    ! fg12 = - matmul(g1, matmul(tau12, g2))
    ! fg21 = - matmul(g2, matmul(tau21, g1))
    ! blockg = bmat(fg11, fg12, fg21, fg22)
    ! ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! !                 "g_level_", level,"_", j, ".dat"
    ! ! open(10, file = fname)
    ! ! call matprinttofile(10, blockg)
    ! ! j = j + 1
    ! ! close(10)
    ! call matprint(blockg)
    ! deallocate(g1)
    ! deallocate(g2)
    ! deallocate(tau12)
    ! deallocate(tau21)
    ! deallocate(fg11)
    ! deallocate(fg12)
    ! deallocate(fg21)
    ! deallocate(fg22)
    ! deallocate(blockg)

    ! write(*,*) "*********TESTING CODE *****************"

    ! call setstate(t_s1, 3, 2, 0, 1)
    ! call setstate(t_s2, 3, 2, 1, 1)
    ! write(*,*) getlstate(t_s1), getlstate(t_s2)
    ! allocate(testmatrix(3,3))
    ! write (*,*) size(testmatrix, 1), size(testmatrix, 2)
    ! testmatrix = g(connectedblock(t_s1, t_s2), 0.0)
    ! call matprint(testmatrix)

    ! i = 1
    ! call setstate(t_s1, smax, ns, 0, 0)
    ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! "g_level_", level-1,"_", i, ".dat"
    ! open(10, file = fname)
    ! allocate(g1(getstatesize(t_s1), getstatesize(t_s1)))
    ! read(10, *) g1
    ! close(10)
    ! i = i + 1
    ! ! read (11,*) ts, tns, ta, tb, tsize
    ! call setstate(t_s2, smax, ns, 1, 0)
    ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! "g_level_", level-1,"_", i, ".dat"
    ! open(10, file = fname)
    ! allocate(g2(getstatesize(t_s2), getstatesize(t_s2)))
    ! read(10, *) g2
    ! close(10)
    ! i = i + 1
    ! allocate(tau12(getstatesize(t_s1), getstatesize(t_s2)))
    ! allocate(tau21(getstatesize(t_s2), getstatesize(t_s1)))
    ! allocate(fg11(getstatesize(t_s1), getstatesize(t_s1)))
    ! allocate(fg12(getstatesize(t_s1), getstatesize(t_s2)))
    ! allocate(fg21(getstatesize(t_s2), getstatesize(t_s1)))
    ! allocate(fg22(getstatesize(t_s2), getstatesize(t_s2)))
    ! allocate(blockg(getstatesize(t_s1) + getstatesize(t_s2), &
    !                 getstatesize(t_s1) + getstatesize(t_s2)))
    ! tau12 = - connecting_tau(t_s1, t_s2)
    ! tau21 = transpose(tau12)
    ! fg11 = inv( inv(g1) - matmul(tau12, matmul(g2, tau21)) )
    ! fg22 = inv( inv(g2) - matmul(tau21, matmul(g1, tau12)) )
    ! fg12 = - matmul(g1, matmul(tau12, g2))
    ! fg21 = - matmul(g2, matmul(tau21, g1))
    ! blockg = bmat(fg11, fg12, fg21, fg22)
    ! ! write(fname, '(a, i2.2,a,i2.2  a)') &
    ! !                 "g_level_", level,"_", j, ".dat"
    ! ! open(10, file = fname)
    ! ! call matprinttofile(10, blockg)
    ! ! j = j + 1
    ! ! close(10)
    ! ! call matprint(blockg)
    ! write(*,*) trace(blockg)
    ! deallocate(g1)
    ! deallocate(g2)
    ! deallocate(tau12)
    ! deallocate(tau21)
    ! deallocate(fg11)
    ! deallocate(fg12)
    ! deallocate(fg21)
    ! deallocate(fg22)
    ! deallocate(blockg)



end program spinless_fermions