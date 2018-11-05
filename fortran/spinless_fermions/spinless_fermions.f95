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

    integer :: i
    integer :: junk

    character(len = 100) :: fname

    ! integer, DIMENSION(3, 3) :: array=reshape( (/ 1, 0, 0, &
    !                                           0, 2, 0, &
    !                                           3, 0, 3 /), &
    !                                        shape(array), order=(/2,1/) )
    integer, allocatable :: tau(:, :)


    
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

    call setstate(s1, 4, 2, 0, 1)
    call setstate(s2, 4, 2, 1, 1)
    allocate(tau(getstatesize(s1), getstatesize(s2)))
    tau = connecting_tau(s1, s2)
    call matprint(tau)




end program spinless_fermions