program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites = 6
    integer :: num_particles = 3
    integer :: level = 2

    integer :: smax
    integer :: ns

    integer :: nlines
    ! integer, allocatable :: rawdata(:, :)
    integer, allocatable :: sizedata(:)
    integer, allocatable :: ordinatedata(:)
    ! integer :: sizedata

    type (state) :: s1
    type (state) :: s2
    type (state) :: s3
    type (state) :: s4

    integer :: i
    integer :: junk

    character(len = 100) :: fname

    smax = num_sites - 1

    ns = num_particles

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

    fname = 'splitstatesatlevel.dat'
    open(10, file = trim(fname))

    call splitstatesatlevel(10, s1, level)
    call splitstatesatlevel(10, s2, level)
    call splitstatesatlevel(10, s3, level)
    call splitstatesatlevel(10, s4, level)

    close(10)

    fname = 'splitstatesraw.dat'
    open(10, file = trim(fname))

    call splitstatesraw(10, s1)
    call splitstatesraw(10, s2)
    call splitstatesraw(10, s3)
    call splitstatesraw(10, s4)

    close(10)

    fname = 'splitstatesatlevelraw.dat'
    open(10, file = trim(fname))

    call splitstatesatlevelraw(10, s1, level)
    call splitstatesatlevelraw(10, s2, level)
    call splitstatesatlevelraw(10, s3, level)
    call splitstatesatlevelraw(10, s4, level)

    close(10)

    fname = 'findmels.dat'
    open(10, file = trim(fname))

    call findmels(10, s1, s1)
    call findmels(10, s1, s2)
    call findmels(10, s2, s2)
    call findmels(10, s2, s3)
    call findmels(10, s3, s3)
    call findmels(10, s3, s4)
    call findmels(10, s4, s4)

    close(10)

    fname = 'findmelsatlevel.dat'
    open(10, file = trim(fname))

    call findmelsatlevel(10, s1, s1, level)
    call findmelsatlevel(10, s1, s2, level)
    call findmelsatlevel(10, s2, s2, level)
    call findmelsatlevel(10, s2, s3, level)
    call findmelsatlevel(10, s3, s3, level)
    call findmelsatlevel(10, s3, s4, level)
    call findmelsatlevel(10, s4, s4, level)

    close(10)

    nlines = 0

    open(11, file = 'splitstatesatlevelraw.dat')
    do
        read(11, *, end = 110)
        nlines = nlines + 1
    end do
    110 close(11)

    open(10, file = 'stateordinates.dat')
    open(11, file = 'splitstatesatlevelraw.dat')

    ! print *, nlines
    
    ! allocate(rawdata(nlines, 5))
    allocate(sizedata(nlines))

    
    ! rawdata = transpose(rawdata)
    
    do i = 1, nlines
        read(11, *) junk, junk, junk, junk, sizedata(i)
        if ( i .eq. 1) then
            write (*,*) sizedata(i), 0
        else
            write (*,*) sizedata(i), sum(sizedata(1:i-1))
        end if
    end do
    close(10)
    close(11)
    deallocate(sizedata)

end program spinless_fermions