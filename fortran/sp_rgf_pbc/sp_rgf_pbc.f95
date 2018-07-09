program sp_rgf_pbc
    use constants
    implicit none

    integer, parameter :: n = 100
    real, parameter :: eta1 = eta * 0

    complex, dimension(:), allocatable :: GF
    real, dimension(:), allocatable :: onsite
    real, dimension(:), allocatable :: w_list 

    integer :: j, iter ! counter variable
    real :: w

    allocate(GF(n))
    allocate(w_list(numpts))
    allocate(onsite(n))

    ! fill the energy list
    do j = 1, numpts
        w_list(j) = minw + (j - 1) * dw
    enddo

    ! fill the on-site energy array
    do j = 1, n
        onsite(j) = epsilon0
    enddo

    open(1, file = 'fgf.dat')
    open(2, file = 'dos.dat')

    !main loop
    do j = 1, numpts
        w = w_list(j)

        GF(1) = 1 / ( w + i * eta - onsite(1) )

        !sweep left once
        do iter = 2, n
            GF(iter) = 1/( w + i * eta1 - onsite(iter) - t * GF(iter-1) * t )
        enddo

        !connect last site to first to get full GF of site 1
        GF(1) = 1/( w + i * eta1 - onsite(1) - t * GF(n) * t )

        !sweep again through the remaining sites to get full GF of each site
        do iter = 2, n
            GF(iter) = 1/( w + i * eta1 - onsite(iter) - t * GF(iter-1) * t )
        enddo

        write (1, *) w, aimag(GF) * (-1/pi)
        write (2, *) w, sum(aimag(GF)) * (-1/pi) / float(n)

    enddo

end program sp_rgf_pbc
