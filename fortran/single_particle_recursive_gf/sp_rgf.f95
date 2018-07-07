program sp_rgf

    implicit none

    10 format ("(", f8.4, ",", x,  f8.4, ")" )

    ! declarations
    integer, parameter :: n = 20                    ! number of sites
    
    complex, parameter :: i = (0, 1)                    ! imaginary i
    real, parameter :: pi = 4 * atan(1.0)
    
    real, parameter :: epsilon0 = 1                     ! on-site energy
    real, parameter :: eta = 0.1                       ! regulator
    real, parameter :: t = 1.0                          ! hopping constant
    integer, parameter :: numpts = 10000                 !number of points
    real, parameter :: minw = -3.0 * t + epsilon0                 !lower bound of energy
    real, parameter :: maxw = 3.0 * t + epsilon0                  !upper bound of energy
    real, parameter :: dw = (maxw - minw) / numpts      !interval length
               
    complex, dimension(:), allocatable :: LGF           ! left connected GF
    complex, dimension(:), allocatable :: RGF           ! right connected GF
    complex, dimension(:), allocatable :: FGF           ! fully connected GF
    real, dimension(:),allocatable :: onsite            ! on-site energy array
    real, dimension(:), allocatable :: w_list           ! energy list

    integer :: j, iterl, iterr, iterf                   !counter variables

    real :: w
    
    ! allocation of arrays
    allocate(LGF(n))
    allocate(RGF(n))
    allocate(FGF(n))
    allocate(onsite(n))
    allocate(w_list(numpts))

    ! fill the energy list
    do j = 1, numpts
        w_list(j) = minw + (j - 1) * dw
    enddo

    ! fill the on-site energy array
    do j = 1, n
        onsite(j) = epsilon0
    enddo

    open(1, file = 'fgf.dat')
    
    !main loop
    do j = 1, numpts
        w = w_list(j)
        
        LGF(1) = 1 / ( w + i * eta - onsite(1) )
        do iterl = 2, n
            LGF(iterl) = 1/( w - onsite(iterl) - t * LGF(iterl-1) * t )
        enddo

        RGF(n) = 1 / ( w + i * eta - onsite(n) )
        do iterr = n-1, 1, -1
            RGF(iterr) = 1/( w - onsite(iterr) - t * LGF(iterr+1) * t )
        enddo

        do iterf = 1, n
            if ( iterf == 1 ) then
                FGF(iterf) = 1 / ( w - onsite(iterf) &
                - t * RGF(iterf+1) * t )
            end if
            if ( iterf == n ) then
                FGF(iterf) = 1 / ( w - onsite(iterf) &
                - t * LGF(iterf-1) * t )
            end if
            if ( iterf >1 .and. iterf < n ) then
                FGF(iterf) = 1 / ( w - onsite(iterf) &
                - t * LGF(iterf-1) * t  - t * RGF(iterf+1) * t )
            end if
        enddo

        write (1, *) w, sum(aimag(FGF)) * (-1/pi) / float(n)

    enddo

    close(1)

end program sp_rgf