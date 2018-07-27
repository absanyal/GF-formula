program sp_rgf_serial

    use constants
    implicit none

    integer, parameter :: n = 4

    real, allocatable :: w_list(:)
    complex, allocatable :: GF(:)
    integer :: iter, mainloop
    real :: w

    complex :: tempGF

    allocate(w_list(numpts))
    allocate(GF(n))

    do iter = 1, numpts
        w_list(iter) = minw + (iter-1) * dw
    enddo

    open(1, file = 'fgf.dat') 
    
    do mainloop = 1, numpts
        w = w_list(mainloop)

        GF(1) = 1/(w + i*eta - epsilon0)

        ! sweep right once
        do iter = 2, n
            GF(iter) = 1/( w + i*eta - epsilon0 - t * GF(iter-1) * t )
        enddo

        ! now sweep left
        do iter = n-1, 1, -1
            ! write (*,*) w, iter
            tempGF = 1/( 1/GF(iter) - t*GF(iter+1)*t )
            GF(iter) = tempGF
        enddo

        write(1, *) w, aimag(GF) * (-1.0/pi)

    enddo

end program sp_rgf_serial
