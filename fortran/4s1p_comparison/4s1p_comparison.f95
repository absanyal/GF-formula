program four_s_one_p_comparison
    use constants
    implicit none
    
    real :: w
    complex :: GF
    complex :: rGF11
    complex :: rGF22
    complex :: rGF33

    complex :: g1
    complex :: g2
    complex :: g3

    open(1, file = 'exact.dat')
    open(2, file = 'rgf.dat')

    w = minw
    do while ( w .le. maxw)
        
        g1 = w + i*eta - epsilon0
        g2 = w + i*eta - epsilon0
        g3 = w + i*eta - epsilon0

        GF = 1/( (1/g1 - t*g3*t) - ( t + t*g3*t ) * (t + t*g3*t) / &
            ( 1/g2 - t*g3*t ) )
        
        write(1, *) w, aimag(GF)

        rGF22 = 1/( 1/g2 - t*t/( 1/g1 - t*t/( 1/g3 - t*t*g2 ) ) )
        rGF33 = 1/( 1/g3 - t*t/( 1/g1 - t*t/( 1/g2 - t*t*g3 ) ) )

        rGF11 = 1/( 1/g1 - t*rGF22*t - t*rGF33*t )

        write(2, *) w, aimag(rGF11)
        
        w = w + dw
    enddo

end program four_s_one_p_comparison
