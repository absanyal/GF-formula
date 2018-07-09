module constants

    implicit none
    
    complex, parameter :: i = (0, 1)                    ! imaginary i
    real, parameter :: pi = 4 * atan(1.0)
    
    real, parameter :: epsilon0 = 0                     ! on-site energy
    real, parameter :: eta = 0.1                        ! regulator
    real, parameter :: t = 1.0                          ! hopping constant
    integer, parameter :: numpts = 10000                !number of points
    real, parameter :: minw = -3.0 * t + epsilon0       !lower bound of energy
    real, parameter :: maxw = 3.0 * t + epsilon0        !upper bound of energy
    real, parameter :: dw = (maxw - minw) / numpts      !interval length

end module constants
