module gf

    implicit none

contains

    function inv(A) result(Ainv)
        real, dimension(:,:), intent(in) :: A
        real, dimension(size(A,1),size(A,2)) :: Ainv

        integer, dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external DGETRF
        external DGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call SGETRF(n, n, Ainv, n, ipiv, info)
        ! print *, info

        if (info /= 0) then
            stop !'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call SGETRI(n, Ainv, n, ipiv, work, n, info)
        !  write (*,*) info

        if (info /= 0) then
            stop !'Matrix inversion failed!'
        end if
    end function inv

    ! function hello(a)
    !     integer , intent(in) :: a
    !     integer :: hello

    !     hello = 2 * a
    ! end function hello

end module gf