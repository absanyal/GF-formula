module statemanip

implicit none

    real(8), parameter :: eta = 0.1

    type state
        integer :: s = 0
        integer :: ns = 0
        integer :: alphas = 0
        integer :: betas = 0
    end type state

contains

    function kdelta(x, y)
        integer, intent(in) :: x
        integer, intent(in) :: y
        integer :: kdelta
        if (x .eq. y) then
            kdelta = 1
        else
            kdelta = 0
        end if
        return
    end function kdelta

    function factorial(n)
        integer, intent(in) :: n
        integer :: factorial
        integer :: i = 1
        factorial = 1
        do i = 2, n, 1
            factorial = factorial * i
        end do
    end function factorial

    subroutine setstate(astate, vs, vn_s, valpha_s, vbeta_s)
        type(state), intent(inout) :: astate
        integer, intent(in) :: vs
        integer, intent(in) :: vn_s
        integer, intent(in) :: valpha_s
        integer, intent(in) :: vbeta_s
        astate%s = vs
        astate%ns = vn_s
        astate%alphas = valpha_s
        astate%betas = vbeta_s
    end subroutine setstate

    function getstate(astate)
        type(state), intent(inout) :: astate
        integer, dimension(4) :: getstate
        getstate(1) = astate%s
        getstate(2) = astate%ns
        getstate(3) = astate%alphas
        getstate(4) = astate%betas
    end function getstate

    function checkvalidity(s1)
        type(state), intent(inout) :: s1
        integer :: checkvalidity
        checkvalidity = 1
        if (s1%s .lt. 1 .or. s1%ns .lt. 0) then
            checkvalidity = 0
        end if
        if (s1%s .eq. 0 .and. s1%ns .eq. 0 .and. &
            s1%alphas .eq. 0 .and. s1%betas .eq. 0) then
            checkvalidity = 0
        end if
        if (s1%s - 1 .lt. s1%ns - s1%alphas) then
            checkvalidity = 0
        end if
        if (s1%alphas .gt. s1%ns) then
            checkvalidity = 0
        end if
    end function

    function getstatesize(s1)
        type(state), intent(inout) :: s1
        integer :: getstatesize
        integer :: vs, vns, va, vb
        vs = s1%s
        vns = s1%ns
        va = s1%alphas
        vb = s1%betas
        getstatesize = 0
        if (checkvalidity(s1) .eq. 1) then
            if (va .eq. 0 .and. vb .eq. 0) then
                getstatesize = &
                factorial(vs-1)/(factorial(vns)*factorial(vs-1-vns))
            end if
            if (va .eq. 0 .and. vb .eq. 1) then
                getstatesize = &
                factorial(vs-1)/(factorial(vns)*factorial(vs-1-vns))
            end if
            if (va .eq. 1 .and. vb .eq. 0) then
                getstatesize = &
                factorial(vs-1)/(factorial(vns-1)*factorial(vs-vns))
            end if
            if (va .eq. 1 .and. vb .eq. 1) then
                getstatesize = &
                factorial(vs-1)/(factorial(vns-1)*factorial(vs-vns))
            end if
        else
            getstatesize = 0
        end if
    end function

    function relegate0(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate0
        integer :: sm1, nsm1, alphasm1, betasm1
        sm1 = s1%s - 1
        betasm1 = s1%alphas
        alphasm1 = 0
        nsm1 = s1%ns - betasm1
        call setstate(relegate0, sm1, nsm1, alphasm1, betasm1)
    end function

    function relegate1(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate1
        integer :: sm1, nsm1, alphasm1, betasm1
        sm1 = s1%s - 1
        betasm1 = s1%alphas
        alphasm1 = 1
        nsm1 = s1%ns - betasm1
        call setstate(relegate1, sm1, nsm1, alphasm1, betasm1)
    end function

    function getlstate(s1)
        type(state), intent(inout) :: s1
        character(len=10) :: getlstate
        character(len=2) :: the_s
        character(len=2) :: k
        character(len=2) :: the_ns

        write(the_s, "(A2)") s1%s
        write(the_ns, "(A2)") s1%ns

        if (s1%alphas .eq. 0 .and. s1%betas .eq. 0) then
            write(k, "(A)") " a"
        end if

        if (s1%alphas .eq. 0 .and. s1%betas .eq. 1) then
            write(k, "(A)") " b"
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 0) then
            write(k, "(A)") " c"
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 1) then
            write(k, "(A)") " d"
        end if

        write(getlstate, "(I2, A2, I2)") s1%s, k, s1%ns
    
    end function

    recursive subroutine splitstatesraw(writeto, s1)
        type(state), intent(inout) :: s1
        integer , intent(in) :: writeto

        type(state) :: s1s0
        type(state) :: s1s1

        if (checkvalidity(s1) .eq. 1) then
            write (writeto,*) getstate(s1), getstatesize(s1)
            s1s0 = relegate0(s1)
            call splitstatesraw(writeto, s1s0)
            s1s1 = relegate1(s1)
            call splitstatesraw(writeto, s1s1)
        end if

    end subroutine

    recursive subroutine splitstatesatlevelraw(writeto, s1, level)
        type(state), intent(inout) :: s1
        integer , intent(in) :: writeto
        integer, intent(in) :: level

        type(state) :: s1s0
        type(state) :: s1s1

        if (checkvalidity(s1) .eq. 1) then
            if (s1%s .eq. level) then
                write (writeto,*) getstate(s1),getstatesize(s1)
            end if
            if (s1%s .gt. level) then
                s1s0 = relegate0(s1)
                call splitstatesatlevelraw(writeto, s1s0, level)
                s1s1 = relegate1(s1)
                call splitstatesatlevelraw(writeto, s1s1, level)
            end if
        end if

    end subroutine

    recursive subroutine splitstates(writeto, s1)
        type(state), intent(inout) :: s1
        integer , intent(in) :: writeto

        type(state) :: s1s0
        type(state) :: s1s1

        if (checkvalidity(s1) .eq. 1) then
            write (writeto,*) getlstate(s1), char(9), getstatesize(s1)
            s1s0 = relegate0(s1)
            call splitstates(writeto, s1s0)
            s1s1 = relegate1(s1)
            call splitstates(writeto, s1s1)
        end if

    end subroutine

    recursive subroutine splitstatesatlevel(writeto, s1, level)
        type(state), intent(inout) :: s1
        integer , intent(in) :: writeto
        integer, intent(in) :: level

        type(state) :: s1s0
        type(state) :: s1s1

        if (checkvalidity(s1) .eq. 1) then
            if (s1%s .eq. level) then
                write (writeto,*) getlstate(s1), char(9), getstatesize(s1)
            end if
            if (s1%s .gt. level) then
                s1s0 = relegate0(s1)
                call splitstatesatlevel(writeto, s1s0, level)
                s1s1 = relegate1(s1)
                call splitstatesatlevel(writeto, s1s1, level)
            end if
        end if

    end subroutine

    function getletter(s1)
        type(state), intent(inout) :: s1
        integer :: getletter
        if (s1%alphas .eq. 0 .and. s1%betas .eq. 0) then
            getletter = 1
        end if

        if (s1%alphas .eq. 0 .and. s1%betas .eq. 1) then
            getletter = 2
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 0) then
            getletter = 3
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 1) then
            getletter = 4
        end if
    end function getletter

    function complex_identity(m)
        integer, intent(in) :: m
        integer :: j
        complex(8), allocatable :: complex_identity(:,:)
        allocate(complex_identity(m, m))
        complex_identity = complex(0.0, 0.0)
        do j = 1, m
            complex_identity(j, j) = complex(1.0, 0.0)
        end do
    end function complex_identity

    function complex_zeros(m, n)
        integer, intent(in) :: m, n
        integer :: j
        complex(8), allocatable :: complex_zeros(:,:)
        allocate(complex_zeros(m, n))
        complex_zeros = complex(0.0, 0.0)
    end function complex_zeros

    subroutine matprint(h)
        integer :: m
        complex(8), dimension(:, :) :: h
        integer :: j
        m = size(h, 1)
        do j = 1, m, 1
            write(*,*)  h(j, :)
        end do
    end subroutine matprint

    subroutine matprinttofile(writeto, h)
        integer :: m
        integer, intent(in) :: writeto
        complex(8), dimension(:, :) :: h
        integer :: j
        m = size(h, 1)
        do j = 1, m, 1
            write(writeto, *) h(j, :)
        end do
    end subroutine matprinttofile

    function connecting_u(s1, s2)
        type(state), intent(inout) :: s1
        type(state) :: s1s0
        type(state) :: s1s1
        type(state), intent(inout) :: s2
        type(state) :: s2s0
        type(state) :: s2s1
        integer :: l1
        integer :: l2
        integer :: l11
        integer :: l12
        integer :: l21
        integer :: l22
        complex(8), allocatable :: connecting_u(:, :)
        integer :: startindex
        integer :: i, j

        allocate(connecting_u(getstatesize(s1), getstatesize(s2)))
        connecting_u = complex(0.0, 0.0)

        if (checkvalidity(s1) .eq. 1 .and. checkvalidity(s2) .eq. 1) then
            l1 = getletter(s1)
            l2 = getletter(s2)
            if (l1 .eq. 3 .and. l2 .eq. 2) then
                connecting_u = complex_zeros(getstatesize(s1), getstatesize(s2))
            else
                s1s0 = relegate0(s1)
                s1s1 = relegate1(s1)
                s2s0 = relegate0(s2)
                s2s1 = relegate1(s2)
                l11 = getletter(s1s0)
                l12 = getletter(s1s1)
                l21 = getletter(s2s0)
                l22 = getletter(s2s1)
                startindex = getstatesize(s1s0) + 1
                j = 1
                do i = startindex, getstatesize(s1), 1
                    connecting_u(i, j) = 1
                    j = j + 1
                end do
            end if
            
        end if
    end function connecting_u

    function hcat( A, B ) result( X )
        complex(8), dimension(:,:) :: A, B
        complex(8) :: X( size(A,1), size(A,2)+size(B,2) )

        X = reshape( [ A, B], shape( X ) )
    end function hcat

    function vcat( A, B ) result( X )
        complex(8), dimension(:,:) :: A, B
        complex(8) :: X( size(A,1)+size(B,1), size(A,2) )

        X = transpose( reshape( &
                [ transpose(A), transpose(B) ], &
                [ size(X,2), size(X,1) ] ) )
    end function vcat

    function bmat(a, b, c, d) result (x)
        complex(8), dimension(:,:) :: a, b, c, d
        complex(8) :: x(size(a, 2) + size(b, 2), size(a, 1) + size(c, 1))
        X = vcat( hcat(a, b), hcat(c, d) )
    end function bmat

    function onebyone()
        complex(8), DIMENSION(1, 1) :: onebyone
        onebyone(1, 1) = complex(0.0, 0.0)
    end function onebyone

    function twobytwo()
        complex(8), DIMENSION(2, 2) :: twobytwo
        twobytwo(1, 1) = complex(0.0, 0.0)
        twobytwo(1, 2) = complex(1.0, 0.0)
        twobytwo(2, 1) = complex(1.0, 0.0)
        twobytwo(2, 2) = complex(0.0, 0.0)
    end function twobytwo

    function diagonalblock(s1)
        type (state), intent(inout) :: s1
        complex(8), allocatable :: diagonalblock(:, :)
        if (getstatesize(s1) .eq. 1) then
            allocate(diagonalblock(1, 1))
            diagonalblock = onebyone()
        end if
        if (getstatesize(s1) .eq. 2) then
            allocate(diagonalblock(2, 2))
            diagonalblock = twobytwo()
        end if
    end function diagonalblock

    function inv(A) result(Ainv)
        complex(8), dimension(:,:), intent(in) :: A
        complex(8), dimension(size(A,1),size(A,2)) :: Ainv

        complex(8), dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external zGETRF
        external zGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call zGETRF(n, n, Ainv, n, ipiv, info)
        ! print *, info

        if (info /= 0) then
            stop !'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call zGETRI(n, Ainv, n, ipiv, work, n, info)
        !  write (*,*) info

        if (info /= 0) then
            stop !'Matrix inversion failed!'
        end if
    end function inv

    function g(h, w)
        complex(8), allocatable :: g(:,:)
        complex(8), allocatable, intent(in) :: h(:,:)
        real(8), intent(in) :: w
        g = inv((w + complex(0, eta)) * complex_identity(size(h, 1)) &
            - h)
    end function g

    function trace(h)
        complex(8), allocatable, intent(in) :: h(:,:)
        complex(8) :: trace
        integer :: j
        j = 1
        trace = 0
        do j = 1, size(h, 1)
            trace = trace + h(j, j)
        end do
    end function

    subroutine setfiletomatrix(inputname, outputmatrix)
        character(len = *) :: inputname
        integer :: nlines
        integer :: i

        complex(8), allocatable :: outputmatrix(:,:)

        nlines = 0

        open(15, file = inputname)
        do
            read(15, *, end = 110)
            nlines = nlines + 1
        end do
        110 close(15)

        open(20, file = inputname)
        allocate(outputmatrix(nlines, nlines))
        outputmatrix = (0.0, 0.0)
        do i = 1, nlines
            read(20,*) outputmatrix(i, :)
        end do
        close(20)

    end subroutine

    function matprod(A, B)
        complex(8), allocatable :: A(:,:)
        complex(8), allocatable :: B(:,:)
        complex(8), allocatable :: matprod(:,:)

        integer :: i, j, n

        allocate(matprod(size(A, 1), size(B, 2)))

        do i = 1, size(matprod, 1)
            do j = 1, size(matprod, 2)
                do n = 1, size(matprod, 2)
                    matprod(i, j) = complex(0.0, 0.0)
                end do
            end do
        end do

        do i = 1, size(matprod, 1)
            do j = 1, size(matprod, 2)
                do n = 1, size(matprod, 2)
                    matprod(i, j) = matprod(i, j) + &
                    A(i, n) * B(n, j)
                end do
            end do
        end do

        ! deallocate(matprod)
    end function

    function tmm(A, B, C)
        complex(8), allocatable :: A(:,:)
        complex(8), allocatable :: B(:,:)
        complex(8), allocatable :: C(:,:)
        complex(8), allocatable :: tmm(:,:)

        tmm = matprod(A, matprod(B, C))
    end function

    function cmplxfilterbelow(z)
        complex(8), intent(in) :: z
        complex(8) :: cmplxfilterbelow
        real(8) :: modsqz
        modsqz = z * conjg(z)
        if ( modsqz .le. 1.0E-20) then
            cmplxfilterbelow = complex(0.0, 0.0)
        end if
    end function

end module statemanip