module statemanip

implicit none

    real, parameter :: eta = 0.1

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

    function mel(s1, s2)
        implicit none

        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2

        ! integer :: isdiagonal

        integer :: Hu, Hl, Hs, Hsp1, Hsm1
        integer :: mel

        integer :: vns1, valphas1, vbetas1
        integer :: vns2, valphas2, vbetas2

        vns1 = s1%ns
        valphas1 = s1%alphas
        vbetas1 = s1%betas

        vns2 = s2%ns
        valphas2 = s2%alphas
        vbetas2 = s2%betas

        Hs = 0
        Hsp1 = 0
        Hsm1 = 0

        if (s1%s .eq. s2%s) then
            ! if (vns1 .ne. vns2 .or. &
            !     valphas1 .ne. valphas2 .or. &
            !     vbetas1 .ne. vbetas2) then
                if (s1%s .ge. 1) then
                    Hsp1 = kdelta(vns1, vns2) &
                    * kdelta(valphas1, valphas2) * &
                            (kdelta(vbetas1, vbetas2 + 1) &
                            + kdelta(vbetas1, vbetas2 - 1))
                    
                    Hsm1 = kdelta(vns1, vns2) &
                    * kdelta(vbetas1, vbetas2) * &
                            (kdelta(valphas1, valphas2 + 1) &
                            + kdelta(valphas1, valphas2 - 1))
                    
                    Hs = kdelta(vns1, vns2-1) &
                    * kdelta(valphas1, valphas2-1) * &
                            kdelta(vbetas1, vbetas2 + 1) + &
                        kdelta(vns1, vns2+1) &
                        * kdelta(valphas1, valphas2+1) * &
                            kdelta(vbetas1, vbetas2 - 1)
                end if
                ! if (s1%s .eq. 1) then
                !     ! Hsm1 = 0
                !     Hs = 0
                !     Hsp1 = 0
                ! end if
            ! end if

            Hl = 0

            if (vns1 .eq. vns2 .and. &
                valphas1 .eq. valphas2 .and. &
                vbetas1 .eq. vbetas2 .and. s1%s .gt. 1) then
                Hl = 1
            end if

            mel = Hs + Hsp1 + Hsm1 + Hl

            ! if (s1%s .eq. 1 .and. mel .eq. 0) then
            !     if (vns1 .eq. vns2 .and. &
            !         valphas1 .eq. valphas2 .and. &
            !         vbetas1 .eq. vbetas2) then
            !         mel = 1
            !     end if
            ! end if

            ! if (s1%s .eq. 1 .and. mel .ne. 1) then
            !     if (vns1 .ne. vns2 .or. &
            !         valphas1 .ne. valphas2 .or. &
            !         vbetas1 .ne. vbetas2) then
            !         mel = 1
            !     end if
            ! end if
            
        else
            mel = 0
        end if

    end function mel

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


    recursive subroutine findmels(writeto, s1, s2, isdiagonal)
        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2
        integer, intent(in) :: isdiagonal
        
        integer :: setisdiag

        type(state) :: s1s0
        type(state) :: s1s1
        type(state) :: s2s0
        type(state) :: s2s1

        integer, intent(in) :: writeto
        character(len=100) :: fmt

        integer :: meltest = 0
        setisdiag = isdiagonal

        if (s1%alphas .ne. s2%alphas .or. s1%betas .ne. s2%betas .or. &
            s1%ns .ne. s2%ns) then
            setisdiag = 0
        end if

        fmt = "(4A, I2, A, 2I2)"

        if (checkvalidity(s1) .eq. 1 .and. &
        checkvalidity(s2) .eq. 1) then
            meltest = mel(s1, s2)
            if (setisdiag .ne. 1 .and. s1%s .eq. 1) then
                meltest = kdelta(s1%ns, s2%ns) &
                * kdelta(s1%alphas, s2%alphas) &
                * kdelta(s1%betas, s2%betas)
            end if
            if (s1%s .eq. 2 .and. setisdiag .eq. 1) then
                meltest = 0
            end if
            ! if (s1%s .eq. 1 .and. setisdiag .eq. 0) then
            !     meltest = mel(s1, s2)
            ! end if
            write (writeto,fmt) getlstate(s1) , char(9), getlstate(s2),&
             char(9), meltest, char(9), getstatesize(s1), getstatesize(s2)
        end if

        
        s1s0 = relegate0(s1)
        s1s1 = relegate1(s1)
        s2s0 = relegate0(s2)
        s2s1 = relegate1(s2)
        

        if (meltest .ne. 0) then
            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmels(writeto, s1s0, s2s0, setisdiag)
            end if

            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmels(writeto, s1s0, s2s1, setisdiag)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmels(writeto, s1s1, s2s0, setisdiag)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmels(writeto, s1s1, s2s1, setisdiag)
            end if

        end if

    end subroutine findmels

    recursive subroutine findmelsatlevel(writeto, s1, s2, level)
        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2

        type(state) :: s1s0
        type(state) :: s1s1
        type(state) :: s2s0
        type(state) :: s2s1

        integer, intent(in) :: writeto
        integer, intent(in) :: level
        integer :: meltest = 0

        character(len=100) :: fmt

        fmt = "(4A, I2, A, 2I2)"

        if (checkvalidity(s1) .eq. 1 .and. &
        checkvalidity(s2) .eq. 1) then
            meltest = mel(s1, s2)
            if (s1%s .eq. level) then
                write (writeto,fmt) getlstate(s1) , char(9), getlstate(s2),&
                char(9), meltest, char(9), getstatesize(s1), getstatesize(s2)
            end if
            if (s1%s .le. level) then
                meltest = 0
            end if
        end if

        s1s0 = relegate0(s1)
        s1s1 = relegate1(s1)
        s2s0 = relegate0(s2)
        s2s1 = relegate1(s2)

        if (meltest .ne. 0) then
            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmelsatlevel(writeto, s1s0, s2s0, level)
            end if

            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmelsatlevel(writeto, s1s0, s2s1, level)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmelsatlevel(writeto, s1s1, s2s0, level)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmelsatlevel(writeto, s1s1, s2s1, level)
            end if

        end if

    end subroutine findmelsatlevel

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
        ! character(len = 1) :: getletter
        integer :: getletter
        if (s1%alphas .eq. 0 .and. s1%betas .eq. 0) then
            ! write(getletter, "(A)") "a"
            getletter = 1
        end if

        if (s1%alphas .eq. 0 .and. s1%betas .eq. 1) then
            ! write(getletter, "(A)") "b"
            getletter = 2
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 0) then
            ! write(getletter, "(A)") "c"
            getletter = 3
        end if

        if (s1%alphas .eq. 1 .and. s1%betas .eq. 1) then
            ! write(getletter, "(A)") "d"
            getletter = 4
        end if
    end function getletter

    function complex_identity(m)
        integer, intent(in) :: m
        integer :: j
        complex, allocatable :: complex_identity(:,:)
        allocate(complex_identity(m, m))
        complex_identity = complex(0.0, 0.0)
        do j = 1, m
            complex_identity(j, j) = complex(1.0, 0.0)
        end do
    end function complex_identity

    function complex_zeros(m, n)
        integer, intent(in) :: m, n
        integer :: j
        complex, allocatable :: complex_zeros(:,:)
        allocate(complex_zeros(m, n))
        complex_zeros = complex(0.0, 0.0)
    end function complex_zeros

    subroutine matprint(h)
        integer :: m
        complex, dimension(:, :) :: h
        integer :: j
        m = size(h, 1)
        do j = 1, m, 1
            write(*,*) h(j, :)
        end do
    end subroutine matprint

    subroutine matprinttofile(writeto, h)
        integer :: m
        integer, intent(in) :: writeto
        complex, dimension(:, :) :: h
        integer :: j
        m = size(h, 1)
        do j = 1, m, 1
            write(writeto, *) h(j, :)
        end do
    end subroutine matprinttofile

    function connecting_tau(s1, s2)
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
        complex, allocatable :: connecting_tau(:, :)
        integer :: startindex
        integer :: i, j

        allocate(connecting_tau(getstatesize(s1), getstatesize(s2)))
        connecting_tau = complex(0, 0)

        if (checkvalidity(s1) .eq. 1 .and. checkvalidity(s2) .eq. 1) then
            l1 = getletter(s1)
            l2 = getletter(s2)
            if (l1 .eq. 3 .and. l2 .eq. 2) then
                connecting_tau = complex_zeros(getstatesize(s1), getstatesize(s2))
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
                    connecting_tau(i, j) = 1
                    j = j + 1
                end do
            end if
            
        end if
    end function connecting_tau

    function hcat( A, B ) result( X )
        complex, dimension(:,:) :: A, B
        complex :: X( size(A,1), size(A,2)+size(B,2) )

        X = reshape( [ A, B], shape( X ) )
    end function hcat

    function vcat( A, B ) result( X )
        complex, dimension(:,:) :: A, B
        complex :: X( size(A,1)+size(B,1), size(A,2) )

        X = transpose( reshape( &
                [ transpose(A), transpose(B) ], &
                [ size(X,2), size(X,1) ] ) )
    end function vcat

    function bmat(a, b, c, d) result (x)
        complex, dimension(:,:) :: a, b, c, d
        complex :: x(size(a, 2) + size(b, 2), size(a, 1) + size(c, 1))
        X = vcat( hcat(a, b), hcat(c, d) )
    end function bmat

    function onebyone()
        complex, DIMENSION(1, 1) :: onebyone
        ! onebyone = reshape( (/ 0 /), &
        !                                    shape(onebyone), order=(/2,1/) )
        onebyone(1, 1) = complex(0.0, 0.0)
    end function onebyone

    function twobytwo()
        complex, DIMENSION(2, 2) :: twobytwo
        ! twobytwo =reshape( (/ 0, 1, 1, 0 /), &
        !                                    shape(twobytwo), order=(/2,1/) )
        twobytwo(1, 1) = complex(0.0, 0.0)
        twobytwo(1, 2) = complex(1.0, 0.0)
        twobytwo(2, 1) = complex(1.0, 0.0)
        twobytwo(2, 2) = complex(0.0, 0.0)
    end function twobytwo

    function diagonalblock(s1)
        type (state), intent(inout) :: s1
        complex, allocatable :: diagonalblock(:, :)
        if (getstatesize(s1) .eq. 1) then
            allocate(diagonalblock(1, 1))
            diagonalblock = onebyone()
        end if
        if (getstatesize(s1) .eq. 2) then
            allocate(diagonalblock(2, 2))
            diagonalblock = twobytwo()
        end if
    end function diagonalblock

    function connectedblock(s1, s2)
        type (state), intent(inout) :: s1
        type (state), intent(inout) :: s2
        complex, allocatable :: connectedblock(:,:)
        complex, allocatable :: h1(:,:)
        complex, allocatable :: h2(:,:)
        complex, allocatable :: htau12(:,:)
        complex, allocatable :: htau21(:,:)
        allocate( h1 ( getstatesize(s1), getstatesize(s1) ) )
        allocate(h2(getstatesize(s2), getstatesize(s2)))
        allocate(htau12(getstatesize(s1), getstatesize(s2)))
        allocate(htau21(getstatesize(s2), getstatesize(s1)))
        allocate(connectedblock( getstatesize(s1) + getstatesize(s2), &
                                 getstatesize(s1) + getstatesize(s2) ) )
        h1 = diagonalblock(s1)
        h2 = diagonalblock(s2)
        htau12 = connecting_tau(s1, s2)
        htau21 = transpose(htau12)
        connectedblock = bmat(h1, htau12, htau21, h2)
        ! connectedblock = reshape([h1, htau12, htau21, h2], &
        !                 shape(connectedblock), order = [2, 1] )
        ! connectedblock = vcat(hcat(h1, htau12), hcat(htau21, h2))
        deallocate(h1)
        deallocate(h2)
        deallocate(htau12)
        deallocate(htau21)
    end function connectedblock

    function inv(A) result(Ainv)
        complex, dimension(:,:), intent(in) :: A
        complex, dimension(size(A,1),size(A,2)) :: Ainv

        complex, dimension(size(A,1)) :: work  ! work array for LAPACK
        integer, dimension(size(A,1)) :: ipiv   ! pivot indices
        integer :: n, info

        ! External procedures defined in LAPACK
        external cGETRF
        external cGETRI

        ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)

        ! DGETRF computes an LU factorization of a general M-by-N matrix A
        ! using partial pivoting with row interchanges.
        call cGETRF(n, n, Ainv, n, ipiv, info)
        ! print *, info

        if (info /= 0) then
            stop !'Matrix is numerically singular!'
        end if

        ! DGETRI computes the inverse of a matrix using the LU factorization
        ! computed by DGETRF.
        call cGETRI(n, Ainv, n, ipiv, work, n, info)
        !  write (*,*) info

        if (info /= 0) then
            stop !'Matrix inversion failed!'
        end if
    end function inv

    function g(h, w)
        complex, allocatable :: g(:,:)
        complex, allocatable, intent(in) :: h(:,:)
        real, intent(in) :: w
        g = inv((w + complex(0, eta)) * complex_identity(size(h, 1)) &
            - h)
    end function g

    function trace(h)
        complex, allocatable, intent(in) :: h(:,:)
        complex :: trace
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

        complex, allocatable :: outputmatrix(:,:)

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
        complex, allocatable :: A(:,:)
        complex, allocatable :: B(:,:)
        complex, allocatable :: matprod(:,:)

        integer :: i, j, n

        allocate(matprod(size(A, 1), size(B, 2)))

        matprod = 0

        do i = 1, size(matprod, 1)
            do j = 1, size(matprod, 2)
                do n = 1, size(matprod, 2)
                    matprod(i, j) = matprod(i, j) + A(i, n) * B(n, j)
                end do
            end do
        end do
    end function

end module statemanip