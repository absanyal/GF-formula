module statemanip

implicit none

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
        if (s1%s .lt. 2 .or. s1%ns .lt. 0) then
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

        ! isdiagonal = 1

        if (s1%s .eq. s2%s) then
            if (vns1 .ne. vns2 .or. &
                valphas1 .ne. valphas2 .or. &
                vbetas1 .ne. vbetas2) then
                ! isdiagonal = 0
                Hsp1 = kdelta(vns1, vns2) * kdelta(valphas1, valphas2) * &
                        (kdelta(vbetas1, vbetas2 + 1) &
                        + kdelta(vbetas1, vbetas2 - 1))
                
                Hsm1 = kdelta(vns1, vns2) * kdelta(vbetas1, vbetas2) * &
                        (kdelta(valphas1, valphas2 + 1) &
                        + kdelta(valphas1, valphas2 - 1))
                
                Hs = kdelta(vns1, vns2-1) * kdelta(valphas1, valphas2-1) * &
                        kdelta(vbetas1, vbetas2 + 1) + &
                    kdelta(vns1, vns2+1) * kdelta(valphas1, valphas2+1) * &
                        kdelta(vbetas1, vbetas2 - 1)
            else
                Hs = 0
                Hsp1 = 0
                Hsm1 = 0
            end if

            if (vns1 .eq. vns2 .and. &
                valphas1 .eq. valphas2 .and. &
                vbetas1 .eq. vbetas2) then
                Hl = 1
            else
                Hl = 0
            end if
            
            mel = Hs + Hsp1 + Hsm1 + Hl
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
            if (s1%s .eq. 2 .and. setisdiag .eq. 1) then
                meltest = 0
            end if
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

end module statemanip