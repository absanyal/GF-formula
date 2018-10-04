module statemanip

implicit none

    type state
        integer :: s = 0
        integer :: n_s = 0
        integer :: alpha_s = 0
        integer :: beta_s = 0
    end type state

contains

    integer function kdelta(x, y) result(output)
        integer, intent(in) :: x
        integer, intent(in) :: y
        if (x .eq. y) then
            output = 1
        else
            output = 0
        end if
        return
    end function kdelta

    subroutine setstate(astate, vs, vn_s, valpha_s, vbeta_s)
        type(state), intent(inout) :: astate
        integer, intent(in) :: vs
        integer, intent(in) :: vn_s
        integer, intent(in) :: valpha_s
        integer, intent(in) :: vbeta_s
        astate%s = vs
        astate%n_s = vn_s
        astate%alpha_s = valpha_s
        astate%beta_s = vbeta_s
    end subroutine setstate

    function getstate(astate)
        type(state), intent(inout) :: astate
        integer, dimension(4) :: getstate
        getstate(1) = astate%s
        getstate(2) = astate%n_s
        getstate(3) = astate%alpha_s
        getstate(4) = astate%beta_s
    end function getstate

    function mel(s1, s2)
        implicit none

        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2

        integer :: Hu, Hl, Hs, Hsp1, Hsm1

        integer :: mel

        integer :: vns1, valphas1, vbetas1
        integer :: vns2, valphas2, vbetas2

        vns1 = s1%n_s
        valphas1 = s1%alpha_s
        vbetas1 = s1%beta_s

        vns2 = s2%n_s
        valphas2 = s2%alpha_s
        vbetas2 = s2%beta_s

        if (s1%s .eq. s2%s) then
            if (vns1 .ne. vns2 .or. &
                valphas1 .ne. valphas2 .or. &
                vbetas1 .ne. vbetas2) then
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
                vbetas1 .eq. vbetas2 .and. s1%s .gt. 2) then
                Hl = 1
            else
                Hl = 0
            end if
            
            mel = Hs + Hsp1 + Hsm1 + Hl
        else
            mel = 0
        end if

    end function mel

    function checkvalidity(s1)
        type(state), intent(inout) :: s1
        integer :: checkvalidity
        checkvalidity = 1
        if (s1%s .lt. 2 .or. s1%n_s .lt. 0) then
            checkvalidity = 0
        end if
        if (s1%s .eq. 0 .and. s1%n_s .eq. 0 .and. &
            s1%alpha_s .eq. 0 .and. s1%beta_s .eq. 0) then
            checkvalidity = 0
        end if
        if (s1%s .le. s1%n_s - s1%alpha_s) then
            checkvalidity = 0
        end if
    end function

    function relegate0(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate0

        if (s1%s - 1 .gt. s1%n_s - (0+s1%alpha_s)) then
            call setstate(relegate0, &
                (s1%s)-1, s1%n_s - (1+s1%alpha_s), 0, s1%alpha_s)
        end if
    end function

    function relegate1(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate1
        if (s1%s .gt. s1%n_s - (s1%alpha_s)) then
            call setstate(relegate1, &
                (s1%s)-1, s1%n_s - (s1%alpha_s), 1, s1%alpha_s)
        end if
    
    end function

    function getlstate(s1)
        type(state), intent(inout) :: s1
        character(len=10) :: getlstate
        character(len=2) :: the_s
        character(len=2) :: k
        character(len=2) :: the_ns

        write(the_s, "(A2)") s1%s
        write(the_ns, "(A2)") s1%n_s

        if (s1%alpha_s .eq. 0 .and. s1%beta_s .eq. 0) then
            write(k, "(A)") " a"
        end if

        if (s1%alpha_s .eq. 0 .and. s1%beta_s .eq. 1) then
            write(k, "(A)") " b"
        end if

        if (s1%alpha_s .eq. 1 .and. s1%beta_s .eq. 0) then
            write(k, "(A)") " c"
        end if

        if (s1%alpha_s .eq. 1 .and. s1%beta_s .eq. 1) then
            write(k, "(A)") " d"
        end if

        write(getlstate, "(I2, A2, I2)") s1%s, k, s1%n_s
    
    end function


    recursive subroutine findmels(s1, s2)
        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2

        type(state) :: s1s0
        type(state) :: s1s1
        type(state) :: s2s0
        type(state) :: s2s1

        integer :: meltest = 0

        character(len=100) :: fmt

        fmt = "(4A, I2)"

        meltest = mel(s1, s2)

        if (checkvalidity(s1) .eq. 1 .and. &
        checkvalidity(s2) .eq. 1) then
            write (*,fmt) getlstate(s1) , char(9), getlstate(s2),&
             char(9), meltest
        end if

        s1s0 = relegate0(s1)
        s1s1 = relegate1(s1)
        s2s0 = relegate0(s2)
        s2s1 = relegate1(s2)

        if (meltest .ne. 0) then

            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmels(s1s0, s2s0)
            end if

            if (checkvalidity(s1s0) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmels(s1s0, s2s1)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s0) .eq. 1) then
                call findmels(s1s1, s2s0)
            end if

            if (checkvalidity(s1s1) .eq. 1 .and. &
            checkvalidity(s2s1) .eq. 1) then
                call findmels(s1s1, s2s1)
            end if

        end if

    end subroutine findmels

end module statemanip