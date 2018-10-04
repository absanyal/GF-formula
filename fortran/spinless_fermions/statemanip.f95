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

        real :: Hu, Hl, Hs, Hsp1, Hsm1

        real :: mel

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
        logical :: checkvalidity
        if (s1%s .eq. 0 .and. s1%n_s .eq. 0 .and. &
            s1%alpha_s .eq. 0 .and. s1%beta .eq. 0) then
            checkvalidity = .false.
        end if
        if (s1%s .ge. 2 .and. s1%n_s .ge. 0) then
            checkvalidity = .true.
        else
            checkvalidity = .false.
        end if
    end function

    function relegate0(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate0

        if (s1%n_s - (0+s1%alpha_s) .ge. 0) then
            call setstate(relegate0, &
                (s1%s)-1, s1%n_s - (0+s1%alpha_s), 0, s1%alpha_s)
        end if
    end function

    function relegate1(s1)
        type(state), intent(inout) :: s1
        type(state) :: relegate1
        if (s1%n_s - (1+s1%alpha_s) .ge. 1) then
            call setstate(relegate1, &
                (s1%s)-1, s1%n_s - (1+s1%alpha_s), 1, s1%alpha_s)
        end if
    
    end function

    ! subroutine findmels
    !     type(state), intent(inout) :: s1
    !     type(state), intent(inout) :: s2


        
    ! end subroutine findmels

end module statemanip