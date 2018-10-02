module statemanip

implicit none

    type state
        integer :: s
        integer :: n_s
        integer :: alpha_s
        integer :: beta_s
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

    subroutine getstate(astate)
        type(state), intent(inout) :: astate
        integer, dimension(4) :: output
        output(1) = astate%s
        output(2) = astate%n_s
        output(3) = astate%alpha_s
        output(4) = astate%beta_s
        write (*,*) output
    end subroutine getstate

    recursive real function mel(s1, s2) result(fullmel)
        implicit none

        type(state), intent(inout) :: s1
        type(state), intent(inout) :: s2

        real :: Hu, Hl, Hs, Hsp1, Hsm1

        integer :: vns1, valphas1, vbetas1
        integer :: vns2, valphas2, vbetas2

        vns1 = s1%n_s
        valphas1 = s1%alpha_s
        vbetas1 = s1%beta_s

        vns2 = s2%n_s
        valphas2 = s2%alpha_s
        vbetas2 = s2%beta_s

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
        
        fullmel = Hs + Hsp1 + Hsm1

    end function mel

end module statemanip