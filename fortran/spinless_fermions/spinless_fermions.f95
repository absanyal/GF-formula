program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites = 4
    integer :: num_particles = 2

    integer :: smax
    integer :: ns

    type (state) :: s1
    type (state) :: s2
    type (state) :: s3
    type (state) :: s4

    smax = num_sites - 1

    ns = num_particles

    call setstate(s1, smax, ns, 0, 0)
    call setstate(s2, smax, ns, 1, 0)
    call setstate(s3, smax, ns-1, 0, 1)
    call setstate(s4, smax, ns-1, 1, 1)
    ! print *, getlstate(s1)
    ! print *, getlstate(s2)
    ! ss1 = relegate0(s1)
    ! print *, getlstate(ss1)
    ! ss2 = relegate1(s1)
    ! print *, getlstate(ss2)
    ! ss3 = relegate0(ss1)
    ! print *, getlstate(ss3)

    ! print *, checkvalidity(s1)

    call findmels(s1, s1)
    write (*,*) "**************************"
    call findmels(s1, s2)
    write (*,*) "**************************"
    call findmels(s2, s2)
    write (*,*) "**************************"
    call findmels(s2, s3)
    write (*,*) "**************************"
    call findmels(s3, s3)
    write (*,*) "**************************"
    call findmels(s3, s4)
    write (*,*) "**************************"
    call findmels(s4, s4)
    write (*,*) "**************************"

end program spinless_fermions