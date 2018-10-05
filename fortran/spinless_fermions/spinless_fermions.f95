program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites
    integer :: num_particles
    integer :: smax
    integer :: ns

    integer :: sigma, ns1, ns2, alpha1, alpha2, beta1, beta2

    type (state) :: s1
    type (state) :: s2
    type (state) :: s3
    type (state) :: s4
    type (state) :: ss1
    type (state) :: ss2
    type (state) :: ss3

    smax = num_sites - 1

    ns = num_particles

    call setstate(s1, 3, 2, 0, 0)
    call setstate(s2, 3, 2, 1, 0)
    call setstate(s3, 3, 1, 0, 1)
    call setstate(s4, 3, 1, 1, 1)
    ! print *, getlstate(s1)
    ! print *, getlstate(s2)
    ! ss1 = relegate0(s1)
    ! print *, getlstate(ss1)
    ! ss2 = relegate1(s1)
    ! print *, getlstate(ss2)
    ! ss3 = relegate0(ss1)
    ! print *, getlstate(ss3)

    ! print *, checkvalidity(s1)

    call findmels(s2, s3)

end program spinless_fermions