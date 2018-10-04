program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites
    integer :: num_particles
    integer :: smax
    integer :: ns

    integer :: sigma, ns1, ns2, alpha1, alpha2, beta1, beta2

    type (state) :: s1
    type (state) :: ss1
    type (state) :: ss2
    type (state) :: ss3

    smax = num_sites - 1

    ns = num_particles

     call setstate(s1, 3, 2, 1, 0)
     print *, getstate(s1)
     ss1 = relegate0(s1)
     print *, getstate(ss1)
     ss2 = relegate1(s1)
     print *, getstate(ss2)
     ss3 = relegate0(ss1)
     print *, getstate(ss3)

end program spinless_fermions