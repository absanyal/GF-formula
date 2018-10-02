program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites
    integer :: num_particles
    integer :: sigma

    type (state) :: s1
    type (state) :: s2

    call setstate(s1, 3, 2, 1, 0)
    call setstate(s2, 3, 1, 0, 1)

    print *, mel(s1, s2)

end program spinless_fermions