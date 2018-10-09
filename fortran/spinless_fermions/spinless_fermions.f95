program spinless_fermions
    
    use statemanip
    implicit none

    integer :: num_sites = 6
    integer :: num_particles = 3
    integer :: level = 4

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! call splitstates(s1)
    ! write (*,*) "********************************************"
    ! call splitstates(s2)
    ! write (*,*) "********************************************"
    ! call splitstates(s3)
    ! write (*,*) "********************************************"
    ! call splitstates(s4)
    ! write (*,*) "********************************************"

    ! call splitstatesatlevel(s1, level)
    ! write (*,*) "********************************************"
    ! call splitstatesatlevel(s2, level)
    ! write (*,*) "********************************************"
    ! call splitstatesatlevel(s3, level)
    ! write (*,*) "********************************************"
    ! call splitstatesatlevel(s4, level)
    ! write (*,*) "********************************************"

    call findmels(s1, s1)
    write (*,*) "********************************************"
    call findmels(s1, s2)
    write (*,*) "********************************************"
    call findmels(s2, s2)
    write (*,*) "********************************************"
    call findmels(s2, s3)
    write (*,*) "********************************************"
    call findmels(s3, s3)
    write (*,*) "********************************************"
    call findmels(s3, s4)
    write (*,*) "********************************************"
    call findmels(s4, s4)
    ! write (*,*) "********************************************"

    ! call findmelsatlevel(s1, s1, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s1, s2, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s2, s2, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s2, s3, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s3, s3, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s3, s4, level)
    ! write (*,*) "********************************************"
    ! call findmelsatlevel(s4, s4, level)
    ! write (*,*) "********************************************"

end program spinless_fermions