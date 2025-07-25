    module mod_create_fcc
    use constants
    implicit none
    real(wp), dimension(3,4)::fcc_coordinates = reshape([0.0_wp, 0.0_wp, 0.0_wp,  &
        0.5_wp, 0.5_wp, 0.0_wp,  &
        0.5_wp, 0.0_wp, 0.5_wp,  &
        0.0_wp, 0.5_wp, 0.5_wp], [3,4])

    contains
    !fcc‚Ì’PˆÊ–E‚ğM–M–MŒÂ‚¾‚¯•À‚×‚½À•W‚ğ‹‚ß‚éŠÖ”
    function get_fcc_lattice_coordinates_cubic(M) result(coordinates)
    integer(i4b),intent(in)::M
    real(wp),dimension(3,4*M**3) :: coordinates

    integer(i4b) :: i, j, k, n

    n = 1
    do i = 1,M
        do j = 1,M
            do k = 1,M
                coordinates(:,n:n+3) = fcc_coordinates + spread([i-1,j-1,k-1], dim=2, ncopies=4)
                n = n + 4
            end do
        end do
    end do

    end function get_fcc_lattice_coordinates_cubic

    function get_fcc_coordinates(N)
    integer(i4b),intent(in)::N
    real(wp),dimension(3,N):: get_fcc_coordinates
    integer(i4b)::M
    real(wp),allocatable :: coordinates(:,:)

    M=int((N/4._wp)**(1/3._wp))+1

    allocate(coordinates(3,4*M**3))
    coordinates=get_fcc_lattice_coordinates_cubic(M)
    get_fcc_coordinates= coordinates(:,1:N)
    end  function get_fcc_coordinates

    end module mod_create_fcc