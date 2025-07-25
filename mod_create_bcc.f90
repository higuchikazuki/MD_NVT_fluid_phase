module mod_create_bcc
  use constants
  implicit none
  ! BCC単位格子：角点（0,0,0）と体心点（0.5,0.5,0.5）
  real(wp), dimension(3,2) :: bcc_coordinates = reshape([ &
       0.0_wp, 0.0_wp, 0.0_wp,  &
       0.5_wp, 0.5_wp, 0.5_wp], [3,2])
  
contains

  ! 立方体状に M×M×M 個のセルを生成し、
  ! 各セルにBCC点2個（角点と体心点）を配置する関数
  function get_bcc_lattice_coordinates_cubic(M) result(coordinates)
      integer(i4b), intent(in) :: M
      real(wp), dimension(3,2*M**3) :: coordinates
      integer(i4b) :: i, j, k, n

      n = 1
      do i = 1, M
          do j = 1, M
              do k = 1, M
                  coordinates(:, n:n+1) = bcc_coordinates + &
                      spread([i-1, j-1, k-1], dim = 2, ncopies = 2)
                  n = n + 2
              end do
          end do
      end do

  end function get_bcc_lattice_coordinates_cubic

  ! N個のBCC格子点を得るための関数
  function get_bcc_coordinates(N)
      integer(i4b), intent(in) :: N
      real(wp), dimension(3,N) :: get_bcc_coordinates
      integer(i4b) :: M
      real(wp), allocatable :: coordinates(:,:)

      ! 各セルに2点があるので、必要なセル数は int((N/2)**(1/3)) に1を足す
      M = int((N/2._wp)**(1/3._wp)) + 1

      allocate(coordinates(3,2*M**3))
      coordinates = get_bcc_lattice_coordinates_cubic(M)
      get_bcc_coordinates = coordinates(:,1:N)
  end function get_bcc_coordinates

end module mod_create_bcc