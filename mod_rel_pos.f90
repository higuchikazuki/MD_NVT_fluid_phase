    module mod_rel_pos
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)
    contains

    function calc_relative_position(L, x1, x2) result(rel_pos)
    real(dp), intent(in) :: L
    real(dp), dimension(3), intent(in) :: x1, x2
    real(dp), dimension(3) :: rel_pos

    rel_pos = x2 - x1
    where (rel_pos > L / 2._dp)
        rel_pos = rel_pos - L
    elsewhere (rel_pos <= -L / 2._dp)
        rel_pos = rel_pos + L
    end where
    end function calc_relative_position

    function calc_relative_position_Ein(L, position, position_Ein, n) result(rel_pos)
      real(dp), intent(in) :: L
      real(dp), dimension(3,n), intent(in) :: position, position_Ein
      real(dp), dimension(3,N) :: rel_pos
      integer, intent(in) :: N
  
      rel_pos = position - position_Ein
      where (rel_pos > L / 2._dp)
          rel_pos = rel_pos - L
      elsewhere (rel_pos <= -L / 2._dp)
          rel_pos = rel_pos + L
      end where
      end function calc_relative_position_Ein

    end module mod_rel_pos