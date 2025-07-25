module mod_DLVOpotential
  use, intrinsic :: iso_fortran_env, only : wp => real64, int64
  implicit none
  real(wp), parameter :: nan = transfer(-1_int64, 0._wp)
  private
  public  g_of_x, get_vx_maximal_minimal_position, get_vx_maximal_minimal_value, &
     & get_DLVO_potential_maxmin, get_DLVO_force_maximal!, v_DLVO, deriv_v_DLVO

contains
  !パラメータ4つを持つDLVO ポテンシャル
  ! pure function v_DLVO(r, theta)
  !   real(wp), intent(in) :: r, theta(3)
  !   real(wp) :: v_DLVO
  !   if (r > theta(4)) then
  !     v_DLVO = - theta(1) / (r / theta(4) - 1) + theta(2) * exp( - (r / theta(4) - 1) / theta(3))
  !   else
  !     v_DLVO = nan
  !   end if
  ! end function v_DLVO
  
  ! !パラメータ4つを持つDLVO ポテンシャルの微分
  ! pure function deriv_v_DLVO(r, theta)
  !   real(wp), intent(in) :: r, theta(4)
  !   real(wp) :: deriv_v_DLVO
  !   if (r > theta(4)) then
  !     deriv_v_DLVO = theta(1) / (theta(4) * (r / theta(4) - 1)**2) - theta(2) * exp(- (r / theta(4) - 1) / theta(3)) / (theta(3) * theta(4))
  !   else
  !     deriv_v_DLVO = nan
  !   end if
  ! end function deriv_v_DLVO
  
  ! DLVO ポテンシャルのピークを求める際に利用する関数。
  ! g_of_x(x) = c の解を求めることで、DLVOポテンシャルの極値の位置を求めることができる。
  ! 右辺の定数 c は DLVO ポテンシャルのパラメータによって決まる定数で、function_c で求められる。
  ! g_of_x(x) = d の解を求めることで、DLVO力の極値の位置を求めることができる。
  ! 右辺の定数 d は DLVO ポテンシャルのパラメータによって決まる定数で、function_d で求められる。
  pure function g_of_x(x)
    real(wp), intent(in) :: x
    real(wp) :: g_of_x
    g_of_x = x * exp(-x)
  end function g_of_x

  ! g_of_x の微分。ニュートンラフソン法に用いる。
  pure function deriv_g_of_x(x)
    real(wp), intent(in) :: x
    real(wp) :: deriv_g_of_x
    deriv_g_of_x = - (x - 1) * exp(-x)
  end function deriv_g_of_x

  ! DLVOポテンシャルの極値の位置の計算に用いる定数を与える関数
  ! g_of_x のコメント参照
  pure function function_c(theta)
    real(wp), intent(in) :: theta(4)
    real(wp) :: function_c
    function_c = (theta(1) / (4 * theta(2) * theta(3)))**0.5_wp
  end function function_c

  ! DLVO力の極値の位置の計算に用いる定数を与える関数
  ! g_of_x のコメント参照
  pure function function_d(theta)
    real(wp), intent(in) :: theta(4)
    real(wp) :: function_d
    function_d = (2 * theta(1) / (27 * theta(2) * theta(3)))**(1._wp / 3)
  end function function_d

  ! DLVOポテンシャルの極大値および極小値を求めるための関数。別文書を参照。
  pure function function_u(x, c)
    real(wp), intent(in) :: x, c
    real(wp) :: function_u
    function_u = (1 - 2 * x) * (c / x)**2
  end function function_u

  ! ニュートンラフソン法で g_of_x(x) = const の１つの解の位置を求める関数
  pure function Newton_Raphson_for_g(const, xinit, threshold) result(solution)
    real(wp), intent(in) :: const, xinit, threshold
    real(wp) :: solution, x, delta_g
    integer :: i1
    x = xinit
    do i1 = 1, 1000
      delta_g = g_of_x(x) - const 
      if (abs(delta_g) < threshold) then
        solution = x
        exit
      else
        x = x - delta_g / deriv_g_of_x(x)
      end if
    end do
    if (i1 >= 1000) solution = nan
  end function Newton_Raphson_for_g

  ! ニュートンラフソン法でDLVOポテンシャルの極大の位置およびその値を求める関数
  ! 戻り値の(1, 1)成分および(2, 1)成分は、それぞれ極大の位置および極大値
  ! 戻り値の(1, 2)成分および(2, 2)成分は、それぞれ極小の位置および極小値
  pure function get_DLVO_potential_maxmin(theta) result(maxmin)
    real(wp), intent(in) :: theta(4)
    real(wp) :: maxmin(2, 2), xMin, xMax, const, threshold
    threshold = 1.e-8_wp
    const = function_c(theta)
    xMax = Newton_Raphson_for_g(const, 0._wp, threshold)
    maxmin(1, 1) = (2 * theta(3) * xMax + 1) * theta(4)
    maxmin(2, 1) = theta(2) * function_u(xMax, const)
    xMin = Newton_Raphson_for_g(const, 2._wp, threshold)
    maxmin(1, 2) = (2 * theta(3) * xMin + 1) * theta(4)
    maxmin(2, 2) = theta(2) * function_u(xMin, const)
  end function get_DLVO_potential_maxmin
  
  ! ニュートンラフソン法でDLVO力の値の極大の位置およびその値を求める関数
  ! 戻り値の第1成分および2成分は、それぞれ極大の位置および極大値
  pure function get_DLVO_force_maximal(theta) result(maximal)
    real(wp), intent(in) :: theta(4)
    real(wp) :: yMax, maximal(2), const, threshold
    threshold = 1.e-8_wp
    const = function_d(theta)
    yMax = Newton_Raphson_for_g(const, 0._wp, threshold)
    maximal(1) = (3 * theta(3) * yMax + 1) * theta(4)
    maximal(2) = (theta(2) / (theta(3) * theta(4))) * ((const / yMax)**3) * (1 - 1.5_wp * yMax)
    ! maximal(2) = (theta(2) / theta(3) ) * ((const / yMax)**3) * (1 - 1.5_wp * yMax)
  end function get_DLVO_force_maximal
  
  ! ニュートンラフソン法で g_of_x(x) = const の２つの解の位置を求める関数
  pure function get_vx_maximal_minimal_position(const) result(maxmin)
    real(wp), intent(in) :: const
    real(wp) :: maxmin(2), threshold
    threshold = 1.e-8_wp
    maxmin(1) = Newton_Raphson_for_g(const, 0._wp, threshold)
    maxmin(2) = Newton_Raphson_for_g(const, 2._wp, threshold)
  end function get_vx_maximal_minimal_position

  ! g_of_x(x) = const の２つの解の位置 maxmin_position から、DLVO ポテンシャルの極大および極小の値を求める関数
  pure function get_vx_maximal_minimal_value(const, theta2, maxmin_position) result(maxminval)
    real(wp), intent(in) :: const, theta2, maxmin_position(2)
    integer :: i1
    real(wp) :: maxminval(2)
    do i1 = 1, 2
      maxminval(i1) = theta2 * function_u(maxmin_position(i1), const)
    end do
  end function get_vx_maximal_minimal_value
end module mod_DLVOpotential
