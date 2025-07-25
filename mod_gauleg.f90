module gaussLegendre
  implicit none
  integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), wp = selected_real_kind(15)
  integer(i4b), parameter :: npar_arth=16, npar2_arth=8
  real(wp), parameter ::  pi = 4*atan(1._wp)
  private
  public :: gauleg

contains
  subroutine gauleg(N1, lambda1, lambda2, x, w, err)
    implicit none
    integer(i4b), intent(in) :: N1
    real(wp), intent(in) :: lambda1, lambda2
    real(wp), intent(out) :: x(N1), w(N1)
    integer(i4b), intent(out) :: err
    real(wp), parameter :: eps = 3.0e-14_wp
    integer(i4b), parameter :: maxit = 1000
    integer(i4b) :: its, j, m
    real(wp) :: xl, xm
    real(wp) :: p1((N1 + 1)/2), p2((N1 + 1)/2), p3((N1 + 1)/2), pp((N1 + 1)/2), z((N1 + 1)/2), z1((N1 + 1)/2)
    logical(kind(.true.)) :: unfinished((N1 + 1)/2)
    err = 0
    m = (N1 + 1)/2
    xm = 0.5_wp*(lambda2 + lambda1)
    xl = 0.5_wp*(lambda2 - lambda1)
    z = cos(pi*(arth(1, 1, m) - 0.25_wp)/(N1 + 0.5_wp))
    unfinished = .true.
    do its = 1, maxit
      where (unfinished)
        p1 = 1.0
        p2 = 0.0
      end where
      do j = 1, N1
        where (unfinished)
          p3 = p2
          p2 = p1
          p1 = ((2.0_wp*j-1.0_wp)*z*p2-(j-1.0_wp)*p3)/j
        end where
      end do
      where (unfinished)
        pp = N1*(z*p1-p2)/(z*z-1.0_wp)
        z1 = z
        z = z1-p1/pp
        unfinished = (abs(z-z1) > eps)
      end where
      if (.not. any(unfinished)) exit
    end do
    if (its == maxit+1) then
      err = 1
      return
    end if
    x(1:m) = xm-xl*z
    x(N1:N1-m+1:-1) = xm+xl*z
    w(1:m) = 2.0_wp*xl/((1.0_wp-z**2)*pp**2)
    w(N1:N1-m+1:-1) = w(1:m)
  end subroutine gauleg

  pure function arth(first, increment, N1)
    implicit none
    integer(i4b), intent(in) :: first, increment, N1
    integer(i4b) :: arth(N1)
    integer(i4b) :: k, k2, temp
    if (N1 > 0) arth(1) = first
    if (N1 <=  npar_arth) then
      do k = 2, N1
        arth(k) = arth(k-1)+increment
      end do
    else
      do k = 2, npar2_arth
        arth(k) = arth(k-1) + increment
      end do
      temp = increment*npar2_arth
      k = npar2_arth
      do
        if (k >=  N1) exit
        k2 = 2*k
        arth(k+1:min(k2, N1)) = temp + arth(1:min(k, N1-k))
        temp = 2*temp
        k = k2
      end do
    end if
  end function arth
end module gaussLegendre
