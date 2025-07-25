module mod_force

  use mod_rel_pos
  implicit none
contains

  function Fpair(theta,param,L,x1, x2)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    real(wp), intent(in) :: theta(4), L, x1(3), x2(3),param(7)
    real(wp) :: durkl, distance, rel_pos(3), Fpair(3)
    rel_pos=calc_relative_position(L, x1, x2)
    distance = sqrt(dot_product(rel_pos, rel_pos))
    if (distance > param(5)) then
      Fpair =  ((theta(1)/(theta(4)*(distance-1)**2)-(theta(2)/theta(3)*theta(4))*exp(-(distance/theta(4) -1)/theta(3)))*(-rel_pos/distance))
    else
      Fpair = param(6)*(rel_pos/distance)
    endif
  end function Fpair


  function Fpair_Ein(theta,param,lambda,L,x1,x2)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    real(wp), intent(in) :: theta(4),param(7), lambda, L, x1(3), x2(3)
    real(wp) :: durkl, distance, rel_pos(3), Fpair_Ein(3)
    rel_pos=calc_relative_position(L, x1, x2)
    distance = sqrt(dot_product(rel_pos, rel_pos))
    ! if (distance > 1.06611160787638) then
    if (distance > param(5)) then
      !Fpair_Ein =  (1-lambda)*((theta(1)/((distance-1)**2)-(theta(2)/theta(3))*exp(-(distance-1)/theta(3)))*(-rel_pos/distance)) !&
      !& + 2*lambda*spring_constant*distance*(-rel_pos_ein/distance_ein)
      Fpair_Ein =  (1-lambda)*((theta(1)/(theta(4)*(distance-1)**2)-(theta(2)/theta(3)*theta(4))*exp(-(distance/theta(4) -1)/theta(3)))*(-rel_pos/distance))
    else
      Fpair_Ein = param(6)*(1-lambda)*(rel_pos/distance)
    endif
  end function Fpair_Ein

  function Elastic_Fpair(lambda,L,x1,x2,x3)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    real(wp), intent(in) :: L, x1(3), x2(3), x3
    real(wp) :: distance_ein, rel_pos_ein(3), Elastic_Fpair(3),lambda
    rel_pos_ein=calc_relative_position(L, x1, x2)
    distance_ein = sqrt(dot_product(rel_pos_ein, rel_pos_ein))
    Elastic_Fpair = 2*lambda*x3*distance_ein*(-rel_pos_ein/distance_ein)
  end function Elastic_Fpair

  !弾性ばねの弾性力,
  function Elastic_Force(lambda,L,position,position_Ein,spring_constant,N)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    real(wp), intent(in) :: L, position(3,n), position_Ein(3,n), spring_constant
    real(wp) :: distance_ein, rel_pos_ein(3,N), Elastic_Force(3,n),lambda
    integer, intent(in) :: N
    rel_pos_ein=calc_relative_position_Ein(L, position, position_Ein, n)
    Elastic_Force = 2*lambda*spring_constant*(-rel_pos_ein)
  end function Elastic_Force

  function force(theta, param, L, N, position)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, position(3, N), param(7)
    integer(wp) :: i, j
    real(wp) :: force(3,N), Fij(3), theta(4)
    force = 0
    do j = 1, N
      do i = j + 1, N
        Fij = Fpair(theta, param, L, position(:, j),position(:, i) )
        force(:, i) = force(:, i) + Fij
        force(:, j) = force(:, j) - Fij
      enddo
    enddo
  end function force

  function total_force(theta, param, L, N, lambda, position)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, position(3, N), param(7)
    integer(wp) :: i, j
    real(wp) :: total_force(3,N), Fij_Ein(3), theta(4), lambda, spring_constant
    total_force = 0
    do j = 1, N
      do i = j + 1, N
        Fij_Ein = Fpair_Ein(theta, param, lambda, L, position(:, j),position(:, i) )
        total_force(:, i) = total_force(:, i) + Fij_Ein
        total_force(:, j) = total_force(:, j) - Fij_Ein
      enddo
    enddo
    ! do i = 1,N
    !   Fij_Ein = Elastic_Fpair(lambda, L, position_Ein(:,i), position(:,i), spring_constant)
    !   total_force(:,i) = total_force(:,i) + Fij_Ein
    ! enddo
  end function total_force

  function virial_internal(theta, param, N, L, position)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, position(3, N)
    integer(i4b) :: i, j
    real(wp) :: virial_internal, rij(3), Fij(3),theta(4), param(7)
    virial_internal = 0
    do j = 1, N
      do i = j + 1, N
        rij = calc_relative_position(L, position(:, i),position(:, j))
        Fij = Fpair(theta, param, L, position(:, i),position(:, j) )
        virial_internal = virial_internal + dot_product(rij, Fij)
      enddo
    enddo
  end function virial_internal

  function instantaneous_P(theta, param, N, L, V, position, momentum)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, V, position(3, N), momentum(3, N), theta(4), param(7)
    real(wp) :: instantaneous_P
    instantaneous_P = (sum(momentum**2) + virial_internal(theta, param, N, L, position))/(3 * V)
  end function instantaneous_P

  function virial_internal_Ein(theta,lambda,param, N, L, position,position_Ein,spring_constant)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, position(3, N),position_Ein(3,N), param(7)
    integer(i4b) :: i, j
    real(wp) :: virial_internal_Ein, rij(3), Fij_Ein(3),theta(4),lambda,spring_constant
    virial_internal_Ein = 0
    do j = 1, N
      do i = j + 1, N
        rij = calc_relative_position(L, position(:, j),position(:, i))
        Fij_Ein = Fpair_Ein(theta,param,lambda,L, position(:, j),position(:, i) ) + Elastic_Fpair(lambda,L,position_Ein(:,j),&
          & position(:,j),spring_constant)
        virial_internal_Ein = virial_internal_Ein + dot_product(rij, Fij_Ein)
      enddo
    enddo
  end function virial_internal_Ein

  function instantaneous_P_Ein(theta,lambda,param, N, L, V, position,position_Ein,spring_constant, momentum)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b), intent(in) :: N
    real(wp), intent(in) :: L, V, position(3, N), momentum(3, N), param(7), theta(4), position_Ein(3,N),spring_constant,lambda
    real(wp) :: instantaneous_P_Ein
    instantaneous_P_Ein = (sum(momentum**2) + virial_internal_Ein(theta,lambda,param, N, L, position,position_Ein,&
    & spring_constant))/(3 * V)
  end function instantaneous_P_Ein


  ! ���q i ���󂯂�͂̑��a��߂��֐��B
  ! �߂�l�� (3, N) �z��ŁA���q i ���󂯂�͂̑��a��(:, i)�����B
  ! function force(N,position,L)
  !   implicit none
  !   integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !   & wp = selected_real_kind(15)
  !   integer(i4b) :: N,i1,i2
  !   real(wp) :: position(3,N),force(3,N), durkl, e,s, relative, rel_pos(3), &
  !   & Fpair_tmp(3),L
  !   e=1.0
  !   s=1.0
  !   force =0
  !   do i2=2,N
  !     do i1=1,i2-1
  !       Fpair_tmp=Fpair(L, position(:,i2),position(:,i1))
  !       force(:,i2)=force(:,i2)-Fpair_tmp
  !       force(:,i1)=force(:,i1)+Fpair_tmp
  !     enddo
  !   enddo
  ! end function force

end module mod_force
