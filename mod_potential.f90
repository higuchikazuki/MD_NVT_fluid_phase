module mod_potential

  use mod_rel_pos
  implicit none
contains

!レナードジョーンズポテンシャル(2粒子間)
  ! function uij(L,x1, x2)

  !   integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !   & wp = selected_real_kind(15)
  !   integer(i4b) :: N
  !   real(wp) ::  relative_pos, rel_pos(3) ,Upair(3),L, x1(3), x2(3),uij

  !   rel_pos=calc_relative_position(L, x1, x2)
  !   relative_pos=norm2(rel_pos)
  !   uij=4*((relative_pos**(-12))-(relative_pos**(-6)))

  ! end function uij


!レナードジョーンズポテンシャル(total)
  ! function total_potential(N,position,L)
  !   implicit none
  !   integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !   & wp = selected_real_kind(15)
  !   integer(i4b) :: N,i2,j2
  !   real(wp) :: position(3,N),  relative_pos, rel_pos(3), L, total_potential

  !   total_potential =0.0
  !   do i2=1,N
  !     do j2=i2+1,N
  !       total_potential=+total_potential+uij(L, position(:,j2),position(:,i2))
  !     enddo
  !   enddo
  ! end function total_potential



  ! function PotentialEnergy_iselect(L,N,position,j2,dammyposition)
  !   integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !   & wp = selected_real_kind(15)
  !   integer(i4b),intent(in):: N,j2
  !   integer(i4b)::i2
  !   real(wp),intent(in)::position(3,N),dammyposition(3),L
  !   real(wp)::PotentialEnergy_iselect
  !   PotentialEnergy_iselect = 0._wp
  !   do i2 = 1, N
  !     if(i2/=j2)then
  !       PotentialEnergy_iselect=PotentialEnergy_iselect+uij(L,dammyposition, position(:, i2))
  !     end if
  !   end do
  ! end function PotentialEnergy_iselect

  function u(theta,L,x1,x2)
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N
    real(wp):: theta(4), relative_pos1, rel_pos1(3) ,L, x1(3), x2(3),u
    rel_pos1=calc_relative_position(L, x1, x2)
    relative_pos1=norm2(rel_pos1)
    ! u = -theta(1)/(relative_pos1-1) + theta(2)*exp(-(relative_pos1-1)/theta(3))
    u = -theta(1)/(relative_pos1/theta(4) -1) + theta(2)*exp(-(relative_pos1/theta(4) -1)/theta(3))
  end function u




  !DLVOポテンシャル(2粒子間)
  function uij(theta, param, L, x1, x2)

    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N
    real(wp):: theta(4), relative_pos, rel_pos(3) ,Upair(3),L, x1(3), x2(3),uij,param(7)

    rel_pos=calc_relative_position(L, x1, x2)
    relative_pos=norm2(rel_pos)

    ! if (relative_pos > 1.06611160787638) then
    if (relative_pos > param(5)) then
      uij = -theta(1)/(relative_pos/theta(4) -1) + theta(2)*exp(-(relative_pos/theta(4) -1)/theta(3))
    else
      ! uij = -224.161864884258*(relative_pos-1.06611160787638) + 22.2505333
      uij = -param(6)*(relative_pos-param(5)) + param(7)
    end if

  end function uij

  ! function uij(theta, L, x1, x2)

  !   integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !   & wp = selected_real_kind(15)
  !   integer(i4b) :: N
  !   real(wp):: theta(3), relative_pos, rel_pos(3) ,Upair(3),L, x1(3), x2(3),uij

  !   rel_pos=calc_relative_position(L, x1, x2)
  !   relative_pos=norm2(rel_pos)
  !   if (relative_pos > 1.03391542648754) then
  !     uij = -theta(1)/(relative_pos-1) + theta(2)*exp(-(relative_pos-1)/theta(3))
  !   else
  !     uij = 27.7040073533468
  !   end if

  ! end function uij


  function total_potential(theta,param,N,position,L)
    implicit none
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N,i2,j2
    real(wp) :: position(3,N),  relative_pos, rel_pos(3), L, total_potential, theta(4), param(7)

    total_potential =0.0
    do i2=1,N
      do j2=i2+1,N
        total_potential=+total_potential+uij(theta,param, L, position(:,j2),position(:,i2))
        ! write(14,*) i2, j2, uij(theta, L, position(:,j2),position(:,i2)), total_potential
        ! write(20,*) 'j2',j2 ,position(:,j2)
        ! write(20,*) 'i2',i2, position(:,i2)
      enddo
    enddo
  end function total_potential

  function elastic_energy(N,L,position, position_Ein, spring_constant)
    implicit none
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N,i5
    real(wp) :: position_Ein(3,N), position(3,N), displacement(3,N), L, elastic_energy, spring_constant

    elastic_energy = 0.0
    do i5 =1,N
      displacement(:,i5) = position(:,i5) - position_Ein(:,i5)
      where (displacement > L / 2._dp)
        displacement = displacement - L
      elsewhere (displacement <= -L / 2._dp)
        displacement = displacement + L
      end where
      elastic_energy = elastic_energy + spring_constant*dot_product(displacement(:,i5),displacement(:,i5))
    end do
  end function elastic_energy


  function free_energy_density(fitting_param_p, kbT0, rho)
    implicit none
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    real(wp), intent(in) :: fitting_param_p(8), kbT0, rho
    integer(i4b) :: i
    real(wp) :: free_energy_density

    free_energy_density = kbT0*rho*(log(rho)-1)
    do i = 1, 8
      free_energy_density = free_energy_density + (fitting_param_p(i)/i)*(rho**(i+1))
    end do
  end function free_energy_density


  function total_potential_Ein(theta,param,N,position_Ein,L)
    implicit none
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N,i2,j2
    real(wp) :: position_Ein(3,N),  relative_pos, rel_pos(3), L, total_potential_Ein, theta(4), param(7)

    total_potential_Ein =0.0
    do i2=1,N
      do j2=i2+1,N
        total_potential_Ein=+total_potential_Ein+uij(theta,param, L, position_Ein(:,j2),position_Ein(:,i2))
      enddo
    enddo
  end function total_potential_Ein


  function switching_potential(theta, param, N, position, lambda, position_Ein, L, spring_constant)
    implicit none
    integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
    & wp = selected_real_kind(15)
    integer(i4b) :: N,i2,j2
    real(wp) :: position_Ein(3,N),  relative_pos, rel_pos(3), L, theta(4), switching_potential, lambda, position(3,N),&
    & spring_constant, param(7)

    switching_potential = total_potential_Ein(theta,param,N,position_Ein,L) + (1-lambda)*(total_potential(theta,param,N,position,L) &
    & - total_potential_Ein(theta,param,N,position_Ein,L) ) + lambda*elastic_energy(N,L,position,position_Ein,spring_constant)
  end function switching_potential



end module mod_potential
