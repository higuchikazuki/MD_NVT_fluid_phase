module mod_time_evolution
  use mod_force,only:force
  use mod_rel_pos
  implicit none
  integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  & wp = selected_real_kind(15)
  private
  public :: RDF
contains

  !real(wp)::position(3,N)
  !position=get_fcc_coordinates(N)*sqrt(2._wp)


  !動径分布関数を求める関数
  function RDF(L, position, delta_r, N, V, ndiv)
    integer(i4b):: N, i1, i2, ir, ndiv
    real(wp) :: RDF(ndiv), ng(ndiv), r, rel_pos(3), delta_r, position(3, N), L, V, &
    & r_max, relative_r(3)
    real(wp), parameter :: pi = 4*atan(1._wp)
    ng = 0
    rel_pos = 0
    r_max = delta_r*ndiv
    do i1 =1, N-1
      do i2 = i1+1, N
        rel_pos =calc_relative_position(L, position(:, i1), position(:, i2))
        r=norm2(rel_pos)
        ir = int(r/delta_r)+1
        if(ir<=ndiv) ng(ir) = ng(ir)+1
      end do
    end do
    do i1 = 1, ndiv
      r = i1*delta_r
      RDF(i1) = ng(i1)/(r**2+r*delta_r+(delta_r**2/3))
      !write(*,*) L, position, deltar, N, V, ndiv
    end do
    RDF = RDF/((4*pi)*delta_r*(N/V)*(N*0.5))
  end function RDF





  ! function time_evolution1(position,momentum_p,dt,s,L)
  ! integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !     & wp = selected_real_kind(15)
  ! integer(i4b),intent(in):: N
  ! real(wp),intent(in)::position(3,N),momentum_p(3,N),dt,s,L
  ! real(wp)::time_evolution1(3,N)
  ! time_evolution1=momentum_p+s*totalforce(L,N,position)*dt/2!(10.28)
  ! end function time_evolution1

  ! function time_evolution2(s_momentum,dt,Q)
  ! integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !     & wp = selected_real_kind(15)
  ! real(wp),intent(in)::s_momentum,dt,Q
  ! real(wp)::time_evolution2
  ! time_evolution2=s_momentum/(1+(s_momentum*dt)/(4*Q))!(10.27)
  ! end function time_evolution2

  ! function time_evolution3(s_momentum,s,dt,Q)
  ! integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  !     & wp = selected_real_kind(15)
  ! real(wp),intent(in)::s_momentum,s,dt,Q
  ! real(wp)::time_evolution3
  ! time_evolution3=s*(1+(s_momentum*dt)/(4*Q))**2!(10.26)
  ! end function time_evolution3

end module mod_time_evolution
