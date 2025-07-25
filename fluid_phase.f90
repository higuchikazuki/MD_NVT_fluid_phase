program fluid_phase
  use mod_create_fcc,only:get_fcc_coordinates
  use mod_rel_pos
  use mod_force
  use mod_potential
  use mod_DLVOpotential
  use gaussLegendre,only:gauleg
  use mod_time_evolution, only:RDF
  implicit none
  integer(selected_int_kind(9)), parameter :: i4b = selected_int_kind(9), &
  & wp = selected_real_kind(15)
  integer(i4b), parameter :: N=500, N1=35
  integer(i4b)::i,i3,j3,j,ndiv,i5,j5,i2,seedsize,i4,err
  integer,allocatable:: seed(:)
  real(wp) :: position(3,N), momentum(3,N), zeta, dt, Q, g, kbT0,kbT, t, L, av_momentum(3),rho,c,V,&
  & position_Ein(3,N), theta(4), spring_constant, square_displacement_ave(N), displacement(3,N),delta_r,gr_dash(100),&
  & lambda1,lambda2,x(N1),w(N1),f(N1), lambda, free_energy, free_energy_, spring_constant_array(N), pressure,&
  & sigma_ln, delta_total_potential_ave, delta_total_potential, integration_potential, delta_lambda,swiching_potential,&
  & elastic_energy_ave,virial_int,virial_int_EIn,param(7), momentum0(3,N), relative_pos, data_potential_force(3,2),&
  & free_energy_fluid, free_energy_ideal, pressure_ave, delta_rho, integrand, L_before, fitting_param_p(8), delta_kbT0,&
  & potential_average
  real(wp), parameter :: pi = 4*atan(1._wp)
  logical, parameter :: test_run = .false.

  open(11,file='data.dat')
  open(12,file='position.dat')
  open(13,file='momentum.dat')
  open(14,file='position_Ein.dat')
  open(15,file='RDF.dat')
  open(16,file='spring_constant.dat')
  open(17,file='data2.dat')
  open(18,file='data_free_energy.dat')
  open(19,file='data3.dat')
  open(20,file='RDF_first.dat')
  open(21,file='force.dat')
  open(22,file='force_first.dat')
  open(23,file='paramater.dat')
  open(24,file='DLVO_potential.dat')
  open(25,file='integrand.dat')
  open(26,file='data_DLVO_ptential_maxmin.dat')
  open(27,file='data_DLVO_force_max.dat')
  open(28,file='position_first.dat')

  rho=0.05_wp
  V = N/rho
  L = V**(1._wp/3)
  Q=1.0
  dt=0.001
  g=3*N
  kbT0=5.0
  ndiv=100
  delta_r=0.0375
  theta = (/0.270833, 203.2973301, 0.0735297, 1.00/)
  param(1:4) = theta(1:4)
  param(5:6) = get_DLVO_force_maximal(theta)
  param(7) = u(theta,L, [0._wp , 0._wp, 0._wp] , [param(5),  0._wp, 0._wp])
  fitting_param_p = (/-2433.43, -37028.6, 83367.8, -71102.9, 31893.7, -7731.52, 1181.56, -80.1666/)
  data_potential_force(1:2,1:2) = get_DLVO_potential_maxmin(theta)
  data_potential_force(3,1:2) = get_DLVO_force_maximal(theta)
  free_energy_ideal = kbT0*rho*(log(rho)-1)
  free_energy = free_energy_ideal
  position_Ein = get_fcc_coordinates(N)*sqrt(2._wp)*1.836
  position = position_Ein
  pressure_ave=0.0
  delta_rho=0.025_wp
  delta_kbT0=0.01_wp

  write(23,'(7g18.8)') param
  write(26, *) get_DLVO_potential_maxmin(theta)
  write(27, *) get_DLVO_force_maximal(theta)

  do i5 = 1, 1000
    relative_pos = param(4) + (i5*0.001_wp)!*1.0E-7
    write(24, *) relative_pos/param(4) , u(theta, L, [0._wp , 0._wp, 0._wp] , [relative_pos,  0._wp, 0._wp])
  end do
  do i5 = 1, 120
    relative_pos = param(4) + (1.0 + i5*0.05_wp)!*1.0E-7
    write(24, *) relative_pos/param(4) , u(theta, L, [0._wp , 0._wp, 0._wp] , [relative_pos,  0._wp, 0._wp])
  end do
  write(24, '(/)')
  do i5 = 1, 1000
    relative_pos = param(4) + (i5*0.001_wp)!*1.0E-7
    write(24, *) relative_pos/param(4) , uij(theta,param,L, [0._wp , 0._wp, 0._wp] , [relative_pos,  0._wp, 0._wp])
  end do
  do i5 = 1, 120
    relative_pos = param(4) + (1.0 + i5*0.05_wp)!*1.0E-7
    write(24, *) relative_pos/param(4) , uij(theta,param,L, [0._wp , 0._wp, 0._wp] , [relative_pos,  0._wp, 0._wp])
  end do


  call RANDOM_SEED(size=seedsize)
  allocate(seed(seedsize))
  if (test_run) then
    write (*, *) 'TEST RUN'
    ! �v�Z�`�F�b�N���̎��͗����V�[�h���Œ肷��
    seed(1) = 3
    do i = 2, seedsize
      seed(i) = seed(i - 1) * (- 5)
    end do
    call random_seed(put = seed)
  else
    do i = 1, seedsize
      call system_clock(count=seed(i)) !���Ԃ����擾
    end do
    call RANDOM_SEED(put=seed)
  end if
  write(14,*) seedsize
  write(14,*) seed
  call RANDOM_NUMBER(momentum0)
  momentum0 = momentum0*2-1
  av_momentum=sum(momentum0,dim=2)/N
  do i=1,N
    momentum0(:,i)=momentum0(:,i)-av_momentum
  enddo
  kbT = (1/g)*(sum(momentum0**2))
  momentum = momentum0
  ! c=(0.5/kbT)**0.5
  ! momentum0=c*momentum0

  write(14,*) position_Ein

  zeta=0.5
  t=0.0

  do i =0,25000
    momentum = momentum*exp(-zeta*dt/2)
    momentum = momentum+force(theta,param,L,N,position)*dt/2
    ! write(22,*) force(theta,param,L,N,position)
    position = position +momentum*dt
    where(position<0)
      position=position+L
    else where(position>=L)
      position=position-L
    endwhere
    zeta=zeta+ (1/Q)*(sum(momentum**2)-g*kbT0)*dt
    momentum = momentum+ force(theta,param,L,N,position)*dt/2
    ! write(22,*) force(theta,param,L,N,position)
    momentum = momentum*exp(-zeta*dt/2)
    t=i*dt

    if(mod(i,250)==0)then
      kbT = (1/g)*(sum(momentum**2))
      write (11,'(5g18.8)') t, kbT,V,total_potential(theta,param,N,position,L),instantaneous_P(theta,param,N,L,V,position,momentum)
    endif

    if(mod(i,1000)==0)then
      gr_dash = RDF(L, position, delta_r, N, V, ndiv)
      do j = 1 , ndiv
        write (20,*) j*delta_r , gr_dash(j)
      end do
    endif

  enddo

  write (28,*) position
  write (12,*) position
  write (13,*) momentum

  kbT0 = 5.0
  do i4 = 1,N1

    pressure_ave = 0.0
    L_before = L


    do i =0,50000
      momentum = momentum*exp(-zeta*dt/2)
      momentum = momentum+force(theta,param,L,N,position)*dt/2
      position = position +momentum*dt
      where(position<0)
        position=position+L
      else where(position>=L)
        position=position-L
      endwhere
      zeta=zeta+ (1/Q)*(sum(momentum**2)-g*kbT0)*dt
      momentum = momentum+ force(theta,param,L,N,position)*dt/2
      momentum = momentum*exp(-zeta*dt/2)
      t=i*dt

      if(mod(i,500)==0)then
        kbT = (1/g)*(sum(momentum**2))
        write (17,'(5g18.8)') t, kbT,V,total_potential(theta,param,N,position,L),instantaneous_P(theta,param,N,L,V,position,momentum)
      endif

      if(mod(i,5000)==0)then
        gr_dash = RDF(L, position, delta_r, N, V, ndiv)
        do j = 1 , ndiv
          write (15,*) j*delta_r , gr_dash(j)
        end do
        write (15,'(/)')
      end if

      if ( i>25000 ) then
        if (mod(i,50) == 0) then
          pressure_ave = pressure_ave + instantaneous_P(theta,param,N,L,V,position,momentum)
        end if
      endif

    enddo

    pressure_ave = pressure_ave/500
    integrand = (pressure_ave - kbT0*rho)/(rho**2)
    free_energy = free_energy + integrand*delta_rho

    write(25,'(5g18.8)') rho, free_energy, integrand, pressure_ave, i4

    write (12,*) position
    write (13,*) momentum
    write (12,'(/)')
    write (13,'(/)')
    write (15,'(/)')
    write (17,'(/)')

    rho = rho + delta_rho
    V = N/rho
    L = V**(1._wp/3)
    position = position*(L/L_before)

  end do

  ! do i4 = 1,N1

  !     potential_ave = 0.0

  !   do i =0,50000
  !     momentum = momentum*exp(-zeta*dt/2)
  !     momentum = momentum+force(theta,param,L,N,position)*dt/2
  !     position = position +momentum*dt
  !     where(position<0)
  !       position=position+L
  !     else where(position>=L)
  !       position=position-L
  !     endwhere
  !     zeta=zeta+ (1/Q)*(sum(momentum**2)-g*kbT0)*dt
  !     momentum = momentum+ force(theta,param,L,N,position)*dt/2
  !     momentum = momentum*exp(-zeta*dt/2)
  !     t=i*dt

  !     if(mod(i,500)==0)then
  !       kbT = (1/g)*(sum(momentum**2))
  !       write (17,'(5g18.8)') t, kbT,V,total_potential(theta,param,N,position,L),instantaneous_P(theta,param,N,L,V,position,momentum)
  !     endif

  !     if(mod(i,5000)==0)then
  !       gr_dash = RDF(L, position, delta_r, N, V, ndiv)
  !       do j = 1 , ndiv
  !         write (15,*) j*delta_r , gr_dash(j)
  !       end do
  !       write (15,'(/)')
  !     end if

  !     if ( i>25000 ) then
  !       if (mod(i,50) == 0) then
  !         potential_ave = potential_ave + total_potential(theta,param,N,position,L)
  !       end if
  !     endif

  !   enddo

  !   potential_ave = potential_ave/500
  !   ! integrand = (pressure_ave - kbT0*rho)/(rho**2)
  !   ! free_energy = free_energy + integrand*delta_rho

  !   write(25,'(2g18.8)') kbT0, potential_ave

  !   write (12,*) position
  !   write (13,*) momentum
  !   write (12,'(/)')
  !   write (13,'(/)')
  !   write (15,'(/)')
  !   write (17,'(/)')

  !   kbT0 = kbT0 - delta_kbT0

  ! end do

end program fluid_phase

