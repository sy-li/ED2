Module mend_rk4
  
  implicit none

  real, parameter :: min_pool_size = 1.0e-30


Contains

  subroutine mend_rk4_scale(y, dy, scale, htry)
    use mend_state_vars, only: mend_model, nwood
    implicit none
    type(mend_model) :: y
    type(mend_model) :: dy
    type(mend_model) :: scale
    real, intent(in) :: htry
!    integer :: iwood

    call mend_rk4_scale_type(y%som, dy%som, scale%som, htry)
!    call mend_rk4_scale_type(y%litt, dy%litt, scale%litt, htry)
!    do iwood = 1, nwood
!       call mend_rk4_scale_type(y%wood(iwood), dy%wood(iwood), scale%wood(iwood), htry)
!    enddo

    return
  end subroutine mend_rk4_scale

  subroutine mend_rk4_scale_type(vars, d_vars, scale_vars, htry)
    use mend_state_vars, only: mend_vars
    implicit none
    
    type(mend_vars) :: vars
    type(mend_vars) :: d_vars
    type(mend_vars) :: scale_vars
    
    real, parameter :: tiny=1.0e-20
    real, intent(in) :: htry

    call mend_rk4_scale_type_org(vars%cvars, d_vars%cvars, scale_vars%cvars, htry, tiny)
    call mend_rk4_scale_type_org(vars%nvars, d_vars%nvars, scale_vars%nvars, htry, tiny)
    call mend_rk4_scale_type_org(vars%pvars, d_vars%pvars, scale_vars%pvars, htry, tiny)

    scale_vars%invars%nh4 = abs(vars%invars%nh4) +   &
         abs(d_vars%invars%nh4 * htry) + tiny
    scale_vars%invars%no3 = abs(vars%invars%no3) +   &
         abs(d_vars%invars%no3 * htry) + tiny
    scale_vars%invars%psol = abs(vars%invars%psol) +  &
         abs(d_vars%invars%psol * htry) + tiny
    scale_vars%invars%plab = abs(vars%invars%plab) +  &
         abs(d_vars%invars%plab * htry) + tiny
    scale_vars%invars%pocc = abs(vars%invars%pocc) +  &
         abs(d_vars%invars%pocc * htry) + tiny
    scale_vars%invars%ppar = abs(vars%invars%ppar) +  &
         abs(d_vars%invars%ppar * htry) + tiny

    return
  end subroutine mend_rk4_scale_type
    
  subroutine mend_rk4_scale_type_org(vars, d_vars, scale_vars, htry, tiny)
    use mend_state_vars, only: org_vars, npom
    implicit none
    
    type(org_vars) :: vars
    type(org_vars) :: d_vars
    type(org_vars) :: scale_vars
    
    integer :: ipom
    real, intent(in) :: htry
    real, intent(in) :: tiny

    do ipom = 1, npom
       scale_vars%pom(ipom,:) = abs(vars%pom(ipom,:)) +   &
            abs(d_vars%pom(ipom,:) * htry) + tiny
       scale_vars%enz_pom(ipom,:) = abs(vars%enz_pom(ipom,:)) +   &
            abs(d_vars%enz_pom(ipom,:) * htry) + tiny
    enddo
    
    scale_vars%dom = abs(vars%dom) + abs(d_vars%dom * htry) + tiny
    scale_vars%mom = abs(vars%mom) + abs(d_vars%mom * htry) + tiny
    scale_vars%enz_mom = abs(vars%enz_mom) + abs(d_vars%enz_mom * htry) + tiny
    scale_vars%qom = abs(vars%qom) + abs(d_vars%qom * htry) + tiny
    scale_vars%amb = abs(vars%amb) + abs(d_vars%amb * htry) + tiny
    scale_vars%dmb = abs(vars%dmb) + abs(d_vars%dmb * htry) + tiny
    scale_vars%enz_ptase = abs(vars%enz_ptase) +  &
         abs(d_vars%enz_ptase * htry) + tiny

    return
  end subroutine mend_rk4_scale_type_org

  subroutine mend_rk4_inc(y, dy, step, ip1, ip2)
    use mend_state_vars, only: mend_model, nwood
    implicit none
    type(mend_model) :: y
    type(mend_model) :: dy
    real, intent(in) :: step
!    integer :: iwood
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2

    call mend_rk4_inc_type(y%som, dy%som, step, ip1, ip2)
!    call mend_rk4_inc_type(y%litt, dy%litt, step, ip1, ip2)
!    do iwood = 1, nwood
!       call mend_rk4_inc_type(y%wood(iwood), dy%wood(iwood), step, ip1, ip2)
!    enddo

    return
  end subroutine mend_rk4_inc

  subroutine mend_rk4_inc_type(vars, d_vars, step, ip1, ip2)
    use mend_state_vars, only: mend_vars
    use ed_max_dims, only: n_pft
    implicit none

    type(mend_vars) :: vars
    type(mend_vars) :: d_vars
    real, intent(in) :: step
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2
    integer :: ipft

    call mend_rk4_inc_type_org(vars%cvars, d_vars%cvars, step, ip1, ip2)
    call mend_rk4_inc_type_org(vars%nvars, d_vars%nvars, step, ip1, ip2)
    call mend_rk4_inc_type_org(vars%pvars, d_vars%pvars, step, ip1, ip2)

    vars%fluxes%co2_lost(ip1) = vars%fluxes%co2_lost(ip1) +  &
         d_vars%fluxes%co2_lost(ip2) * step
    vars%fluxes%nmin(ip1) = vars%fluxes%nmin(ip1) +  &
         d_vars%fluxes%nmin(ip2) * step
    vars%fluxes%nitr(ip1) = vars%fluxes%nitr(ip1) +  &
         d_vars%fluxes%nitr(ip2) * step
    vars%fluxes%ngas_lost(ip1) = vars%fluxes%ngas_lost(ip1) +   &
         d_vars%fluxes%ngas_lost(ip2) * step
    vars%invars%nh4(ip1) = vars%invars%nh4(ip1) + d_vars%invars%nh4(ip2) * step
    vars%invars%no3(ip1) = vars%invars%no3(ip1) + d_vars%invars%no3(ip2) * step
    vars%invars%psol(ip1) = vars%invars%psol(ip1) +  &
         d_vars%invars%psol(ip2) * step
    vars%invars%plab(ip1) = vars%invars%plab(ip1) +  &
         d_vars%invars%plab(ip2) * step
    vars%invars%pocc(ip1) = vars%invars%pocc(ip1) +  &
         d_vars%invars%pocc(ip2) * step
    vars%invars%psec(ip1) = vars%invars%psec(ip1) +  &
         d_vars%invars%psec(ip2) * step
    vars%invars%ppar(ip1) = vars%invars%ppar(ip1) +   &
         d_vars%invars%ppar(ip2) * step

    do ipft = 1, n_pft
       vars%fluxes%nh4_plant(ipft,ip1) = vars%fluxes%nh4_plant(ipft,ip1) +  &
            d_vars%fluxes%nh4_plant(ipft,ip2) * step
       vars%fluxes%no3_plant(ipft,ip1) = vars%fluxes%no3_plant(ipft,ip1) +  &
            d_vars%fluxes%no3_plant(ipft,ip2) * step
       vars%fluxes%p_plant(ipft,ip1) = vars%fluxes%p_plant(ipft,ip1) +   &
            d_vars%fluxes%p_plant(ipft,ip2) * step
    enddo

    vars%fluxes%nh4_bnf(ip1) = vars%fluxes%nh4_bnf(ip1) +   &
         d_vars%fluxes%nh4_bnf(ip2) * step
    vars%fluxes%c_leach(ip1) = vars%fluxes%c_leach(ip1) +   &
         d_vars%fluxes%c_leach(ip2) * step
    vars%fluxes%n_leach(ip1) = vars%fluxes%n_leach(ip1) +   &
         d_vars%fluxes%n_leach(ip2) * step
    vars%fluxes%p_leach(ip1) = vars%fluxes%p_leach(ip1) +   &
         d_vars%fluxes%p_leach(ip2) * step
    vars%fluxes%nh4_dep(ip1) = vars%fluxes%nh4_dep(ip1) +   &
         d_vars%fluxes%nh4_dep(ip2) * step
    vars%fluxes%no3_dep(ip1) = vars%fluxes%no3_dep(ip1) +   &
         d_vars%fluxes%no3_dep(ip2) * step
    vars%fluxes%ppar_dep(ip1) = vars%fluxes%ppar_dep(ip1) +   &
         d_vars%fluxes%ppar_dep(ip2) * step

    vars%plvars%plant_input_C_pom(:,ip1) =  &
         vars%plvars%plant_input_C_pom(:,ip1) +  &
         d_vars%plvars%plant_input_C_pom(:,ip2) * step
    vars%plvars%plant_input_C_dom(ip1) =   &
         vars%plvars%plant_input_C_dom(ip1) +  &
         d_vars%plvars%plant_input_C_dom(ip2) * step
    vars%plvars%plant_input_N_pom(:,ip1) =   &
         vars%plvars%plant_input_N_pom(:,ip1) +  &
         d_vars%plvars%plant_input_N_pom(:,ip2) * step
    vars%plvars%plant_input_N_dom(ip1) =   &
         vars%plvars%plant_input_N_dom(ip1) +  &
         d_vars%plvars%plant_input_N_dom(ip2) * step
    vars%plvars%plant_input_P_pom(:,ip1) =   &
         vars%plvars%plant_input_P_pom(:,ip1) +  &
         d_vars%plvars%plant_input_P_pom(:,ip2) * step
    vars%plvars%plant_input_P_dom(ip1) =   &
         vars%plvars%plant_input_P_dom(ip1) +  &
         d_vars%plvars%plant_input_P_dom(ip2) * step

    return
  end subroutine mend_rk4_inc_type

  subroutine mend_rk4_inc_type_org(vars, dvars, step, ip1, ip2)
    use mend_state_vars, only: org_vars, npom
    implicit none

    type(org_vars) :: vars, dvars
    real, intent(in) :: step
    integer :: ipom
    integer, intent(in) :: ip1
    integer, intent(in) :: ip2

    do ipom = 1, npom
       vars%pom(ipom,ip1) = vars%pom(ipom,ip1) + dvars%pom(ipom,ip2) * step
       vars%enz_pom(ipom,ip1) = vars%enz_pom(ipom,ip1) +   &
            dvars%enz_pom(ipom,ip2) * step
    enddo

    vars%dom(ip1) = vars%dom(ip1) + dvars%dom(ip2) * step
    vars%mom(ip1) = vars%mom(ip1) + dvars%mom(ip2) * step
    vars%enz_mom(ip1) = vars%enz_mom(ip1) + dvars%enz_mom(ip2) * step
    vars%qom(ip1) = vars%qom(ip1) + dvars%qom(ip2) * step
    vars%amb(ip1) = vars%amb(ip1) + dvars%amb(ip2) * step
    vars%dmb(ip1) = vars%dmb(ip1) + dvars%dmb(ip2) * step
    vars%enz_ptase(ip1) = vars%enz_ptase(ip1) + dvars%enz_ptase(ip2) * step

    return
  end subroutine mend_rk4_inc_type_org

  subroutine mend_rk4_errmax(scalev, errv, errmax)
    use mend_state_vars, only: mend_model, nwood
    implicit none
    type(mend_model) :: errv
    type(mend_model) :: scalev
    real(kind=8), intent(inout) :: errmax
!    integer :: iwood

    call mend_rk4_errmax_type(scalev%som, errv%som, errmax)
!    call mend_rk4_errmax_type(scalev%litt, errv%litt, errmax)
!    do iwood = 1, nwood
!       call mend_rk4_errmax_type(scalev%wood(iwood), errv%wood(iwood), errmax)
!    enddo

    return
  end subroutine mend_rk4_errmax

  subroutine mend_rk4_errmax_type(scale_detr, err_detr, errmax)
    use mend_state_vars, only: mend_vars
    implicit none

    real :: err
    real(kind=8), intent(inout) :: errmax
    type(mend_vars) :: scale_detr
    type(mend_vars) :: err_detr

    call mend_rk4_errmax_type_org(scale_detr%cvars, err_detr%cvars, errmax)
    call mend_rk4_errmax_type_org(scale_detr%nvars, err_detr%nvars, errmax)
    call mend_rk4_errmax_type_org(scale_detr%pvars, err_detr%pvars, errmax)

    err = abs(err_detr%invars%nh4(1) / scale_detr%invars%nh4(1))
    errmax = max(errmax, err)
    err = abs(err_detr%invars%no3(1) / scale_detr%invars%no3(1))
    errmax = max(errmax, err)
    err = abs(err_detr%invars%psol(1) / scale_detr%invars%psol(1))
    errmax = max(errmax, err)
    err = abs(err_detr%invars%plab(1) / scale_detr%invars%plab(1))
    errmax = max(errmax, err)
    err = abs(err_detr%invars%pocc(1) / scale_detr%invars%pocc(1))
    errmax = max(errmax, err)
    err = abs(err_detr%invars%ppar(1) / scale_detr%invars%ppar(1))
    errmax = max(errmax, err)

    return
  end subroutine mend_rk4_errmax_type

  subroutine mend_rk4_errmax_type_org(scale_vars, err_vars, errmax)
    use mend_state_vars, only: org_vars, npom
    implicit none

    type(org_vars) :: scale_vars, err_vars
    real(kind=8), intent(inout) :: errmax
    integer :: ipom
    real :: err

    do ipom = 1, npom
       err = abs(err_vars%pom(ipom,1) / scale_vars%pom(ipom,1))
       errmax = max(errmax, err)
       err = abs(err_vars%enz_pom(ipom,1) / scale_vars%enz_pom(ipom,1))
       errmax = max(errmax, err)
    enddo

    err = abs(err_vars%dom(1) / scale_vars%dom(1))
    errmax = max(errmax, err)
    err = abs(err_vars%mom(1) / scale_vars%mom(1))
    errmax = max(errmax, err)
    err = abs(err_vars%enz_mom(1) / scale_vars%enz_mom(1))
    errmax = max(errmax, err)
    err = abs(err_vars%qom(1) / scale_vars%qom(1))
    errmax = max(errmax, err)
    err = abs(err_vars%amb(1) / scale_vars%amb(1))
    errmax = max(errmax, err)
    err = abs(err_vars%dmb(1) / scale_vars%dmb(1))
    errmax = max(errmax, err)
    err = abs(err_vars%enz_ptase(1) / scale_vars%enz_ptase(1))
    errmax = max(errmax, err)

    return
  end subroutine mend_rk4_errmax_type_org

  subroutine mend_rk4_sanity(mend, reject_step)
    use mend_state_vars, only: mend_model, nwood
    implicit none
    
    type(mend_model) :: mend
!    integer :: iwood
    logical, intent(inout) :: reject_step

    call mend_rk4_sanity_type(mend%som, reject_step)
!    call mend_rk4_sanity_type(mend%litt, reject_step)
!    do iwood = 1, nwood
!       call mend_rk4_sanity_type(mend%wood(iwood), reject_step)
!    enddo

    return
  end subroutine mend_rk4_sanity

  subroutine mend_rk4_sanity_type(vars, reject_step)
    use mend_state_vars, only: mend_vars
    implicit none
    type(mend_vars) :: vars
    logical, intent(inout) :: reject_step

!if(reject_step)print*,1,reject_step
    call mend_rk4_sanity_type_org(vars%cvars, reject_step)
!if(reject_step)print*,11,reject_step
    call mend_rk4_sanity_type_org(vars%nvars, reject_step)
!if(reject_step)print*,12,reject_step
    call mend_rk4_sanity_type_org(vars%pvars, reject_step)
!if(reject_step)print*,2,reject_step
    if(vars%fluxes%nitr(1) < 0.)reject_step = .true.
!if(reject_step)print*,21,reject_step
    if(vars%fluxes%co2_lost(1) < 0.)reject_step = .true.
    if(vars%fluxes%ngas_lost(1) < 0.)reject_step = .true.
!if(reject_step)print*,22,reject_step
    if(vars%invars%nh4(1) < min_pool_size)reject_step = .true.
    if(vars%invars%no3(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,23,reject_step,vars%invars%nh4(1),vars%invars%no3(1)
    if(vars%invars%psol(1) < min_pool_size)reject_step = .true.
    if(vars%invars%plab(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,24,reject_step
 !   if(vars%invars%pocc(1) < min_pool_size)reject_step = .true.
!    if(vars%invars%ppar(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,3,reject_step
    if(any(vars%fluxes%nh4_plant(:,1) < 0.))reject_step = .true.
    if(vars%fluxes%nh4_bnf(1) < 0.)reject_step = .true.
    if(any(vars%fluxes%no3_plant(:,1) < 0.))reject_step = .true.
    if(vars%fluxes%c_leach(1) < 0.)reject_step = .true.
    if(vars%fluxes%n_leach(1) < 0.)reject_step = .true.
    if(vars%fluxes%p_leach(1) < 0.)reject_step = .true.
    if(any(vars%fluxes%p_plant(:,1) < 0.))reject_step = .true.
!if(reject_step)print*,4,reject_step

!    if(reject_step)then
!       print*,vars%cvars%pom
!       print*,vars%nvars%pom
!       print*,vars%pvars%pom
!       print*,vars%cvars%pom/vars%nvars%pom
!       print*,vars%cvars%pom/vars%pvars%pom
!       print*,vars%cvars%enz_pom
!       print*,vars%nvars%enz_pom
!       print*,vars%pvars%enz_pom
!       print*,vars%cvars%enz_pom/vars%nvars%enz_pom
!       print*,vars%cvars%enz_pom/vars%pvars%enz_pom
!    endif

    return
  end subroutine mend_rk4_sanity_type

  subroutine mend_rk4_sanity_type_org(vars, reject_step)
    use mend_state_vars, only: org_vars, npom
    implicit none

    type(org_vars) :: vars
    logical, intent(inout) :: reject_step
    integer :: ipom

    do ipom = 1, npom
       if(vars%pom(ipom,1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,111,reject_step,ipom,vars%pom(ipom,1)
       if(vars%enz_pom(ipom,1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,112,reject_step,ipom,vars%enz_pom(ipom,1)
    enddo

    if(vars%dom(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,113,reject_step,vars%dom(1)
    if(vars%mom(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,114,reject_step,vars%mom(1)
    if(vars%enz_mom(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,115,reject_step,vars%enz_mom(1),vars%dom(1),vars%amb(1)
    if(vars%qom(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,116,reject_step,ipom
    if(vars%amb(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,117,reject_step,ipom
    if(vars%dmb(1) < min_pool_size)reject_step = .true.
!if(reject_step)print*,118,reject_step,ipom
    if(vars%enz_ptase(1) < min_pool_size)reject_step = .true.

    return
  end subroutine mend_rk4_sanity_type_org

  subroutine mend_zero_fluxes(mend)
    use mend_state_vars, only: mend_model
    implicit none
    type(mend_model) :: mend

    mend%som%fluxes%co2_lost = 0.
    mend%som%fluxes%ngas_lost = 0.
    mend%som%fluxes%c_leach = 0.
    mend%som%fluxes%n_leach = 0.
    mend%som%fluxes%p_leach = 0.
    mend%som%fluxes%nh4_bnf = 0.
    mend%som%fluxes%nh4_plant = 0.
    mend%som%fluxes%no3_plant = 0.
    mend%som%fluxes%p_plant = 0.
    mend%som%fluxes%nmin = 0.
    mend%som%fluxes%nitr = 0.
    mend%som%fluxes%nh4_dep = 0.
    mend%som%fluxes%no3_dep = 0.
    mend%som%fluxes%ppar_dep = 0.
    mend%som%plvars%plant_input_C_pom = 0.
    mend%som%plvars%plant_input_C_dom = 0.
    mend%som%plvars%plant_input_N_pom = 0.
    mend%som%plvars%plant_input_N_dom = 0.
    mend%som%plvars%plant_input_P_pom = 0.
    mend%som%plvars%plant_input_P_dom = 0.

    return
  end subroutine mend_zero_fluxes


end Module mend_rk4
