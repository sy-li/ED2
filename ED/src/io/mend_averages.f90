Module mend_averages
  implicit none

Contains

  subroutine mend_normalize(y, step)
    use mend_state_vars, only: mend_model, nwood
    implicit none
    type(mend_model) :: y
    real, intent(inout) :: step
    integer :: iwood

    call mend_normalize_type(y%som, step)
    call mend_normalize_type(y%litt, step)
    do iwood = 1, nwood
       call mend_normalize_type(y%wood(iwood), step)
    enddo

    step = 0.

    return
  end subroutine mend_normalize

  subroutine mend_normalize_type(vars, step)
    use mend_state_vars, only: mend_vars
    use ed_misc_coms, only: dtlsm
    implicit none

    type(mend_vars) :: vars
    real, intent(in) :: step

    call mend_normalize_type_org(vars%cvars, step)
    call mend_normalize_type_org(vars%nvars, step)
    call mend_normalize_type_org(vars%pvars, step)

    vars%fluxes%co2_lost = vars%fluxes%co2_lost / step / dtlsm
    vars%fluxes%nmin = vars%fluxes%nmin / step / dtlsm
    vars%fluxes%nitr = vars%fluxes%nitr / step / dtlsm
    vars%fluxes%ngas_lost = vars%fluxes%ngas_lost / step / dtlsm
    vars%invars%nh4 = vars%invars%nh4 / step
    vars%invars%no3 = vars%invars%no3 / step
    vars%invars%psol = vars%invars%psol / step
    vars%invars%plab = vars%invars%plab / step
    vars%invars%pocc = vars%invars%pocc / step
    vars%invars%psec = vars%invars%psec / step
    vars%invars%ppar = vars%invars%ppar /step

    vars%fluxes%nh4_plant = vars%fluxes%nh4_plant / step / dtlsm
    vars%fluxes%nh4_bnf = vars%fluxes%nh4_bnf / step / dtlsm
    vars%fluxes%no3_plant = vars%fluxes%no3_plant / step / dtlsm
    vars%fluxes%c_leach = vars%fluxes%c_leach / step / dtlsm
    vars%fluxes%n_leach = vars%fluxes%n_leach / step / dtlsm
    vars%fluxes%p_leach = vars%fluxes%p_leach / step / dtlsm
    vars%fluxes%p_plant = vars%fluxes%p_plant / step / dtlsm
    vars%fluxes%nh4_dep = vars%fluxes%nh4_dep / step
    vars%fluxes%no3_dep = vars%fluxes%no3_dep / step
    vars%fluxes%ppar_dep = vars%fluxes%ppar_dep / step
    
    vars%plvars%plant_input_C_pom = vars%plvars%plant_input_C_pom / step / dtlsm
    vars%plvars%plant_input_C_dom = vars%plvars%plant_input_C_dom / step / dtlsm
    vars%plvars%plant_input_N_pom = vars%plvars%plant_input_N_pom / step / dtlsm
    vars%plvars%plant_input_N_dom = vars%plvars%plant_input_N_dom / step / dtlsm
    vars%plvars%plant_input_P_pom = vars%plvars%plant_input_P_pom / step / dtlsm
    vars%plvars%plant_input_P_dom = vars%plvars%plant_input_P_dom / step / dtlsm

    return
  end subroutine mend_normalize_type

  subroutine mend_normalize_type_org(vars, step)
    use mend_state_vars, only: org_vars, npom
    implicit none

    type(org_vars) :: vars
    real, intent(in) :: step
    integer :: ipom

    do ipom = 1, npom
       vars%pom(ipom,:) = vars%pom(ipom,:) / step
       vars%enz_pom(ipom,:) = vars%enz_pom(ipom,:) / step
    enddo

    vars%dom = vars%dom / step
    vars%mom = vars%mom / step
    vars%enz_mom = vars%enz_mom / step
    vars%qom = vars%qom / step
    vars%amb = vars%amb / step
    vars%dmb = vars%dmb / step
    vars%enz_ptase = vars%enz_ptase / step

    return
  end subroutine mend_normalize_type_org

end Module mend_averages
