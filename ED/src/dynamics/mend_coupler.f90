Module mend_coupler
  implicit none

Contains

  subroutine mend_update_parameters_coupler(soil_tempk, soil_water,  &
       ntext_soil, pH)
    use mend_consts_coms, only: som_consts, som_consts_base
    use grid_coms, only: nzg
    use soil_coms, only: soil
    use consts_coms, only: wdns, grav
    use mend_diagnose, only: mend_update_parameters
    implicit none

    real, dimension(nzg), intent(in) :: soil_tempk
    real, dimension(nzg), intent(in) :: soil_water
    integer, dimension(nzg), intent(in) :: ntext_soil
    
    real :: tp
    real :: wp
    real, intent(in) :: pH !=6.24  ! SROAK: 5.82; SRTDF: 6.24; PVTDF: 6.83
    real :: wfp
    integer :: ntxt

    ntxt = ntext_soil(nzg)
    tp = soil_tempk(nzg) - 273.15
    wfp = soil_water(nzg) / soil(ntxt)%slmsts

!tp = tp + 2.
!wfp = 0.8

!    root_avail_water = 0.
!    do k = krdepth, nzg
!       root_avail_water = root_avail_water + soil_water(k)*dslz(k) !m3/m2
!    enddo
!    root_avail_water = root_avail_water / (-slz(krdepth)) !m3/m3
!    wgpfrac = root_avail_water / slmsts
    wp = wdns * grav * soil(ntxt)%slpots / wfp**soil(ntxt)%slbs * 1.0e-6

    call mend_update_parameters(tp, wp, pH, wfp, som_consts, som_consts_base, &
         soil(ntxt)%slmsts, soil(ntxt)%sfldcap)

    return
  end subroutine mend_update_parameters_coupler

  subroutine mend_init(sens_params)
    use ed_state_vars, only: edgrid_g, edtype, polygontype, sitetype
    use nutrient_constants, only: soil_bulk_den, soil_ph
    use mend_consts_coms, only: mend_init_consts, litt_consts, wood_consts, &
         som_consts
    use mend_som, only: som_init
    use mend_state_vars, only: npom, nwood, mend_zero_vars, mend_mm_time
    use nutrient_constants, only: soil_cpct, soil_som_c2n, soil_totp,   &
         soil_extrp
    implicit none

    type(edtype), pointer :: cgrid
    type(polygontype), pointer :: cpoly
    type(sitetype), pointer :: csite
    integer :: ipy
    integer :: isi
    integer :: ipa
    integer :: iwood
    integer, intent(in) :: sens_params

    write(*,*)'sens_params = ',sens_params
    call mend_init_consts(sens_params)

    mend_mm_time = 0.

    cgrid => edgrid_g(1)
    do ipy = 1, cgrid%npolygons
       cpoly => cgrid%polygon(ipy)
       do isi = 1, cpoly%nsites
          csite => cpoly%site(isi)
          do ipa = 1, csite%npatches

             call mend_zero_vars(csite%mend_mm, ipa, ipa)

             csite%mend%bulk_den(ipa) = soil_bulk_den
             csite%mend%pH(ipa) = soil_ph
             csite%mend_mm%bulk_den(ipa) = csite%mend%bulk_den(ipa)
             csite%mend_mm%pH(ipa) = csite%mend%pH(ipa)

             call som_init(npom,  &
                  csite%mend%som%cvars%pom(:,ipa), csite%mend%som%cvars%dom(ipa), &
                  csite%mend%som%cvars%enz_pom(:,ipa), csite%mend%som%cvars%mom(ipa),  &
                  csite%mend%som%cvars%qom(ipa), csite%mend%som%cvars%enz_mom(ipa),  &
                  csite%mend%som%cvars%amb(ipa), csite%mend%som%cvars%dmb(ipa),   &
                  csite%mend%som%nvars%pom(:,ipa), csite%mend%som%nvars%dom(ipa),  &
                  csite%mend%som%nvars%enz_pom(:,ipa), csite%mend%som%nvars%mom(ipa),  &
                  csite%mend%som%nvars%qom(ipa), csite%mend%som%nvars%enz_mom(ipa),  &
                  csite%mend%som%nvars%amb(ipa), csite%mend%som%nvars%dmb(ipa),   &
                  csite%mend%som%pvars%pom(:,ipa), csite%mend%som%pvars%dom(ipa), &
                  csite%mend%som%pvars%enz_pom(:,ipa), csite%mend%som%pvars%mom(ipa),  &
                  csite%mend%som%pvars%qom(ipa), csite%mend%som%pvars%enz_mom(ipa),  &
                  csite%mend%som%pvars%amb(ipa), csite%mend%som%pvars%dmb(ipa),   &
                  csite%mend%som%fluxes%co2_lost(ipa),  &
                  csite%mend%som%fluxes%nmin(ipa), &
                  csite%mend%som%fluxes%nitr(ipa), &
                  csite%mend%som%invars%nh4(ipa), &
                  csite%mend%som%invars%no3(ipa), csite%mend%som%invars%psol(ipa), &
                  csite%mend%som%invars%plab(ipa), csite%mend%som%fluxes%ngas_lost(ipa), &
                  csite%mend%som%cvars%enz_ptase(ipa), csite%mend%som%nvars%enz_ptase(ipa), &
                  csite%mend%som%pvars%enz_ptase(ipa), csite%mend%som%fluxes%nh4_dep(ipa),  &
                  csite%mend%som%fluxes%no3_dep(ipa), csite%mend%som%fluxes%ppar_dep(ipa), &
                  csite%mend%som%fluxes%nh4_plant(:,ipa), csite%mend%som%fluxes%nh4_bnf(ipa), &
                  csite%mend%som%fluxes%no3_plant(:,ipa), csite%mend%som%fluxes%c_leach(ipa), &
                  csite%mend%som%fluxes%n_leach(ipa), csite%mend%som%fluxes%p_leach(ipa),  &
                  csite%mend%som%fluxes%p_plant(:,ipa), csite%mend%som%invars%pocc(ipa),  &
                  csite%mend%som%invars%psec(ipa), &
                  csite%mend%som%invars%ppar(ipa), csite%mend%som%plvars%enz_plant_n(:,ipa),  &
                  csite%mend%som%plvars%enz_plant_p(:,ipa), csite%mend%som%plvars%vnh4up_plant(:,ipa),  &
                  csite%mend%som%plvars%vno3up_plant(:,ipa),  &
                  csite%mend%som%plvars%vpup_plant(:,ipa),   &
                  csite%mend%som%plvars%plant_input_C_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_N_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_P_pom(:,ipa), &
                  csite%mend%som%plvars%plant_input_C_dom(ipa), &
                  csite%mend%som%plvars%plant_input_N_dom(ipa), &
                  csite%mend%som%plvars%plant_input_P_dom(ipa), &
                  som_consts, csite%mend%bulk_den(ipa),   &
                  soil_cpct, soil_som_c2n, soil_totp, soil_extrp)
          enddo
       enddo
    enddo
    
    return
  end subroutine mend_init

  subroutine mend_extern_forcing(mend, ipa, ncohorts, broot, nplant, &
       pft, krdepth, slden, nstorage, pstorage,   &
       nstorage_max, pstorage_max, water_supply_nl, lai)
    use mend_state_vars, only: mend_model, nwood
    use mend_som, only: som_extern_forcing
    use mend_consts_coms, only: som_consts
    use nutrient_constants, only: ndep_rate, pdep_rate, ndep_appl, pdep_appl, nlsl
    use ed_misc_coms, only: current_time
    use mend_plant, only: som_plant_enzymes, litt_plant_enzymes, wood_plant_enzymes
    use soil_coms, only: nzg
    implicit none
    type(mend_model) :: mend
    integer, intent(in) :: ipa
    integer :: iwood
    integer, intent(in) :: ncohorts
    integer, intent(in), dimension(ncohorts) :: pft
    integer, intent(in), dimension(ncohorts) :: krdepth
    real, intent(in), dimension(ncohorts) :: broot
    real, intent(in), dimension(ncohorts) :: nplant
    real, intent(in), dimension(ncohorts) :: nstorage
    real, intent(in), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: nstorage_max
    real, intent(in), dimension(ncohorts) :: pstorage_max
    real, intent(in), dimension(ncohorts) :: lai
    real, intent(in), dimension(ncohorts) :: water_supply_nl
    real :: broot_total
    integer :: ico
    real, intent(in) :: slden

print*,water_supply_nl;stop
    call som_extern_forcing(ndep_rate, som_consts, slden, &
         mend%som%fluxes%nh4_dep(ipa), mend%som%fluxes%no3_dep(ipa), &
         pdep_rate, mend%som%fluxes%ppar_dep(ipa), current_time%year, &
         ndep_appl, pdep_appl)
    
    call som_plant_enzymes(ncohorts, broot, nplant, pft,  &
         krdepth, slden, mend%som%plvars%enz_plant_n(:,ipa),  &
         mend%som%plvars%enz_plant_p(:,ipa), &
         mend%som%plvars%vnh4up_plant(:,ipa),  &
         mend%som%plvars%vno3up_plant(:,ipa), &
         mend%som%plvars%vpup_plant(:,ipa), som_consts, nstorage, pstorage, &
         nstorage_max, pstorage_max, water_supply_nl, lai)

!    call litt_extern_forcing(mend%litt%fluxes%nh4_dep(ipa), &
!         mend%litt%fluxes%no3_dep(ipa), mend%litt%fluxes%ppar_dep(ipa))
    
!    call litt_plant_enzymes(mend%litt%plvars%enz_plant_n(:,ipa), &
!         mend%litt%plvars%enz_plant_p(:,ipa), &
!         mend%litt%plvars%vnh4up_plant(:,ipa),  &
!         mend%litt%plvars%vno3up_plant(:,ipa), &
!         mend%litt%plvars%vpup_plant(:,ipa))
    
!    do iwood = 1, nwood
!       call wood_extern_forcing(iwood, mend%wood(iwood)%fluxes%nh4_dep(ipa), &
!            mend%wood(iwood)%fluxes%no3_dep(ipa),  &
!            mend%wood(iwood)%fluxes%ppar_dep(ipa))
!       call wood_plant_enzymes(iwood, &
!            mend%wood(iwood)%plvars%enz_plant_n(:,ipa), &
!            mend%wood(iwood)%plvars%enz_plant_p(:,ipa), &
!            mend%wood(iwood)%plvars%vnh4up_plant(:,ipa), &
!            mend%wood(iwood)%plvars%vno3up_plant(:,ipa), &
!            mend%wood(iwood)%plvars%vpup_plant(:,ipa))
!    enddo

    return
  end subroutine mend_extern_forcing

  subroutine mend_derivs_coupler(som, d_som, litt, d_litt, wood, d_wood, &
       csite, ipa, som_water_drainage, soil_water, d_can_co2, &
       d_co2budget_storage,ccapcani, ntext_soil)
    use mend_exchange, only: litt2som_exchange, plant2litt_exchange, &
         mend_plant2som_exchange, plant2wood_exchange, wood2litt_exchange, &
         wood2som_exchange, zero_exchange_vars, inc_exchange_vars, &
         wood2wood_exchange, som2canopy_exchange
    use ed_state_vars, only: sitetype
    use mend_derivs, only: mend_derivs_layer
    use mend_consts_coms, only: som_consts, litt_consts, wood_consts
    use mend_state_vars, only: mend_vars, exchange_vars, npom, nwood
    use grid_coms, only: nzg
    use soil_coms, only: dslz, soil
    use consts_coms, only: wdns, pi1
    use nutrient_constants, only: nlsl

    implicit none

    real(kind=8), intent(in) :: som_water_drainage
    real(kind=8), dimension(nzg), intent(in) :: soil_water
    real :: litt_water_drainage
    real :: wood_water_drainage
    integer, intent(in) :: ipa
    type(sitetype), target :: csite
    type(mend_vars) :: som
    type(mend_vars) :: d_som
    type(mend_vars) :: litt
    type(mend_vars) :: d_litt
    type(mend_vars), dimension(nwood) :: wood
    type(mend_vars), dimension(nwood) :: d_wood
    real(kind=8), intent(in) :: ccapcani
    type(exchange_vars) :: litt2som
    type(exchange_vars) :: plant2litt
    type(exchange_vars) :: plant2som
    type(exchange_vars), dimension(nwood) :: plant2wood
    type(exchange_vars), dimension(nwood) :: wood2litt
    type(exchange_vars), dimension(nwood) :: wood2som
    type(exchange_vars) :: wood2litt_sum
    type(exchange_vars) :: wood2som_sum
    type(exchange_vars), dimension(nwood, nwood) :: wood2wood

    real, dimension(npom) :: input_pom_c_net
    real, dimension(npom) :: input_pom_n_net
    real, dimension(npom) :: input_pom_p_net
    real :: input_dom_c_net
    real :: input_dom_n_net
    real :: input_dom_p_net
    real :: input_nh4_net
    real :: input_no3_net
    real :: input_psol_net
    real :: input_ppar_net
    integer :: iwood
    integer :: jwood
    real :: gm2_mgg
    real :: total_water
    integer :: k
    real :: som_water_drainage_ps  ! units: 1/s
    real(kind=8), intent(inout) :: d_can_co2
    real(kind=8), intent(inout) :: d_co2budget_storage
    real :: wfp
    integer, dimension(nzg) :: ntext_soil

    wfp = soil_water(nzg) / soil(ntext_soil(nzg))%slmsts
!    wfp = soil_water(nzg) / soil(csite%ntext_soil(nzg,ipa))%slmsts
    gm2_mgg = 1. / (som_consts%eff_soil_depth * csite%mend%bulk_den(ipa))
    total_water = 0.
    do k = nlsl, nzg
       total_water = total_water + soil_water(k) * dslz(k)
    enddo
    ! No negative drainage.
    som_water_drainage_ps = max(0., som_water_drainage) / (total_water * wdns)

    call mend_plant2som_exchange(npom, &
         plant2som%pom_c, plant2som%pom_n, plant2som%pom_p, &
         plant2som%dom_c, plant2som%dom_n, plant2som%dom_p, &
         plant2som%nh4, plant2som%no3, plant2som%psol,   &
         csite%plant_input_C(:,ipa), csite%plant_input_N(:,ipa),  &
         csite%plant_input_P(:,ipa))

    d_som%plvars%plant_input_C_pom(:,1) = plant2som%pom_c * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_N_pom(:,1) = plant2som%pom_n * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_P_pom(:,1) = plant2som%pom_p * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_C_dom = plant2som%dom_c * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_N_dom = plant2som%dom_n * gm2_mgg * 1000./86400.
    d_som%plvars%plant_input_P_dom = plant2som%dom_p * gm2_mgg * 1000./86400.

    input_pom_c_net = (plant2som%pom_c+litt2som%pom_c+wood2som_sum%pom_c) * &
         gm2_mgg * 1000./86400.
    input_dom_c_net = (plant2som%dom_c+litt2som%dom_c+wood2som_sum%dom_c) * &
         gm2_mgg * 1000./86400.
    input_pom_n_net = (plant2som%pom_n+litt2som%pom_n+wood2som_sum%pom_n) * &
         gm2_mgg * 1000./86400.
    input_dom_n_net = (plant2som%dom_n+litt2som%dom_n+wood2som_sum%dom_n) * &
         gm2_mgg * 1000./86400.
    input_pom_p_net = (plant2som%pom_p+litt2som%pom_p+wood2som_sum%pom_p) * &
         gm2_mgg * 1000./86400.
    input_dom_p_net = (plant2som%dom_p+litt2som%dom_p+wood2som_sum%dom_p) * &
         gm2_mgg * 1000./86400.
    input_nh4_net = (plant2som%nh4 + litt2som%nh4 + wood2som_sum%nh4) * &
         gm2_mgg  * 1000./86400. + som%fluxes%nh4_dep(1)
    input_no3_net = (plant2som%no3 + litt2som%no3 + wood2som_sum%no3) * &
         gm2_mgg  * 1000./86400. + som%fluxes%no3_dep(1)
    input_psol_net = (plant2som%psol + litt2som%psol + wood2som_sum%psol) * &
         gm2_mgg  * 1000./86400. + som%fluxes%ppar_dep(1)
    input_ppar_net = 0.
    call mend_derivs_layer(npom, som_consts,  &
         som%cvars%pom(:,1), input_pom_c_net, d_som%cvars%pom(:,1),  &
         som%cvars%dom(1), input_dom_c_net, d_som%cvars%dom(1),  &
         som%cvars%enz_pom(:,1), d_som%cvars%enz_pom(:,1),  &
         som%cvars%mom(1), d_som%cvars%mom(1),  &
         som%cvars%qom(1), d_som%cvars%qom(1),  &
         som%cvars%enz_mom(1), d_som%cvars%enz_mom(1),  &
         som%cvars%amb(1), d_som%cvars%amb(1),   &
         som%cvars%dmb(1), d_som%cvars%dmb(1),   &
         som%fluxes%co2_lost(1), d_som%fluxes%co2_lost(1),  &
         som%fluxes%nmin(1), d_som%fluxes%nmin(1),  &
         som%fluxes%nitr(1), d_som%fluxes%nitr(1),  &
         som%nvars%pom(:,1), input_pom_n_net, d_som%nvars%pom(:,1),  &
         som%nvars%dom(1), input_dom_n_net, d_som%nvars%dom(1),  &
         som%nvars%enz_pom(:,1), d_som%nvars%enz_pom(:,1),  &
         som%nvars%mom(1), d_som%nvars%mom(1),   &
         som%nvars%qom(1), d_som%nvars%qom(1),   &
         som%nvars%enz_mom(1), d_som%nvars%enz_mom(1),   &
         som%nvars%amb(1), d_som%nvars%amb(1),  &
         som%nvars%dmb(1), d_som%nvars%dmb(1),  &
         som%invars%nh4(1), input_nh4_net, d_som%invars%nh4(1), &
         som%invars%no3(1), input_no3_net, d_som%invars%no3(1),  &
         som%pvars%pom(:,1), input_pom_p_net, d_som%pvars%pom(:,1),  &
         som%pvars%dom(1), input_dom_p_net, d_som%pvars%dom(1),  &
         som%pvars%enz_pom(:,1), d_som%pvars%enz_pom(:,1),  &
         som%pvars%mom(1), d_som%pvars%mom(1),   &
         som%pvars%qom(1), d_som%pvars%qom(1),  &
         som%pvars%enz_mom(1), d_som%pvars%enz_mom(1),   &
         som%pvars%amb(1), d_som%pvars%amb(1),  &
         som%pvars%dmb(1), d_som%pvars%dmb(1),   &
         som%invars%psol(1), input_psol_net, d_som%invars%psol(1),  &
         som%invars%plab(1), d_som%invars%plab(1), &
         som%fluxes%ngas_lost(1), d_som%fluxes%ngas_lost(1), &
         som%cvars%enz_ptase(1), d_som%cvars%enz_ptase(1), &
         som%nvars%enz_ptase(1), d_som%nvars%enz_ptase(1), &
         som%pvars%enz_ptase(1), d_som%pvars%enz_ptase(1),  &
         d_som%fluxes%nh4_plant(:,1), &
         som%fluxes%nh4_bnf(1), d_som%fluxes%nh4_bnf(1), &
         d_som%fluxes%no3_plant(:,1), &
         som%fluxes%c_leach(1), d_som%fluxes%c_leach(1), &
         som%fluxes%n_leach(1), d_som%fluxes%n_leach(1),  &
         som%fluxes%p_leach(1), d_som%fluxes%p_leach(1), &
         d_som%fluxes%p_plant(:,1), &
         som%invars%pocc(1), d_som%invars%pocc(1), &
         som%invars%ppar(1), d_som%invars%ppar(1), input_ppar_net,  &
         som%plvars%enz_plant_n(:,1), som%plvars%enz_plant_p(:,1), &
         som%plvars%vnh4up_plant(:,1), som%plvars%vno3up_plant(:,1),  &
         som%plvars%vpup_plant(:,1), som_water_drainage_ps, &
         csite%mend%bulk_den(ipa), pi1, wfp)

    call som2canopy_exchange(d_som%fluxes%co2_lost(1),  &
         csite%mend%bulk_den(ipa), som_consts, &
         d_can_co2, d_co2budget_storage, ccapcani)

    return
  end subroutine mend_derivs_coupler

  subroutine mend_update_diag(mend)
    use mend_state_vars, only: mend_model, nwood
    use mend_consts_coms, only: som_consts, litt_consts, wood_consts
    use mend_diagnose, only: mend_update_diag_layer
    implicit none
    type(mend_model) :: mend
    integer :: iwood

    call mend_update_diag_layer(mend%som, som_consts, 1, mend%bulk_den(1))
!    call mend_update_diag_layer(mend%litt, litt_consts, 1, mend%bulk_den(1))
!    do iwood = 1, nwood
!       call mend_update_diag_layer(mend%wood(iwood), wood_consts(iwood), 1, &
!            mend%bulk_den(1))
!    enddo

    return
  end subroutine mend_update_diag

  subroutine mend_slow_P(mend, ipa)
    use mend_state_vars, only: mend_model, nwood
    use mend_derivs, only: mend_slow_P_layer
    use mend_consts_coms, only: som_consts, litt_consts, wood_consts
    use mend_diagnose, only: mend_update_diag_layer
    implicit none
    type(mend_model) :: mend
    integer, intent(in) :: ipa
    integer :: iwood

    call mend_slow_P_layer(som_consts, mend%som%invars%plab(ipa), &
         mend%som%invars%pocc(ipa), mend%som%invars%psec(ipa), &
         mend%som%invars%ppar(ipa), &
         mend%som%invars%psol(ipa), mend%bulk_den(ipa))
    call mend_update_diag_layer(mend%som, som_consts, ipa, mend%bulk_den(ipa))

    return
  end subroutine mend_slow_P

end Module mend_coupler
