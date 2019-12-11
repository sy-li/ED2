Module mend_plant
  implicit none

  real, parameter :: water_supply_scale = 0.01

Contains

  subroutine som_plant_enzymes(ncohorts, broot, nplant, pft, &
       krdepth, slden, enz_plant_n, &
       enz_plant_p, vnh4up_plant, vno3up_plant, vpup_plant, consts, &
       nstorage, pstorage, nstorage_min, pstorage_min, water_supply, lai)
    use mend_consts_coms, only: decomp_consts
    use pft_coms, only: root_beta
    use soil_coms, only: slz
    use nutrient_constants, only: nstorage_max_factor, pstorage_max_factor, nlsl
    use ed_max_dims, only: n_pft
    implicit none

    integer, intent(in) :: ncohorts
    integer, dimension(ncohorts), intent(in) :: pft
    integer, dimension(ncohorts), intent(in) :: krdepth
    real, dimension(ncohorts), intent(in) :: broot
    real, dimension(ncohorts), intent(in) :: nplant
    real, dimension(ncohorts), intent(in) :: water_supply
    real, dimension(ncohorts), intent(in) :: lai
    real, intent(in) :: slden
    real, intent(out), dimension(n_pft) :: enz_plant_n
    real, intent(out), dimension(n_pft) :: enz_plant_p
    real, intent(out), dimension(n_pft) :: vnh4up_plant
    real, intent(out), dimension(n_pft) :: vno3up_plant
    real, intent(out), dimension(n_pft) :: vpup_plant
    real, intent(in), dimension(ncohorts) :: nstorage
    real, intent(in), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: nstorage_min
    real, intent(in), dimension(ncohorts) :: pstorage_min
    type(decomp_consts) :: consts
    integer :: ico
    real :: broot_nl
    real :: n_limit_factor
    real :: p_limit_factor
    real :: nstorage_total
    real :: nstorage_min_total
    real :: pstorage_total
    real :: pstorage_min_total
    real :: nplant_total
    integer, parameter :: nutr_limit_scheme = 1
    real :: transp_fact
    real, save :: transp_fact_max = 0.

    enz_plant_n = 0.
    enz_plant_p = 0.
    vnh4up_plant = 0.
    vno3up_plant = 0.
    vpup_plant = 0.

!    print*,sum(water_supply/lai)/600*86400.*5.
! Aim for 3
! Normalization factor = 3/5/86400*600*lai = 0.004 * lai

    ! Plant enzyme carrier abundance.
    ! (1) In a typical tropical forest, Zhu et al. assume that the number of plant
    ! carrier enzymes is equal to the number of microbial carrier
    ! enzymes. But this statement is quite ambiguous. Does it apply to NH4
    ! carriers, PO4 carriers, sum of all carriers? Does it mean that the number of 
    ! enzymes is the same, or the amount of C in enzymes is the same? 
    ! Say that it applies to the number of enzymes, in moles. Further assume
    ! that it applies to each enzyme type individually.
    ! (2) Even more confusingly, they say that the scaling factor between
    ! microbial biomass and enzyme abundance is 0.05. Does this mean 0.05 g C enzyme
    ! per g C microbe? Does it apply to all types of transport enzymes - is it a sum
    ! over types? Are transport enzymes generic? Further, Tang and Riley suggest a huge
    ! possible range for this number. Should transport enzymes rather scale with surface
    ! area? I will assume a generic transport enzyme. 
    ! "Abundance" is 0.05 amb_c [mg enz C / g soil].
    ! They assume amb_c = 0.1 gC/m2 or, say, 7e-4 mgC/gSoil.
    ! Enzyme amount is then 4e-5 mgC/gSoil.
    ! They estimate plant root biomass at 400 gC/m2, or 3 mgC/gSoil.
    ! Then, enzyme C per fine root C should be: 1e-5. This agrees with their 0.0000125.
    ! ANOTHER APPROACH.
    ! Say that it applies to number counts (moles).
    ! Number of enzymes is 0.05 * amb_c * (g soil) * (# microbes/ mg micr biomass C)
    ! Plants have same number.
    ! 0.05 * amb_c * (g soil) * (# microbes / mg micr biomass C) = (# plant enz / mg C fine root) (mg C fine root)
    ! 0.05 * amb_c * (g soil) * (# microbes / mg micr biomass C) = (# plant enz / mg C fine root) (3 mgC/gSoil) (g soil)
    ! (# plant enz / mg C fine root) = 1e-5 * (# microbes / mg micr biomass C) = K1
    ! This, I guess, would represent the sum of all transporter enzymes.
    ! 5e7 cells/gram soil, 5e10 cell/kgSoil, 4e13 cell/m3 soil, 8e12 cell/m2, 8e13 cell/g micr C, 8e10 cell/mg micr C
    ! (3) The MM constants are listed in grams/m2. Presumably this is in gN/m2, but I am
    ! not 100% sure. If so, I convert to gN/kgSoil. So conversion of the plant enzymes is:
    ! [E_plant_sum] = (8e5) * (1e6 mg C fine root / kg C fine root) * broot * nplant / soildepth / bulkden / AvoNum * (14 gN/mol)  
    
    ! (4) Going back to the first alternative,
    ! (g plant enz C /g plant fine root C) = 1e-5
    ! (mol plant enz C /g plant fine root C) = 1e-5 / 12 [g C / mol C]
    ! (mol plant enz N /g plant fine root C) = 1e-5 / 12 / C:N_plant_enz
    ! (g plant enz N /g plant fine root C) = 1e-5 / 12 / C:N_plant_enz * 14
    ! (g plant enz N /kg plant fine root C) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000
    ! (g plant enz N /plant) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot
    ! (g plant enz N /m2) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant
    ! (g plant enz N /m3) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant / soildepth
    ! (g plant enz N /kg soil) = 1e-5 / 12 / C:N_plant_enz * 14 * 1000 * broot * nplant / soildepth / bulkden
    ! Still, this doesn't make sense. We need the number of moles of transporter enzyme, and then
    ! multiply that by the MM. 
    ! moles enzyme = (# enzyme) / AvoNum = (total mg C in enzyme) / (mg C / enzyme) / AvoNum
    ! moles enzyme = (0.05 * amb_c) * (g soil) / (mgC/enzyme) / AvoNum
    ! (moles enzyme / gSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum
    ! (gN / gSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14
    ! (gN / kgSoil) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000
    ! (gN / m3) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000 * 750
    ! (gN / m2) = (0.05 * amb_c) / (mgC/enzyme) / AvoNum * 14 * 1000 * 750 * 0.2
    ! amb_c = 0.1 gC/m2 / (0.2 m) / (750 kg/m3) = 7e-4 mgC/gSoil
    ! (gN / kgSoil) = 1e-22 / (mgC/enzyme)

    ! Scaling factor (gCenz / gCrootbiomass) * 
    ! fine root biomass (kgCfineroot/m2) *  ! kgCenz / m2
    ! 1000gCfineroot / kgCfineroot  /       ! gCenz / m2
    ! characteristic soil depth (m) /       ! gCenz / m3
    ! soil bulk density (g/m3) *            ! gCenz / gsoil
    ! 1000mgCenz / 1 gCenz     /            ! mgCenz / gsoil
    ! enzC:Nratio                           ! mgNenz / gsoil

    ! Ecosystem-level nutrient limitation
    if(nutr_limit_scheme == 2)then
       nstorage_total = 0.
       nstorage_min_total = 0.
       pstorage_total = 0.
       pstorage_min_total = 0.
       nplant_total = 0.
       do ico = 1, ncohorts
          nstorage_total = nstorage_total + nstorage(ico) * nplant(ico) 
          nstorage_min_total = nstorage_min_total + nstorage_min(ico) *  &
               nplant(ico) 
          pstorage_total = pstorage_total + pstorage(ico) * nplant(ico) 
          pstorage_min_total = pstorage_min_total + pstorage_min(ico) *  &
               nplant(ico) 
          nplant_total = nplant_total + nplant(ico)
       enddo
       nstorage_total = nstorage_total / nplant_total
       nstorage_min_total = nstorage_min_total / nplant_total
       pstorage_total = pstorage_total / nplant_total
       pstorage_min_total = pstorage_min_total / nplant_total
       n_limit_factor = 1.0 - ((nstorage_total - nstorage_min_total) / &
            (nstorage_min_total * (nstorage_max_factor - 1.)))**consts%wexp
       p_limit_factor = 1.0 - ((pstorage_total - pstorage_min_total) / &
            (pstorage_min_total * (pstorage_max_factor - 1.)))**consts%wexp
       n_limit_factor = max(min(1., n_limit_factor), 0.)
       p_limit_factor = max(min(1., p_limit_factor), 0.)

       do ico = 1, ncohorts
          broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
               (-slz(nlsl)/(-slz(krdepth(ico)))))
          enz_plant_n(pft(ico)) = enz_plant_n(pft(ico)) + &
               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   &
               consts%eff_soil_depth / slden / consts%enz_plant_c2n * &
               n_limit_factor
          enz_plant_p(pft(ico)) = enz_plant_p(pft(ico)) + &
               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   & 
               consts%eff_soil_depth / slden / consts%enz_plant_c2p * &
               p_limit_factor
       enddo
    elseif(nutr_limit_scheme == 1)then
       do ico = 1, ncohorts
             broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
                  (-slz(nlsl)/(-slz(krdepth(ico)))))

!          if(lai(ico) > 0.)then
!             broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
!                  (-slz(nlsl)/(-slz(krdepth(ico))))) * min(1., sum(water_supply) / (0.004 * lai(ico)))
!             print*,1,lai(ico),sum(water_supply),sum(water_supply)/(0.004*lai(ico))
!          else
!             broot_nl = 0.
!          endif
          n_limit_factor = 1.0 - ((nstorage(ico) - nstorage_min(ico)) / &
               (nstorage_min(ico)*(nstorage_max_factor - 1.)))**consts%wexp
          n_limit_factor = max(min(1., n_limit_factor), 0.)
          p_limit_factor = 1.0 - ((pstorage(ico) - pstorage_min(ico)) / &
               (pstorage_min(ico)*(pstorage_max_factor - 1.)))**consts%wexp
          p_limit_factor = max(min(1., p_limit_factor), 0.)
          
!          enz_plant_n(pft(ico)) = enz_plant_n(pft(ico)) + &
!               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   &
!               consts%eff_soil_depth / slden / consts%enz_plant_c2n !* &
!!!!!               n_limit_factor
!          enz_plant_p(pft(ico)) = enz_plant_p(pft(ico)) + &
!               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   & 
!               consts%eff_soil_depth / slden / consts%enz_plant_c2p !* &
!!!!!               p_limit_factor

          enz_plant_n(pft(ico)) = enz_plant_n(pft(ico)) + &
               consts%enz2biomass_plant * nplant(ico) * broot_nl /   &
               consts%eff_soil_depth / slden * 14.
          enz_plant_p(pft(ico)) = enz_plant_p(pft(ico)) + &
               consts%enz2biomass_plant * nplant(ico) * broot_nl /   & 
               consts%eff_soil_depth / slden * 31.


!          if(lai(ico) > 0.)then
!             print*,1,ico,water_supply(ico),lai(ico),water_supply(ico)/(water_supply_scale*lai(ico))
!             transp_fact = water_supply(ico)/(water_supply_scale*lai(ico))
!             if(transp_fact_max < transp_fact)print*,'new max transp_fact',transp_fact
!             transp_fact_max = max(transp_fact_max,transp_fact)
!             transp_fact = min(1., transp_fact)
             transp_fact = 1.
             vnh4up_plant(pft(ico)) = vnh4up_plant(pft(ico)) + consts%vnh4up_plant_base *  &
                  consts%enz2biomass_plant * nplant(ico) * broot_nl /   & 
                  consts%eff_soil_depth / slden * 14. * &
                  n_limit_factor * transp_fact
             vno3up_plant(pft(ico)) = vno3up_plant(pft(ico)) + consts%vno3up_plant_base *  &
                  consts%enz2biomass_plant * nplant(ico) * broot_nl /   & 
                  consts%eff_soil_depth / slden * 14. * &
                  n_limit_factor * transp_fact
             vpup_plant(pft(ico)) = vpup_plant(pft(ico)) + consts%vpup_plant_base *  &
                  consts%enz2biomass_plant * nplant(ico) * broot_nl /   & 
                  consts%eff_soil_depth / slden * 31. * &
                  p_limit_factor * transp_fact
!          endif

!          vnh4up_plant(pft(ico)) = vnh4up_plant(pft(ico)) + consts%vnh4up_plant_base *  &
!               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   & 
!               consts%eff_soil_depth / slden / consts%enz_plant_c2n * &
!               n_limit_factor
!          vno3up_plant(pft(ico)) = vno3up_plant(pft(ico)) + consts%vno3up_plant_base *  &
!               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   & 
!               consts%eff_soil_depth / slden / consts%enz_plant_c2n * &
!               n_limit_factor
!          vpup_plant(pft(ico)) = vpup_plant(pft(ico)) + consts%vpup_plant_base *  &
!               consts%enz2biomass_plant * nplant(ico) * broot_nl * 1000. /   & 
!               consts%eff_soil_depth / slden / consts%enz_plant_c2p * &
!               p_limit_factor
          
       enddo
    endif

    ! These rates from Zhu et al. (2016) Biogeosciences, Table 2 footnotes.
    ! Units: 1/s
!    vnh4up_plant = consts%vnh4up_plant_base
!    vno3up_plant = consts%vno3up_plant_base
!    vpup_plant = consts%vpup_plant_base

!    vnh4up_plant(24) = vnh4up_plant(24) * 1000. / 120.
!    vno3up_plant(24) = vno3up_plant(24) * 1000. / 2.
!    vpup_plant(24) = vpup_plant(24) * 400. / 12.

!    vnh4up_plant(36) = vnh4up_plant(36) * 1000. / 120.
!    vno3up_plant(36) = vno3up_plant(36) * 1000. / 2.
!    vpup_plant(36) = vpup_plant(36) * 400. / 12.

    return
  end subroutine som_plant_enzymes

  subroutine som_plant_feedback(nh4_plant, no3_plant, p_plant, slden,  &
       consts, ncohorts, nstorage, pstorage, nstorage_min, pstorage_min, &
       nplant, broot, rh, co2_lost, pft, krdepth, water_supply_layer_frac, lai)
    use ed_misc_coms, only: dtlsm
    use mend_consts_coms, only: decomp_consts
    use nutrient_constants, only: nstorage_max_factor, pstorage_max_factor, nlsl
    use ed_max_dims, only: n_pft
    use pft_coms, only: root_beta
    use soil_coms, only: slz, nzg
    implicit none

    real :: broot_nl
    integer, dimension(ncohorts), intent(in) :: pft
    integer, dimension(ncohorts), intent(in) :: krdepth
    integer, intent(in) :: ncohorts
    real, intent(inout), dimension(n_pft) :: nh4_plant
    real, intent(inout), dimension(n_pft) :: no3_plant
    real, intent(inout), dimension(n_pft) :: p_plant
    real, intent(in) :: slden
    type(decomp_consts) :: consts
    real, dimension(n_pft) :: plant_n_uptake
    real, dimension(n_pft) :: plant_p_uptake
    real, dimension(n_pft) :: total_n_activity
    real, dimension(n_pft) :: total_n_activity_nolimit
    real, dimension(n_pft) :: total_p_activity
    real, dimension(n_pft) :: total_p_activity_nolimit
    integer :: ico
    real, dimension(ncohorts) :: n_limit_factor
    real, dimension(ncohorts) :: p_limit_factor
    real, intent(inout), dimension(ncohorts) :: nstorage
    real, intent(inout), dimension(ncohorts) :: pstorage
    real, intent(in), dimension(ncohorts) :: lai
    real, intent(in), dimension(nzg,ncohorts) :: water_supply_layer_frac
    real, dimension(ncohorts) :: water_supply
    real, intent(in), dimension(ncohorts) :: nstorage_min
    real, intent(in), dimension(ncohorts) :: pstorage_min
    real, intent(in), dimension(ncohorts) :: nplant
    real, intent(in), dimension(ncohorts) :: broot
    real, intent(inout) :: rh
    real :: co2_lost_units
    real, intent(in) :: co2_lost
    real :: nstorage_total
    real :: nstorage_min_total
    real :: pstorage_total
    real :: pstorage_min_total
    real :: nplant_total
    integer, parameter :: nutr_limit_scheme = 1
    real :: total_nlim
    real :: total_plim
    real :: transp_fact

    do ico =1, ncohorts
       water_supply(ico) = sum(water_supply_layer_frac(nlsl:nzg,ico))
    enddo

    ! kgN/m2
    plant_n_uptake = (nh4_plant + no3_plant) * 1.0e-6  * &
         1000. * slden * consts%eff_soil_depth  
    ! kgP/m2
    plant_p_uptake = p_plant * 1.0e-6  * &
         1000. * slden * consts%eff_soil_depth  

!print*
!print*,plant_p_uptake(24:29)

    if(nutr_limit_scheme == 1)then
       total_n_activity = 0.
       total_p_activity = 0.
       total_n_activity_nolimit = 0.
       total_p_activity_nolimit = 0.
       do ico = 1, ncohorts
          broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
               (-slz(nlsl)/(-slz(krdepth(ico)))))

          n_limit_factor(ico) = 1.0 - ((nstorage(ico) - nstorage_min(ico)) / &
               (nstorage_min(ico)*(nstorage_max_factor - 1.)))**consts%wexp
          n_limit_factor(ico) = max(min(1., n_limit_factor(ico)), 0.)
!BUG          n_limit_factor = max(min(1., n_limit_factor), 0.)
          p_limit_factor(ico) = 1.0 - ((pstorage(ico) - pstorage_min(ico)) / &
               (pstorage_min(ico)*(pstorage_max_factor - 1.)))**consts%wexp
          p_limit_factor(ico) = max(min(1., p_limit_factor(ico)), 0.)
!BUG          p_limit_factor = max(min(1., p_limit_factor), 0.)

!          if(lai(ico) > 0.)then
!             transp_fact = min(1.,water_supply(ico) / (water_supply_scale * lai(ico)))
!          else
!             transp_fact = 0.
!          endif
          transp_fact = 1.

!          print*,2,ico,water_supply(ico),lai(ico),water_supply(ico)/(water_supply_scale*lai(ico))

          total_n_activity(pft(ico)) = total_n_activity(pft(ico)) + nplant(ico) *   &
               broot_nl * n_limit_factor(ico) * transp_fact
          total_p_activity(pft(ico)) = total_p_activity(pft(ico)) + nplant(ico) *   &
               broot_nl * p_limit_factor(ico) * transp_fact
          total_n_activity_nolimit(pft(ico)) = total_n_activity_nolimit(pft(ico)) + nplant(ico) *   &
               broot_nl
          total_p_activity_nolimit(pft(ico)) = total_p_activity_nolimit(pft(ico)) + nplant(ico) *   &
               broot_nl
       enddo
!print*,total_p_activity(24:29)       
       do ico = 1, ncohorts
          broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
               (-slz(nlsl)/(-slz(krdepth(ico)))))

!          if(lai(ico) > 0.)then
!             transp_fact = min(1.,water_supply(ico) / (water_supply_scale * lai(ico)))
!          else
!             transp_fact = 0.
!          endif
          transp_fact = 1.

          if(total_n_activity(pft(ico)) > 1.e-30)then
             nstorage(ico) = nstorage(ico) + plant_n_uptake(pft(ico)) * broot_nl /  &
                  total_n_activity(pft(ico)) * n_limit_factor(ico) * transp_fact
!          else
!             nstorage(ico) = nstorage(ico) + plant_n_uptake(pft(ico)) * broot_nl /  &
!                  total_n_activity_nolimit(pft(ico))
          endif

          if(total_p_activity(pft(ico)) > 1.e-30)then
             pstorage(ico) = pstorage(ico) + plant_p_uptake(pft(ico)) * broot_nl /  &
                  total_p_activity(pft(ico)) * p_limit_factor(ico) * transp_fact
!          else
!             pstorage(ico) = pstorage(ico) + plant_p_uptake(pft(ico)) * broot_nl /  &
!                  total_p_activity_nolimit(pft(ico))
          endif
       enddo
    elseif(nutr_limit_scheme == 2)then
       total_n_activity = 0.
       total_p_activity = 0.
       do ico = 1, ncohorts
          broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
               (-slz(nlsl)/(-slz(krdepth(ico)))))
          total_n_activity(pft(ico)) = total_n_activity(pft(ico)) + nplant(ico) *   &
               broot_nl
          total_p_activity(pft(ico)) = total_p_activity(pft(ico)) + nplant(ico) *   &
               broot_nl
       enddo
       
       do ico = 1, ncohorts
          broot_nl = broot(ico) * (1. - root_beta(pft(ico))**  &
               (-slz(nlsl)/(-slz(krdepth(ico)))))
          if(total_n_activity(pft(ico)) > 1.e-30)  &
               nstorage(ico) = nstorage(ico) + plant_n_uptake(pft(ico)) * broot_nl /  &
               total_n_activity(pft(ico))
          if(total_p_activity(pft(ico)) > 1.e-30)  &
               pstorage(ico) = pstorage(ico) + plant_p_uptake(pft(ico)) * broot_nl /  &
               total_p_activity(pft(ico))
       enddo
    endif

    ! gC/kgSoil
    co2_lost_units = co2_lost
    ! gC/m3Soil
    co2_lost_units = co2_lost_units * slden
    ! gC/m2Soil
    co2_lost_units = co2_lost_units * consts%eff_soil_depth
    ! molC/m2Soil
    co2_lost_units = co2_lost_units / 12.
    ! umolC/m2Soil
    co2_lost_units = co2_lost_units * 1.0e6
    ! umolC/m2Soil/s
    co2_lost_units = co2_lost_units / dtlsm

    ! Averaging this over the day.
    rh = rh + co2_lost_units * dtlsm / 86400.

!    nh4_plant = 0.
!    no3_plant = 0.
!    p_plant = 0.
!    co2_lost = 0.

    return
  end subroutine som_plant_feedback

  subroutine litt_plant_enzymes(enz_plant_n, enz_plant_p, vnh4up_plant,  &
       vno3up_plant, vpup_plant)
    use ed_max_dims, only: n_pft
    implicit none
    
    real, intent(out), dimension(n_pft) :: enz_plant_n
    real, intent(out), dimension(n_pft) :: enz_plant_p
    real, intent(out), dimension(n_pft) :: vnh4up_plant
    real, intent(out), dimension(n_pft) :: vno3up_plant
    real, intent(out), dimension(n_pft) :: vpup_plant

    enz_plant_n = 0.
    enz_plant_p = 0.

    vnh4up_plant = 0.
    vno3up_plant = 0.
    vpup_plant = 0.

    return
  end subroutine litt_plant_enzymes

  subroutine wood_plant_enzymes(iwood, enz_plant_n, enz_plant_p,  &
       vnh4up_plant, vno3up_plant, vpup_plant)
    use ed_max_dims, only: n_pft
    implicit none
    
    integer, intent(in) :: iwood
    real, intent(out), dimension(n_pft) :: enz_plant_n
    real, intent(out), dimension(n_pft) :: enz_plant_p
    real, intent(out), dimension(n_pft) :: vnh4up_plant
    real, intent(out), dimension(n_pft) :: vno3up_plant
    real, intent(out), dimension(n_pft) :: vpup_plant

    enz_plant_n = 0.
    enz_plant_p = 0.

    vnh4up_plant = 0.
    vno3up_plant = 0.
    vpup_plant = 0.

    return
  end subroutine wood_plant_enzymes


end Module mend_plant
