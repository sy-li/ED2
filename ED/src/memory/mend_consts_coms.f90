Module mend_consts_coms

  use mend_state_vars, only: nwood, npom
  use ed_max_dims, only: n_pft

  implicit none

  type decomp_consts

     character(len=8) :: type
     real, dimension(npom) :: vmax_decomp_pom
     real :: vmax_decomp_mom
     real, dimension(npom) :: kmm_decomp_pom
     real :: kmm_decomp_mom
     real :: decomp_dissolved_frac
     real :: binding_affin
     real :: desorp_rate
     real :: sorp_cap
     real :: growth_yield
     real :: kmm_amb
     real :: vmax_amb
     real :: spec_maint_rate
     real :: amb_size_scale
     real :: dorm_resp_frac
     real :: prod_frac_enz_pom
     real :: prod_frac_enz_mom
     real :: prod_frac_enz_ptase
     real :: turnover_frac_dom
     real, dimension(npom) :: turnover_frac_pom
     real :: mb_cn_max
     real :: mb_cn_min
     real :: enz_turnover_rate
     real, dimension(npom) :: enz_pom_c2n
     real :: enz_mom_c2n
     real :: enz_ptase_c2n
     real, dimension(npom) :: enz_pom_c2p
     real :: enz_mom_c2p
     real :: enz_ptase_c2p
     real :: mb_cp_max
     real :: mb_cp_min
     real :: v_nitr
     real :: v_denitr
     real :: kmm_decomp_ptase_dom
     real, dimension(npom) :: kmm_decomp_ptase_pom
     real :: vmax_decomp_ptase_dom
     real, dimension(npom) :: vmax_decomp_ptase_pom
     real :: wexp
     real :: ygn_base
     real :: ygp_base
     real :: vnup_amb_base
     real :: vpup_amb_base
     real :: ksnh4_amb
     real :: ksno3_amb
     real, dimension(n_pft) :: ksnh4_plant
     real, dimension(n_pft) :: ksno3_plant
     real :: kspsol_amb
     real, dimension(n_pft) :: kspsol_plant
     real :: stoich_adj_rate
     real :: plab_max
     real :: kmm_plang
     real :: occlu_rate
     real :: weath_rate
     real :: enz2biomass_plant
     real :: eff_soil_depth
     real :: enz_plant_c2n
     real :: enz_plant_c2p
     real :: vnh4up_plant_base
     real :: vno3up_plant_base
     real :: vpup_plant_base
     real :: dom_leach_effic
     real :: no3_leach_effic
     real :: psol_leach_effic
     real :: amb_turnover_rate
     real :: sorp_rate
     real :: swp_a2d
     real :: wdorm
     real :: swp_d2a
     real :: vmax_a2d
     real :: vmax_d2a
     real :: p_adsorp_rate
     real :: p_desorp_rate
     real :: micr_bio2enz

     ! Not currently used
     real :: vnup_nit
     real :: vnup_den
     real :: ksnh4_nit
     real :: enz_nit_n
     real :: ksno3_den
     real :: enz_den_n
     real :: kspsol_lab
     real :: vmax_surf
     real :: frag_rate

!!!!!!!!!!!!
     

  end type decomp_consts

  type(decomp_consts) :: som_consts
  type(decomp_consts) :: litt_consts
  type(decomp_consts), dimension(nwood) :: wood_consts

  type(decomp_consts) :: som_consts_base
  type(decomp_consts) :: litt_consts_base
  type(decomp_consts), dimension(nwood) :: wood_consts_base

Contains

  subroutine mend_init_consts(sens_params)

    implicit none

    integer :: iwood
    integer, intent(in) :: sens_params
    integer, parameter :: sens_type = 2

    if(npom /= 2)then
       print*,'npom not set to 2',npom
       print*, 'npom=2 assumed in consts_coms.f90 initialization'
       stop
    endif

    ! SOM

    som_consts%type = 'som'   ! nt

    ! Wang et al. (2015), ISME J, eyeballed from graphs

    som_consts%growth_yield = 0.45 ! dimensionless; t
    som_consts%kmm_amb = 0.45  ! mgC/gSOIL; t
    som_consts%vmax_amb = 0.025 / 3600. ! 1/s; t
    som_consts%spec_maint_rate = 0.01 / 3600. ! 1/s; t
    som_consts%prod_frac_enz_pom = 0.07 ! unitless; nt
    som_consts%prod_frac_enz_mom = 0.06  ! unitless; nt
    som_consts%enz_turnover_rate = 0.0025 / 3600.  ! 1/s; nt

    som_consts%vmax_decomp_pom(1) = 1.2 / 3600. / 1.2 * 0.7 ! 1/s; t
    som_consts%vmax_decomp_pom(2) = 40. / 3600. /40. *47. ! 1/s; t
    som_consts%vmax_decomp_mom = 1.0 / 3600.  * 0.3 ! 1/s; t

! This set swapped after run n2 on 11/22/2017.
    som_consts%kmm_decomp_pom(1) = 80. ! mgC/gsoil; t 
    som_consts%kmm_decomp_pom(2) = 1. ! mgC/gsoil; t 
    som_consts%kmm_decomp_mom = 780. ! mgC/gsoil; t
!    som_consts%kmm_decomp_pom(1) = 180. ! mgC/gsoil; t 
!    som_consts%kmm_decomp_pom(2) = 30. ! mgC/gsoil; t 
!    som_consts%kmm_decomp_mom = 1500. ! mgC/gsoil; t

    som_consts%binding_affin = 8.  ! mgC/gSOIL; nt
    som_consts%sorp_cap = 4.2  ! mgC/gSOIL; nt
    som_consts%desorp_rate = 0.003 / 3600.  ! 1/s; t
    som_consts%dorm_resp_frac = 0.001  ! unitless; nt
    som_consts%decomp_dissolved_frac = 0.9  ! dimensionless; nt
    som_consts%turnover_frac_dom = 0.5   ! unitless; nt
    som_consts%turnover_frac_pom(1) = 0.5   ! unitless; nt
    som_consts%turnover_frac_pom(2) = 0.0  ! unitless; nt

    som_consts%sorp_rate = som_consts%binding_affin *  &
         som_consts%desorp_rate  ! 1/s; t
    som_consts%amb_turnover_rate = som_consts%spec_maint_rate * &
         (1. - som_consts%prod_frac_enz_pom -   &
         som_consts%prod_frac_enz_mom)  !  1/s; t
    som_consts%vmax_a2d = som_consts%spec_maint_rate !  1/s; t
    som_consts%vmax_d2a = som_consts%spec_maint_rate !  1/s; t

    ! Sinsabaugh et al. (2009), Nature 462, 795-798. Ecoezymatic stoichiometry
    ! of microbial organic nutrient acquisition in soil and sediment.
    som_consts%mb_cp_max = 120. ! mgC/mgP; nt
    som_consts%mb_cp_min = 30.  ! mgC/mgP; nt

    ! Numbers that I estimated.
    som_consts%amb_size_scale = 0.05  ! mgC/gSOIL; nt
    som_consts%enz_pom_c2p(1) = 60.  ! mgC/mgP. Typical microbial C:P.; nt
    som_consts%enz_pom_c2p(2) = 60.  ! mgC/mgP; nt
    som_consts%enz_mom_c2p = 60.  ! mgC/mgP; nt
    som_consts%ygn_base = 0.9  ! dimensionless; nt; G has different number
    som_consts%ygp_base = 1. ! dimensionless; nt
    som_consts%stoich_adj_rate = 1.0 / 3600. ! 1/s; nt
!    som_consts%stoich_adj_rate = 1.0e-6 / 3600. ! 1/s; nt
!standard    som_consts%stoich_adj_rate = 0.00001 * 0.1 / 3600. ! 1/s; nt
    som_consts%eff_soil_depth = 0.2  ! m; nt
    som_consts%enz_plant_c2n = 3.  ! gC/gN  ! assumed similar to microbes; nt
    som_consts%enz_plant_c2p = 60.  ! gC/gP ! assumed similar to microbes; nt
    som_consts%dom_leach_effic = 0.   ! dimensionless; nt
    som_consts%no3_leach_effic = 1.   ! dimensionless; nt
    som_consts%psol_leach_effic = 1.   ! dimensionless; nt
   ! The following deal with phosphatase.  Turned off.
    som_consts%prod_frac_enz_ptase = 0.035  ! unitless; nt
    som_consts%enz_ptase_c2n = 3.  ! mgC/mgN; nt
    som_consts%enz_ptase_c2p = 60.  ! mgC/mgP; nt
    som_consts%kmm_decomp_ptase_pom(1) = 80. ! mgP/gSOIL; ADD T
    som_consts%kmm_decomp_ptase_pom(2) = 1. ! mgP/gSOIL; ADD T
    som_consts%kmm_decomp_ptase_dom = 780. ! mgP/gSOIL; ADD T
!    som_consts%kmm_decomp_ptase_pom(1) = 180.  ! mgP/gSOIL; ADD T
!    som_consts%kmm_decomp_ptase_pom(2) = 30.  ! mgP/gSOIL; ADD T
!    som_consts%kmm_decomp_ptase_dom = 1500.  ! mgP/gSOIL; ADD T

!    som_consts%vmax_decomp_ptase_pom(1) = 0.7 / 3600. ! 1/s; ADD T
!    som_consts%vmax_decomp_ptase_pom(2) = 47. / 3600. ! 1/s; ADD T
!    som_consts%vmax_decomp_ptase_dom = 0.  ! 1/s
    som_consts%vmax_decomp_ptase_pom(1) = 1. * 0.7 / 3600. ! 1/s; ADD T
    som_consts%vmax_decomp_ptase_pom(2) = 1. * 47. / 3600. ! 1/s; ADD T
    som_consts%vmax_decomp_ptase_dom = 0.  ! 1/s
!    som_consts%vmax_decomp_ptase_pom(1) = 0.0 * 0.7 / 3600. ! 1/s; ADD T
!    som_consts%vmax_decomp_ptase_pom(2) = 0.0 * 47. / 3600. ! 1/s; ADD T
!    som_consts%vmax_decomp_ptase_dom = 0.  ! 1/s

    ! Taken from Gangsheng's code.
    som_consts%mb_cn_max = 14.  ! mgC/mgN; nt
    som_consts%mb_cn_min = 3.  ! mgC/mgN; nt
    som_consts%enz_pom_c2n(1) = 3.  ! mgC/mgN; nt
    som_consts%enz_pom_c2n(2) = 3.  ! mgC/mgN; nt
    som_consts%enz_mom_c2n = 3.  ! mgC/mgN; nt
    som_consts%v_nitr = 3.0e-3 / 3600.  ! 1/s; t
    som_consts%v_denitr = 1.0e-2 / 3600.  * 0. ! 1/s; t
    som_consts%wexp=4.  ! dimensionless; nt
    som_consts%swp_a2d = -0.4; !  MPa; nt
    som_consts%wdorm = 4. !  dimensionless; nt
    som_consts%swp_d2a = -0.4*0.25  ! MPa; nt

    ! Taken from Zhu et al. (2016), Biogeosciences
    som_consts%vnup_amb_base = 1000. / 86400. ! 1/s; nt
    ! Gangsheng's estimate, below, is orders of magnitude smaller.
!    som_consts%vnup_amb_base = 0.01 / 3600. ! 1/s
    som_consts%vpup_amb_base = 400. / 86400. ! 1/s; t
    som_consts%vnh4up_plant_base = 120. / 86400.  ! 1/s; t
    som_consts%vno3up_plant_base = 2. / 86400.  ! 1/s; t
    som_consts%vpup_plant_base = 12. / 86400.  ! 1/s; t
! NCOM Experiment
!    som_consts%ksnh4_amb = 0.173 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%ksno3_amb = 0.173 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%ksno3_plant = 0.173 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%kspsol_amb = 0.11 / (som_consts%eff_soil_depth) ! mgP/gsoil; t
!    som_consts%ksnh4_nit = 0.173 !  gN/m2SOIL
!    som_consts%ksno3_den = 0.173 !  gN/m2SOIL
! Base values
    som_consts%ksnh4_amb = 0.071 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
    som_consts%ksno3_amb = 0.096 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
    som_consts%kspsol_amb = 0.037 / (som_consts%eff_soil_depth) ! mgP/gsoil; t
    som_consts%ksnh4_nit = 0.082 !  gN/m2SOIL
    som_consts%ksno3_den = 0.022 !  gN/m2SOIL
    som_consts%kspsol_lab = 200. / som_consts%eff_soil_depth !  gP/m2SOIL to gP/m3SOIL
    som_consts%kmm_plang = 200. / (som_consts%eff_soil_depth) !gP/m2 converted to mgP/g; nt
! End base values

    ! Gangsheng's numbers, below, are very different.
!    som_consts%ksnh4_amb = 0.005
!    som_consts%ksno3_amb = 0.009

    ! Density now assigned elsewhere.

    ! KSNH4_PLANT
    ! Zhu et al. (2016) says 0.173 gN/m2. Need to divide by soil depth (m) and density (kg/m3).
    ! Depth division is here; density division is site dependent and is done online.
    ! Actually, I should be dividing by Zhu's bulk density, I think. But they didn't list their
    ! bulk density. Ignore this for now.
    som_consts%ksnh4_plant = 0.173 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
    ! Similar concerns apply for ksno3_plant and kspsol_plant.
    som_consts%ksno3_plant =  0.085 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
    som_consts%kspsol_plant = 0.11 / (som_consts%eff_soil_depth) ! mgN/gsoil; t

    ! Mystery parameter. Should really be subject to sensitivity.
    som_consts%enz2biomass_plant = 0.0000125*1000.
    som_consts%micr_bio2enz = 1. * 0.001

! Override for n7
!    som_consts%ksnh4_plant(24) = 0.071 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%ksno3_plant(24) = 0.096 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%kspsol_plant(24) = 0.037 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%ksnh4_plant(36) = 0.071 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%ksno3_plant(36) = 0.096 / (som_consts%eff_soil_depth) ! mgN/gsoil; t
!    som_consts%kspsol_plant(36) = 0.037 / (som_consts%eff_soil_depth) ! mgN/gsoil; t

    ! Taken from Yang et al. (2014), Glob Change Biol
!    som_consts%plab_max = 10. / som_consts%eff_soil_depth !  gP/m3

    ! Density now assigned elsewhere.
!    som_consts%kmm_plang = 0.001 / (som_consts%eff_soil_depth) !gP/m2 converted to mgP/g; nt
    som_consts%occlu_rate = 1.0e-6 / 30.4 / 24. / 3600.  !  1/month to 1/s; nt
!    som_consts%weath_rate = 1.0e-3 / 30.4 / 24. / 3600. !  1/month to 1/s; nt
    som_consts%weath_rate = 1.0e-4 / 30.4 / 24. / 3600. !  1/month to 1/s; nt
!    som_consts%weath_rate = 0. * 1.0e-3 / 30.4 / 24. / 3600. !  1/month to 1/s; nt
    som_consts%p_adsorp_rate = 0.003 / 30.4 / 24/ 3600.
    som_consts%p_desorp_rate = 0.00022 / 30.4 / 24/ 3600.

    som_consts%weath_rate = 200. * 1.0e-4 / 30.4 ! gP/m3/day  

    if(sens_type == 2)then
       if(sens_params == 1)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 0.01
       elseif(sens_params == 2)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 0.03
       elseif(sens_params == 3)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 0.1
       elseif(sens_params == 4)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 0.3
       elseif(sens_params == 5)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 3.0
       elseif(sens_params == 6)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 10.0
       elseif(sens_params == 7)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 30.0
       elseif(sens_params == 8)then
          som_consts%enz2biomass_plant = 0.0000125*1000. * 100.
       elseif(sens_params == 9)then
          som_consts%micr_bio2enz = 0.01 * 0.001
       elseif(sens_params == 10)then
          som_consts%micr_bio2enz = 0.03 * 0.001
       elseif(sens_params == 11)then
          som_consts%micr_bio2enz = 0.1 * 0.001
       elseif(sens_params == 12)then
          som_consts%micr_bio2enz = 0.3 * 0.001
       elseif(sens_params == 13)then
          som_consts%micr_bio2enz = 3. * 0.001
       elseif(sens_params == 14)then
          som_consts%micr_bio2enz = 10. * 0.001
       elseif(sens_params == 15)then
          som_consts%micr_bio2enz = 30. * 0.001
       elseif(sens_params == 16)then
          som_consts%micr_bio2enz = 100. * 0.001
       else
          ! Remain at baseline.
       endif
    endif

    if(sens_type == 1)then
       if(sens_params == 0)then
          ! Keep things at baseline values
          !       som_consts%weath_rate = 200. * 1.0e-4 / 30.4 ! gP/m3/day  
          !       som_consts%vmax_decomp_ptase_pom(1) = 1. * 0.7 / 3600. ! 1/s; ADD T
          !       som_consts%vmax_decomp_ptase_pom(2) = 1. * 47. / 3600. ! 1/s; ADD T
          !       som_consts%dom_leach_effic = 0.   ! dimensionless; nt
          !       som_consts%psol_leach_effic = 1.   ! dimensionless; nt
       elseif(sens_params > 0 .and. sens_params <= 10)then
          som_consts%weath_rate = 1.0e-6 * 1000. / 0.5 / 365.25 *  &  ! gP/m3/day; min value
               real(sens_params-1) * 10. 
       elseif(sens_params > 10 .and. sens_params <= 20)then
          som_consts%vmax_decomp_ptase_pom(1) = 0.7 / 3600. * & ! 1/s
               real(sens_params-11) * 0.1
          som_consts%vmax_decomp_ptase_pom(2) = 47. / 3600. * & ! 1/s
               real(sens_params-11) * 0.1
       elseif(sens_params > 20 .and. sens_params <= 30)then
          if(sens_params == 21)som_consts%dom_leach_effic = 0.07
          if(sens_params == 22)som_consts%dom_leach_effic = 0.13
          if(sens_params == 23)som_consts%dom_leach_effic = 0.21
          if(sens_params == 24)som_consts%dom_leach_effic = 0.29
          if(sens_params == 25)som_consts%dom_leach_effic = 0.38
          if(sens_params == 26)som_consts%dom_leach_effic = 0.47
          if(sens_params == 27)som_consts%dom_leach_effic = 0.58
          if(sens_params == 28)som_consts%dom_leach_effic = 0.70
          if(sens_params == 29)som_consts%dom_leach_effic = 0.84
          if(sens_params == 30)som_consts%dom_leach_effic = 1.
       elseif(sens_params > 30 .and. sens_params <= 40)then
          if(sens_params == 31)som_consts%psol_leach_effic = 0. 
          if(sens_params == 32)som_consts%psol_leach_effic = 0.07
          if(sens_params == 33)som_consts%psol_leach_effic = 0.13
          if(sens_params == 34)som_consts%psol_leach_effic = 0.21
          if(sens_params == 35)som_consts%psol_leach_effic = 0.29
          if(sens_params == 36)som_consts%psol_leach_effic = 0.38
          if(sens_params == 37)som_consts%psol_leach_effic = 0.47
          if(sens_params == 38)som_consts%psol_leach_effic = 0.58
          if(sens_params == 39)som_consts%psol_leach_effic = 0.70
          if(sens_params == 40)som_consts%psol_leach_effic = 0.84
       endif
    endif
    
    ! Not used
    som_consts%vnup_den = 0.1/24./3600.  ! 1/s
    som_consts%vnup_nit = 0.1/24./3600.  ! 1/s
    som_consts%enz_nit_n = 1.2e-3 !  gN/m2SOIL
    som_consts%enz_den_n = 1.2e-3 !  gN/m2SOIL
! Base values
! End Base values
    som_consts%plab_max = 182. / som_consts%eff_soil_depth !  gP/m3
    som_consts%vmax_surf = 182. / som_consts%eff_soil_depth !  gP/m2SOIL to gP/m3SOIL
    som_consts%frag_rate = 0.

!!!!
     

    ! LITTER

    litt_consts%type = 'litter'
    litt_consts%vmax_decomp_pom(1) = 1.2 / 3600. ! 1/s
    litt_consts%vmax_decomp_pom(2) = 38. / 3600. ! 1/s
    litt_consts%vmax_decomp_mom = 0.  ! 1/s
    litt_consts%kmm_decomp_pom(1) = 5.0e3 ! gC/m2
    litt_consts%kmm_decomp_pom(2) = 5.0e3 ! gC/m2 
    litt_consts%kmm_decomp_mom = 500.  ! mgC/gSOIL
    litt_consts%decomp_dissolved_frac = 1.0
    litt_consts%binding_affin = 7.  ! mgC/gSOIL
    litt_consts%desorp_rate = 0.0  ! 1/s
    litt_consts%sorp_cap = 3.8  ! mgC/gSOIL
    litt_consts%growth_yield = 0.45  ! dimensionless
    litt_consts%kmm_amb = 0.4  * 160. ! mg/g converted to g/m2
    litt_consts%vmax_amb = 0.02 / 3600. ! 1/s
    litt_consts%spec_maint_rate = 0.01 / 3600.  ! 1/s
    litt_consts%amb_size_scale = 1.  ! gC/m2
    litt_consts%dorm_resp_frac = 0.001  ! dimensionless
    litt_consts%prod_frac_enz_pom = 0.13  ! dimensionless
    litt_consts%prod_frac_enz_mom = 0.0  ! dimensionless
    litt_consts%prod_frac_enz_ptase = 0.  ! dimensionless
    litt_consts%turnover_frac_dom = 0.5  ! dimensionless
    litt_consts%turnover_frac_pom(1) = 0.5  ! dimensionless
    litt_consts%turnover_frac_pom(2) = 0.0  ! dimensionless
    litt_consts%mb_cn_max = 25.  ! mgC/mgN
    litt_consts%mb_cn_min = 3.  ! mgC/mgN
    litt_consts%enz_turnover_rate = 0.0025 / 3600.  ! 1/s
    litt_consts%enz_pom_c2n(1) = 3. ! mgC/mgN
    litt_consts%enz_pom_c2n(2) = 3. ! mgC/mgN
    litt_consts%enz_mom_c2n = 3. ! mgC/mgN
    litt_consts%enz_ptase_c2n = 10. ! mgC/mgN
    litt_consts%enz_pom_c2p(1) = 300. ! mgC/mgP
    litt_consts%enz_pom_c2p(2) = 300. ! mgC/mgP
    litt_consts%enz_mom_c2p = 300. ! mgC/mgP
    litt_consts%enz_ptase_c2p = 1000. ! mgC/mgP
    litt_consts%mb_cp_max = 100. ! mgC/mgP
    litt_consts%mb_cp_min = 6. ! mgC/mgP
    litt_consts%v_nitr = 3.0e-3 /3600. ! 1/s
    litt_consts%v_denitr = 1.0e-2 / 3600. ! 1/s
    litt_consts%kmm_decomp_ptase_pom(1) = 0.06 * 160.  ! mg/g converted to g/m2
    litt_consts%kmm_decomp_ptase_pom(2) = 0.01 * 160.  ! mg/g converted to g/m2
    litt_consts%kmm_decomp_ptase_dom = 0.001 * 160.  ! mg/g converted to g/m2
    litt_consts%vmax_decomp_ptase_pom(1) = 0.  ! 1/s
    litt_consts%vmax_decomp_ptase_pom(2) = 0.  ! 1/s
    litt_consts%vmax_decomp_ptase_dom = 0.  ! 1/s
    litt_consts%wexp=4.  ! dimensionless
    litt_consts%ygn_base = 0.9 ! dimensionless
    litt_consts%ygp_base = 1. ! dimensionless
    litt_consts%vnup_amb_base = 0.01 / 3600. ! 1/s
    litt_consts%vpup_amb_base = 0.01 / 3600. ! 1/s
    litt_consts%ksnh4_amb = 0.005 * 160.  ! mg/g converted to g/m2
    litt_consts%ksno3_amb = 0.009 * 160.  ! mg/g converted to g/m2
    litt_consts%kspsol_amb = 0.005 * 160.  ! mg/g converted to g/m2
    litt_consts%ksnh4_plant = 1.2 * 0.2 !  gN/m3SOIL to g/m2
    litt_consts%ksno3_plant =  1.8 * 0.2 !  gN/m3SOIL to g/m2
    litt_consts%kspsol_plant = 0.11 * 0.2 !  gP/m3SOIL to g/m2
    litt_consts%stoich_adj_rate = 0.1/3600.  ! 1/s
    litt_consts%plab_max = 10.   ! gP/m2
    litt_consts%kmm_plang = 0.001  !gP/m2
    litt_consts%occlu_rate = 1.0e-6 / 30.4 / 24. / 3600.  !  1/month to 1/s
    litt_consts%weath_rate = 1.0e-3 / 30.4 / 24. / 3600. !  1/month to 1/s
    litt_consts%enz2biomass_plant = 0.0000125 ! gCenz/gCfineroot
    litt_consts%eff_soil_depth = 0.2  ! m
    litt_consts%enz_plant_c2n = 10.  ! gC/gN
    litt_consts%enz_plant_c2p = 300.  ! gC/gP
    litt_consts%vnh4up_plant_base = 120. / 86400.  ! 1/s
    litt_consts%vno3up_plant_base = 2. / 86400.  ! 1/s
    litt_consts%vpup_plant_base = 12. / 86400.  ! 1/s
    litt_consts%dom_leach_effic = 1.  ! dimensionless
    litt_consts%no3_leach_effic = 1.  ! dimensionless
    litt_consts%psol_leach_effic = 1.  ! dimensionless
    litt_consts%p_adsorp_rate = 0.003 / 30.4 / 24/ 3600.
    litt_consts%p_desorp_rate = 0.00022 / 30.4 / 24/ 3600.
    
    litt_consts%vnup_den = 0.1/24./3600. ! 1/s
    litt_consts%vnup_nit = 0.1/24./3600. ! 1/s
    litt_consts%enz_nit_n = 1.2e-3  !  gN/m2SOIL
    litt_consts%ksnh4_nit = 0.082  !  gN/m2SOIL
    litt_consts%enz_den_n = 1.2e-3  !  gN/m2SOIL
    litt_consts%ksno3_den = 0.022  !  gN/m2SOIL
    litt_consts%frag_rate = 0.

!!!!!!!!!!!!!!!!!
    
    do iwood = 1, nwood
       wood_consts(iwood) = litt_consts
       write(wood_consts(iwood)%type,'(a4,i1.1)')'wood',iwood
    enddo

    wood_consts(1)%frag_rate = 1. / (10. * 365.25 * 24. * 3600.)
    wood_consts(2)%frag_rate = 1. / (2. * 365.25 * 24. * 3600.)
    wood_consts(3)%frag_rate = 1. / (5. * 365.25 * 24. * 3600.)
    wood_consts(4)%frag_rate = 1. / (1. * 365.25 * 24. * 3600.)

    som_consts_base = som_consts
    litt_consts_base = litt_consts
    do iwood = 1, nwood
       wood_consts_base(iwood) = wood_consts(iwood)
    enddo


    return
  end subroutine mend_init_consts

end Module mend_consts_coms
