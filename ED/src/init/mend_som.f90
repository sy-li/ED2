Module mend_som
  implicit none

Contains

  subroutine som_init(npom,  &
       pom_c, dom_c, enz_pom_c, mom_c, qom_c, enz_mom_c, amb_c, dmb_c,  &
       pom_n, dom_n, enz_pom_n, mom_n, qom_n, enz_mom_n, amb_n, dmb_n,  &
       pom_p, dom_p, enz_pom_p, mom_p, qom_p, enz_mom_p, amb_p, dmb_p,  &
       co2_lost, nmin, nitr, nh4, no3, psol, plab, ngas_lost, enz_ptase_c, enz_ptase_n, &
       enz_ptase_p, nh4_dep, no3_dep, ppar_dep, nh4_plant, nh4_bnf, &
       no3_plant, c_leach, n_leach, p_leach, p_plant, pocc, psec, ppar,  &
       enz_plant_n, enz_plant_p, vnh4up_plant, vno3up_plant, &
       vpup_plant, &
       plant_input_C_pom, plant_input_N_pom, plant_input_P_pom, &
       plant_input_C_dom, plant_input_N_dom, plant_input_P_dom, &
       consts, bulk_den, soil_cpct, soil_som_c2n, soil_totp, soil_extrp)
    use ed_max_dims, only: n_pft
    use mend_consts_coms, only: decomp_consts
    implicit none

!    character(len=1), intent(in) :: fertex
    real, intent(in) :: soil_cpct
    real, intent(in) :: soil_som_c2n
    real, intent(in) :: soil_totp
    real, intent(in) :: soil_extrp
    type(decomp_consts) :: consts
    integer, intent(in) :: npom
    real, intent(out) :: dom_c
    real, intent(out) :: amb_c
    real, intent(out) :: dmb_c
    real, intent(out) :: mom_c
    real, intent(out) :: qom_c
    real, dimension(npom), intent(out) :: pom_c
    real, dimension(npom), intent(out) :: enz_pom_c
    real, intent(out) :: enz_mom_c
    real, intent(out) :: enz_ptase_c

    real, intent(out) :: dom_n
    real, intent(out) :: amb_n
    real, intent(out) :: dmb_n
    real, intent(out) :: mom_n
    real, intent(out) :: qom_n
    real, dimension(npom), intent(out) :: pom_n
    real, dimension(npom), intent(out) :: enz_pom_n
    real, intent(out) :: enz_mom_n
    real, intent(out) :: enz_ptase_n

    real, intent(out) :: dom_p
    real, intent(out) :: amb_p
    real, intent(out) :: dmb_p
    real, intent(out) :: mom_p
    real, intent(out) :: qom_p
    real, dimension(npom), intent(out) :: pom_p
    real, dimension(npom), intent(out) :: enz_pom_p
    real, intent(out) :: enz_mom_p
    real, intent(out) :: enz_ptase_p
    real, intent(out) :: no3_dep
    real, intent(out) :: nh4_dep

    real, intent(out) :: co2_lost
    real, intent(out) :: ngas_lost
    real, intent(out) :: nh4
    real, intent(out) :: no3
    real, intent(out) :: psol
    real, intent(out) :: plab
    real, intent(out) :: pocc
    real, intent(out) :: psec
    real, intent(out) :: ppar
    real, intent(out) :: ppar_dep

    real, intent(out), dimension(n_pft) :: nh4_plant
    real, intent(out) :: nh4_bnf
    real, intent(out), dimension(n_pft) :: no3_plant
    real, intent(out) :: nmin
    real, intent(out) :: nitr
    real, intent(out) :: c_leach
    real, intent(out) :: n_leach
    real, intent(out) :: p_leach
    real, intent(out), dimension(n_pft) :: p_plant
    real, intent(out), dimension(n_pft) :: enz_plant_n
    real, intent(out), dimension(n_pft) :: enz_plant_p
    real, intent(out), dimension(n_pft) :: vnh4up_plant
    real, intent(out), dimension(n_pft) :: vno3up_plant
    real, intent(out), dimension(n_pft) :: vpup_plant
    real, intent(out), dimension(npom) :: plant_input_C_pom
    real, intent(out), dimension(npom) :: plant_input_N_pom
    real, intent(out), dimension(npom) :: plant_input_P_pom
    real :: my_total, my_b, my_c
    real, intent(out) :: plant_input_C_dom
    real, intent(out) :: plant_input_N_dom
    real, intent(out) :: plant_input_P_dom
    real, intent(in) :: bulk_den
    real :: cpct
    real :: c2n
    real :: totp
    real :: extrp
    real :: c2p

    cpct = soil_cpct
    c2n = soil_som_c2n
    totp = soil_totp
    extrp = soil_extrp

!    if(fertex == 'A')then
!       cpct = 50.4
    !    c2n = 11.2
    !    totp = 395.33
    !    extrp = 7.872
    ! elseif(fertex == 'B')then
    !    cpct = 59.6
    !    c2n = 10.78
    !    totp = 550.93
    !    extrp = 13.02
    ! elseif(fertex == 'C')then
    !    cpct = 51.8
    !    c2n = 10.8
    !    totp = 614.83
    !    extrp = 19.474
    ! elseif(fertex == 'D')then
    !    cpct = 43.3
    !    c2n = 10.78
    !    totp = 452.5
    !    extrp = 7.159
    ! elseif(fertex == 'E')then
    !    cpct = 48.1
    !    c2n = 10.66
    !    totp = 441.03
    !    extrp = 6.12
    ! elseif(fertex == 'F')then
    !    cpct = 58.7
    !    c2n = 11.11
    !    totp = 607.33
    !    extrp = 10.011
    ! elseif(fertex == 'G')then
    !    cpct = 60.9
    !    c2n = 11.38
    !    totp = 638.59
    !    extrp = 13.32
    ! elseif(fertex == 'H')then
    !    cpct = 39.3
    !    c2n = 11.3
    !    totp = 540.37
    !    extrp = 6.087
    ! elseif(fertex == 'I')then
    !    cpct = 49.7
    !    c2n = 10.73
    !    totp = 459.95
    !    extrp = 8.439
    ! elseif(fertex == 'J')then
    !    cpct = 40.1
    !    c2n = 10.61
    !    totp = 442.44
    !    extrp = 6.507
    ! elseif(fertex == 'K')then
    !    cpct = 38.9
    !    c2n = 10.77
    !    totp = 596.51
    !    extrp = 6.746
    ! elseif(fertex == 'L')then
    !    cpct = 42.2
    !    c2n = 10.9
    !    totp = 646.95
    !    extrp = 9.963
    ! elseif(fertex == 'M')then
    !    cpct = 43.7
    !    c2n = 10.64
    !    totp = 402.22
    !    extrp = 6.38
    ! elseif(fertex == 'N')then
    !    cpct = 43.3
    !    c2n = 10.72
    !    totp = 551.04
    !    extrp = 7.126
    ! elseif(fertex == 'O')then
    !    cpct = 36.5
    !    c2n = 10.67
    !    totp = 525.29
    !    extrp = 7.047
    ! elseif(fertex == 'P')then
    !    cpct = 35.8
    !    c2n = 10.77
    !    totp = 518.76
    !    extrp = 6.818
    ! endif

!    cpct = 37.9
!    c2n = 10.8
!    extrp = 0.0455*1000.
!    totp = 150. * 2. / bulk_den * 1000. + extrp

!    c2p = 60.
!    c2p = consts%mb_cp_max
! Tipping et al. (2016) Biogeochemistry
    c2p = (cpct*0.1) / (0.012*(cpct*0.1)**0.57)
!    c2p = 800.

! density (g/cm3):
! SROAK:0.87  SRTDF:0.90  PVTDF:0.77
! %C: 2.65  3.79  4.30
! So, for SRTDF, SOM_C=37.9
! %N: 0.22  0.35  0.44
    pom_c(1) = 7./31. * cpct
    pom_n(1) = pom_c(1) / c2n
    pom_p(1) = pom_c(1) / c2p

    pom_c(2) = 2.5/31. * cpct
    pom_n(2) = pom_c(2) / c2n
    pom_p(2) = pom_c(2) / c2p

    dom_c = 0.15/31. * cpct
    dom_n = dom_c / c2n
    dom_p = dom_c / c2p

! a1
!    enz_pom_c(1) = 0.001 / 31. * cpct
! a2-4
!    enz_pom_c(1) = 0.025
! a5
!    enz_pom_c(1) = 0.0042
! a6
!    enz_pom_c(1) = 0.022
! a7
!    enz_pom_c(1) = 0.098
! a8
!    enz_pom_c(1) = 0.0026
! a9
    enz_pom_c(1) = 0.014
! a10
!    enz_pom_c(1) = 0.066

    enz_pom_n(1) = enz_pom_c(1) / consts%enz_pom_c2n(1)
    enz_pom_p(1) = enz_pom_c(1) / consts%enz_pom_c2p(1)

! a1
!    enz_pom_c(2) = 0.0004 / 31. * cpct
! a2-4
!    enz_pom_c(2) = 0.025
! a5
!    enz_pom_c(2) = 0.00013
! a6
!    enz_pom_c(2) = 0.00083
! a7
!    enz_pom_c(2) = 0.0019
! a8
!    enz_pom_c(2) = 0.000071
! a9
    enz_pom_c(2) = 0.00044
! a10
!    enz_pom_c(2) = 0.0015

    enz_pom_n(2) = enz_pom_c(2) / consts%enz_pom_c2n(2)
    enz_pom_p(2) = enz_pom_c(2) / consts%enz_pom_c2p(2)

    mom_c = 21./31. * cpct
    mom_n = mom_c / c2n
    mom_p = mom_c / c2p

    qom_c = 0.8 / 31. * cpct
    qom_n = qom_c / c2n
    qom_p = qom_c / c2p

! a1
!    enz_mom_c = 0.001 / 31. * cpct
! a2-4
!    enz_mom_c = 0.025
! a5
!    enz_mom_c = 0.0037
! a6
!    enz_mom_c = 0.020
! a7
!    enz_mom_c = 0.086
! a8
!    enz_mom_c = 0.0023
! a9
    enz_mom_c = 0.013
! a10
!    enz_mom_c = 0.058

    enz_mom_n = enz_mom_c / consts%enz_mom_c2n
    enz_mom_p = enz_mom_c / consts%enz_mom_c2p

! a1
!    amb_c = 0.25 / 31. * cpct
! a2-4
!    amb_c = 1.
! a5
!    amb_c = 0.0021
! a6
!    amb_c = 0.043
! a7
!    amb_c = 0.31
! a8
!    amb_c = 0.0014
! a9
    amb_c = 0.023
! a10
!    amb_c = 0.21

! a1-4
!    amb_n = amb_c / (0.5 * (consts%mb_cn_min + consts%mb_cn_max))
! a5
!    amb_n = 0.00027
! a6
!    amb_n = 0.0085
! a7
!    amb_n = 0.059
! a8
!    amb_n = 0.00017
! a9
    amb_n = 0.0046
! a10
!    amb_n = 0.038

! a1-4
!    amb_p = amb_c / (0.5 * (consts%mb_cp_min + consts%mb_cp_max))
! a5
!    amb_p = 3.5e-5
! a6
!    amb_p = 0.00098
! a7
!    amb_p = 0.0068
! a8
!    amb_p = 2.2e-5
! a9
    amb_p = 0.00055
! a10
!    amb_p = 0.0044

! a1
!    dmb_c = 0.25 / 31. * cpct
! a2-4
!    dmb_c = 0.001
! a5
!    dmb_c = 0.36
! a6
!    dmb_c = 0.34
! a7
!    dmb_c = 0.13
! a8
!    dmb_c = 0.30
! a9
    dmb_c = 0.27
! a10
!    dmb_c = 0.12

! a1-4
!    dmb_n = dmb_c / (0.5 * (consts%mb_cn_min + consts%mb_cn_max))
! a5
!    dmb_n = 0.072
! a6
!    dmb_n = 0.070
! a7
!    dmb_n = 0.024
! a8
!    dmb_n = 0.052
! a9
    dmb_n = 0.046
! a10
!    dmb_n = 0.021

! a1-4
!    dmb_p = dmb_c / (0.5 * (consts%mb_cp_min + consts%mb_cp_max))
! a5
!    dmb_p = 0.0083
! a6
!    dmb_p = 0.0078
! a7
!    dmb_p = 0.0028
! a8
!    dmb_p = 0.0061
! a9
    dmb_p = 0.0055
! a10
!    dmb_p = 0.0024

! a1
!    enz_ptase_c = 0.0001 / 31. * cpct
! a2-4
!    enz_ptase_c = 0.025
! a5
!    enz_ptase_c = 0.0021
! a6
!    enz_ptase_c = 0.012
! a7
!    enz_ptase_c = 0.050
! a8
!    enz_ptase_c = 0.0013
! a9
    enz_ptase_c = 0.0075
! a10
!    enz_ptase_c = 0.034

    enz_ptase_n = enz_ptase_c / consts%enz_ptase_c2n
    enz_ptase_p = enz_ptase_c / consts%enz_ptase_c2p

    co2_lost = 0.
    ngas_lost = 0.
    nh4 = 1.0e-5
    no3 = 1.0e-5

    ! See Yang et al. 2013, Biogeosciences, 10, 2525-2537
    !  the total is from measurements, but I am very rougly assuming
    !  a 3:1 ratio of occluded:secondary P, based on Yang et al. (2013)
    !  results for inceptisols.
    pocc = (totp - extrp) * 0.001 * 0.75
    psec = (totp - extrp) * 0.001 * 0.25
    ppar = 100. * 2. / bulk_den ! gP/m2 * 1000mg/g / 0.5 m / bulk_den / 1000kg/g
!    ppar = 0.

    ! SROAK:84 SRTDF:358 PV:527  mgP/kg
    my_total = extrp * 0.001
    my_b = consts%kmm_plang / bulk_den + consts%plab_max / bulk_den - my_total
    my_c = - consts%kmm_plang / bulk_den * my_total
    psol = (-my_b + sqrt(my_b**2-4.*my_c))*0.5
    plab = my_total - psol
!print*,plab,psol,consts%plab_max/bulk_den
!    plab = extrp * 0.001 ! mgP/gSOIL
    if(plab > consts%plab_max / bulk_den)then
       print*,'initialization problem.'
       print*,'initial plab > plab_max'
       stop
    endif
!    psol = plab * consts%kmm_plang / bulk_den / (consts%plab_max / bulk_den - plab) 
!print*,plab,psol,consts%plab_max/bulk_den;stop

    nh4_plant = 0.
    nh4_bnf = 0.
    no3_plant = 0.
    nmin = 0.
    nitr = 0.
    c_leach = 0.
    n_leach = 0.
    p_leach = 0.
    p_plant = 0.

    nh4_dep = 0.
    no3_dep = 0.
    ppar_dep = 0.

    enz_plant_n = 0.
    enz_plant_p = 0.

    vnh4up_plant = 0.
    vno3up_plant = 0.
    vpup_plant = 0.

    plant_input_C_pom = 0.
    plant_input_C_dom = 0.
    plant_input_N_pom = 0.
    plant_input_N_dom = 0.
    plant_input_P_pom = 0.
    plant_input_P_dom = 0.

    return
  end subroutine som_init

  subroutine som_extern_forcing(ndep_rate, consts, slden, input_nh4, &
       input_no3, pdep_rate, input_ppar, year, ndep_appl, pdep_appl)

    use mend_consts_coms, only: decomp_consts

    implicit none

    type(decomp_consts) :: consts
    real, intent(in) :: slden
    real, intent(in) :: ndep_rate
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(in) :: pdep_rate
    real, intent(out) :: input_ppar
    real, intent(in) :: ndep_appl
    real, intent(in) :: pdep_appl
    integer, intent(in) :: year

    input_nh4 = ndep_rate / (consts%eff_soil_depth * slden * 1000.)

    input_no3 = 0.

    input_ppar = pdep_rate / (consts%eff_soil_depth * slden * 1000.)

    if(year >= 2015)then
       input_nh4 = input_nh4 + ndep_appl / (consts%eff_soil_depth * slden * 1000.)
       input_ppar = input_ppar + pdep_appl / (consts%eff_soil_depth * slden * 1000.)
    endif

    return
  end subroutine som_extern_forcing

end Module mend_som
