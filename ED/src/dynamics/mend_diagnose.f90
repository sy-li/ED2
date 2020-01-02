Module mend_diagnose
  implicit none

Contains

  real function fTArh0(T, Tref, Ea)
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Ea [KJ/mol]: Arrhenius activation energy
    !temp (C): temperature
    !!ARGUMENTS:
    use consts_coms, only: rmol
    implicit none
    real, intent(in) :: T
    real, intent(in) :: Tref
    real, intent(in) :: Ea
    real :: TKref
    real :: TK

    TKref = Tref + 273.15
    TK = T + 273.15

    fTArh0 = exp(Ea*1.e3/rmol * (1.e0/TKref - 1.e0/TK))

    return
  END function fTArh0 !!function fTArh0


  real function fTArh(sCase, T, Tref)
    implicit none

    character(len=*), intent(in) :: sCase
    real :: Ea

    real, intent(in) :: T
    real, intent(in) :: Tref

    SELECT CASE (trim(sCase))
    CASE ("BG")  !! Beta-glucosidase
       Ea = 42.2
    CASE ("CBH") !! Cellobiohydrolase
       Ea = 32.2
    CASE ("EG")  !! Endo-glucanase
       Ea = 34.4
    CASE ("PER") !! Peroxidase
       Ea = 52.9
    CASE ("POX") !! Phenol oxidase
       Ea = 53.1
    CASE ("LIG") !! Ligninases
       Ea = 53.0
    CASE ("CEL") !! Cellulases
       Ea = 36.3
    CASE ("NAG") !! N-acetylglutamate synthase
       Ea = 52.9
    CASE ("PAC") !! Acid phosphatases
       Ea = 36.3
    CASE ("PAL") !! Alkaline phosphatases
       Ea = 23.6
    CASE ("PHO") !! PHOSPHATASES
       Ea = 34.4
    CASE ("Km")  !! half-saturation constant
       Ea = 30.
    CASE ("MR")  !! Microbial Maintenance
       Ea = 20.
    CASE ("Kads")!! adsorption
       Ea = 5.
    CASE ("Kdes")!! desorption
       Ea = 20.
    CASE ("DOM") !! DOM uptake
       Ea = 47.
    CASE DEFAULT
       Ea = 47.
    END SELECT

    fTArh = fTArh0(T, Tref, Ea)

    return
  end function fTArh

  real function fSWP_OPT(SWP)
    implicit none
    !SWP Scalar for SOM (lignin) decomposition
    !Hansen et al (1990), DAISY Model, page 105, Eq (6-16)
    
    REAL, PARAMETER :: SWPmin = -exp(2.5*log(10.))   !![MPa]
    REAL, PARAMETER :: SWPlow = -exp(-1.5*log(10.))
    REAL, PARAMETER :: SWPhigh= -exp(-2.5*log(10.))
    REAL, PARAMETER :: SWPmax = -exp(-4.0*log(10.))
    !!ARGUMENTS:
    REAL, intent(in) :: SWP ![MPa]
    
    if (SWP.le.SWPmin) then
       fSWP_OPT = 0.0  
    else if (SWP.le.SWPlow) then
       fSWP_OPT = 0.625-0.25*log10(abs(SWP))
    else if (SWP.le.SWPhigh) then
       fSWP_OPT = 1.0
    else if (SWP.le.SWPmax) then
       fSWP_OPT = (2.5+0.4*log10(abs(SWP)))/1.5
    else
       fSWP_OPT = 0.6
    end if
    return
  end function fSWP_OPT

  real function fpH(sCase,pH)
    implicit none
    character(len=*), intent(in) :: sCase
    real, intent(in) :: pH !pH value

    !!LOCAL VARIABLES
    real :: pHopt !optimum pH
    real :: pHsen !pH sensitivity
    
    SELECT CASE (trim(sCase))
    CASE ("BG")  !! Beta-glucosidase
       pHopt = 5.6
       pHsen = 1.7
    CASE ("CBH") !! Cellobiohydrolase
       pHopt = 5.1
       pHsen = 2.1
    CASE ("EG")  !! Endo-glucanase
       pHopt = 5.1
       pHsen = 1.6
    CASE ("PER") !! Peroxidase
       pHopt = 4.5
       pHsen = 1.5
    CASE ("POX") !! Phenol oxidase
       pHopt = 4.1
       pHsen = 1.4
    CASE ("LIG") !! Ligninases
       pHopt = 4.2
       pHsen = 1.4
    CASE ("CEL") !! Cellulases
       pHopt = 5.3
       pHsen = 1.7
    CASE ("PAC") !! Acid phosphatases
       pHopt = 5.2
       pHsen = 1.8
    CASE ("PAL") !! Alkaline phosphatases
       pHopt = 9.5
       pHsen = 2.6
    CASE ("PHO") !! PHOSPHATASES
       pHopt = 6.0
       pHsen = 2.0
    CASE ("ENZ") !! ENZ for MOM (Mineral-Associated Organic Matter)
       pHopt = 4.8
       pHsen = 1.6
    CASE DEFAULT
       pHopt = 6.0 !!mean pH of 763 soil samples
       pHsen = 2.0
    END SELECT
    
    fpH = fpH0(pH,pHopt,pHsen)
    return
  end function fpH
!------------------------------------------------------------------
  real function fpH0(pH, pHopt, pHsen)
    implicit none
    real, intent(in) :: pH !pH value
    real, intent(in) :: pHopt !optimum pH
    real, intent(in) :: pHsen !pH sensitivity

!    fpH0 = exp(-((pH - pHopt)/pHsen)**2)
    fpH0 = 1.

    return
  end function fpH0

  real function fSWP(SWP,BIOME,SOM)
    !SWP Scalar for SOM (Cellulose) decomposition
    !Manzoni et al (2012), Ecology, 93: 930-938
    !!fSWP increases with increasing SWP (wetter condition)
    implicit none
    !REAL(8), PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 

    !!ARGUMENTS:
    real, intent(in) :: SWP ![MPa]
    CHARACTER(len=3), intent(in) :: BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3), intent(in) :: SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
    !!LOCAL VARIABLES:
    real :: SWPmin, w ![MPa]
    
    if (trim(SOM).eq."SOD") then
       SWPmin = -1.71e0
       w      = 1.43e0
    else if (trim(SOM).eq."SOL") then
       SWPmin = -13.86e0
       w      = 1.20e0
    else if (trim(SOM).eq."LIT") then
       SWPmin = -36.49e0
       w      = 1.04e0
    end if
    
    SELECT CASE (trim(BIOME)) 
    CASE ("ASM")  !!Arid, Semiarid, & Mediterranean
       if (trim(SOM).eq."SOL") then
          SWPmin = -10.95e0
          w      = 1.26e0
       end if
    CASE ("MGC")  !!Mesic Grassland & Cropland
       if (trim(SOM).eq."SOD") then
          SWPmin = -1.71e0
          w      = 1.43e0
       else if (trim(SOM).eq."SOL") then
          SWPmin = -22.61e0
          w      = 1.11e0
       else if (trim(SOM).eq."LIT") then
          SWPmin = -39.73e0
          w      = 0.89e0
       end if
    CASE ("MDF")    !!Mesic Deciduous Forest
       if (trim(SOM).eq."SOL") then
          SWPmin = -4.97e0
          w      = 1.07e0
       else if (trim(SOM).eq."LIT") then
          SWPmin = -29.00e0
          w      = 1.27e0
       end if
    CASE ("MCF")    !!MEsic Conifer Forest
       if (trim(SOM).eq."SOL") then
          SWPmin = -8.24e0
          w      = 1.40e0
       else if (trim(SOM).eq."LIT") then
          SWPmin = -39.85e0
          w      = 1.06e0
       end if
    CASE DEFAULT    !!All Biome Average
       if (trim(SOM).eq."SOD") then
          SWPmin = -1.71e0
          w      = 1.43e0
       else if (trim(SOM).eq."SOL") then
          SWPmin = -13.86e0
          w      = 1.20e0
       else if (trim(SOM).eq."LIT") then
          SWPmin = -36.49e0
          w      = 1.04e0
       end if
    END SELECT
    fSWP = fSWP0(SWP,SWPmin,w)
    return
  end function fSWP

  real function fSWP0(SWP,SWPmin,w)
    implicit none
    !Rate Scalar for Soil Water Potential (SWP)
    !Manzoni et al (2012), Ecology, 93: 930-938
!    USE MOD_MEND
!    IMPLICIT NONE
    
    REAL, PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 
    !!ARGUMENTS:
    REAL :: SWP, SWPmin ![MPa]
    REAL :: w
    
    if (SWP.lt.SWPmin) then
       fSWP0 = 0.0e0  
    else if (SWP.lt.SWP_FC) then
       fSWP0 = 1-(log(SWP/SWP_FC)/log(SWPmin/SWP_FC))**w
    else 
       fSWP0 = 1.0
    end if
    return
  end function fSWP0

  REAL function fT_Linear(T, Tref, slope, intercept)
    implicit none
    !fKmT: [mg/m3], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [mg/m3/C]
    !intercept: [mg/m3]
    !!ARGUMENTS:
    REAL, intent(in) :: T, Tref, slope, intercept
    
    fT_Linear = slope * (T - Tref) + intercept
    return
  END function fT_Linear
  
  real function fT_CUE(T,Tref,slope,CUEref)
    !fT_CUE: [-], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [1/degree C]
    !intercept: [-]
    !! parameter values: Wang et al. (2015), ISMEJ 
    !!ARGUMENTS:
    implicit none
    REAL, intent(in):: T, Tref,slope,CUEref
!    REAL(8), PARAMETER :: Tref      = 0    !![degree C]
!    REAL(8), PARAMETER :: slope     = -0.01
!    REAL(8), PARAMETER :: intercept = 0.56  !! []
    
    !!LOCAL VARIABLES:
    REAL, PARAMETER:: CUEmax = 0.9D0
    REAL, PARAMETER:: CUEmin = 0.1D-1
    
    fT_CUE = fT_Linear(T, Tref, slope, CUEref)
    if(fT_CUE.gt.CUEmax) then
       fT_CUE = CUEmax
    elseif(fT_CUE.lt.CUEmin) then
       fT_CUE = CUEmin
    end if
    
    return
  END function fT_CUE

  real function fSWP_A2D(SWP, SWP_A2D, w)
    implicit none
    !Soil Water Scalar for Microbial Dormancy
    !! Manzoni (2014) SBB, 73: 69-83
    
    !!ARGUMENTS:
    REAL, intent(in):: SWP_A2D !! = -0.4 [MPa], 
    REAL, intent(in):: SWP
    REAL, intent(in):: w       !! = 4 
    
    fSWP_A2D = abs(SWP)**w/(abs(SWP)**w + abs(SWP_A2D)**w)
!    if(ISNAN(fSWP_A2D)) then
!        print*,"wp_scalar_low=",fSWP_A2D
!    end if
    return
  END function fSWP_A2D !!function fSWP_A2D

  real function fSWP_D2A(SWP, SWP_D2A, w)
    implicit none
    !!Soil Water Scalar for Microbial reactivation
    !!ARGUMENTS:
    REAL, intent(in):: SWP_D2A !! = 1/4*SWP_A2D [MPa]
    REAL, intent(in):: SWP
    REAL, intent(in):: w       !! = 4 
    
    fSWP_D2A = abs(SWP_D2A)**w/(abs(SWP)**w + abs(SWP_D2A)**w)
    return
  END function fSWP_D2A !!function fSWP_D2A
  
  REAL function fWFP_PieceWise0(WFP,WFPcr,Slope,Intercept)
    implicit none
    !SWP Scalar for Nitrification/Denitrification
    !Muller (1999)
    
    !!ARGUMENTS:
    REAL, intent(in) :: WFP ![-], water-filled porosity
    REAL, intent(in), dimension(4) :: WFPcr
    REAL, intent(in), dimension(2) :: Slope
    REAL, intent(in), dimension(2) :: Intercept
    
    if (WFP.le.WFPcr(1)) then
       fWFP_PieceWise0 = 0.0d0  
    else if (WFP.le.WFPcr(2)) then
       fWFP_PieceWise0 = Intercept(1) + Slope(1)*WFP
    else if (WFP.le.WFPcr(3)) then
       fWFP_PieceWise0 = 1.0d0
    else if (WFP.le.WFPcr(4)) then
       fWFP_PieceWise0 = Intercept(2) + Slope(2)*WFP
    else
       fWFP_PieceWise0 = 0.0d0  
    end if
    fWFP_PieceWise0 = max(0.d0, min(1.d0, fWFP_PieceWise0))
    return
  END function fWFP_PieceWise0 !!fWFP_PieceWise0
!------------------------------------------------------------------
  REAL function fWFP_PieceWise(sCase,WFP)
    implicit none
    !SWP Scalar for Nitrification/Denitrification
    !Muller (1999)
    
    !!ARGUMENTS:
    CHARACTER(len=*), intent(in) :: sCase
    REAL, intent(in) :: WFP ![-], water-filled porosity
    
    !!LOCAL VARIABLES:
    REAL :: WFPcr(4), Slope(2), Intercept(2)
    
    SELECT CASE (trim(sCase))
    CASE ("NITRIF")  !!
       WFPcr = (/0.09d0, 0.54d0,0.69d0,1.00d0/)
       Slope       = (/ 2.20d0, -3.23d0/)
       Intercept   = (/-0.19d0,  3.23d0 /)
    CASE ("DENITRIF_NO3")
       WFPcr       = (/0.36d0, 1.d0, 1.d0, 1.01d0 /)
       Slope       = (/ 1.56d0, 0.d0 /)
       Intercept   = (/-0.56d0, 1.d0/)
    CASE ("DENITRIF_NO2")
       WFPcr       = (/0.4d0, 0.6d0, 0.66d0, 0.7d0/)
       Slope       = (/5.d0, -20d0/)
       Intercept   = (/-2.d0, 14d0/)    
    CASE ("DENITRIF_NO")
       WFPcr       = (/ 0.1d0, 0.8d0, 0.9d0, 1.d0/)
       Slope       = (/1.43d0, -10d0/)
       Intercept   = (/-0.14d0, 10.d0/)
    CASE ("DENITRIF_N2O")
       WFPcr       = (/0.40d0, 0.85d0, 1.d0, 1.01d0/)
       Slope       = (/2.22d0, 0.d0/)
       Intercept   = (/-0.89d0,1.d0/)
    CASE DEFAULT
       
    END SELECT
    
    fWFP_PieceWise = fWFP_PieceWise0(WFP,WFPcr,Slope,Intercept)

    return
  END function fWFP_PieceWise !!fWFP_PieceWise

  subroutine mend_update_parameters(tp, wp, pH, wfp, consts, consts_base, &
       slmsts, sfldcap)
    use mend_consts_coms, only: decomp_consts
    implicit none
    type(decomp_consts) :: consts
    type(decomp_consts) :: consts_base
    real, parameter :: const_Tref = 20.
    real, intent(in) :: tp ! temperature in C
    real, intent(in) :: wp ! soil water potential, MPa
    real, intent(in) :: pH
    character(len=16) :: sCase
    real :: tp_scalar
    real :: wp_scalar_opt
    real :: pH_scalar
    character(len=3), parameter :: biome='tdf'
    character(len=3), parameter :: som='SOL'
    real :: wp_scalar
    real :: cue_slope
    real :: wp_scalar_low
    real, intent(in) :: wfp
    real, intent(in) :: slmsts
    real, intent(in) :: sfldcap
    real :: Potter_p
    real :: Potter_y

    !![1] DECOMPOSITION of POM1 & POM2
    sCase = "LIG"
    tp_scalar       = fTArh(sCase, tp, const_Tref)
    wp_scalar_opt   = fSWP_OPT(wp)
    pH_scalar       = fpH(sCase,pH)
    consts%vmax_decomp_pom(1) = consts_base%vmax_decomp_pom(1) * &
         tp_scalar * wp_scalar_opt * pH_scalar
    
    sCase = "CEL"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, biome, som)
    pH_scalar = fpH(sCase,pH)
    consts%vmax_decomp_pom(2) = consts_base%vmax_decomp_pom(2) *  &
         tp_scalar * wp_scalar * pH_scalar
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%kmm_decomp_pom(1) = consts_base%kmm_decomp_pom(1) * tp_scalar ![mg POM/cm3], half-saturation constant for conversion of POM by ENZP (/xx(6), xx(7)/)
    consts%kmm_decomp_pom(2) = consts_base%kmm_decomp_pom(2) * tp_scalar ![mg POM/cm3], half-saturation constant for conversion of POM by ENZP (/xx(6), xx(7)/)

    !![2] DECOMPOSITION of MOM
    sCase = "ENZ"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
    consts%vmax_decomp_mom = consts_base%vmax_decomp_mom * tp_scalar * wp_scalar * pH_scalar ![mg MOM/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC xx(5)
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%kmm_decomp_mom = consts_base%kmm_decomp_mom * tp_scalar ![mg MOM/cm3], half-saturation constant for conversion of MAOC by ENZMAOC xx(8)

    !![3] ADSORPTION-DESORPTION
    sCase = "Kdes"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%desorp_rate = consts_base%desorp_rate * tp_scalar ![mg DOM/h],desorption rate constant  xx(11)
    
    sCase = "Kads"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%sorp_rate = consts_base%sorp_rate * tp_scalar ![mg DOM/mg DOM/h], adsorption rate constant = Kba*Kdes
    
    !![6] MICROBIAL UPTAKE of DOM
    sCase = "DOM"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)  !!same as Cellulose-decomposition
    consts%vmax_amb = consts_base%vmax_amb * tp_scalar !!*wp_scalar !![mg DOM/mg MB/h], maximum uptake rate of DOM by MB  xx(17)
    
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%spec_maint_rate = consts_base%spec_maint_rate * tp_scalar !!*wp_scalar ![1/h], specific microbial maintenance rate
     
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)  
    consts%kmm_amb = consts_base%kmm_amb * tp_scalar ![mg DOM/cm3], half-saturation constant for uptake of DOM by MB  xx(19)
    
    sCase = "CUE"  !!label, no-use
    cue_slope = -0.012
    consts%growth_yield = fT_CUE(tp, const_Tref, CUE_slope, consts_base%growth_yield) ![-], carbon use efficiency in uptake of DOM by MB 
    
    !![7] MICROBIAL MORTALITY
    wp_scalar_low = fSWP_A2D(wp, consts%SWP_A2D, consts%wdorm)  !!negative effect  !!wdie
    consts%amb_turnover_rate = consts_base%amb_turnover_rate * wp_scalar_low
       
    !![8] MICROBIAL DORMANCY & RESUSCITATION

    
    wp_scalar     = fSWP_D2A(wp, consts%SWP_D2A, consts%wdorm)  !!positive effect, increase with increasing SWC or SWP
    wp_scalar_low = fSWP_A2D(wp, consts%SWP_A2D, consts%wdorm)  !!negative effect
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%vmax_a2d = consts_base%vmax_a2d * tp_scalar * wp_scalar_low 
    consts%vmax_d2a = consts_base%vmax_d2a * tp_scalar * wp_scalar
    
    !![9] MICROBIAL UPTAKE of NH4 & NO3
    sCase = "NH4"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%vnup_amb_base = consts_base%vnup_amb_base * tp_scalar !![mg N/mg MBN/h], maximum uptake rate of NH4 or NO3 by MB  xx(28)
    consts%vpup_amb_base = consts_base%vpup_amb_base * tp_scalar !![mg P/mg MBP/h], maximum uptake rate of P by MB
    consts%vnh4up_plant_base = consts_base%vnh4up_plant_base * tp_scalar !![1 / s], maximum uptake rate of NH4 by plants
    consts%vno3up_plant_base = consts_base%vno3up_plant_base * tp_scalar !![1/s], maximum uptake rate of NO3 by plants
    consts%vpup_plant_base = consts_base%vpup_plant_base * tp_scalar !![1/s], maximum uptake rate of P by plants
   
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    consts%ksnh4_amb = consts_base%ksnh4_amb * tp_scalar ![mg N/cm3], half-saturation constant xx(29)
    consts%ksno3_amb = consts_base%ksno3_amb * tp_scalar ![mg N/cm3], half-saturation constant xx(30)
    consts%kspsol_amb = consts_base%kspsol_amb * tp_scalar ![mg N/cm3], half-saturation constant
    consts%ksnh4_plant = consts_base%ksnh4_plant * tp_scalar ![mg N/cm3], half-saturation constant xx(29)
    consts%ksno3_plant = consts_base%ksno3_plant * tp_scalar ![mg N/cm3], half-saturation constant xx(30)
    consts%kspsol_plant = consts_base%kspsol_plant * tp_scalar ![mg N/cm3], half-saturation constant
    
    !![10] NITRIFICATION & DENITRIFICATION
    sCase = "NITRIF"
    tp_scalar       = fTArh(sCase, tp, const_Tref)
    wp_scalar_opt   = fWFP_PieceWise(sCase,wfp)
    consts%v_nitr = consts_base%v_nitr * tp_scalar * wp_scalar_opt !!xx(31)
    
    sCase = "DENITRIF_NO3"
    tp_scalar       = fTArh(sCase, tp, const_Tref)
    wp_scalar_opt   = fWFP_PieceWise(sCase,wfp)
    consts%v_denitr = consts_base%v_denitr * tp_scalar*  wp_scalar_opt  !!xx(32)

!   Potter et al. 1996
    Potter_P = slmsts - sfldcap
    Potter_y = 0.477 * Potter_P**3 - 0.596 * Potter_P**2 +   &
         0.437 * Potter_P + 0.564
    consts%v_denitr = Potter_P**Potter_y

    return
  end subroutine mend_update_parameters

  subroutine mend_update_diag_layer(mvars, consts, ipa, bulk_den)
    use mend_state_vars, only: mend_vars, npom
    use mend_consts_coms, only: decomp_consts
    implicit none
    real, intent(in) :: bulk_den
    integer, intent(in) :: ipa
    type(mend_vars) :: mvars
    type(decomp_consts) :: consts
    real :: my_total, my_k, my_v
    real :: my_b
    real :: my_c
    integer :: ipom
    real :: tot_enz_plant_p

    my_total = mvars%invars%plab(ipa) + mvars%invars%psol(ipa)
    tot_enz_plant_p = sum(mvars%plvars%enz_plant_p(:,ipa)/consts%kspsol_plant(:)) * &
         bulk_den

    my_k = consts%kspsol_lab / bulk_den * (1.0 +  &
         tot_enz_plant_p  +   &
         mvars%pvars%amb(ipa) / (consts%kspsol_amb / bulk_den) + &
         (consts%vmax_surf/bulk_den) / (consts%kspsol_lab / bulk_den))

    my_v = consts%vmax_surf / bulk_den

    mvars%invars%plab(ipa) = my_v * my_total / (my_k + my_total)
    mvars%invars%psol(ipa) = my_total - mvars%invars%plab(ipa)

    do ipom = 1, npom
       mvars%nvars%enz_pom(ipom,ipa) = mvars%cvars%enz_pom(ipom,ipa) / consts%enz_pom_c2n(ipom)
       mvars%pvars%enz_pom(ipom,ipa) = mvars%cvars%enz_pom(ipom,ipa) / consts%enz_pom_c2p(ipom)
    enddo
    mvars%nvars%enz_mom(ipa) = mvars%cvars%enz_mom(ipa) / consts%enz_mom_c2n
    mvars%pvars%enz_mom(ipa) = mvars%cvars%enz_mom(ipa) / consts%enz_mom_c2p
    mvars%nvars%enz_ptase(ipa) = mvars%cvars%enz_ptase(ipa) / consts%enz_ptase_c2n
    mvars%pvars%enz_ptase(ipa) = mvars%cvars%enz_ptase(ipa) / consts%enz_ptase_c2p

    return
  end subroutine mend_update_diag_layer

end Module mend_diagnose
