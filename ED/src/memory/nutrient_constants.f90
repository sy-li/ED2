! Defines constants used in nitrogen cycle 

module nutrient_constants

  !! Constants for Disturbance settings
  integer, parameter :: dist_start = 50 ! year disturbance starts
  integer, parameter :: dist_end = 100  ! year disturbance reaches max_treefall_distrubance_rate
  real, parameter :: max_treefall_disturbance_rate = 0.008 !  1/years

  !!Constants for N inputs and losses
  real, parameter :: N_deposition = 0.000383/365! 0.0009/365          !(kg N/m2/day) (Hedin et al.2009)
  ! N deposition for boreal forests take from Piiraine et al 1998 Biogeochemical Investigations at Watershed, Landscape, and Regional Scales 
  real, parameter :: leaching_efficiency_factor = 0.0 !0.01  !fraction of N that can be lost by leaching from mineralized_soil_N 
  real, parameter :: soil_depth = 1.5                   !used to calcualte leaching (meters) 1.5 from (Jackson 1996)
  real, parameter :: bstorage_max_factor = 2.0          !maximum storage capacity of N = amount in leaves * nstorge_max_factor
  real, parameter :: nstorage_max_factor = 2.0          !maximum storage capacity of N = amount in leaves * nstorge_max_factor
  real, parameter :: pstorage_max_factor = 2.0          !maximum storage capacity of N = amount in leaves * nstorge_max_factor
  real, parameter :: gas_loss_rate = 0.2                !gas loss is approximately 20% of total N loss (Houlton et al 2006 PNAS) if the site has                                                                         !2700 mm rain/yr or less but is about 50% gas loss if > 2700 mm 
  !NOTE: 1.784159e-05 is the mean observed rate at Sarah Batterman's field site (Agua Salud)
  real, parameter :: Cost_BNF = 9.12                    !( 9.12 g C g N -1), (gutschick, 1981) metabolic and nodule maintence cost
  real, parameter :: Cost_fixer = 0                     !( g C gN -1), additional cost for being a fixer
  !!Constants for N inputs and losses
!  real, parameter :: BNF_rate= 0.0003941912             !kg N fixed/total tree biomass(above and belowground)/ day
  real, parameter :: BNF_rate= 1. * 0.0003941912             !kg N fixed/total tree biomass(above and belowground)/ day

  !! Constants for forest organic layer
  real, parameter :: organic_matter_density = 200.      !(kgC/m^3). Used to calculate litter layer depth from fast_soil_C + slow_soil_C
  real, parameter :: crit_olayer_min =  0.02            !(cm) minimum organic layer depth where seedling mortality starts increasing
  real, parameter :: crit_olayer_max =  0.04            !(cm) organic layer depth where seedling mortality stops increasing with increased depth
  real, parameter :: blackspruce_seedling_mort = 0.95   !factor that determins the seedling mortality at organic_layers_depth > crit_olayer_max
  real, parameter :: aspen_seedling_mort =  1.0         !factor that determins the seedling mortality at organic_layers_depth > crit_olayer_max
  integer, parameter :: organic_soil_texture =  12      !soil type in ed_params.f90 that corresponds to the organic layer texture

  real :: Ndep_rate, Pdep_rate, Ccwd_init, cwd_resp_frac
  real, dimension(4):: Csom_init !AMT 1/20/17
  real, dimension(3):: Clitt_init !AMT 1/20/17
  real :: Cmyco_init, Nmin_init, Psol_init, Ndep_appl,Pdep_appl

  real :: soil_bulk_den, soil_ph, soil_cpct, soil_som_c2n, soil_totp, soil_extrp
  integer :: nlsl

! when creating a new file, add to bin/rules.mk and bin/objects.mk

end module nutrient_constants
