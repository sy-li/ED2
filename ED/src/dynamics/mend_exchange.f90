Module mend_exchange
  implicit none

  real, parameter :: pom_burial_rate = 1./(365.*2.*24)  ! 1/hr
  real, parameter :: dom_burial_rate = 1./240.  ! 1/hr
  real, parameter :: som_depth = 0.2 ! m

Contains

  subroutine plant2litt_exchange(npom, input_pom_c, input_pom_n, input_pom_p, &
       input_dom_c, input_dom_n, input_dom_p, &
       input_nh4, input_no3, input_psol)
    implicit none
    
    integer, intent(in) :: npom
    real, dimension(npom), intent(out) :: input_pom_c
    real, intent(out) :: input_dom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_psol

    input_pom_c(1) = 0.04  ! gC/m2/hr
    input_pom_n(1) = input_pom_c(1) / 100.
    input_pom_p(1) = input_pom_c(1) / 100.

    input_pom_c(2) = 0.04  ! gC/m2/hr
    input_pom_n(2) = input_pom_c(2) / 25.
    input_pom_p(2) = input_pom_c(2) / 100.

    input_dom_c = 0.04  ! gC/m2/hr
    input_dom_n = input_dom_c / 15.
    input_dom_p = input_dom_c / 100.

    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.

    return
  end subroutine plant2litt_exchange

  subroutine plant2som_exchange(npom, input_pom_c, input_pom_n, input_pom_p, &
       input_dom_c, input_dom_n, input_dom_p, &
       input_nh4, input_no3, input_psol, plant_input_C, plant_input_N, &
       plant_input_P)
    implicit none
    
    integer, intent(in) :: npom
    real, dimension(npom), intent(out) :: input_pom_c
    real, intent(out) :: input_dom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_psol
    real, dimension(4), intent(in) :: plant_input_C
    real, dimension(4), intent(in) :: plant_input_N
    real, dimension(4), intent(in) :: plant_input_P

!    input_pom_c(1) = (plant_input_C(2) + plant_input_C(3) +   &
!         plant_input_C(4))*1000./86400. ! gC/m2/s
!    input_pom_n(1) = (plant_input_N(2) + plant_input_N(3) +   &
!         plant_input_N(4))*1000./86400. ! gN/m2/s
!    input_pom_p(1) = (plant_input_P(2) + plant_input_P(3) +   &
!         plant_input_P(4))*1000./86400. ! gP/m2/s

!    input_pom_c(2) = (plant_input_C(1) +   &
!         0. * plant_input_C(4))*1000./86400. ! gC/m2/s
!    input_pom_n(2) = (plant_input_N(1) +   &
!         0. * plant_input_N(4))*1000./86400. ! gN/m2/s
!    input_pom_p(2) = (plant_input_P(1) +   &
!         0. * plant_input_P(4))*1000./86400. ! gP/m2/s

!    input_dom_c = 0. * plant_input_C(1) * 1000. / 86400. ! gC/m2/s
!    input_dom_n = 0. * plant_input_N(1) * 1000. / 86400. ! gC/m2/s
!    input_dom_p = 0. * plant_input_P(1) * 1000. / 86400. ! gC/m2/s

    input_pom_c(1) = (plant_input_C(3) +   &
         0.75 * plant_input_C(4))*1000./86400. ! gC/m2/s
    input_pom_n(1) = (plant_input_N(3) +   &
         0.75 * plant_input_N(4))*1000./86400. ! gN/m2/s
    input_pom_p(1) = (plant_input_P(3) +   &
         0.75 * plant_input_P(4))*1000./86400. ! gP/m2/s

    input_pom_c(2) = (plant_input_C(2) +   &
         0.25 * plant_input_C(4))*1000./86400. ! gC/m2/s
    input_pom_n(2) = (plant_input_N(2) +   &
         0.25 * plant_input_N(4))*1000./86400. ! gN/m2/s
    input_pom_p(2) = (plant_input_P(2) +   &
         0.25 * plant_input_P(4))*1000./86400. ! gP/m2/s

    input_dom_c = plant_input_C(1) * 1000. / 86400. ! gC/m2/s
    input_dom_n = plant_input_N(1) * 1000. / 86400. ! gC/m2/s
    input_dom_p = plant_input_P(1) * 1000. / 86400. ! gC/m2/s
!if(plant_input_C(4) > 0.)then
!print*
!print*,plant_input_C
!print*,plant_input_C/plant_input_N
!print*,plant_input_C/plant_input_P
!endif
    input_nh4 = 0. 
    input_no3 = 0.
    input_psol = 0.

    return
  end subroutine plant2som_exchange

  subroutine litt2som_exchange(npom,  &
       litt_pom_c, litt_pom_n, litt_pom_p, &
       litt_dom_c, litt_dom_n, litt_dom_p, &
       litt_nh4, litt_no3, litt_psol, &
       input_pom_c, input_pom_n, input_pom_p, &
       input_dom_c, input_dom_n, input_dom_p, &
       input_nh4, input_no3, input_psol)
    implicit none
    
    integer, intent(in) :: npom

    real, dimension(npom), intent(in) :: litt_pom_c
    real, intent(in) :: litt_dom_c
    real, dimension(npom), intent(in) :: litt_pom_n
    real, intent(in) :: litt_dom_n
    real, intent(in) :: litt_nh4
    real, intent(in) :: litt_no3
    real, dimension(npom), intent(in) :: litt_pom_p
    real, intent(in) :: litt_dom_p
    real, intent(in) :: litt_psol

    real, dimension(npom), intent(out) :: input_pom_c
    real, intent(out) :: input_dom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_psol

    input_pom_c(1) = 0. * pom_burial_rate * litt_pom_c(1)
    input_pom_n(1) = 0. * pom_burial_rate * litt_pom_n(1)
    input_pom_p(1) = 0. * pom_burial_rate * litt_pom_p(1)

    input_pom_c(2) = 0. * pom_burial_rate * litt_pom_c(2)
    input_pom_n(2) = 0. * pom_burial_rate * litt_pom_n(2)
    input_pom_p(2) = 0. * pom_burial_rate * litt_pom_p(2)

    input_dom_c = 0. * dom_burial_rate * litt_dom_c
    input_dom_n = 0. * dom_burial_rate * litt_dom_n
    input_dom_p = 0. * dom_burial_rate * litt_dom_p

    input_nh4 = 0. * dom_burial_rate * litt_nh4
    input_no3 = 0. * dom_burial_rate * litt_no3
    input_psol = 0. * dom_burial_rate * litt_psol

    return
  end subroutine litt2som_exchange

  subroutine plant2wood_exchange(iwood, npom, &
       input_pom_c, input_pom_n,   &
       input_pom_p, input_dom_c,   &
       input_dom_n, input_dom_p, &
       input_nh4, input_no3,   &
       input_psol)
    implicit none

    integer, intent(in) :: iwood
    integer, intent(in) :: npom
    real, intent(out), dimension(npom) :: input_pom_c
    real, intent(out), dimension(npom) :: input_pom_n
    real, intent(out), dimension(npom) :: input_pom_p
    real, intent(out) :: input_dom_c
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(out) :: input_psol

    ! g/m2/hr
    input_pom_c = 0.
    input_pom_n = 0.
    input_pom_p = 0.
    input_dom_c = 0.
    input_dom_n = 0.
    input_dom_p = 0.
    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.
    return

    if(iwood == 1)then
       input_pom_c(1) = 100. / 365.25 / 24.
       input_pom_n(1) = input_pom_c(1) / 150.
       input_pom_p(1) = input_pom_c(1) / 300.
    elseif(iwood == 2)then
       input_pom_c(1) = 100. / 365.25 / 24.
       input_pom_n(1) = input_pom_c(1) / 150.
       input_pom_p(1) = input_pom_c(1) / 300.
    elseif(iwood == 3)then
       input_pom_c(1) = 20. / 365.25 / 24.
       input_pom_n(1) = input_pom_c(1) / 150.
       input_pom_p(1) = input_pom_c(1) / 300.
    elseif(iwood == 4)then
       input_pom_c(1) = 500. / 365.25 / 24.
       input_pom_n(1) = input_pom_c(1) / 40.
       input_pom_p(1) = input_pom_c(1) / 100.
    else
       input_pom_c(1) = 0.
       input_pom_n(1) = 0.
       input_pom_p(1) = 0.
    endif

    if(iwood == 4)then
       input_pom_c(2) = 500. / 365.25 / 24.
       input_pom_n(2) = input_pom_c(1) / 40.
       input_pom_p(2) = input_pom_c(1) / 100.
    else
       input_pom_c(2) = 0.
       input_pom_n(2) = 0.
       input_pom_p(2) = 0.
    endif

    input_dom_c = 0.
    input_dom_n = 0.
    input_dom_p = 0.

    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.

    return
  end subroutine plant2wood_exchange

  subroutine wood2litt_exchange(iwood, npom, &
       pom_c, pom_n, pom_p, &
       dom_c, dom_n, dom_p, &
       nh4, no3, psol, &
       input_pom_c, input_pom_n,   &
       input_pom_p, input_dom_c,  &
       input_dom_n, input_dom_p,   &
       input_nh4, input_no3, input_psol, consts)
    
    use mend_consts_coms, only: decomp_consts

    implicit none

    integer, intent(in) :: iwood
    integer, intent(in) :: npom

    real, dimension(npom), intent(in) :: pom_c
    real, dimension(npom), intent(in) :: pom_n
    real, dimension(npom), intent(in) :: pom_p
    real, intent(in) :: dom_c
    real, intent(in) :: dom_n
    real, intent(in) :: dom_p
    real, intent(in) :: nh4
    real, intent(in) :: no3
    real, intent(in) :: psol

    real, dimension(npom), intent(out) :: input_pom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_c
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(out) :: input_psol

    integer :: ipom
    type(decomp_consts) :: consts

    input_pom_c = 0.
    input_pom_n = 0.
    input_pom_p = 0.
    input_dom_c = 0.
    input_dom_n = 0.
    input_dom_p = 0.
    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.
    return

    if(iwood /= 2)then
       input_pom_c(:) = 0.
       input_pom_n(:) = 0.
       input_pom_p(:) = 0.
       input_dom_c = 0.
       input_dom_n = 0.
       input_dom_p = 0.
       input_nh4 = 0.
       input_no3 = 0.
       input_psol = 0.
       return
    endif

    do ipom = 1, npom
       input_pom_c(ipom) = consts%frag_rate * pom_c(ipom)
       input_pom_n(ipom) = consts%frag_rate * pom_n(ipom)
       input_pom_p(ipom) = consts%frag_rate * pom_p(ipom)
    enddo

    input_dom_c = consts%frag_rate * dom_c
    input_dom_n = consts%frag_rate * dom_n
    input_dom_p = consts%frag_rate * dom_p

    input_nh4 = consts%frag_rate * nh4
    input_no3 = consts%frag_rate * no3
    input_psol = consts%frag_rate * psol

    return
  end subroutine wood2litt_exchange

  subroutine wood2som_exchange(iwood, npom, &
       pom_c, pom_n, pom_p, &
       dom_c, dom_n, dom_p, &
       nh4, no3, psol, &
       input_pom_c, input_pom_n,   &
       input_pom_p, input_dom_c,  &
       input_dom_n, input_dom_p,   &
       input_nh4, input_no3, input_psol, consts)
    use mend_consts_coms, only: decomp_consts
    implicit none

    integer, intent(in) :: iwood
    integer, intent(in) :: npom

    real, dimension(npom), intent(in) :: pom_c
    real, dimension(npom), intent(in) :: pom_n
    real, dimension(npom), intent(in) :: pom_p
    real, intent(in) :: dom_c
    real, intent(in) :: dom_n
    real, intent(in) :: dom_p
    real, intent(in) :: nh4
    real, intent(in) :: no3
    real, intent(in) :: psol

    real, dimension(npom), intent(out) :: input_pom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_c
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(out) :: input_psol
    integer :: ipom
    type(decomp_consts) :: consts

    input_pom_c = 0.
    input_pom_n = 0.
    input_pom_p = 0.
    input_dom_c = 0.
    input_dom_n = 0.
    input_dom_p = 0.
    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.

    return

    if(iwood /= 4)then
       input_pom_c(:) = 0.
       input_pom_n(:) = 0.
       input_pom_p(:) = 0.
       input_dom_c = 0.
       input_dom_n = 0.
       input_dom_p = 0.
       input_nh4 = 0.
       input_no3 = 0.
       input_psol = 0.
       return
    endif

    do ipom = 1, npom
       input_pom_c(ipom) = consts%frag_rate * pom_c(ipom)
       input_pom_n(ipom) = consts%frag_rate * pom_n(ipom)
       input_pom_p(ipom) = consts%frag_rate * pom_p(ipom)
    enddo

    input_dom_c = consts%frag_rate * dom_c
    input_dom_n = consts%frag_rate * dom_n
    input_dom_p = consts%frag_rate * dom_p

    input_nh4 = consts%frag_rate * nh4
    input_no3 = consts%frag_rate * no3
    input_psol = consts%frag_rate * psol

    return
  end subroutine wood2som_exchange

  subroutine wood2wood_exchange(iw1, iw2, npom, pom_c, pom_n, pom_p, &
       dom_c, dom_n, dom_p, nh4, no3, psol, input_pom_c, input_pom_n,  &
       input_pom_p, input_dom_c, input_dom_n, input_dom_p, input_nh4,   &
       input_no3, input_psol, consts)

    use mend_consts_coms, only: decomp_consts

    implicit none

    integer, intent(in) :: iw1
    integer, intent(in) :: iw2
    integer, intent(in) :: npom

    real, dimension(npom), intent(in) :: pom_c
    real, dimension(npom), intent(in) :: pom_n
    real, dimension(npom), intent(in) :: pom_p
    real, intent(in) :: dom_c
    real, intent(in) :: dom_n
    real, intent(in) :: dom_p
    real, intent(in) :: nh4
    real, intent(in) :: no3
    real, intent(in) :: psol

    real, dimension(npom), intent(out) :: input_pom_c
    real, dimension(npom), intent(out) :: input_pom_n
    real, dimension(npom), intent(out) :: input_pom_p
    real, intent(out) :: input_dom_c
    real, intent(out) :: input_dom_n
    real, intent(out) :: input_dom_p
    real, intent(out) :: input_nh4
    real, intent(out) :: input_no3
    real, intent(out) :: input_psol
    integer :: ipom
    type(decomp_consts) :: consts

    input_pom_c(:) = 0.
    input_pom_n(:) = 0.
    input_pom_p(:) = 0.
    input_dom_c = 0.
    input_dom_n = 0.
    input_dom_p = 0.
    input_nh4 = 0.
    input_no3 = 0.
    input_psol = 0.
    return

    if( (iw1 == 1 .and. iw2 == 2) .or. &
         (iw1 == 3 .and. iw2 == 4))then
       do ipom = 1, npom
          input_pom_c(ipom) = consts%frag_rate * pom_c(ipom)
          input_pom_n(ipom) = consts%frag_rate * pom_n(ipom)
          input_pom_p(ipom) = consts%frag_rate * pom_p(ipom)
       enddo
       
       input_dom_c = consts%frag_rate * dom_c
       input_dom_n = consts%frag_rate * dom_n
       input_dom_p = consts%frag_rate * dom_p
       
       input_nh4 = consts%frag_rate * nh4
       input_no3 = consts%frag_rate * no3
       input_psol = consts%frag_rate * psol
    endif

    return
  end subroutine wood2wood_exchange

  subroutine zero_exchange_vars(evars)
    use mend_state_vars, only: exchange_vars, npom
    implicit none
    type(exchange_vars) :: evars
    integer :: ipom

    do ipom = 1, npom
       evars%pom_c(ipom) = 0.
       evars%pom_n(ipom) = 0.
       evars%pom_p(ipom) = 0.
    enddo

    evars%dom_c = 0.
    evars%dom_n = 0.
    evars%dom_p = 0.
    
    evars%nh4 = 0.
    evars%no3 = 0.
    evars%psol = 0.

    return
  end subroutine zero_exchange_vars

  subroutine inc_exchange_vars(esum, einc)
    use mend_state_vars, only: exchange_vars, npom
    implicit none
    type(exchange_vars) :: esum
    type(exchange_vars) :: einc
    integer :: ipom

    do ipom = 1, npom
       esum%pom_c(ipom) = esum%pom_c(ipom) + einc%pom_c(ipom)
       esum%pom_n(ipom) = esum%pom_n(ipom) + einc%pom_n(ipom)
       esum%pom_p(ipom) = esum%pom_p(ipom) + einc%pom_p(ipom)
    enddo

    esum%dom_c = esum%dom_c + einc%dom_c
    esum%dom_n = esum%dom_n + einc%dom_n
    esum%dom_p = esum%dom_p + einc%dom_p
    
    esum%nh4 = esum%nh4 + einc%nh4
    esum%no3 = esum%no3 + einc%no3
    esum%psol = esum%psol + einc%psol

    return
  end subroutine inc_exchange_vars

  subroutine som2canopy_exchange(d_co2_lost, slden, consts, d_can_co2, &
       d_co2budget_storage,ccapcani)
    use mend_consts_coms, only: decomp_consts
    implicit none

    type(decomp_consts) :: consts
    real, intent(in) :: slden
    real, intent(in) :: d_co2_lost
    real :: d_co2_lost_units
    real(kind=8), intent(inout) :: d_co2budget_storage
    real(kind=8), intent(inout) :: d_can_co2
    real(kind=8), intent(in) :: ccapcani

    ! gC/kgSoil/s
    d_co2_lost_units = d_co2_lost
    ! gC/m3Soil/s
    d_co2_lost_units = d_co2_lost_units * slden
    ! gC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units * consts%eff_soil_depth
    ! molC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units / 12.
    ! umolC/m2Soil/s
    d_co2_lost_units = d_co2_lost_units * 1.0e6

    d_can_co2 = d_can_co2 + d_co2_lost_units * ccapcani
    d_co2budget_storage = d_co2budget_storage + d_co2_lost_units

    return
  end subroutine som2canopy_exchange

end Module mend_exchange
