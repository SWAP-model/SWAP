!> Main meteorological processing coordinator module
!!
!! This module serves as the primary coordinator for daily meteorological processing
!! in the SWAP model. It orchestrates the calculation of evapotranspiration components,
!! interception, and partitioning between soil evaporation and transpiration.
!!
!! The module integrates functionality from several specialized modules:
!! - meteo_process_mod: Reading and resetting meteorological data
!! - interception_mod: Various interception calculation methods
!! - et_mod: Evapotranspiration calculations (Penman-Monteith)
!! - runoff_mod: Surface runoff calculations (SCS Curve Number method)
!!
!! ## Main Functionality
!! The ProcessMeteoDay subroutine handles:
!! - Daily and sub-daily meteorological data processing
!! - Interception calculations using multiple methods (Von Hoyningen-Hune & Braden, Gash, Rutter)
!! - Reference evapotranspiration via Penman-Monteith or user input
!! - Partitioning of atmospheric demand into soil evaporation and transpiration
!! - Wet fraction calculations for crop canopy
!! - Atmospheric CO2 corrections for transpiration
!!
!! @author Original SWAP development team
!! @date Last modified February 2014, refactored February 2026
module meteo_mod

  use meteo_process_mod, only: ReadMeteoDay, ResetMetFlx
  use interception_mod, only: VonHHBraden, Gash, DivIntercep, ruttervw
  use et_mod, only: PenMon, reduceva_daily, pm_inputs_t, pm_outputs_t
  use runoff_mod, only: cn_step
  use swap_state_mod, only: swap_state_t
  use swap_config_mod, only: swap_config_t

  implicit none

  private
  public :: ProcessMeteoDay

contains

  !> Process daily meteorological data and calculate evapotranspiration components
  !!
  !! Main orchestrator for daily meteorological processing. Performs:
  !! - Interception calculations (Von Hoyningen-Hune & Braden, Gash, or Rutter methods)
  !! - Reference evapotranspiration calculation or input processing
  !! - Potential soil evaporation and transpiration partitioning
  !! - Wet fraction calculations for crop canopy
  !! - Sub-daily timestep handling for detailed meteorology
  !!
  !! ## Processing Steps
  !! 1. Calculate interception and net rainfall/irrigation
  !! 2. Loop over day parts (1 for daily, nmetdetail for sub-daily)
  !! 3. Calculate ET0, EW0, ES0 via Penman-Monteith or use specified values
  !! 4. Apply interception for adapted Rutter model (if selected)
  !! 5. Calculate wet fraction of canopy
  !! 6. Partition into potential soil evaporation and transpiration
  !! 7. Apply atmospheric CO2 corrections if enabled
  !! 8. Aggregate results for detailed meteorology
  !!
  !! @note
  !! **Original documentation:**
  !! Last modified: February 2014
  !! Supports both daily (swmetdetail=0) and detailed (swmetdetail=1) meteorology
  !! Uses module variables from Variables (MeteoVars retired in GR-ATM-CLEAN Phase D)
  !! @endnote
  subroutine ProcessMeteoDay(state, config)
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    if (config%meteo%swmetdetail == 0) then
      call process_meteo_day_daily(state, config)
    else
      call process_meteo_day_subdaily(state, config)
    end if

  end subroutine ProcessMeteoDay

  !> Private orchestrator: daily-meteo path (swmetdetail==0).
  !! Single ET record per day. Body is the daily-only slice of the former
  !! ProcessMeteoDay: Sections 3, 4, 5(daily branch), 6, 7, 9.
  subroutine process_meteo_day_daily(state, config)
    use swap_constants, only: nihil, small
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), parameter :: rcs = 0.15d0
    type(pm_outputs_t) :: pmo

    real(8) :: aintc, interc, eintc, etr, hum, win
    real(8) :: wfrac, netrainflux, rainflux
    real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond
    real(8) :: gctp

    associate (time => state%timecontrol,   &
               atmo => state%atmosphere)

    ! Daily-mode meteo scratch values were previously
    ! supplied via module-level MeteoVars from ReadMeteoDay. Now read
    ! directly from state%atmosphere arrays — same index as meteo_io.f90.
    hum = atmo%ahum(time%daymeteo+1-atmo%daynrfirst)
    win = atmo%awin(time%daymeteo+1-atmo%daynrfirst)
    etr = atmo%aetr(time%daymeteo+1-atmo%daynrfirst)

    ! === Section 3: Interception calculations ===
    call apply_interception_step(state, config, aintc)

    ! Daily mode: single iteration (irecord=1).

    ! === Section 4: Calculate evapotranspiration (et0, ew0, es0) ===
    call compute_reference_et(state, config, 1, etr, hum, win, rcs, pmo)
    Edirect     = pmo%Edirect
    Tdirect     = pmo%Tdirect
    Tdirectwet  = pmo%Tdirectwet
    Edirectpond = pmo%Edirectpond

    ! === Section 5: Interception option NHI (adapted Rutter model) ===
    ! Daily branch only — siccapact is set in the crop runtime for daily mode.
    if (state%crop%common%swinter .eq. 3) then
      if (state%crop%common%croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
        gctp = state%crop%common%gc
      else
        gctp = 1.0d0 - exp(-1.0d0*state%crop%kdir*state%crop%kdif*state%exchange%crop_water%lai)
        if (gctp .lt. 1.0d-5) atmo%siccapact = 0.0d0
      endif

      ! Calculate interception, adapted Rutter method
      call ruttervw(gctp, time%dt, atmo%siccapact, &
                    atmo%fimin, state%atmosphere%ew0, atmo%grai, &
                    atmo%sicact, aintc, eintc)

      ! Divide interception into rain and irrigation parts; compute net rain
      ! (nraida) and net sprinkling irrigation (nird).
      call DivIntercep(aintc, state)
    endif

    ! === Section 6: Fraction of the day the crop is wet ===
    call compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, 1, wfrac)

    ! === Section 7: Potential soil evaporation & transpiration ===
    call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

    ! === Section 9: Actual daily rain/snow fluxes and soil evaporation for Daily Meteo ===

    ! Finterception: ratio net / gross rain flux; net rainflux = gross - interception
    if (atmo%grai.gt.1.d-5) then
      ! finterception is exclusively meant for dividing rain flux into interception part
      ! and net rain part; not for sprinkler irrigation!
      atmo%finterception = atmo%nraida / atmo%grai
      if (aintc.lt.1.0d-5) atmo%finterception = 1.0d0
    else
      atmo%finterception = 1.0d0
    endif

    ! In case of daily precipitation sum: set actual gross and net rainflux,
    ! and interception on TIMESTEP basis
    if (config%meteo%swrain.eq.0) then
      rainflux    = atmo%fprecnosnow * atmo%grai
      netrainflux = atmo%finterception * rainflux
      atmo%graidt  = rainflux
      atmo%nraidt  = netrainflux
      if (atmo%swusecn == 1) then
        call cn_step(state)
        atmo%nraidt = atmo%nraidt - atmo%Runoff_CN
      end if
      atmo%aintcdt = rainflux - atmo%nraidt  ! aintcdt involves ONLY interception of RAIN
    endif

    ! Soil evaporation rate of today
    if (.not. time%fletsine) then
      call reduceva_daily(atmo%nraida, state)
    endif

    ! Save daily potential values for use in ETSine
    atmo%ptraday = atmo%ptra
    atmo%pevaday = atmo%peva

    ! Calculate atmospheric demand [cm]
    atmo%atmdem = state%atmosphere%et0*0.1d0

    end associate

  end subroutine process_meteo_day_daily

  !> Private orchestrator: sub-daily (detailed) meteo path (swmetdetail==1).
  !! Loops over nmetdetail records. Body is the sub-daily-only slice of the former
  !! ProcessMeteoDay: Section 3 (once); per-record Sections 4, 5(subdaily branch),
  !! 6, 7, 8; then Section 10 (daily totals).
  subroutine process_meteo_day_subdaily(state, config)
    use swap_array_dimensions, only: magrs
    use array_utils, only: afgen
    use swap_constants, only: nihil
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), parameter :: rcs = 0.15d0
    type(pm_outputs_t) :: pmo

    integer :: i, irecord, count, first, last
    real(8) :: aintc, interc, eintc, hum, win, svp
    real(8) :: wfrac, sumtav
    real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond

    associate (time => state%timecontrol,   &
               atmo => state%atmosphere)

    ! Sub-daily mode: hum/win/etr come from the per-record refresh inside
    ! compute_reference_et; etr is unused for swmetdetail==1.
    hum = 0.0d0
    win = 0.0d0

    ! === Section 3: Interception calculations ===
    call apply_interception_step(state, config, aintc)

    do irecord = 1, config%meteo%nmetdetail
      ! Refresh caller-local hum/win so that the value left after the loop
      ! reflects the last record (consumed by Section 10's rh computation).
      hum = atmo%ahum(irecord)
      win = atmo%awind_subdaily(irecord)

      ! === Section 4: Calculate evapotranspiration (et0, ew0, es0) ===
      call compute_reference_et(state, config, irecord, 0.0d0, hum, win, rcs, pmo)
      Edirect     = pmo%Edirect
      Tdirect     = pmo%Tdirect
      Tdirectwet  = pmo%Tdirectwet
      Edirectpond = pmo%Edirectpond

      ! === Section 5: Interception option NHI (adapted Rutter model) ===
      ! swinter=3 (adapted-Rutter / MetaSWAP msw1eic) is rejected by the TOML
      ! config validators, so this branch is unreachable on the TOML path; the
      ! ruttervw/msw1eic kernel was removed with the MetaSWAP drop.

      ! === Section 6: Fraction of the period the crop is wet ===
      call compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, irecord, wfrac)

      ! === Section 7: Potential soil evaporation & transpiration ===
      call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

      ! === Section 8: Results for detailed weather records ===
      atmo%tpot(irecord) = atmo%ptra
      atmo%epot(irecord) = atmo%peva
      if (atmo%grai .lt. 1.0d-12) then
        atmo%grain(irecord) = 0.0d0
        atmo%nrain(irecord) = 0.0d0
      else
        atmo%grain(irecord) = atmo%arain_subdaily(irecord) / time%metperiod
        atmo%nrain(irecord) = atmo%arain_subdaily(irecord) / time%metperiod * atmo%nraida / atmo%grai
      endif
    enddo

    ! === Section 10: Set daily weather values for Detailed Meteo ===

    ! Average temperature of today
    sumtav = 0.d0
    do i = 1, config%meteo%nmetdetail
      sumtav = sumtav + atmo%atav(i)
    enddo
    atmo%Tav = sumtav * time%metperiod

    ! Minimum and maximum temperature of today
    atmo%tmx = -50.d0
    atmo%tmn = 99.d0
    do i = 1, config%meteo%nmetdetail
      atmo%tmx = max(atmo%tmx, atmo%atav(i))
      atmo%tmn = min(atmo%tmn, atmo%atav(i))
    enddo

    ! Calculate saturated vapour pressure [kpa]
    svp = 0.3055d0*(exp(17.27d0*atmo%tmn/(atmo%tmn+237.3d0)) + &
                    exp(17.27d0*atmo%tmx/(atmo%tmx+237.3d0)))
    ! Calculate relative humidity [fraction]
    atmo%rh = min(hum/svp, 1.0d0)

    ! Average temperature between 6 and 18 hour
    sumtav = 0.d0
    count = 0
    first = int(0.25/time%metperiod) + 1
    last = int(0.75/time%metperiod)
    do i = first, last
      sumtav = sumtav + atmo%atav(i)
      count = count + 1
    enddo
    atmo%tavd = sumtav / count

    ! Daily radiation (J/m2/d) and atmospheric demand (cm/d)
    atmo%rad = 0.d0
    atmo%atmdem = 0.d0
    do i = 1, config%meteo%nmetdetail
      atmo%rad = atmo%rad + atmo%arad(i)
      atmo%atmdem            = atmo%atmdem            + atmo%tpot(i)
    enddo

    ! Fluxes of current time step (start of the day)
    atmo%ptra = atmo%tpot(1)
    atmo%peva = atmo%epot(1)
    atmo%graidt  = atmo%grain(1)
    atmo%nraidt  = atmo%nrain(1)
    atmo%aintcdt = atmo%graidt - atmo%nraidt  ! aintcdt involves ONLY interception of RAIN

    end associate

  end subroutine process_meteo_day_subdaily

  !> Private helper: Section 3 interception calculation (VonHHBraden / Gash + DivIntercep).
  !! Returns aintc and updates state%atmosphere%nraida via DivIntercep.
  subroutine apply_interception_step(state, config, aintc)
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), intent(out) :: aintc

    ! Calculation of interception and net rain & net irrigation depth [cm]
    if ((state%exchange%crop_water%lai .lt. 1.d-3) .or. (state%atmosphere%grai+state%exchange%crop_water%interc_demand .lt. 1.d-5) .or. &
        (state%crop%common%swinter.eq.0) .or. (state%atmosphere%gsnow.gt.0.0d0) .or.(state%atmosphere%ssnow.gt.0.0d0)) then

      ! No vegetation, rainfall/irrigation or interception calculation
      aintc = 0.d0

    else if (state%crop%common%swinter .eq. 1) then
      ! Calculate interception, method Von Hoyningen-Hune and Braden
      aintc = VonHHBraden(state%atmosphere%grai, state%exchange%crop_water%interc_demand, &
                          state%atmosphere%isua, state%crop%kdif, &
                          state%crop%kdir, state%exchange%crop_water%lai, state%crop%cofab)
    else if (state%crop%common%swinter .eq. 2) then
      ! Calculate interception, method Gash (1995). Evaluate the 5 AFGEN
      ! tables at the current time and pass the raw scalars to the pure
      ! Gash function.
      block
        use array_utils, only: afgen
        use swap_array_dimensions, only: magrs
        real(8) :: pfree, pstem, scanopy_raw, avprec_raw, avevap_raw
        pfree       = afgen(state%atmosphere%pfreetb,   (2*magrs), state%timecontrol%t)
        pstem       = afgen(state%atmosphere%pstemtb,   (2*magrs), state%timecontrol%t)
        scanopy_raw = afgen(state%atmosphere%scanopytb, (2*magrs), state%timecontrol%t)
        avprec_raw  = afgen(state%atmosphere%avprectb,  (2*magrs), state%timecontrol%t)
        avevap_raw  = afgen(state%atmosphere%avevaptb,  (2*magrs), state%timecontrol%t)
        aintc = Gash(state%atmosphere%grai, state%exchange%crop_water%interc_demand, state%atmosphere%isua, &
                     pfree, pstem, scanopy_raw, avprec_raw, avevap_raw)
      end block
    end if

    ! Divide interception into rain part and irrigation part and
    ! calculate net rain (nraida) and net sprinkling irrigation (nird)
    if (state%crop%common%swinter.ne.3) &
      call DivIntercep (aintc, state)

  end subroutine apply_interception_step

  !> Private helper: Section 4 reference ET (etr-direct or PenMon + crop-factor adjustments).
  !! For swmetdetail==0 .and. swetr==1: uses supplied etr directly.
  !! Otherwise: PenMon pack/call/unpack + post-call swcf/swcfbs adjustments.
  !! Writes state%atmosphere%es0/et0/ew0 and returns pmo for caller to unpack
  !! Edirect/Tdirect/Tdirectwet/Edirectpond.
  subroutine compute_reference_et(state, config, irecord, etr, hum_in, win_in, rcs, pmo)
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    integer, intent(in) :: irecord
    real(8), intent(in) :: etr, hum_in, win_in, rcs
    type(pm_outputs_t), intent(out) :: pmo

    type(pm_inputs_t) :: pmi
    real(8) :: rad_loc

    associate (time => state%timecontrol)

    ! Reference evapotranspiration has been specified
    if (config%meteo%swmetdetail.eq.0 .and. config%meteo%swetr.eq.1) then
      if (.not. state%crop%flCropEmergence) then
        ! no crop
        state%atmosphere%et0 = 0.0d0
        state%atmosphere%ew0 = 0.0d0
        state%atmosphere%es0 = etr
        if (state%crop%swcfbs.eq.1) state%atmosphere%es0 = state%crop%cfbs*etr
      else
        ! crop is present
        if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
          state%atmosphere%et0 = state%crop%common%cf*etr
          if (state%crop%swcf .eq. 1) then
            state%atmosphere%ew0 = state%crop%common%cf*etr
          else
            state%atmosphere%ew0 = state%crop%fixed%cfeic*etr
          endif
        endif
        state%atmosphere%es0 = etr
        if (state%crop%swcfbs.eq.1) state%atmosphere%es0 = state%crop%cfbs*etr
      endif

    ! Reference evapotranspiration must be calculated
    elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0) then

      if (config%meteo%swmetdetail.eq.1) then
        ! Define weather variables of current record. hum/win come from the
        ! caller (sub-daily orchestrator already indexed by irecord); rad and
        ! Tav are refreshed here because they are written into state.
        rad_loc = state%atmosphere%arad(irecord) / time%metperiod     ! from j/m2/period to j/m2/d
        state%atmosphere%Tav = state%atmosphere%atav(irecord)
      else
        rad_loc = state%atmosphere%rad
      endif

      ! Calculate evapotranspiration using Penman-Monteith: et0, ew0, es0 (mm/d)
      ! in case of daily meteo (swmetdetail = 0) irecord is always 1
      ! Pack PM inputs.
      pmi%daynr           = time%daynr
      pmi%irecord         = irecord
      pmi%nmetdetail      = config%meteo%nmetdetail
      pmi%flmetdetail     = time%flmetdetail
      pmi%flCropEmergence = state%crop%flCropEmergence
      pmi%swcf            = state%crop%swcf
      pmi%swdivide        = config%meteo%swdivide

      pmi%lat  = config%meteo%lat
      pmi%alt  = config%meteo%alt
      pmi%altw = config%meteo%altw
      pmi%a    = state%atmosphere%angstroma
      pmi%b    = state%atmosphere%angstromb
      pmi%rcs  = rcs

      pmi%rad    = rad_loc
      pmi%tav    = state%atmosphere%Tav
      pmi%tmn    = state%atmosphere%tmn
      pmi%tmx    = state%atmosphere%tmx
      pmi%hum    = hum_in
      pmi%win    = win_in
      pmi%atmtr  = state%atmosphere%atmtr
      pmi%difpp  = state%atmosphere%difpp
      pmi%dsinbe = state%atmosphere%dsinbe
      pmi%daylp  = state%atmosphere%daylp

      pmi%rsc    = state%crop%common%rsc
      pmi%rsw    = state%crop%common%rsw
      pmi%ch     = state%exchange%crop_water%crop_height
      pmi%albedo = state%crop%common%albedo
      pmi%kdif   = state%crop%kdif
      pmi%kdir   = state%crop%kdir
      pmi%lai    = state%exchange%crop_water%lai

      pmi%rsoil  = state%atmosphere%rsoil

      call PenMon(pmi, pmo)

      ! Unpack PM outputs to existing state.
      state%atmosphere%es0 = pmo%es0
      state%atmosphere%et0 = pmo%et0
      state%atmosphere%ew0 = pmo%ew0

      if (.not. state%crop%flCropEmergence) then
        ! no crop
        if (state%crop%swcfbs .eq. 1) then
          if (state%crop%swcf .eq. 1) then
            state%atmosphere%es0 = state%crop%cfbs*state%atmosphere%et0
          else
            state%atmosphere%es0 = state%crop%cfbs*state%atmosphere%es0
          endif
        endif
        state%atmosphere%et0 = 0.0d0
        if (config%meteo%swmetdetail.eq.1 .and. (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3)) then
          if (state%crop%swcf.eq.1) then
            state%atmosphere%ew0 = state%crop%common%cf*state%atmosphere%ew0
          else
            state%atmosphere%ew0 = state%crop%fixed%cfeic*state%atmosphere%ew0
          endif
        endif
      else
        ! crop is present
        if (state%crop%swcfbs .eq. 1) then
          if (state%crop%swcf .eq. 1) then
            state%atmosphere%es0 = state%crop%cfbs*state%atmosphere%et0
          else
            state%atmosphere%es0 = state%crop%cfbs*state%atmosphere%es0
          endif
        endif
        if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
          state%atmosphere%et0 = state%crop%common%cf*state%atmosphere%et0
          if (state%crop%swcf.eq.1) then
            state%atmosphere%ew0 = state%crop%common%cf*state%atmosphere%ew0
          else
            state%atmosphere%ew0 = state%crop%fixed%cfeic*state%atmosphere%ew0
          endif
        endif
      endif

    endif

    end associate

  end subroutine compute_reference_et

  !> Private helper: Section 6 wet-fraction calculation.
  !! For swmetdetail==0: simple ratio based on aintc/eintc and ew0.
  !! For swmetdetail==1: accumulates interc across sub-daily records and
  !! updates state%atmosphere%restint at end of each call.
  subroutine compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, irecord, wfrac)
    use swap_constants, only: nihil
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    real(8), intent(in)    :: aintc, eintc, Tdirectwet
    real(8), intent(inout) :: interc
    integer, intent(in)    :: irecord
    real(8), intent(out)   :: wfrac

    associate (time => state%timecontrol)

    ! Calculate fraction of the day or period the crop is wet
    if (config%meteo%swmetdetail.eq.0) then
      ! Fraction of the day the crop is wet
      if (state%atmosphere%ew0.lt.0.0001d0) then
        wfrac = 0.0d0
      else
        if (state%atmosphere%ew0.lt.0.0001d0) then
          wfrac = 0.0d0
        else
          if (state%crop%common%swinter .ne. 3) then
            if (config%meteo%swdivide .eq. 0) then
              wfrac = max(min(aintc*10.0d0/state%atmosphere%ew0,1.0d0),0.0d0)
            else
              if(tdirectwet.gt.nihil) then
                wfrac = max(min(aintc*10.0d0/tdirectwet,1.0d0),0.0d0)
              else
                wfrac = 0.0d0
              endif
            endif
          else
            wfrac = max(min(eintc*10.0d0/state%atmosphere%ew0,1.0d0),0.0d0)
          endif
        endif
      endif
    ! Fraction of the period the crop is wet
    elseif (config%meteo%swmetdetail.eq.1) then
      if (state%atmosphere%grai .lt. 1.0d-12) then
        interc = 0.0d0
        wfrac  = 0.0d0
      else
        interc = state%atmosphere%restint + aintc * state%atmosphere%arain_subdaily(irecord) / state%atmosphere%grai
        if (state%atmosphere%ew0.lt.0.0001d0) then
          wfrac = 0.0d0
        else
          if (state%crop%swcf.ne.3) then
            wfrac = max(min(interc*10.0d0/state%atmosphere%ew0/time%metperiod,1.d0),0.d0)
          else
            wfrac = max(min(eintc/state%atmosphere%ew0,1.0d0),0.0d0)
          endif
        endif
      endif
      ! Remaining amount of interception for swmetdetail = 1
      state%atmosphere%restint = max(interc - wfrac * time%metperiod * state%atmosphere%ew0 * 0.1d0, 0.d0)
    endif

    end associate

  end subroutine compute_wet_fraction

  !> Private helper: Section 7 partition into peva/ptra.
  !! Applies cover-fraction correction, ponding correction, PMdirect override,
  !! and CO2 correction. Writes state%atmosphere%peva and state%atmosphere%ptra.
  subroutine partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)
    use swap_constants, only: nihil, small
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    real(8), intent(in) :: wfrac, Edirect, Tdirect, Edirectpond

    associate (atmo => state%atmosphere)

    ! Potential soil evaporation (peva) [cm/d]
    atmo%peva = max(0.0d0, (state%atmosphere%es0*exp(-1.0d0*state%crop%kdir*state%crop%kdif*state%exchange%crop_water%lai)*0.1d0))
    if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
      atmo%peva = max(0.0d0,(1.0d0-wfrac)*atmo%peva)
    end if

    ! Alternative for peva (simple model, soil cover fraction specified)
    if (state%crop%common%flCropCalendar .and. .not. state%crop%common%flCropHarvest) then
      if (state%crop%common%croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
        atmo%peva = (1.0d0 - state%crop%common%gc)*state%atmosphere%es0*0.1d0
        if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
          atmo%peva = (1.0d0-wfrac)*atmo%peva
        end if
      endif
    endif

    ! Adapt peva in case of ponding — [SS-SWC S-2.12B] state%exchange%atmos_soil%pond
    if (state%exchange%atmos_soil%pond .gt. 1.0d-10) then
      if (config%meteo%swetr.eq.0 .and. state%atmosphere%es0.gt.1.0d-8) then
        atmo%peva = state%atmosphere%ew0/state%atmosphere%es0 * atmo%peva
      elseif (state%atmosphere%es0.gt.1.0d-8) then
        if (state%crop%swcfbs .eq. 1 .and. state%crop%cfbs .gt. small) then
          atmo%peva = atmo%cfevappond * atmo%peva / state%crop%cfbs
        else
          atmo%peva = atmo%cfevappond * atmo%peva
        endif
      endif
    endif

    ! Potential soil evaporation [cm/d] according to PMdirect
    if (config%meteo%swdivide .eq. 1) then
      if (state%exchange%atmos_soil%pond .gt. 1.0d-10) then
        atmo%peva = Edirectpond*0.1d0
      else
        atmo%peva = Edirect*0.1d0
      endif
    endif

    ! Potential transpiration (ptra) [cm/d]
    if (state%crop%swcf .ne. 3) then
      atmo%ptra = ((1.0d0-wfrac)*state%atmosphere%et0-atmo%peva*10.0d0)*0.1d0
    else
      atmo%ptra = (1.0d0-wfrac)*state%atmosphere%et0*0.1d0
    endif
    atmo%ptra = max(atmo%ptra,(1.01d0*nihil))

    ! Potential transpiration [cm/d] according to PMdirect
    if (config%meteo%swdivide .eq. 1) then
      atmo%ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
      atmo%ptra = max(atmo%ptra,(1.01d0*nihil))
    endif

    ! Correction of potential transpiration as a function of atmospheric CO2 concentration
    if (atmo%flco2 .and. state%crop%flCropEmergence) then
      atmo%ptra = state%exchange%crop_water%co2_transp_fac * atmo%ptra
    endif

    end associate

  end subroutine partition_peva_ptra

end module meteo_mod










