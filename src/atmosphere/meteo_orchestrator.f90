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
  use interception_mod, only: VonHHBraden, Gash, ruttervw, DivIntercep
  use et_mod, only: PenMon, reduceva_daily, pm_inputs_t, pm_outputs_t
  use runoff_mod, only: cn_step
  use swap_state_mod, only: swap_state_t
  use swap_config_mod, only: swap_config_t   ! [SS-GR-ATM B24] config added for meteo switches

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
    use variables, only: rad, croptype, gc, flCropHarvest, finterception, swusecn
    use swap_constants, only: nihil, small
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), parameter :: rcs = 0.15d0
    type(pm_outputs_t) :: pmo

    real(8) :: aintc, interc, eintc, dttp, gctp, etr, hum, win
    real(8) :: wfrac, netrainflux, rainflux
    real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond

    associate( &
       at_peva     => state%atmosphere%peva,     &
       at_ptra     => state%atmosphere%ptra,     &
       at_atmdem   => state%atmosphere%atmdem,   &
       at_pevaday  => state%atmosphere%pevaday,  &
       at_ptraday  => state%atmosphere%ptraday,  &
       tc_fletsine => state%timecontrol%fletsine, &
       tc_daymeteo => state%timecontrol%daymeteo  )

    ! [GR-ATM-CLEAN D.2] Daily-mode meteo scratch values were previously
    ! supplied via module-level MeteoVars from ReadMeteoDay. Now read
    ! directly from state%atmosphere arrays — same index as meteo_io.f90.
    hum = state%atmosphere%ahum(tc_daymeteo+1-state%atmosphere%daynrfirst)
    win = state%atmosphere%awin(tc_daymeteo+1-state%atmosphere%daynrfirst)
    etr = state%atmosphere%aetr(tc_daymeteo+1-state%atmosphere%daynrfirst)

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
    ! Daily branch only — siccapact is set in cropgrowth module for daily mode.
    if (state%crop%common%swinter .eq. 3) then
      if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then   ! [GR-CROP C3]
        gctp  = gc
      else
        gctp  = 1.0d0 - exp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)
        if (gctp .lt. 1.0d-5) then
          state%atmosphere%siccapact = 0.0d0   ! [SS-GR-ATM B24] direct state write
        endif
      endif
      dttp = 1.0d0  ! value of 1 d required for the daily meteo option

      ! Calculate interception, method Rutter
      call ruttervw(gctp, state%timecontrol%dt, state%atmosphere%siccapact, &
                    state%atmosphere%fimin, state%crop%ew0, state%atmosphere%grai, &
                    state%atmosphere%sicact, aintc, eintc)

      ! Divide interception into rain part and irrigation part and
      ! calculate net rain (nraida) and net sprinkling irrigation (nird)
      call DivIntercep(aintc, state)
    endif

    ! === Section 6: Fraction of the day the crop is wet ===
    call compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, 1, wfrac)

    ! === Section 7: Potential soil evaporation & transpiration ===
    call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

    ! === Section 9: Actual daily rain/snow fluxes and soil evaporation for Daily Meteo ===

    ! Finterception: ratio net / gross rain flux; net rainflux = gross - interception
    if (state%atmosphere%grai.gt.1.d-5) then
      ! finterception is exclusively meant for dividing rain flux into interception part
      ! and net rain part; not for sprinkler irrigation!
      finterception = state%atmosphere%nraida / state%atmosphere%grai
      if (aintc.lt.1.0d-5) finterception = 1.0d0
    else
      finterception = 1.0d0
    endif

    ! In case of daily precipitation sum: set actual gross and net rainflux,
    ! and interception on TIMESTEP basis
    if (config%meteo%swrain.eq.0) then
      rainflux    = state%atmosphere%fprecnosnow * state%atmosphere%grai
      netrainflux = finterception * rainflux
      state%atmosphere%graidt  = rainflux
      state%atmosphere%nraidt  = netrainflux
      if (swuseCN == 1) then
        ! SS-ATM A-2.6: state added — CNmethod needs nraidt/melt from state%atmosphere
        call cn_step(state)
        state%atmosphere%nraidt = state%atmosphere%nraidt - state%atmosphere%Runoff_CN   ! [SS-GR-ATM B24]
      end if
      state%atmosphere%aintcdt = rainflux - state%atmosphere%nraidt  ! aintcdt involves ONLY interception of RAIN
    endif

    ! Soil evaporation rate of today
    if (.not. tc_fletsine) then
      call reduceva_daily(state%atmosphere%nraida, state)
    endif

    ! Save daily potential values for use in ETSine
    at_ptraday = at_ptra
    at_pevaday = at_peva

    ! Calculate atmospheric demand [cm]
    at_atmdem = state%crop%et0*0.1d0

    end associate

  end subroutine process_meteo_day_daily

  !> Private orchestrator: sub-daily (detailed) meteo path (swmetdetail==1).
  !! Loops over nmetdetail records. Body is the sub-daily-only slice of the former
  !! ProcessMeteoDay: Section 3 (once); per-record Sections 4, 5(subdaily branch),
  !! 6, 7, 8; then Section 10 (daily totals).
  subroutine process_meteo_day_subdaily(state, config)
    use variables, only: rad, tmn, tmx, croptype, gc, siccaptb, tav
    use swap_array_dimensions, only: magrs
    use array_utils, only: afgen
    use swap_constants, only: nihil
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), parameter :: rcs = 0.15d0
    type(pm_outputs_t) :: pmo

    integer :: i, irecord, count, first, last
    real(8) :: aintc, interc, eintc, dttp, gctp, hum, win, svp
    real(8) :: wfrac, sumtav
    real(8) :: Edirect, Tdirect, Tdirectwet, Edirectpond

    associate( &
       at_peva   => state%atmosphere%peva,   &
       at_ptra   => state%atmosphere%ptra,   &
       at_atmdem => state%atmosphere%atmdem, &
       tc_t      => state%timecontrol%t,     &
       tc_dt     => state%timecontrol%dt,    &
       metperiod => state%timecontrol%metperiod )

    ! Sub-daily mode: hum/win/etr come from the per-record refresh inside
    ! compute_reference_et; etr is unused for swmetdetail==1.
    hum = 0.0d0
    win = 0.0d0

    ! === Section 3: Interception calculations ===
    call apply_interception_step(state, config, aintc)

    do irecord = 1, config%meteo%nmetdetail
      ! Refresh caller-local hum/win so that the value left after the loop
      ! reflects the last record (consumed by Section 10's rh computation).
      hum = state%atmosphere%ahum(irecord)
      win = state%atmosphere%awind_subdaily(irecord)

      ! === Section 4: Calculate evapotranspiration (et0, ew0, es0) ===
      call compute_reference_et(state, config, irecord, 0.0d0, hum, win, rcs, pmo)
      Edirect     = pmo%Edirect
      Tdirect     = pmo%Tdirect
      Tdirectwet  = pmo%Tdirectwet
      Edirectpond = pmo%Edirectpond

      ! === Section 5: Interception option NHI (adapted Rutter model) ===
      ! Sub-daily branch only.
      if (state%crop%common%swinter .eq. 3) then
        state%atmosphere%siccapact = afgen(siccaptb,(2*magrs),tc_t)   ! [SS-GR-ATM B24] direct state write
        if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
          gctp  = gc
        elseif (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.1) then
          gctp  = 1.0d0 - exp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)
        endif
        if (gctp .lt. 1.0d-5) then
          state%atmosphere%siccapact = 0.0d0   ! [SS-GR-ATM B24] direct state write
        endif
        dttp = tc_dt

        ! Calculate interception, method Rutter
        call ruttervw(gctp, state%timecontrol%dt, state%atmosphere%siccapact, &
                      state%atmosphere%fimin, state%crop%ew0, state%atmosphere%grai, &
                      state%atmosphere%sicact, aintc, eintc)

        ! Divide interception into rain part and irrigation part and
        ! calculate net rain (nraida) and net sprinkling irrigation (nird)
        call DivIntercep(aintc, state)
      endif

      ! === Section 6: Fraction of the period the crop is wet ===
      call compute_wet_fraction(state, config, aintc, eintc, Tdirectwet, interc, irecord, wfrac)

      ! === Section 7: Potential soil evaporation & transpiration ===
      call partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)

      ! === Section 8: Results for detailed weather records ===
      state%atmosphere%tpot(irecord) = at_ptra
      state%atmosphere%epot(irecord) = at_peva
      if (state%atmosphere%grai .lt. 1.0d-12) then
        state%atmosphere%grain(irecord) = 0.0d0
        state%atmosphere%nrain(irecord) = 0.0d0
      else
        state%atmosphere%grain(irecord) = state%atmosphere%arain_subdaily(irecord) / metperiod
        state%atmosphere%nrain(irecord) = state%atmosphere%arain_subdaily(irecord) / metperiod * state%atmosphere%nraida / state%atmosphere%grai
      endif
    enddo

    ! === Section 10: Set daily weather values for Detailed Meteo ===

    ! Average temperature of today
    sumtav = 0.d0
    do i = 1, config%meteo%nmetdetail
      sumtav = sumtav + state%atmosphere%atav(i)
    enddo
    state%atmosphere%Tav = sumtav * metperiod
    tav = state%atmosphere%Tav   ! [SS-GR-ATM B24] dual-write — legacy tav consumed by snow.f90/swapoutput.f90

    ! Minimum and maximum temperature of today
    tmx = -50.d0
    tmn = 99.d0
    do i = 1, config%meteo%nmetdetail
      tmx = max(tmx,state%atmosphere%atav(i))
      tmn = min(tmn,state%atmosphere%atav(i))
    enddo

    ! Calculate saturated vapour pressure [kpa]
    svp = 0.3055d0*(exp(17.27d0*tmn/(tmn+237.3d0)) + &
                    exp(17.27d0*tmx/(tmx+237.3d0)))
    ! Calculate relative humidity [fraction]
    state%atmosphere%rh = min(hum/svp,1.0d0)   ! [SS-GR-ATM B24] direct state write
    ! [SS-GR-ATM B.5] rh dual-write to legacy global RETIRED: no crop consumers remain

    ! Average temperature between 6 and 18 hour
    sumtav = 0.d0
    count = 0
    first = int(0.25/metperiod) + 1
    last = int(0.75/metperiod)
    do i = first, last
      sumtav = sumtav + state%atmosphere%atav(i)
      count = count + 1
    enddo
    state%atmosphere%tavd = sumtav / count   ! [SS-GR-ATM B24] direct state write
    ! [SS-GR-ATM B.5] tavd dual-write to legacy global RETIRED: no crop consumers remain

    ! Daily radiation (J/m2/d) and atmospheric demand (cm/d)
    rad = 0.d0
    at_atmdem = 0.d0
    do i = 1,config%meteo%nmetdetail
      rad = rad + state%atmosphere%arad(i)
      at_atmdem = at_atmdem + state%atmosphere%tpot(i)
    enddo

    ! Fluxes of current time step (start of the day)
    at_ptra = state%atmosphere%tpot(1)
    at_peva = state%atmosphere%epot(1)
    state%atmosphere%graidt  = state%atmosphere%grain(1)
    state%atmosphere%nraidt  = state%atmosphere%nrain(1)
    state%atmosphere%aintcdt = state%atmosphere%graidt - state%atmosphere%nraidt  ! aintcdt involves ONLY interception of RAIN

    end associate

  end subroutine process_meteo_day_subdaily

  !> Private helper: Section 3 interception calculation (VonHHBraden / Gash + DivIntercep).
  !! Returns aintc and updates state%atmosphere%nraida via DivIntercep.
  subroutine apply_interception_step(state, config, aintc)
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config

    real(8), intent(out) :: aintc

    ! Calculation of interception and net rain & net irrigation depth [cm]
    if ((state%crop%lai .lt. 1.d-3) .or. (state%atmosphere%grai+state%crop%gird .lt. 1.d-5) .or. &
        (state%crop%common%swinter.eq.0) .or. (state%atmosphere%gsnow.gt.0.0d0) .or.(state%atmosphere%ssnow.gt.0.0d0)) then

      ! No vegetation, rainfall/irrigation or interception calculation
      aintc = 0.d0

    else if (state%crop%common%swinter .eq. 1) then
      ! Calculate interception, method Von Hoyningen-Hune and Braden
      aintc = VonHHBraden(state%atmosphere%grai, state%crop%gird, &
                          state%atmosphere%isua, state%crop%kdif, &
                          state%crop%kdir, state%crop%lai, state%crop%cofab)
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
        aintc = Gash(state%atmosphere%grai, state%crop%gird, state%atmosphere%isua, &
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
  !! Writes state%crop%es0/et0/ew0 and returns pmo for caller to unpack
  !! Edirect/Tdirect/Tdirectwet/Edirectpond.
  subroutine compute_reference_et(state, config, irecord, etr, hum_in, win_in, rcs, pmo)
    use variables, only: rad, angstroma, angstromb, daylp, tmn, tmx, difpp, &
                         dsinbe, atmtr, rsoil
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    integer, intent(in) :: irecord
    real(8), intent(in) :: etr, hum_in, win_in, rcs
    type(pm_outputs_t), intent(out) :: pmo

    type(pm_inputs_t) :: pmi
    real(8) :: rad_loc

    associate( &
       tc_daynr       => state%timecontrol%daynr,       &
       tc_flmetdetail => state%timecontrol%flmetdetail, &
       metperiod      => state%timecontrol%metperiod    )

    ! Reference evapotranspiration has been specified
    if (config%meteo%swmetdetail.eq.0 .and. config%meteo%swetr.eq.1) then
      if (.not. state%crop%flCropEmergence) then
        ! no crop
        state%crop%et0 = 0.0d0
        state%crop%ew0 = 0.0d0
        state%crop%es0 = etr
        if (state%crop%swcfbs.eq.1) state%crop%es0 = state%crop%cfbs*etr
      else
        ! crop is present
        if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
          state%crop%et0 = state%crop%common%cf*etr
          if (state%crop%swcf .eq. 1) then
            state%crop%ew0 = state%crop%common%cf*etr
          else
            state%crop%ew0 = state%crop%fixed%cfeic*etr
          endif
        endif
        state%crop%es0 = etr
        if (state%crop%swcfbs.eq.1) state%crop%es0 = state%crop%cfbs*etr
      endif

    ! Reference evapotranspiration must be calculated
    elseif (config%meteo%swmetdetail.eq.1 .or. config%meteo%swetr.eq.0) then

      if (config%meteo%swmetdetail.eq.1) then
        ! Define weather variables of current record. hum/win come from the
        ! caller (sub-daily orchestrator already indexed by irecord); rad and
        ! Tav are refreshed here because they are written into state.
        rad_loc = state%atmosphere%arad(irecord) / metperiod     ! from j/m2/period to j/m2/d
        state%atmosphere%Tav = state%atmosphere%atav(irecord)
      else
        rad_loc = rad
      endif

      ! Calculate evapotranspiration using Penman-Monteith: et0, ew0, es0 (mm/d)
      ! in case of daily meteo (swmetdetail = 0) irecord is always 1
      ! Pack PM inputs.
      pmi%daynr           = tc_daynr
      pmi%irecord         = irecord
      pmi%nmetdetail      = config%meteo%nmetdetail
      pmi%flmetdetail     = tc_flmetdetail
      pmi%flCropEmergence = state%crop%flCropEmergence
      pmi%swcf            = state%crop%swcf
      pmi%swdivide        = state%cfg%meteo%swdivide

      pmi%lat  = state%cfg%meteo%lat
      pmi%alt  = state%cfg%meteo%alt
      pmi%altw = state%cfg%meteo%altw
      pmi%a    = angstroma
      pmi%b    = angstromb
      pmi%rcs  = rcs

      pmi%rad    = rad_loc
      pmi%tav    = state%atmosphere%Tav
      pmi%tmn    = tmn
      pmi%tmx    = tmx
      pmi%hum    = hum_in
      pmi%win    = win_in
      pmi%atmtr  = atmtr
      pmi%difpp  = difpp
      pmi%dsinbe = dsinbe
      pmi%daylp  = daylp

      pmi%rsc    = state%crop%common%rsc
      pmi%rsw    = state%crop%common%rsw
      pmi%ch     = state%crop%common%ch
      pmi%albedo = state%crop%common%albedo
      pmi%kdif   = state%crop%kdif
      pmi%kdir   = state%crop%kdir
      pmi%lai    = state%crop%lai

      pmi%rsoil  = rsoil

      call PenMon(pmi, pmo)

      ! Unpack PM outputs to existing state.
      state%crop%es0 = pmo%es0
      state%crop%et0 = pmo%et0
      state%crop%ew0 = pmo%ew0

      if (.not. state%crop%flCropEmergence) then
        ! no crop
        if (state%crop%swcfbs .eq. 1) then
          if (state%crop%swcf .eq. 1) then
            state%crop%es0 = state%crop%cfbs*state%crop%et0
          else
            state%crop%es0 = state%crop%cfbs*state%crop%es0
          endif
        endif
        state%crop%et0 = 0.0d0
        if (config%meteo%swmetdetail.eq.1 .and. (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3)) then
          if (state%crop%swcf.eq.1) then
            state%crop%ew0 = state%crop%common%cf*state%crop%ew0
          else
            state%crop%ew0 = state%crop%fixed%cfeic*state%crop%ew0
          endif
        endif
      else
        ! crop is present
        if (state%crop%swcfbs .eq. 1) then
          if (state%crop%swcf .eq. 1) then
            state%crop%es0 = state%crop%cfbs*state%crop%et0
          else
            state%crop%es0 = state%crop%cfbs*state%crop%es0
          endif
        endif
        if (state%crop%swcf.eq.1 .or. state%crop%swcf.eq.3) then
          state%crop%et0 = state%crop%common%cf*state%crop%et0
          if (state%crop%swcf.eq.1) then
            state%crop%ew0 = state%crop%common%cf*state%crop%ew0
          else
            state%crop%ew0 = state%crop%fixed%cfeic*state%crop%ew0
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

    associate( metperiod => state%timecontrol%metperiod )

    ! Calculate fraction of the day or period the crop is wet
    if (config%meteo%swmetdetail.eq.0) then
      ! Fraction of the day the crop is wet
      if (state%crop%ew0.lt.0.0001d0) then
        wfrac = 0.0d0
      else
        if (state%crop%ew0.lt.0.0001d0) then
          wfrac = 0.0d0
        else
          if (state%crop%common%swinter .ne. 3) then
            if (state%cfg%meteo%swdivide .eq. 0) then
              wfrac = max(min(aintc*10.0d0/state%crop%ew0,1.0d0),0.0d0)
            else
              if(tdirectwet.gt.nihil) then
                wfrac = max(min(aintc*10.0d0/tdirectwet,1.0d0),0.0d0)
              else
                wfrac = 0.0d0
              endif
            endif
          else
            wfrac = max(min(eintc*10.0d0/state%crop%ew0,1.0d0),0.0d0)
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
        if (state%crop%ew0.lt.0.0001d0) then
          wfrac = 0.0d0
        else
          if (state%crop%swcf.ne.3) then
            wfrac = max(min(interc*10.0d0/state%crop%ew0/metperiod,1.d0),0.d0)
          else
            wfrac = max(min(eintc/state%crop%ew0,1.0d0),0.0d0)
          endif
        endif
      endif
      ! Remaining amount of interception for swmetdetail = 1
      state%atmosphere%restint = max(interc - wfrac * metperiod * state%crop%ew0 * 0.1d0, 0.d0)
    endif

    end associate

  end subroutine compute_wet_fraction

  !> Private helper: Section 7 partition into peva/ptra.
  !! Applies cover-fraction correction, ponding correction, PMdirect override,
  !! and CO2 correction. Writes state%atmosphere%peva and state%atmosphere%ptra.
  subroutine partition_peva_ptra(state, config, wfrac, Edirect, Tdirect, Edirectpond)
    use swap_constants, only: nihil, small
    use variables, only: cfevappond, flco2, croptype, flCropHarvest, gc
    type(swap_state_t),  intent(inout) :: state
    type(swap_config_t), intent(in)    :: config
    real(8), intent(in) :: wfrac, Edirect, Tdirect, Edirectpond

    associate( &
       at_peva => state%atmosphere%peva, &
       at_ptra => state%atmosphere%ptra )

    ! Potential soil evaporation (peva) [cm/d]
    at_peva = max(0.0d0, (state%crop%es0*exp(-1.0d0*state%crop%kdir*state%crop%kdif*state%crop%lai)*0.1d0))
    if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
      at_peva = max(0.0d0,(1.0d0-wfrac)*at_peva)
    end if

    ! Alternative for peva (simple model, soil cover fraction specified)
    if (state%crop%common%flCropCalendar .and. .not.flCropHarvest) then   ! [GR-CROP C3]
      if (croptype(state%crop%common%icrop).eq.1 .and. state%crop%common%swgc.eq.2) then
        at_peva = (1.0d0-gc)*state%crop%es0*0.1d0
        if (state%crop%swcf.ne.3 .or. (config%meteo%swmetdetail.eq.0 .and. state%crop%common%swinter.ne.3)) then
          at_peva = (1.0d0-wfrac)*at_peva
        end if
      endif
    endif

    ! Adapt peva in case of ponding — [SS-SWC S-2.12B] state%soilwater%pond
    if (state%soilwater%pond .gt. 1.0d-10) then
      if (config%meteo%swetr.eq.0 .and. state%crop%es0.gt.1.0d-8) then
        at_peva = state%crop%ew0/state%crop%es0 * at_peva
      elseif (state%crop%es0.gt.1.0d-8) then
        if (state%crop%swcfbs .eq. 1 .and. state%crop%cfbs .gt. small) then
          at_peva = cfevappond * at_peva / state%crop%cfbs
        else
          at_peva = cfevappond * at_peva
        endif
      endif
    endif

    ! Potential soil evaporation [cm/d] according to PMdirect
    if (state%cfg%meteo%swdivide .eq. 1) then
      if (state%soilwater%pond .gt. 1.0d-10) then  ! [SS-SWC S-2.12B]
        at_peva = Edirectpond*0.1d0
      else
        at_peva = Edirect*0.1d0
      endif
    endif

    ! Potential transpiration (ptra) [cm/d]
    if (state%crop%swcf .ne. 3) then
      at_ptra = ((1.0d0-wfrac)*state%crop%et0-at_peva*10.0d0)*0.1d0
    else
      at_ptra = (1.0d0-wfrac)*state%crop%et0*0.1d0
    endif
    at_ptra = max(at_ptra,(1.01d0*nihil))

    ! Potential transpiration [cm/d] according to PMdirect
    if (state%cfg%meteo%swdivide .eq. 1) then
      at_ptra = (1.0d0-wfrac) * Tdirect * 0.1d0
      at_ptra = max(at_ptra,(1.01d0*nihil))
    endif

    ! Correction of potential transpiration as a function of atmospheric CO2 concentration
    if (flCO2 .and. state%crop%flCropEmergence) then
      at_ptra = state%crop%wofost%fco2tra * at_ptra
    endif

    end associate

  end subroutine partition_peva_ptra

end module meteo_mod










