!> Runtime initialization for the type-3 (grass / WOFOST grass) rotation,
!! replacing the per-rotation runtime side-effects of legacy `readgrass`.
!!
!! Operates on `cropgrass_config_t` (already validated and populated at
!! config-load time) and writes to the same `variables`-module globals
!! that legacy `readgrass` writes to.  After the copy, runs the
!! tail-of-readgrass init math (cumdens build, dateharvest conversion,
!! DelayRegrowthTab packing).
!!
!! Note: `swharvest`, `dmharvest`, `daylastharvest`, `dmlastharvest`,
!! `swdmmow`, `maxdaymow`, `swlossmow`, `swlossgrz`, `swdmgrz`,
!! `maxdaygrz`, `dmgrazing`, `lsdb`, `tagprest`, `swhydrlift` are LOCAL
!! save variables in `grass()` (cropgrowth.f90), NOT in the variables
!! module.  They are NOT returned from this subroutine; Task 8 (dispatch
!! wiring in cropgrowth.f90) will read them from the typed config and
!! assign them directly to the `grass()` locals.
!!
!! The `dateharvest` and `seqgrazmow` globals ARE in the variables module
!! and are populated here.
!!
!! Scope (Phase 3, ADR 0015): only the supported subset of switches.
!! Defense-in-depth runtime guards mirror the validator stub-errors.
module cropgrass_init_mod
   use iso_fortran_env, only: real64
   use cropgrass_config_mod, only: cropgrass_config_t
   implicit none
   private

   public :: cropgrass_init_from_config

contains

   subroutine cropgrass_init_from_config(cfg, icrop)
      use variables, only: &
         ! ET-related
         swcf, albedo, rsw, rsc, cftb, chtb,                                &
         ! Interception
         swinter, cofab,                                                      &
         ! Crop state
         tdwi, laiem, rgrlai,                                                 &
         ! Start-of-growth trigger
         swtsum, tsumtemp, tsumtime, tsumdepth,                               &
         ! Green area
         slatb, ssa, span, tbase,                                             &
         ! Assimilation
         kdif, kdir, eff, amaxtb, tmpftb, tmnftb,                            &
         ! Biomass conversion
         cvl, cvr, cvs,                                                       &
         ! Maintenance respiration
         q10, rml, rmr, rms, rfsetb,                                         &
         ! Partitioning
         frtb, fltb, fstb,                                                    &
         ! Death rates
         perdl, rdrrtb, rdrstb,                                               &
         ! Root depth and density
         swrd, swdmi2rd, swrdc, rdi, rri, rdc, rdctb, rdtb, rlwtb, wrtmax,  &
         cumdens,                                                             &
         ! Oxygen stress
         swoxygen, swWrtNonox, aeratecrit,                                    &
         hlim1, hlim2u, hlim2l,                                              &
         q10_microbial, specific_resp_humus, srl, swrootradius,              &
         dry_mat_cont_roots, air_filled_root_por, spec_weight_root_tissue,   &
         var_a, root_radiusO2, swoxygentype,                                  &
         ! Drought stress
         swdrought, hlim3h, hlim3l, hlim4, adcrh, adcrl,                    &
         ! Salinity stress (guarded; set to 0 only)
         swsalinity,                                                          &
         ! Compensation
         swcompensate, swstressor, alphacrit, dcritrtz,                      &
         ! Management
         mowrest, swpotrelmf, relmf, seqgrazmow,                             &
         ! Mowing / harvest
         dateharvest, dmmowtb, DelayRegrowthTab,                             &
         ! CO2 (flCO2 only; swco2 is a local in readgrass, not a global)
         flCO2,                                                               &
         ! Irrigation scheduling (set to 0; schedule=1 stub-errored)
         schedule,                                                             &
         ! Simulation time globals needed for dateharvest sentinel
         tend, yearmeteo
      use array_utils,  only: afgen
      use error_mod,    only: fatalerr_collected
      implicit none

      type(cropgrass_config_t), intent(in) :: cfg
      integer,                  intent(in) :: icrop  ! rotation slot (reserved)

      integer      :: i
      real(real64) :: depth, sum_val
      real(real64) :: rootdis(202)

      ! ----------------------------------------------------------------
      ! Defense-in-depth guards — the validator should have caught these,
      ! but we guard at runtime too per ADR 0015.
      ! ----------------------------------------------------------------
      if (cfg%swdrought == 2) &
         call fatalerr_collected('cropgrass_init', &
            'swdrought=2 not supported on TOML path; validator should have rejected.')
      if (cfg%swsalinity /= 0) &
         call fatalerr_collected('cropgrass_init', &
            'swsalinity/=0 not supported on TOML path; validator should have rejected.')
      if (cfg%swcompensate == 2) &
         call fatalerr_collected('cropgrass_init', &
            'swcompensate=2 not supported on TOML path; validator should have rejected.')
      if (cfg%swco2 == 1) &
         call fatalerr_collected('cropgrass_init', &
            'swco2=1 not supported on TOML path; validator should have rejected.')
      if (cfg%swlossmow == 1) &
         call fatalerr_collected('cropgrass_init', &
            'swlossmow=1 not supported on TOML path; validator should have rejected.')
      if (cfg%swlossgrz == 1) &
         call fatalerr_collected('cropgrass_init', &
            'swlossgrz=1 not supported on TOML path; validator should have rejected.')
      if (cfg%swoxygen == 2 .and. cfg%swoxygentype == 2) &
         call fatalerr_collected('cropgrass_init', &
            'swoxygen=2 swoxygentype=2 not supported on TOML path; validator should have rejected.')
      if (cfg%swrd == 1) &
         call fatalerr_collected('cropgrass_init', &
            'swrd=1 not supported on TOML path; validator should have rejected.')
      if (cfg%swrd == 3) &
         call fatalerr_collected('cropgrass_init', &
            'swrd=3 not supported on TOML path; validator should have rejected.')
      if (cfg%swcf == 3) &
         call fatalerr_collected('cropgrass_init', &
            'swcf=3 not supported on TOML path; validator should have rejected.')
      if (cfg%swrdc == 1) &
         call fatalerr_collected('cropgrass_init', &
            'swrdc=1 not supported on TOML path; validator should have rejected.')
      if (cfg%swinter == 2 .or. cfg%swinter == 3) &
         call fatalerr_collected('cropgrass_init', &
            'swinter=2 or 3 not supported on TOML path; validator should have rejected.')
      if (cfg%swtsum == 2) &
         call fatalerr_collected('cropgrass_init', &
            'swtsum=2 not supported on TOML path; validator should have rejected.')
      if (cfg%schedule%schedule == 1) &
         call fatalerr_collected('cropgrass_init', &
            'schedule=1 not supported on TOML path; validator should have rejected.')
      if (allocated(cfg%seqgrazmow)) then
         block
            integer :: jseq
            do jseq = 1, size(cfg%seqgrazmow)
               if (cfg%seqgrazmow(jseq) /= 2) then
                  call fatalerr_collected('cropgrass_init', &
                     'grazing/dewooling (seqgrazmow/=2) not supported on TOML path; ' // &
                     'validator should have rejected.')
                  exit
               end if
            end do
         end block
      end if

      ! ================================================================
      ! Config → module globals copy (mirror readgrass rd* calls 1:1)
      ! ================================================================

      ! Part 1: crop factor / crop height (readgrass lines 3508-3557)
      ! Legacy: swcf=1 → cftb; swcf=2 → chtb; swcf=3 guarded above.
      swcf = cfg%swcf
      if (cfg%swcf == 1) then
         ! ETref standard defaults for albedo/rsc/rsw
         albedo = 0.23d0
         rsc    = 70.0d0
         rsw    = 0.0d0
         if (allocated(cfg%cftb)) call copy_table(cfg%cftb, cftb)
         chtb = -99.99d0
      else if (cfg%swcf == 2) then
         albedo = cfg%albedo
         rsc    = cfg%rsc
         rsw    = cfg%rsw
         if (allocated(cfg%chtb)) call copy_table(cfg%chtb, chtb)
         cftb = -99.99d0
      end if

      ! Part 2: interception (readgrass lines 3560-3585)
      swinter = cfg%swinter
      if (cfg%swinter == 1) cofab = cfg%cofab

      ! Part 3: initial crop state (readgrass lines 3604-3606)
      tdwi   = cfg%tdwi
      laiem  = cfg%laiem
      rgrlai = cfg%rgrlai

      ! Part 4: start-of-growth trigger (readgrass lines 3609-3614)
      swtsum = cfg%swtsum
      ! swtsum=2 is stub-guarded above; swtsum=0,1 need no extra fields.

      ! Part 5: green area (readgrass lines 3617-3620)
      if (allocated(cfg%slatb)) call copy_table(cfg%slatb, slatb)
      ssa   = cfg%ssa
      span  = cfg%span
      tbase = cfg%tbase

      ! Part 6: assimilation (readgrass lines 3623-3628)
      kdif = cfg%kdif
      kdir = cfg%kdir
      eff  = cfg%eff
      if (allocated(cfg%amaxtb))  call copy_table(cfg%amaxtb,  amaxtb)
      if (allocated(cfg%tmpftb))  call copy_table(cfg%tmpftb,  tmpftb)
      if (allocated(cfg%tmnftb))  call copy_table(cfg%tmnftb,  tmnftb)

      ! Part 7: conversion of assimilates (readgrass lines 3631-3633)
      cvl = cfg%cvl
      cvr = cfg%cvr
      cvs = cfg%cvs

      ! Part 8: maintenance respiration (readgrass lines 3636-3640)
      q10 = cfg%q10
      rml = cfg%rml
      rmr = cfg%rmr
      rms = cfg%rms
      if (allocated(cfg%rfsetb))  call copy_table(cfg%rfsetb,  rfsetb)

      ! Part 9: partitioning (readgrass lines 3643-3645)
      if (allocated(cfg%frtb))    call copy_table(cfg%frtb,    frtb)
      if (allocated(cfg%fltb))    call copy_table(cfg%fltb,    fltb)
      if (allocated(cfg%fstb))    call copy_table(cfg%fstb,    fstb)

      ! Part 10: death rates (readgrass lines 3648-3650)
      perdl = cfg%perdl
      if (allocated(cfg%rdrrtb))  call copy_table(cfg%rdrrtb,  rdrrtb)
      if (allocated(cfg%rdrstb))  call copy_table(cfg%rdrstb,  rdrstb)

      ! Part 11: oxygen stress (readgrass lines 3653-3723)
      ! Legacy default: swoxygen = 1 (readgrass line 3653)
      swoxygen = cfg%swoxygen
      if (cfg%swoxygen == 1) then
         hlim1  = cfg%hlim1
         hlim2u = cfg%hlim2u
         hlim2l = cfg%hlim2l
      else if (cfg%swoxygen == 2) then
         ! swoxygentype=1 physical path (swoxygentype=2 stub-guarded above)
         swoxygentype        = cfg%swoxygentype
         q10_microbial       = cfg%q10_microbial
         specific_resp_humus = cfg%specific_resp_humus
         srl                 = cfg%srl
         swrootradius        = cfg%swrootradius
         if (cfg%swrootradius == 1) then
            dry_mat_cont_roots     = cfg%dry_mat_cont_roots
            air_filled_root_por    = cfg%air_filled_root_por
            spec_weight_root_tissue = cfg%spec_weight_root_tissue
            var_a                  = cfg%var_a
         else
            root_radiusO2 = cfg%root_radiusO2
         end if
      end if
      ! Growth of roots during oxygen stress (readgrass lines 3715-3723)
      swWrtNonox = cfg%swwrtnonox
      aeratecrit = cfg%aeratecrit

      ! Part 12: drought stress (readgrass lines 3726-3754)
      ! swdrought=2 is stub-guarded above.
      swdrought = cfg%swdrought
      if (cfg%swdrought == 1) then
         hlim3h = cfg%hlim3h
         hlim3l = cfg%hlim3l
         hlim4  = cfg%hlim4
         adcrh  = cfg%adcrh
         adcrl  = cfg%adcrl
      end if

      ! Part 13: salt stress (readgrass lines 3757-3778)
      ! swsalinity /= 0 is stub-guarded above; always 0 on TOML path.
      swsalinity = cfg%swsalinity

      ! Part 14: compensation (readgrass lines 3781-3830)
      swcompensate = cfg%swcompensate
      if (cfg%swcompensate > 0) then
         ! swstressor defaults to 1 per legacy; only set when enabled.
         swstressor = cfg%swstressor
      end if
      if (cfg%swcompensate == 1) then
         alphacrit = cfg%alphacrit
      end if
      ! swcompensate=2 (Walsum dcritrtz) is stub-guarded above.

      ! Part 15: rooting (readgrass lines 3835-3873)
      swrdc = cfg%swrdc
      if (allocated(cfg%rdctb)) call copy_table(cfg%rdctb, rdctb)
      swrd = cfg%swrd
      ! swrd=1 and swrd=3 are stub-guarded above; only swrd=2 active.
      if (cfg%swrd == 2) then
         rdi      = cfg%rdi
         rri      = cfg%rri
         rdc      = cfg%rdc
         swdmi2rd = cfg%swdmi2rd
      end if

      ! Part 16: management factors (readgrass lines 3876-3885)
      relmf      = cfg%relmf
      swpotrelmf = cfg%swpotrelmf

      ! Part 17: sequence of mowing / grazing (readgrass lines 3891-3905)
      ! SeqGrazMow is a fixed-size integer array in variables (size 366).
      ! seqgrazmow(i) /= 2 is stub-guarded above.
      if (allocated(cfg%seqgrazmow)) then
         do i = 1, cfg%nseqgrazmow
            seqgrazmow(i) = cfg%seqgrazmow(i)
         end do
      end if

      ! Part 18: mowing settings (readgrass lines 3985-4035)
      ! Only the mowing block is active (SeqGrazMow all-2; grazing block
      ! is entirely guarded). mowrest always set.
      mowrest = cfg%mowrest

      ! swharvest (mowing trigger) is a LOCAL in grass(); handled in Task 8.
      ! dmharvest, daylastharvest, dmlastharvest, swdmmow, maxdaymow
      ! are also grass() locals; handled in Task 8.
      ! swlossmow, swlossgrz are locals and stub-guarded above (= 0 only).

      ! swharvest=1: mowing by DM threshold.
      ! swharvest=2: mowing by fixed dates → populate dateharvest global.
      if (cfg%swharv == 1) then
         ! swdmmow=1: fixed DM threshold (dmharvest, daylastharvest,
         !   dmlastharvest are grass() locals; set in Task 8).
         ! swdmmow=2: flexible DM threshold → dmmowtb global.
         if (cfg%swdmmow == 2) then
            if (allocated(cfg%dmmowtb)) call copy_table(cfg%dmmowtb, dmmowtb)
         end if
      else if (cfg%swharv == 2) then
         ! Populate dateharvest from cfg%mowing_dates (DOY floats).
         ! Legacy rdatim reads calendar dates like '1980-05-06' and
         ! converts them to t1900-relative real(8) timestamps.
         ! TOML stores them as DOY floats spanning the full simulation.
         ! Conversion: walk the dates; when DOY decreases, advance year.
         !
         ! This must be called every crop rotation (not just icrop==1)
         ! because InitializeCrop zeroes the dateharvest array at the start
         ! of each crop period.  populate_dateharvest anchors the DOY→t1900
         ! mapping to tstart (first simulation year), not yearmeteo, so the
         ! same full date sequence is reproduced correctly every time.
         call populate_dateharvest(cfg, tend)
      end if

      ! Regrowth delay table (readgrass lines 4028-4035).
      ! DelayRegrowthTab is a global (variables.f90 line 501).
      ! daydelay is a readgrass local; pack via the interleaved pattern.
      ! Config has dmmowdelay as flat pairs (already interleaved).
      if (allocated(cfg%dmmowdelay)) call copy_table(cfg%dmmowdelay, DelayRegrowthTab)

      ! Part 19: irrigation scheduling (readgrass line 4040)
      ! schedule=1 is stub-guarded above; always 0.
      schedule = cfg%schedule%schedule

      ! Part 20: CO2 correction (readgrass lines 4050-4083)
      ! swco2=1 is stub-guarded above; flCO2 is always .false.
      flCO2 = .false.

      ! ================================================================
      ! Runtime init math — cumdens (readgrass lines 4091-4121)
      ! Only when swdrought=1 (Feddes). Verbatim port of legacy.
      ! swdrought=2 is stub-guarded above.
      ! ================================================================
      if (cfg%swdrought == 1) then
         ! Build rootdis: 101 points at 0.01-spaced depths
         do i = 0, 100
            depth          = 0.01d0 * dble(i)
            rootdis(i*2+1) = depth
            rootdis(i*2+2) = afgen(rdctb, 22, depth)
         end do

         ! Copy depths to odd cumdens indices
         do i = 1, 202, 2
            cumdens(i) = rootdis(i)
         end do

         ! Trapezoidal cumulative integration
         sum_val    = 0.0d0
         cumdens(2) = 0.0d0
         do i = 4, 202, 2
            sum_val    = sum_val + (rootdis(i-2) + rootdis(i)) * 0.5d0 &
                                 * (cumdens(i-1) - cumdens(i-3))
            cumdens(i) = sum_val
         end do

         ! Normalize to 1
         if (sum_val > 0.0d0) then
            do i = 2, 202, 2
               cumdens(i) = cumdens(i) / sum_val
            end do
         end if
      end if

   end subroutine cropgrass_init_from_config

   ! ------------------------------------------------------------------
   ! Populate the dateharvest global from cfg%mowing_dates (DOY floats).
   !
   ! Legacy readgrass uses rdatim which reads calendar-date strings like
   ! '1980-05-06' and converts them to t1900-relative timestamps (days
   ! since 1-Jan-1900, Julian Day basis).  The TOML schema stores these
   ! as DOY floats spanning the full simulation period (e.g., 1980: DOY
   ! 127, 149, ...; 1981: DOY 104, ...).
   !
   ! Conversion:
   !   t1900_of_jan1(Y) = jday(Y, 1, 1) - 2415020    [jd1900 = 2415020]
   !   dateharvest(i)   = t1900_of_jan1(year_i) + mowing_dates(i) - 1.0
   !
   ! The year advances when mowing_dates(i) < mowing_dates(i-1) (year
   ! roll-over).  The starting year is yearmeteo (set by timecontrol.f90
   ! to the year of the current t1900, which at crop init time equals the
   ! first simulation year).
   !
   ! Sentinel: dateharvest(nmow+1) = tend + 1.0  (legacy line 4023).
   ! ------------------------------------------------------------------
   subroutine populate_dateharvest(cfg, tend_val)
      use variables, only: dateharvest, tstart
      implicit none
      type(cropgrass_config_t), intent(in) :: cfg
      real(real64),             intent(in) :: tend_val

      integer      :: i, cur_year, start_year
      real(real64) :: t_jan1

      if (.not. allocated(cfg%mowing_dates)) return
      if (cfg%nmow <= 0) return

      ! Derive the simulation start year from tstart.
      ! This is stable across all crop rotations: mowing_dates span the
      ! entire simulation, so we always anchor the DOY→t1900 mapping to
      ! the first simulation year regardless of which rotation icrop we are.
      start_year = year_from_t1900(int(tstart))
      cur_year   = start_year
      t_jan1     = real(t1900_from_year(cur_year), real64)

      do i = 1, cfg%nmow
         ! Year roll-over: DOY decreases means we crossed Jan 1.
         if (i > 1) then
            if (cfg%mowing_dates(i) < cfg%mowing_dates(i-1)) then
               cur_year = cur_year + 1
               t_jan1   = real(t1900_from_year(cur_year), real64)
            end if
         end if
         dateharvest(i) = t_jan1 + cfg%mowing_dates(i) - 1.0d0
      end do

      ! Sentinel past the last mowing date (legacy readgrass line 4023).
      dateharvest(cfg%nmow + 1) = tend_val + 1.0d0

   end subroutine populate_dateharvest

   ! ------------------------------------------------------------------
   ! Compute the t1900 value for January 1 of year Y.
   ! t1900 = Julian Day Number of Jan-1-Y  minus  jd1900 (= 2415020).
   ! jday formula (Gregorian): same as src/io/readmeteo.f90:765-772.
   ! ------------------------------------------------------------------
   pure function t1900_from_year(y) result(t)
      integer, intent(in) :: y
      integer :: t
      integer :: a, yy, mm, jd
      ! jday(Y, 1, 1):  a=(14-1)/12=1, yy=Y+4800-1=Y+4799, mm=1+12-3=10
      a  = 1
      yy = y + 4799
      mm = 10
      jd = 1 + (153*mm + 2)/5 + 365*yy + yy/4 - yy/100 + yy/400 - 32045
      t  = jd - 2415020
   end function t1900_from_year

   ! ------------------------------------------------------------------
   ! Compute the calendar year containing a given t1900 value.
   ! Uses a bisection / forward-search approach on t1900_from_year.
   ! ------------------------------------------------------------------
   pure function year_from_t1900(t) result(y)
      integer, intent(in) :: t
      integer :: y, lo, hi, mid
      ! Rough estimate: t1900 for Jan-1-1900 ≈ 0; ~365.25 d/yr.
      y  = 1900 + max(0, (t - 1) / 366)
      ! Refine: find year such that t1900_from_year(y) <= t < t1900_from_year(y+1)
      lo = y - 2
      hi = y + 2
      do while (t1900_from_year(hi) <= t)
         lo = hi
         hi = hi + 10
      end do
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (t1900_from_year(mid) <= t) then
            lo = mid
         else
            hi = mid
         end if
      end do
      y = lo
   end function year_from_t1900

   ! ------------------------------------------------------------------
   ! Copy a flat allocatable array into a fixed-size module global.
   ! Copies up to min(size(src), size(dst)) elements.
   ! ------------------------------------------------------------------
   subroutine copy_table(src, dst)
      real(real64), intent(in)    :: src(:)
      real(real64), intent(inout) :: dst(:)
      integer :: n
      n = min(size(src), size(dst))
      dst(1:n) = src(1:n)
   end subroutine copy_table

end module cropgrass_init_mod
