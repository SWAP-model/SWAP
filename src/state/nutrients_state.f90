!> @file nutrients_state.f90
!! SS-GR-CROP: typed nutrients runtime state. Hosts management_soil +
!! wofost_soil_* nutrient pool and flow data.
!!
!! Audit basis: wofost_soil_declarations.f90 + wofost_soil_interface.f90.
!! The WSN (Water-Soil-Nitrogen) model operates on a single representative
!! soil layer (dz_WSN), so all fields are scalars or small fixed-size arrays
!! indexed by FOM fraction (maxfn=8), NOT per-soil-layer arrays.
!! nlay is accepted by init for future extensibility but not currently used.
module nutrients_state_mod
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private
   public :: nutrients_state_t
   public :: seed_nutrients_from_config
   public :: load_nutrients_events

   !> Maximum number of FOM (Fresh Organic Matter) fractions.
   !! Mirrors Wofost_Soil_Declarations maxfn = 8.
   integer, parameter :: maxfn_ns = 8

   type :: nutrients_state_t

      ! -----------------------------------------------------------------------
      ! Organic matter pools — end-of-timestep values (kg/m3 in WSN layer)
      ! -----------------------------------------------------------------------
      real(real64) :: fom_t(maxfn_ns) = 0.0_real64   !! fresh OM fractions at t
      real(real64) :: bio_t           = 0.0_real64   !! microbial biomass at t
      real(real64) :: hum_t           = 0.0_real64   !! humus at t

      ! -----------------------------------------------------------------------
      ! Organic matter pools — start-of-timestep values (kg/m3)
      ! -----------------------------------------------------------------------
      real(real64) :: fom_t0(maxfn_ns) = 0.0_real64  !! FOM fractions at t0
      real(real64) :: bio_t0            = 0.0_real64  !! biomass at t0
      real(real64) :: hum_t0            = 0.0_real64  !! humus at t0

      ! -----------------------------------------------------------------------
      ! Mineral N concentrations (kg/m3 in soil solution)
      ! -----------------------------------------------------------------------
      real(real64) :: cnh4_t  = 0.0_real64  !! ammonium conc at end of timestep
      real(real64) :: cnh4_t0 = 0.0_real64  !! ammonium conc at start of timestep
      real(real64) :: cnh4_av = 0.0_real64  !! time-averaged ammonium conc

      real(real64) :: cno3_t  = 0.0_real64  !! nitrate conc at end of timestep
      real(real64) :: cno3_t0 = 0.0_real64  !! nitrate conc at start of timestep
      real(real64) :: cno3_av = 0.0_real64  !! time-averaged nitrate conc

      ! -----------------------------------------------------------------------
      ! Per-timestep process rates / scratchpads
      ! -----------------------------------------------------------------------
      real(real64) :: nminer      = 0.0_real64  !! net N mineralisation (kg/m3)
      real(real64) :: cdissi      = 0.0_real64  !! C dissimilation (kg/m2)
      real(real64) :: nsupplynh4n = 0.0_real64  !! NH4-N supply rate (kg/m3)
      real(real64) :: nsupplyno3n = 0.0_real64  !! NO3-N supply rate (kg/m3)

      ! -----------------------------------------------------------------------
      ! Crop–soil N coupling (from Wofost_Soil_Interface)
      ! -----------------------------------------------------------------------
      real(real64) :: ndemandsoil  = 0.0_real64  !! total crop N demand (kg/ha)
      real(real64) :: nsupplysoil  = 0.0_real64  !! total mineral N supply (kg/ha)
      real(real64) :: ndemand      = 0.0_real64  !! crop N demand (kg/m2)
      real(real64) :: nsupply      = 0.0_real64  !! mineral N supply (kg/m2)
      real(real64) :: laicritnupt  = 0.0_real64  !! critical LAI for N uptake

      ! -----------------------------------------------------------------------
      ! Balance-check state variables (reset after each output period)
      ! — Organic matter balances (kg/m2)
      ! -----------------------------------------------------------------------
      real(real64) :: fom_old  = 0.0_real64
      real(real64) :: fom_end  = 0.0_real64
      real(real64) :: fom_add  = 0.0_real64   !! OM added by amendments
      real(real64) :: fom_cres = 0.0_real64   !! OM added by crop residues
      real(real64) :: fom2bio  = 0.0_real64
      real(real64) :: fom2hum  = 0.0_real64
      real(real64) :: fom_dis  = 0.0_real64

      real(real64) :: bio_old  = 0.0_real64
      real(real64) :: bio_end  = 0.0_real64
      real(real64) :: bio2bio  = 0.0_real64
      real(real64) :: bio2hum  = 0.0_real64
      real(real64) :: bio_dis  = 0.0_real64

      real(real64) :: hum_old  = 0.0_real64
      real(real64) :: hum_end  = 0.0_real64
      real(real64) :: hum_add  = 0.0_real64
      real(real64) :: hum_cres = 0.0_real64
      real(real64) :: hum2bio  = 0.0_real64
      real(real64) :: hum2hum  = 0.0_real64
      real(real64) :: hum_dis  = 0.0_real64

      ! — Organic N balances (kg/m2)
      real(real64) :: nfom_old  = 0.0_real64
      real(real64) :: nfom_end  = 0.0_real64
      real(real64) :: nfom_add  = 0.0_real64
      real(real64) :: nfom_cres = 0.0_real64
      real(real64) :: nfom2bio  = 0.0_real64
      real(real64) :: nfom2hum  = 0.0_real64
      real(real64) :: nfom_min  = 0.0_real64

      real(real64) :: nbio_old  = 0.0_real64
      real(real64) :: nbio_end  = 0.0_real64
      real(real64) :: nbio2bio  = 0.0_real64
      real(real64) :: nbio2hum  = 0.0_real64
      real(real64) :: nbio_min  = 0.0_real64

      real(real64) :: nhum_old  = 0.0_real64
      real(real64) :: nhum_end  = 0.0_real64
      real(real64) :: nhum_add  = 0.0_real64
      real(real64) :: nhum_cres = 0.0_real64
      real(real64) :: nhum2bio  = 0.0_real64
      real(real64) :: nhum2hum  = 0.0_real64
      real(real64) :: nhum_min  = 0.0_real64

      ! — NH4 balance terms (kg/m2)
      real(real64) :: nh4_old    = 0.0_real64
      real(real64) :: nh4_end    = 0.0_real64
      real(real64) :: nh4_miner  = 0.0_real64
      real(real64) :: nh4_intop  = 0.0_real64
      real(real64) :: nh4_inlat  = 0.0_real64
      real(real64) :: nh4_inbot  = 0.0_real64
      real(real64) :: nh4_upt    = 0.0_real64
      real(real64) :: nh4_out    = 0.0_real64
      real(real64) :: nh4_nitrif = 0.0_real64

      ! — NO3 balance terms (kg/m2)
      real(real64) :: no3_old    = 0.0_real64
      real(real64) :: no3_end    = 0.0_real64
      real(real64) :: no3_intop  = 0.0_real64
      real(real64) :: no3_inlat  = 0.0_real64
      real(real64) :: no3_inbot  = 0.0_real64
      real(real64) :: no3_upt    = 0.0_real64
      real(real64) :: no3_out    = 0.0_real64
      real(real64) :: no3_denitr = 0.0_real64

      ! — Amendment/residue N tracking (per output period, kg/m2)
      real(real64) :: nh4n_amend = 0.0_real64
      real(real64) :: no3n_amend = 0.0_real64
      real(real64) :: nh4n_cres  = 0.0_real64
      real(real64) :: no3n_cres  = 0.0_real64
      real(real64) :: nh4n_volat = 0.0_real64

      ! — Previous-period crop residue DM and N losses (for ANIMO cropext output)
      real(real64) :: idwrt_1   = 0.0_real64  !! root DM at previous harvest
      real(real64) :: idwlv_1   = 0.0_real64  !! leaf DM at previous harvest
      real(real64) :: idwst_1   = 0.0_real64  !! stem DM at previous harvest
      real(real64) :: idwso_1   = 0.0_real64  !! storage organ DM at previous harvest
      real(real64) :: inlossl_1 = 0.0_real64  !! leaf N loss at previous harvest
      real(real64) :: inlossr_1 = 0.0_real64  !! root N loss at previous harvest
      real(real64) :: inlosss_1 = 0.0_real64  !! stem N loss at previous harvest
      real(real64) :: inlosso_1 = 0.0_real64  !! storage organ N loss at previous harvest

      ! -----------------------------------------------------------------------
      ! ANIMO cropext file accumulators (cumulative per output period)
      ! -----------------------------------------------------------------------
      real(real64) :: ntotuptake = 0.0_real64  !! total N uptake (kg/ha)
      real(real64) :: ptotuptake = 0.0_real64  !! total P uptake (kg/ha, P/N ratio applied)
      real(real64) :: dmcressur  = 0.0_real64  !! aboveground crop residue DM (kg/ha)
      real(real64) :: ncressurf  = 0.0_real64  !! aboveground crop residue N (kg/ha)
      real(real64) :: pcressurf  = 0.0_real64  !! aboveground crop residue P (kg/ha)
      real(real64) :: dmcresbott = 0.0_real64  !! belowground crop residue DM (kg/ha)
      real(real64) :: ncresbott  = 0.0_real64  !! belowground crop residue N (kg/ha)
      real(real64) :: pcresbott  = 0.0_real64  !! belowground crop residue P (kg/ha)

   contains
      procedure :: init => nutrients_state_init
   end type nutrients_state_t

contains

   subroutine nutrients_state_init(self, nlay, config_nut, pathwork_in)
      use nutrients_config_mod, only: nutrients_config_t
      class(nutrients_state_t), intent(inout) :: self
      integer,                  intent(in)    :: nlay
      type(nutrients_config_t), intent(in)    :: config_nut
      character(len=*),         intent(in)    :: pathwork_in
      ! All fields carry default initializers (= 0.0_real64 or = 0).
      ! nlay accepted for API consistency; WSN operates on a single
      ! representative layer so no per-layer allocations are needed here.
      ! Explicit zero to be safe at re-init if ever called again.
      self%fom_t   = 0.0_real64
      self%bio_t   = 0.0_real64
      self%hum_t   = 0.0_real64
      self%fom_t0  = 0.0_real64
      self%bio_t0  = 0.0_real64
      self%hum_t0  = 0.0_real64
      self%cnh4_t  = 0.0_real64 ;  self%cnh4_t0 = 0.0_real64 ;  self%cnh4_av = 0.0_real64
      self%cno3_t  = 0.0_real64 ;  self%cno3_t0 = 0.0_real64 ;  self%cno3_av = 0.0_real64
      self%nminer      = 0.0_real64
      self%cdissi      = 0.0_real64
      self%nsupplynh4n = 0.0_real64
      self%nsupplyno3n = 0.0_real64
      self%ndemandsoil = 0.0_real64 ;  self%nsupplysoil = 0.0_real64
      self%ndemand     = 0.0_real64 ;  self%nsupply     = 0.0_real64
      self%laicritnupt = 0.0_real64

      ! Seed legacy WSN compute globals (Wofost_Soil_Declarations) from
      ! typed config. The body lives below; relocated from the deleted
      ! config_to_variables_mod (apply_nutrients).
      call seed_nutrients_from_config(config_nut, pathwork_in)

      ! Mirror the same five config-derived initial pools onto the typed
      ! state fields. Replaces the WSN dual-write block previously in
      ! swap_init_body, which copied the just-seeded globals back onto
      ! state. The legacy WSN compute path still reads the globals; the
      ! state fields are the typed mirror for future readers and the
      ! management_soil.f90 runtime dual-writes.
      self%fom_t(:) = config_nut%initial%fom(:)
      self%bio_t    = config_nut%initial%bio
      self%hum_t    = config_nut%initial%hum
      self%cnh4_t   = config_nut%initial%cnh4
      self%cno3_t   = config_nut%initial%cno3
   end subroutine nutrients_state_init


   !> Seed Wofost_Soil_Declarations globals from typed nutrients config.
   !! Body relocated from src/io/toml/config_to_variables.f90 (apply_nutrients).
   !! The destination is a separate legacy module (WSN); migrating those
   !! globals is out of scope for this arc.
   subroutine seed_nutrients_from_config(cfg, pathwork_in)
      use nutrients_config_mod, only: nutrients_config_t
      use wofost_soil_declarations, only: FOM_t, Bio_t, Hum_t, &
                                           cNH4_t, cNO3_t, SorpCoef
      type(nutrients_config_t), intent(in) :: cfg
      character(len=*),         intent(in) :: pathwork_in
      integer :: i

      SorpCoef = cfg%sorp_coef
      do i = 1, 8
         FOM_t(i) = cfg%initial%fom(i)
      end do
      Bio_t  = cfg%initial%bio
      Hum_t  = cfg%initial%hum
      cNH4_t = cfg%initial%cnh4
      cNO3_t = cfg%initial%cno3

      ! N2b (ADR 0027): stage timed amendments from the CSV companion.
      call load_nutrients_events(cfg, pathwork_in)
   end subroutine seed_nutrients_from_config


   !> Stage timed soil management events from the CSV companion
   !! at cfg%events_file. Sets the legacy globals consumed by
   !! SoilManagement(3) and Wofost_SoilAmendents.
   !!
   !! Default (empty events_file or empty CSV): namend = 0, isme = 1.
   !!
   !! Body relocated from src/io/toml/config_to_variables.f90 (apply_nutrients_events).
   !! See ADR 0027 ([nutrients] N2b).
   subroutine load_nutrients_events(cfg, pathwork_in)
      use, intrinsic :: iso_fortran_env, only: real64
      use error_mod, only: error_collection_t, fatalerr_collected
      use nutrients_config_mod, only: nutrients_config_t
      use nutrients_csv_mod, only: amendment_events_table_t
      use wofost_soil_declarations, only: MatNum, Amend, VolaFrac, &
                                            TimeAmend, NuAmend, iamend, &
                                            namend, isme, maxamn
      type(nutrients_config_t), intent(in) :: cfg
      character(len=*),         intent(in) :: pathwork_in

      type(amendment_events_table_t) :: typed_tbl
      type(error_collection_t)       :: errs
      character(len=300) :: csvpath
      integer :: i, j, n
      real(real64) :: tmp_date, tmp_amount, tmp_volat
      integer :: tmp_mat

      ! Default: no amendments. Reset legacy globals to a known state.
      namend = 0
      isme   = 1

      ! events_file is always allocated by get_optional_string_with_default
      ! (defaults to empty string on missing key), so check len_trim instead.
      if (.not. allocated(cfg%events_file)) return
      if (len_trim(cfg%events_file) == 0)   return

      ! Typed loader: parse + validate (material 1..20, amount [0,5e5], volat [0,1]).
      ! Range/enum validation is owned by the loader; this site only handles
      ! the post-load conversion and global population.
      csvpath = trim(pathwork_in) // trim(cfg%events_file)
      call typed_tbl%load(trim(csvpath), errs)
      call errs%abort_if_fatal()

      n = 0
      if (typed_tbl%is_loaded) n = size(typed_tbl%rows)
      if (n < 1) return     ! Empty CSV: no amendments. Not an error.
      if (n > maxamn) then
         call fatalerr_collected('load_nutrients_events', &
            'CSV row count exceeds maxamn (1000)')
         return
      end if

      ! Defense-in-depth validation (parallel to the loader checks; keeps
      ! the existing fatalerr contract for callers that bypass the typed path).
      do i = 1, n
         if (typed_tbl%rows(i)%material < 1 .or. typed_tbl%rows(i)%material > 20) then
            call fatalerr_collected('load_nutrients_events', &
               'material out of range [1, 20]')
            return
         end if
         if (typed_tbl%rows(i)%amount_kgha < 0.0_real64 .or. &
             typed_tbl%rows(i)%amount_kgha > 500000.0_real64) then
            call fatalerr_collected('load_nutrients_events', &
               'amount_kgha out of range [0, 500000]')
            return
         end if
         if (typed_tbl%rows(i)%volat_fraction < 0.0_real64 .or. &
             typed_tbl%rows(i)%volat_fraction > 1.0_real64) then
            call fatalerr_collected('load_nutrients_events', &
               'volat_fraction out of range [0, 1]')
            return
         end if
      end do

      ! Sort by date (in-place bubble sort, mirrors deleted SoilManagement(1)).
      ! Acceptable O(n^2) given n <= 1000 and this runs once at config-load.
      do i = 1, n - 1
         do j = i + 1, n
            if (typed_tbl%rows(i)%date > typed_tbl%rows(j)%date) then
               tmp_date                    = typed_tbl%rows(i)%date
               typed_tbl%rows(i)%date      = typed_tbl%rows(j)%date
               typed_tbl%rows(j)%date      = tmp_date
               tmp_mat                     = typed_tbl%rows(i)%material
               typed_tbl%rows(i)%material  = typed_tbl%rows(j)%material
               typed_tbl%rows(j)%material  = tmp_mat
               tmp_amount                       = typed_tbl%rows(i)%amount_kgha
               typed_tbl%rows(i)%amount_kgha    = typed_tbl%rows(j)%amount_kgha
               typed_tbl%rows(j)%amount_kgha    = tmp_amount
               tmp_volat                         = typed_tbl%rows(i)%volat_fraction
               typed_tbl%rows(i)%volat_fraction  = typed_tbl%rows(j)%volat_fraction
               typed_tbl%rows(j)%volat_fraction  = tmp_volat
            end if
         end do
      end do

      ! Populate per-event legacy globals.
      ! Unit conversion: kg/ha -> kg/m^2 via ×1e-4 applied at the copy site.
      do i = 1, n
         MatNum(i)   = typed_tbl%rows(i)%material
         Amend(i)    = 1.0e-4_real64 * typed_tbl%rows(i)%amount_kgha   ! kg/ha -> kg/m^2
         VolaFrac(i) = typed_tbl%rows(i)%volat_fraction
      end do

      ! Group dosages per date (mirrors deleted SoilManagement(1)).
      j = 1
      NuAmend(j)   = 1
      TimeAmend(j) = typed_tbl%rows(1)%date
      iamend(1, 1) = 1
      do i = 2, n
         if (abs(typed_tbl%rows(i)%date - typed_tbl%rows(i - 1)%date) < 1.0e-3_real64) then
            NuAmend(j) = NuAmend(j) + 1
         else
            j = j + 1
            NuAmend(j)   = 1
            TimeAmend(j) = typed_tbl%rows(i)%date
         end if
         iamend(j, NuAmend(j)) = i
      end do

      namend = j
      isme   = 1
   end subroutine load_nutrients_events

end module nutrients_state_mod
