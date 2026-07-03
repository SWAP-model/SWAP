!> @file swap_var_registry_mod.f90
!! Single source of truth for the C-ABI variable surface (arc T1-G′ sub-arc 1).
!!
!! A `var_registry_t` maps a variable name to metadata + a pointer into the live
!! `swap_state_t`, so the BMI / CAPI accessors resolve variables by table lookup
!! instead of hand-maintained `select case` switches. Each entry carries a
!! namespace tag (BMI vs CAPI) so that enumerate-style calls (item counts,
!! name lists) return exactly one facade's vocabulary, while name resolution
!! stays namespace-scoped and therefore byte-for-byte compatible with the
!! switches it replaces.
!!
!! Pointers are bound AFTER `swap_init` (arrays are allocated by then) via
!! `build_variable_registry`. They alias components of a state instance that must
!! have the TARGET attribute (the CAPI singleton `capi_state` does), so the
!! associations remain valid after the builder returns.
!!
!! Scope note: like `capi_state`, a registry instance is single-kernel-per-process
!! (it points into one state). Per-instance registries are a Phase-3 / sub-arc-2
!! concern (handle-based multi-instance).
module swap_var_registry_mod
   use iso_fortran_env, only: real64
   use iso_c_binding,   only: c_ptr, c_loc
   use swap_state_mod,  only: swap_state_t
   implicit none
   private

   public :: var_registry_t, build_variable_registry
   public :: NS_BMI, NS_CAPI

   ! Namespace (facade) tags. A variable belongs to exactly one namespace; the
   ! same underlying field may be registered under several namespaces/names
   ! (e.g. `theta` as BMI `soil_water_content` and CAPI `theta`).
   integer, parameter :: NS_BMI  = 1
   integer, parameter :: NS_CAPI = 2

   integer, parameter :: NAME_LEN = 64, UNIT_LEN = 16, MAX_VARS = 32

   type :: var_entry_t
      character(len=NAME_LEN) :: name  = ''
      character(len=UNIT_LEN) :: units = ''
      integer               :: ns       = NS_BMI
      logical               :: is_array = .false.
      logical               :: readable = .true.
      logical               :: settable = .false.
      real(real64)          :: scale    = 1.0_real64   !< read returns scale*value (recharge = -qbot)
      real(real64), pointer :: v1(:)    => null()       !< bound field, rank-1
      real(real64), pointer :: v0       => null()       !< bound field, scalar
   end type var_entry_t

   type :: var_registry_t
      type(var_entry_t) :: e(MAX_VARS)
      integer           :: n = 0
   contains
      procedure :: clear      => reg_clear
      procedure :: add_array  => reg_add_array
      procedure :: add_scalar => reg_add_scalar
      procedure :: find       => reg_find
      procedure :: is_readable => reg_is_readable
      procedure :: is_settable => reg_is_settable
      procedure :: read_into  => reg_read_into
      procedure :: write_from => reg_write_from
      procedure :: c_view     => reg_c_view
      procedure :: nbytes     => reg_nbytes
      procedure :: units      => reg_units
      procedure :: pack_names => reg_pack_names
      procedure :: count_ns   => reg_count_ns
   end type var_registry_t

contains

   subroutine reg_clear(self)
      class(var_registry_t), intent(inout) :: self
      self%n = 0
   end subroutine reg_clear

   subroutine reg_add_array(self, ns, name, units, arr, readable)
      class(var_registry_t), intent(inout)      :: self
      integer,               intent(in)         :: ns
      character(*),          intent(in)         :: name, units
      real(real64),          intent(in), target :: arr(:)
      logical,               intent(in), optional :: readable
      integer :: k
      self%n = self%n + 1
      k = self%n
      self%e(k)%name     = name
      self%e(k)%units    = units
      self%e(k)%ns       = ns
      self%e(k)%is_array = .true.
      self%e(k)%readable = .true.
      if (present(readable)) self%e(k)%readable = readable
      self%e(k)%v1      => arr
   end subroutine reg_add_array

   subroutine reg_add_scalar(self, ns, name, units, field, readable, settable, scale)
      class(var_registry_t), intent(inout)      :: self
      integer,               intent(in)         :: ns
      character(*),          intent(in)         :: name, units
      real(real64),          intent(in), target :: field
      logical,               intent(in), optional :: readable, settable
      real(real64),          intent(in), optional :: scale
      integer :: k
      self%n = self%n + 1
      k = self%n
      self%e(k)%name     = name
      self%e(k)%units    = units
      self%e(k)%ns       = ns
      self%e(k)%is_array = .false.
      self%e(k)%readable = .true.
      self%e(k)%settable = .false.
      self%e(k)%scale    = 1.0_real64
      if (present(readable)) self%e(k)%readable = readable
      if (present(settable)) self%e(k)%settable = settable
      if (present(scale))    self%e(k)%scale    = scale
      self%e(k)%v0      => field
   end subroutine reg_add_scalar

   !> Index of the entry named `name` in namespace `ns`, or 0 if none.
   integer function reg_find(self, ns, name) result(idx)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: ns
      character(*),          intent(in) :: name
      integer :: i
      idx = 0
      do i = 1, self%n
         if (self%e(i)%ns == ns .and. trim(self%e(i)%name) == trim(name)) then
            idx = i
            return
         end if
      end do
   end function reg_find

   logical function reg_is_readable(self, idx) result(ok)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: idx
      ok = self%e(idx)%readable
   end function reg_is_readable

   logical function reg_is_settable(self, idx) result(ok)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: idx
      ok = self%e(idx)%settable
   end function reg_is_settable

   !> Copy entry `idx` into `dest(1:n)` with BMI get_value semantics: arrays fill
   !! the first min(n,size) slots and zero the remainder; scalars fill dest(1) and
   !! zero the rest; `scale` is applied on read.
   subroutine reg_read_into(self, idx, dest, n)
      class(var_registry_t), intent(in)  :: self
      integer,               intent(in)  :: idx, n
      real(real64),          intent(out) :: dest(n)
      integer :: m
      if (self%e(idx)%is_array) then
         m = min(n, size(self%e(idx)%v1))
         dest(1:m) = self%e(idx)%scale * self%e(idx)%v1(1:m)
         if (m < n) dest(m+1:n) = 0.0_real64
      else
         if (n >= 1) dest(1) = self%e(idx)%scale * self%e(idx)%v0
         if (n > 1)  dest(2:n) = 0.0_real64
      end if
   end subroutine reg_read_into

   !> Write src(1) into a settable scalar entry (the BMI imposed-BC path).
   subroutine reg_write_from(self, idx, src, n)
      class(var_registry_t), intent(inout) :: self
      integer,               intent(in)    :: idx, n
      real(real64),          intent(in)    :: src(n)
      if (n >= 1) self%e(idx)%v0 = src(1)
   end subroutine reg_write_from

   !> Zero-copy C pointer + element count for an array entry (CAPI view path).
   subroutine reg_c_view(self, idx, ptr, nsize)
      class(var_registry_t), intent(in)  :: self
      integer,               intent(in)  :: idx
      type(c_ptr),           intent(out) :: ptr
      integer,               intent(out) :: nsize
      ptr   = c_loc(self%e(idx)%v1(1))
      nsize = size(self%e(idx)%v1)
   end subroutine reg_c_view

   integer function reg_nbytes(self, idx) result(nb)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: idx
      if (self%e(idx)%is_array) then
         nb = size(self%e(idx)%v1) * 8
      else
         nb = 8
      end if
   end function reg_nbytes

   function reg_units(self, idx) result(u)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: idx
      character(len=UNIT_LEN) :: u
      u = self%e(idx)%units
   end function reg_units

   !> Count entries in namespace `ns` that are readable (want_settable=.false.)
   !! or settable (want_settable=.true.).
   integer function reg_count_ns(self, ns, want_settable) result(c)
      class(var_registry_t), intent(in) :: self
      integer,               intent(in) :: ns
      logical,               intent(in) :: want_settable
      integer :: i
      c = 0
      do i = 1, self%n
         if (self%e(i)%ns /= ns) cycle
         if (want_settable) then
            if (self%e(i)%settable) c = c + 1
         else
            if (self%e(i)%readable) c = c + 1
         end if
      end do
   end function reg_count_ns

   !> Pack the NUL-delimited name list for namespace `ns` (readable, or settable
   !! when want_settable) into buf_names in registration order. `nbytes` returns
   !! the number of bytes written (names + their NUL terminators), so the caller
   !! copies exactly that many into the C buffer.
   subroutine reg_pack_names(self, ns, want_settable, buf_names, nbytes)
      class(var_registry_t), intent(in)  :: self
      integer,               intent(in)  :: ns
      logical,               intent(in)  :: want_settable
      character(len=*),      intent(out) :: buf_names
      integer,               intent(out) :: nbytes
      integer :: i, k, j
      character(len=1), parameter :: nul = achar(0)
      buf_names = ''
      k = 0
      do i = 1, self%n
         if (self%e(i)%ns /= ns) cycle
         if (want_settable) then
            if (.not. self%e(i)%settable) cycle
         else
            if (.not. self%e(i)%readable) cycle
         end if
         do j = 1, len_trim(self%e(i)%name)
            if (k + 1 > len(buf_names)) then; nbytes = k; return; end if
            k = k + 1
            buf_names(k:k) = self%e(i)%name(j:j)
         end do
         if (k + 1 > len(buf_names)) then; nbytes = k; return; end if
         k = k + 1
         buf_names(k:k) = nul
      end do
      nbytes = k
   end subroutine reg_pack_names

   !> Register the C-ABI variable surface over a (post-init) state instance.
   !! `st` must be a TARGET (the CAPI singleton is) so the bound pointers survive.
   subroutine build_variable_registry(reg, st)
      type(var_registry_t), intent(inout)      :: reg
      type(swap_state_t),   intent(inout), target :: st
      integer :: numnod
      call reg%clear()
      numnod = st%mesh%numnod

      ! --- BMI namespace: 8 readable outputs (profile arrays + scalars) --------
      call reg%add_array (NS_BMI, 'soil_water_content',        'm3 m-3', st%soilwater%theta)
      call reg%add_array (NS_BMI, 'pressure_head',             'cm',     st%soilwater%h)
      call reg%add_array (NS_BMI, 'soil_temperature',          'degC',   st%heat%tsoil)
      call reg%add_scalar(NS_BMI, 'groundwater_level',         'cm',     st%soilwater%gwl)
      call reg%add_scalar(NS_BMI, 'bottom_flux',               'cm d-1', st%soilwater%qbot)
      call reg%add_scalar(NS_BMI, 'actual_evapotranspiration', 'cm d-1', st%soilwater%iqrot)
      call reg%add_scalar(NS_BMI, 'recharge',                  'cm d-1', st%soilwater%qbot, scale=-1.0_real64)
      call reg%add_scalar(NS_BMI, 'surface_runoff',            'cm d-1', st%soilwater%runots)

      ! --- BMI namespace: 2 settable boundary conditions -----------------------
      ! groundwater_level_imposed overwrites the bottom-node pressure head.
      call reg%add_scalar(NS_BMI, 'groundwater_level_imposed', 'cm',     st%soilwater%h(numnod), &
                          readable=.false., settable=.true.)
      call reg%add_scalar(NS_BMI, 'bottom_flux_imposed',       'cm d-1', st%soilwater%qbot, &
                          readable=.false., settable=.true.)

      ! --- CAPI namespace: raw-named zero-copy array views ---------------------
      call reg%add_array (NS_CAPI, 'theta',  'm3 m-3', st%soilwater%theta)
      call reg%add_array (NS_CAPI, 'h',      'cm',     st%soilwater%h)
      call reg%add_array (NS_CAPI, 'tsoil',  'degC',   st%heat%tsoil)
      call reg%add_array (NS_CAPI, 'inqrot', 'cm',     st%soilwater%inqrot)
   end subroutine build_variable_registry

end module swap_var_registry_mod
