!> TOML emitter: writes a swap_config_t to a TOML file.
!!
!! Inverse of load_swap_config. Round-trip invariant: load -> emit -> load
!! produces an identical swap_config_t for all scalar fields (within
!! floating-point tolerance). Array-of-tables (drainage levels, crop
!! rotation) are Phase 4b; this emitter intentionally omits them.
module write_swap_config_mod
   use iso_fortran_env, only: real64
   use tomlf, only: toml_table, toml_error, toml_datetime, &
                    set_value, add_table, new_table, toml_dump
   use swap_config_mod, only: swap_config_t
   use general_config_mod,      only: general_config_t
   use simulation_config_mod,   only: simulation_config_t
   use meteorology_config_mod,  only: meteorology_config_t
   use drainage_config_mod,     only: drainage_config_t
   use soil_config_mod,         only: soil_config_t
   use crop_config_mod,         only: crop_config_t
   use error_mod, only: error_collection_t, ERR_IO_WRITE_FAILED
   implicit none
   private

   public :: write_swap_config

contains

   !> Write config to TOML file at path.
   subroutine write_swap_config(config, path, errors)
      type(swap_config_t),      intent(in)    :: config
      character(len=*),         intent(in)    :: path
      type(error_collection_t), intent(inout) :: errors

      type(toml_table)           :: doc
      type(toml_error), allocatable :: toml_err

      call new_table(doc)
      call emit_general    (doc, config%general)
      call emit_simulation (doc, config%simulation)
      call emit_meteorology(doc, config%meteo)
      call emit_drainage   (doc, config%drain)
      call emit_soil       (doc, config%soil)
      call emit_crop       (doc, config%crop)

      call toml_dump(doc, trim(path), toml_err)
      if (allocated(toml_err)) then
         call errors%append(ERR_IO_WRITE_FAILED, &
                            "could not write " // trim(path) // ": " // toml_err%message, &
                            trim(path))
      end if
   end subroutine write_swap_config

   subroutine emit_general(doc, c)
      type(toml_table),       intent(inout) :: doc
      type(general_config_t), intent(in)    :: c
      type(toml_table), pointer :: sec, paths

      call add_table(doc, 'general', sec)
      if (allocated(c%project))  call set_value(sec, 'project', c%project)
      call set_value(sec, 'swscre',  c%swscre)
      call set_value(sec, 'swerror', c%swerror)

      call add_table(sec, 'paths', paths)
      if (allocated(c%pathwork))  call set_value(paths, 'work',       c%pathwork)
      if (allocated(c%pathatm))   call set_value(paths, 'atmosphere', c%pathatm)
      if (allocated(c%pathcrop))  call set_value(paths, 'crop',       c%pathcrop)
      if (allocated(c%pathdrain)) call set_value(paths, 'drain',      c%pathdrain)
   end subroutine emit_general

   subroutine emit_simulation(doc, c)
      type(toml_table),          intent(inout) :: doc
      type(simulation_config_t), intent(in)    :: c
      type(toml_table), pointer :: sec, output
      type(toml_datetime) :: dtv

      call add_table(doc, 'simulation', sec)

      ! Emit tstart/tend as TOML local-date literals (matching the reader keys
      ! start_date / end_date). This ensures load -> emit -> load round-trips.
      dtv = days1900_to_datetime(c%tstart)
      call set_value(sec, 'start_date', dtv)
      dtv = days1900_to_datetime(c%tend)
      call set_value(sec, 'end_date', dtv)

      call set_value(sec, 'nprintday', c%nprintday)

      call add_table(sec, 'output', output)
      call set_value(output, 'swmonth', c%swmonth)
      call set_value(output, 'period',  c%period)
      call set_value(output, 'swres',   c%swres)
      call set_value(output, 'swodat',  c%swodat)
      call set_value(output, 'swyrvar', c%swyrvar)
   end subroutine emit_simulation

   subroutine emit_meteorology(doc, c)
      type(toml_table),           intent(inout) :: doc
      type(meteorology_config_t), intent(in)    :: c
      type(toml_table), pointer :: sec, et, rain, inter, temporal

      call add_table(doc, 'meteorology', sec)
      if (allocated(c%metfile)) call set_value(sec, 'file', c%metfile)
      call set_value(sec, 'lat',  c%lat)
      call set_value(sec, 'alt',  c%alt)
      call set_value(sec, 'altw', c%altw)

      call add_table(sec, 'evapotranspiration', et)
      call set_value(et, 'swetr',      c%swetr)
      call set_value(et, 'swdivide',   c%swdivide)
      call set_value(et, 'angstrom_a', c%angstroma)
      call set_value(et, 'angstrom_b', c%angstromb)

      call add_table(sec, 'temporal', temporal)
      call set_value(temporal, 'swmetdetail',  c%swmetdetail)
      call set_value(temporal, 'nmetdetail',   c%nmetdetail)
      call set_value(temporal, 'swmetfilall',  c%swmetfilall)

      call add_table(sec, 'rain', rain)
      call set_value(rain, 'swrain',   c%swrain)
      call set_value(rain, 'swetsine', c%swetsine)

      call add_table(sec, 'interception', inter)
      call set_value(inter, 'swinter', c%swinter)
   end subroutine emit_meteorology

   subroutine emit_drainage(doc, c)
      type(toml_table),       intent(inout) :: doc
      type(drainage_config_t), intent(in)   :: c
      type(toml_table), pointer :: sec, basic

      call add_table(doc, 'drainage', sec)
      call set_value(sec, 'swdra',    c%swdra)
      call set_value(sec, 'dramet',   c%dramet)
      call set_value(sec, 'swdivd',   c%swdivd)
      call set_value(sec, 'swdislay', c%swdislay)
      call set_value(sec, 'nrlevs',   c%nrlevs)
      call set_value(sec, 'altcu',    c%altcu)

      call add_table(sec, 'basic', basic)
      call set_value(basic, 'basegw', c%basegw)
      call set_value(basic, 'entres', c%entres)
      call set_value(basic, 'shape',  c%shape)
      ! NOTE: array-of-tables [[drainage.levels]] emit is Phase 4b.
   end subroutine emit_drainage

   subroutine emit_soil(doc, c)
      type(toml_table),   intent(inout) :: doc
      type(soil_config_t), intent(in)   :: c
      type(toml_table), pointer :: sec, initial

      call add_table(doc, 'soil', sec)
      call set_value(sec, 'swsophy', c%swsophy)
      call set_value(sec, 'swhyst',  c%swhyst)
      call set_value(sec, 'swinco',  c%swinco)
      call set_value(sec, 'swmacro', c%swmacro)
      call set_value(sec, 'swscal',  c%swscal)
      call set_value(sec, 'ksatexm', c%ksatexm)
      call set_value(sec, 'rsoil',   c%rsoil)
      call set_value(sec, 'reva_top', c%reva_top)

      call add_table(sec, 'initial', initial)
      call set_value(initial, 'gwli',    c%gwli)
      call set_value(initial, 'pondini', c%pondini)
      call set_value(initial, 'pondmx',  c%pondmx)
   end subroutine emit_soil

   subroutine emit_crop(doc, c)
      type(toml_table),   intent(inout) :: doc
      type(crop_config_t), intent(in)   :: c
      type(toml_table), pointer :: sec

      call add_table(doc, 'crop', sec)
      call set_value(sec, 'swcrop', c%swcrop)
      ! NOTE: array-of-tables [[crop.rotation]] emit is Phase 4b.
   end subroutine emit_crop

   !> Convert days-since-1900 (SWAP time axis) back to a TOML local-date.
   !! Uses the inverse of parse_date_to_days1900 in toml_field_helpers_mod.
   !! JD1900 = 2415021 (Julian day of 1900-01-01).
   !! Richards' algorithm (Wikipedia: Julian day number calculation).
   pure function days1900_to_datetime(days) result(dtv)
      real(real64),  intent(in) :: days
      type(toml_datetime)       :: dtv

      integer :: jd, y, m, d
      integer :: f, e, g, h

      jd = nint(days) + 2415021   ! days-since-1900 -> Julian Day Number

      ! Richards' JDN -> Gregorian calendar algorithm
      f = jd + 1401 + (((4*jd + 274277) / 146097) * 3) / 4 - 38
      e = 4 * f + 3
      g = mod(e, 1461) / 4
      h = 5 * g + 2
      d = mod(h, 153) / 5 + 1
      m = mod(h / 153 + 2, 12) + 1
      y = e / 1461 - 4716 + (14 - m) / 12

      dtv%date%year  = y
      dtv%date%month = m
      dtv%date%day   = d
      ! Leave dtv%time fields at defaults (-1) so serializer outputs a local-date.
   end function days1900_to_datetime

end module write_swap_config_mod
