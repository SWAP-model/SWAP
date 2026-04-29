!> Reader for the [simulation] section of a SWAP TOML.
module read_simulation_toml_mod
   use tomlf, only: toml_table, toml_datetime, get_value
   use simulation_config_mod, only: simulation_config_t
   use toml_field_helpers_mod, only: get_table, &
                                     get_optional_int_with_default, &
                                     get_optional_real_with_default, &
                                     parse_date_to_days1900
   use error_mod, only: error_collection_t, ERR_PARSE_MISSING_REQUIRED, ERR_PARSE_TYPE_MISMATCH
   implicit none
   private

   public :: read_simulation_toml

contains

   subroutine read_simulation_toml(doc, config, errors)
      type(toml_table), pointer,    intent(in)    :: doc
      type(simulation_config_t),    intent(inout) :: config
      type(error_collection_t),     intent(inout) :: errors

      type(toml_table), pointer :: sec, timing, num
      type(toml_datetime)       :: dtv
      integer                   :: stat

      call get_table(doc, 'simulation', sec, 'simulation', errors)
      if (.not. associated(sec)) return

      call get_value(sec, 'start_date', dtv, stat=stat)
      if (stat == 0) then
         config%tstart = parse_date_to_days1900(dtv)
      else
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "simulation.start_date is required", &
                            "simulation.start_date")
      end if

      call get_value(sec, 'end_date', dtv, stat=stat)
      if (stat == 0) then
         config%tend = parse_date_to_days1900(dtv)
      else
         call errors%append(ERR_PARSE_MISSING_REQUIRED, &
                            "simulation.end_date is required", &
                            "simulation.end_date")
      end if

      call get_optional_int_with_default(sec, 'nprintday', config%nprintday, 1, 'simulation.nprintday', errors)

      call get_table(sec, 'output', timing, 'simulation.output', errors)
      if (associated(timing)) then
         call get_optional_int_with_default(timing, 'swmonth', config%swmonth, 0, 'simulation.output.swmonth', errors)
         call get_optional_int_with_default(timing, 'period',  config%period,  1, 'simulation.output.period',  errors)
         call get_optional_int_with_default(timing, 'swres',   config%swres,   0, 'simulation.output.swres',   errors)
         call get_optional_int_with_default(timing, 'swodat',  config%swodat,  0, 'simulation.output.swodat',  errors)
         call get_optional_int_with_default(timing, 'swyrvar', config%swyrvar, 0, 'simulation.output.swyrvar', errors)
      end if

      call get_table(sec, 'numerical', num, 'simulation.numerical', errors)
      if (associated(num)) then
         call get_optional_real_with_default(num, 'dt',     config%numerical%dt, &
                                             config%numerical%dt,     'simulation.numerical.dt',     errors)
         call get_optional_real_with_default(num, 'dtmin',  config%numerical%dtmin, &
                                             config%numerical%dtmin,  'simulation.numerical.dtmin',  errors)
         call get_optional_real_with_default(num, 'dtmax',  config%numerical%dtmax, &
                                             config%numerical%dtmax,  'simulation.numerical.dtmax',  errors)
         call get_optional_int_with_default(num, 'maxit',     config%numerical%MaxIt, &
                                            config%numerical%MaxIt,     'simulation.numerical.MaxIt',     errors)
         call get_optional_int_with_default(num, 'maxbacktr', config%numerical%MaxBackTr, &
                                            config%numerical%MaxBackTr, 'simulation.numerical.MaxBackTr', errors)
         call get_optional_real_with_default(num, 'taccur', config%numerical%taccur, &
                                             config%numerical%taccur, 'simulation.numerical.taccur', errors)
      end if
   end subroutine read_simulation_toml

end module read_simulation_toml_mod
