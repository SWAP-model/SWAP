!> TOML-based SWAP configuration reader.
!!
!! Loads SWAP configuration from TOML files and maps supported keys directly to
!! `swap_state_t`, without using the legacy `variables` module.
!!
!! The module is intentionally state-first and can be reused by future BMI-based
!! in-memory configuration injection.
module readswaptoml_mod
   use tomlf, only: toml_table, toml_array, toml_datetime, toml_error, toml_load, get_value
   use tomlf_type, only: len
   use swap_state_mod, only: swap_state_t
   use readdrainagetoml_mod, only: ReadDrainageToml_state
   use swap_log, only: log_info
   implicit none
   private

   public :: ReadSwapToml_state

contains

!> Read SWAP configuration from TOML into explicit state.
!!
!! @param[inout] state SWAP model state container
!! @param[in]    swp_toml_path Path to main TOML configuration file
   subroutine ReadSwapToml_state(state, swp_toml_path)
      type(swap_state_t), intent(inout) :: state
      character(len=*),    intent(in)    :: swp_toml_path

      type(toml_table), allocatable :: doc
      type(toml_table), pointer     :: tab, subtab
      type(toml_array), pointer     :: arr
      type(toml_error), allocatable :: err
      type(toml_datetime)           :: dtv

      integer :: istat, i, n
      character(len=:), allocatable :: s

      call toml_load(doc, trim(swp_toml_path), error=err)
      if (allocated(err)) then
         call fatalerr('ReadSwapToml_state', trim(err%message))
      end if

      ! ------------------------------------------------------------------
      ! [general]
      ! ------------------------------------------------------------------
      call get_value(doc, 'general', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'project', s, stat=istat)
         if (istat == 0) state%time%project = trim(s)
         call get_value(tab, 'swscre',  state%time%swscre,  default=state%time%swscre)

         call get_value(tab, 'paths', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'work', s, stat=istat)
            if (istat == 0) state%time%pathwork = trim(s)
            call get_value(subtab, 'atmosphere', s, stat=istat)
            if (istat == 0) state%atm%pathatm = trim(s)
            call get_value(subtab, 'crop', s, stat=istat)
            if (istat == 0) state%crop%pathcrop = trim(s)
            call get_value(subtab, 'drain', s, stat=istat)
            if (istat == 0) state%drain%pathdrain = trim(s)
         end if
      end if

      state%time%swpfile = trim(swp_toml_path)

      ! ------------------------------------------------------------------
      ! [simulation]
      ! ------------------------------------------------------------------
      call get_value(doc, 'simulation', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'start_date', dtv, stat=istat)
         if (istat == 0) state%time%tstart = toml_datetime_to_t1900(dtv)

         call get_value(tab, 'end_date', dtv, stat=istat)
         if (istat == 0) state%time%tend = toml_datetime_to_t1900(dtv)
      end if

      ! ------------------------------------------------------------------
      ! [output]
      ! ------------------------------------------------------------------
      call get_value(doc, 'output', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'file_prefix', s, stat=istat)
         if (istat == 0) state%time%outfil = trim(s)
         call get_value(tab, 'nprintday',   state%time%nprintday, default=state%time%nprintday)
         call get_value(tab, 'swheader',    state%time%swheader, default=state%time%swheader)

         call get_value(tab, 'timing', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swmonth', state%time%swmonth, default=state%time%swmonth)
            call get_value(subtab, 'period',  state%time%period,  default=state%time%period)
            call get_value(subtab, 'swres',   state%time%swres,   default=state%time%swres)
            call get_value(subtab, 'swodat',  state%time%swodat,  default=state%time%swodat)
         end if
      end if

      ! ------------------------------------------------------------------
      ! [meteorology]
      ! ------------------------------------------------------------------
      call get_value(doc, 'meteorology', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'file', s, stat=istat)
         if (istat == 0) state%atm%metfil = trim(s)
         call get_value(tab, 'lat',  state%atm%lat,    default=state%atm%lat)

         call get_value(tab, 'evapotranspiration', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swetr',      state%atm%swetr,      default=state%atm%swetr)
            call get_value(subtab, 'alt',        state%atm%alt,        default=state%atm%alt)
            call get_value(subtab, 'altw',       state%atm%altw,       default=state%atm%altw)
            call get_value(subtab, 'angstrom_a', state%atm%angstroma,  default=state%atm%angstroma)
            call get_value(subtab, 'angstrom_b', state%atm%angstromb,  default=state%atm%angstromb)
            call get_value(subtab, 'swdivide',   state%atm%swdivide,   default=state%atm%swdivide)
         end if

         call get_value(tab, 'temporal', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swmetdetail', state%atm%swmetdetail, default=state%atm%swmetdetail)
            call get_value(subtab, 'nmetdetail',  state%atm%nmetdetail,  default=state%atm%nmetdetail)
            call get_value(subtab, 'swetsine',    state%atm%swetsine,    default=state%atm%swetsine)
         end if

         call get_value(tab, 'rainfall', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swrain', state%atm%swrain, default=state%atm%swrain)
            call get_value(subtab, 'rainfall_file', s, stat=istat)
            if (istat == 0) state%atm%rainfil = trim(s)
         end if
      end if

      ! ------------------------------------------------------------------
      ! [crop] and [crop.rotation]
      ! ------------------------------------------------------------------
      call get_value(doc, 'crop', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'rotation', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'crop_file', arr, requested=.false., stat=istat)
            if (istat == 0) then
               n = len(arr)
               if (n > 0) then
                  state%ncrop = n
                  call ensure_crop_arrays(state, n)
                  do i = 1, n
                     call get_value(arr, i, s, stat=istat)
                     if (istat == 0) state%crop%cropfil(i) = trim(s)
                  end do
               end if
            end if

            call get_value(subtab, 'start_date', arr, requested=.false., stat=istat)
            if (istat == 0 .and. allocated(state%crop%cropstart)) then
               n = min(len(arr), size(state%crop%cropstart))
               do i = 1, n
                  call get_value(arr, i, s, stat=istat)
                  if (istat == 0) state%crop%cropstart(i) = iso_date_to_t1900(trim(s))
               end do
            end if

            call get_value(subtab, 'end_date', arr, requested=.false., stat=istat)
            if (istat == 0 .and. allocated(state%crop%cropend)) then
               n = min(len(arr), size(state%crop%cropend))
               do i = 1, n
                  call get_value(arr, i, s, stat=istat)
                  if (istat == 0) state%crop%cropend(i) = iso_date_to_t1900(trim(s))
               end do
            end if
         end if
      end if

      ! ------------------------------------------------------------------
      ! [soil]
      ! ------------------------------------------------------------------
      call get_value(doc, 'soil', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'initial', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'gwli', state%soil%gwli, default=state%soil%gwli)
            state%soil%gwl = state%soil%gwli
         end if

         call get_value(tab, 'surface', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'pondmx', state%soil%pondmx, default=state%soil%pondmx)
         end if

         call get_value(tab, 'evaporation', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'cofredbl', state%soil%cofred, default=state%soil%cofred)
            call get_value(subtab, 'cofredbo', state%soil%cofred, default=state%soil%cofred)
         end if

         call get_value(tab, 'hysteresis', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swhyst', state%soil%swhyst, default=state%soil%swhyst)
            call get_value(subtab, 'tau', state%soil%tau, default=state%soil%tau)
         end if

         call get_value(tab, 'numerical', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'gwlconv', state%soil%gwlconv, default=state%soil%gwlconv)
            call get_value(subtab, 'critdevh1cp', state%soil%CritDevh1Cp, default=state%soil%CritDevh1Cp)
            call get_value(subtab, 'critdevh2cp', state%soil%CritDevh2Cp, default=state%soil%CritDevh2Cp)
            call get_value(subtab, 'maxit', state%soil%msteps, default=state%soil%msteps)
         end if
      end if

      ! ------------------------------------------------------------------
      ! [drainage]
      ! ------------------------------------------------------------------
      call get_value(doc, 'drainage', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'drainage_file', s, stat=istat)
         if (istat == 0) state%drain%drfil = trim(s)
         call read_drainage_toml_if_present(state)
      end if

      ! ------------------------------------------------------------------
      ! [boundary]
      ! ------------------------------------------------------------------
      call get_value(doc, 'boundary', tab, requested=.false., stat=istat)
      if (istat == 0) then
         call get_value(tab, 'bottom', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swbotb', state%boundary%swbotb, default=state%boundary%swbotb)
            call get_value(subtab, 'swqhbot', state%boundary%swqhbot, default=state%boundary%swqhbot)
            call get_value(subtab, 'deepgw', state%boundary%deepgw, default=state%boundary%deepgw)
            call get_value(subtab, 'hbot', state%boundary%hbot, default=state%boundary%hbot)
            call get_value(subtab, 'qbot', state%boundary%qbot, default=state%boundary%qbot)
         end if

         call get_value(tab, 'top', subtab, requested=.false., stat=istat)
         if (istat == 0) then
            call get_value(subtab, 'swpondmx', state%boundary%swpondmx, default=state%boundary%swpondmx)
            call get_value(subtab, 'pondmx', state%boundary%pondmx, default=state%boundary%pondmx)
            call get_value(subtab, 'rsro', state%boundary%rsro, default=state%boundary%rsro)
            call get_value(subtab, 'runon', state%boundary%runon, default=state%boundary%runon)
            call get_value(subtab, 'qtop', state%boundary%qtop, default=state%boundary%qtop)
         end if
      end if

      call log_info('ReadSwapToml', 'Loaded TOML configuration: ' // trim(swp_toml_path))

   end subroutine ReadSwapToml_state


!> Ensure crop arrays are allocated for `n` crops.
   subroutine ensure_crop_arrays(state, n)
      type(swap_state_t), intent(inout) :: state
      integer,            intent(in)    :: n

      if (.not. allocated(state%crop%cropfil)) then
         allocate(state%crop%cropfil(n))
      else if (size(state%crop%cropfil) /= n) then
         deallocate(state%crop%cropfil)
         allocate(state%crop%cropfil(n))
      end if

      if (.not. allocated(state%crop%cropstart)) then
         allocate(state%crop%cropstart(n))
      else if (size(state%crop%cropstart) /= n) then
         deallocate(state%crop%cropstart)
         allocate(state%crop%cropstart(n))
      end if

      if (.not. allocated(state%crop%cropend)) then
         allocate(state%crop%cropend(n))
      else if (size(state%crop%cropend) /= n) then
         deallocate(state%crop%cropend)
         allocate(state%crop%cropend(n))
      end if
   end subroutine ensure_crop_arrays


!> Load detailed drainage TOML file when present.
   subroutine read_drainage_toml_if_present(state)
      type(swap_state_t), intent(inout) :: state

      character(len=:), allocatable :: filename
      character(len=:), allocatable :: drainage_path
      logical :: exists

      if (len_trim(state%drain%drfil) == 0) return

      filename = trim(state%drain%drfil)
      if (index(filename, '.toml') == 0) filename = trim(filename) // '.dra.toml'

      if (index(filename, '/') == 1) then
         drainage_path = trim(filename)
      else
         if (len_trim(state%drain%pathdrain) > 0) then
            if (state%drain%pathdrain(len_trim(state%drain%pathdrain):len_trim(state%drain%pathdrain)) == '/') then
               drainage_path = trim(state%drain%pathdrain) // trim(filename)
            else
               drainage_path = trim(state%drain%pathdrain) // '/' // trim(filename)
            end if
         else
            drainage_path = trim(filename)
         end if
      end if

      inquire(file=drainage_path, exist=exists)
      if (.not. exists) then
         call log_info('ReadSwapToml', 'Drainage TOML not found, skipping detailed read: ' // trim(drainage_path))
         return
      end if

      call ReadDrainageToml_state(state%drain, drainage_path)

   end subroutine read_drainage_toml_if_present


!> Convert TOML datetime to SWAP time (days since 1900-01-01).
   real(8) function toml_datetime_to_t1900(dt)
      type(toml_datetime), intent(in) :: dt
      integer :: y, m, d
      real(8) :: frac

      y = dt%date%year
      m = dt%date%month
      d = dt%date%day

      frac = 0.0d0
      if (dt%time%hour >= 0) then
         frac = (dble(dt%time%hour) + dble(max(0, dt%time%minute))/60.0d0 + &
                 dble(max(0, dt%time%second))/3600.0d0) / 24.0d0
      end if

      toml_datetime_to_t1900 = dble(gregorian_to_jdn(y, m, d) - gregorian_to_jdn(1900, 1, 1)) + frac
   end function toml_datetime_to_t1900


!> Convert ISO date string (YYYY-MM-DD) to SWAP time (days since 1900-01-01).
   real(8) function iso_date_to_t1900(s)
      character(len=*), intent(in) :: s
      integer :: y, m, d

      if (len_trim(s) < 10) then
         iso_date_to_t1900 = 0.0d0
         return
      end if

      read(s(1:4), *, err=100) y
      read(s(6:7), *, err=100) m
      read(s(9:10), *, err=100) d

      iso_date_to_t1900 = dble(gregorian_to_jdn(y, m, d) - gregorian_to_jdn(1900, 1, 1))
      return
100   continue
      iso_date_to_t1900 = 0.0d0
   end function iso_date_to_t1900


!> Gregorian date to Julian Day Number.
   integer function gregorian_to_jdn(y, m, d)
      integer, intent(in) :: y, m, d
      integer :: a, y1, m1

      a = (14 - m) / 12
      y1 = y + 4800 - a
      m1 = m + 12*a - 3
      gregorian_to_jdn = d + (153*m1 + 2)/5 + 365*y1 + y1/4 - y1/100 + y1/400 - 32045
   end function gregorian_to_jdn

end module readswaptoml_mod
