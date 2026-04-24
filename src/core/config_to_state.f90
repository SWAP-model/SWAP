!> TEMPORARY Phase 4a adapter: populate legacy swap_state_t fields from
!! the new swap_config_t. This module disappears when state types fold
!! into config types in Phase 4+.
!!
!! Contains NO business logic. Pure field copy. If a transformation
!! needs logic, it lives in `finalize` on the config type instead.
module config_to_state_mod
   use swap_config_mod, only: swap_config_t
   use swap_state_mod,  only: swap_state_t
   implicit none
   private

   public :: config_to_state

contains

   subroutine config_to_state(config, state)
      type(swap_config_t), intent(in)    :: config
      type(swap_state_t),  intent(inout) :: state

      ! [general] -> state%time
      if (allocated(config%general%project))  state%time%project  = config%general%project
      if (allocated(config%general%pathwork)) state%time%pathwork = config%general%pathwork
      state%time%swscre  = config%general%swscre

      ! [simulation] -> state%time
      state%time%tstart    = config%simulation%tstart
      state%time%tend      = config%simulation%tend
      state%time%nprintday = config%simulation%nprintday
      state%time%swmonth   = config%simulation%swmonth
      state%time%period    = config%simulation%period
      state%time%swres     = config%simulation%swres
      state%time%swodat    = config%simulation%swodat
      ! config%simulation%swyrvar has no corresponding field in time_state_t — omitted.

      ! [meteorology] -> state%atm
      if (allocated(config%meteo%metfile)) state%atm%metfil = config%meteo%metfile
      state%atm%lat         = config%meteo%lat
      state%atm%alt         = config%meteo%alt
      state%atm%altw        = config%meteo%altw
      state%atm%swetr       = config%meteo%swetr
      state%atm%swdivide    = config%meteo%swdivide
      state%atm%swmetdetail = config%meteo%swmetdetail
      state%atm%nmetdetail  = config%meteo%nmetdetail
      state%atm%swrain      = config%meteo%swrain
      state%atm%swetsine    = config%meteo%swetsine
      state%atm%swinter     = config%meteo%swinter
      state%atm%swMetFilAll = config%meteo%swmetfilall
      state%atm%angstroma   = config%meteo%angstroma
      state%atm%angstromb   = config%meteo%angstromb

      ! [drainage] -> state%drain
      state%drain%dramet   = config%drain%dramet
      state%drain%swdivd   = config%drain%swdivd
      state%drain%swdislay = config%drain%swdislay
      state%drain%nrlevs   = config%drain%nrlevs
      ! config%drain%swdra has no corresponding field in drainage_state_t — omitted.
      ! Array fields (swdtyp, zbotdr, drares, etc.) are allocated on
      ! swap_state_init; copying them requires matching sizes.
      ! The parity test (Task 27) drives completion of array copies.
   end subroutine config_to_state

end module config_to_state_mod
