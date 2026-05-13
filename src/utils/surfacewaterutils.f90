!> Module containing surface water level and storage conversion utilities
!!
!! This module provides functions for converting between:
!! - Surface water storage and water level
!! - Water level and discharge
!! - Calculating surface runoff
!!
!! @author Original SWAP team
!! @date February 2026 (modularization)
module surfacewater_utils
   use error_mod, only: fatalerr_collected
   use iso_fortran_env, only: real64
   ! SS-SWST Phase 2 Task 11: imper removed from globals; callers pass it explicitly.
   ! SS-SWC Phase 2 S-2.8: pond removed from use-list; read from state%soilwater%pond in runoff().
   ! SS-TC TC-12: dt retired from only-list; read via state%timecontrol in runoff().
   use variables, only: hqhtab, qqhtab, swdra, pondmx, rsro, rsroexp
   use swap_state_mod, only: swap_state_t
   ! swstlev_from_table is defined in surfacewater_state_mod to avoid a
   ! circular dependency (this module already pulls in swap_state_mod,
   ! which transitively depends on surfacewater_state_mod).
   use surfacewater_state_mod, only: swstlev_from_table

   implicit none

   private
   public :: wlevst, swstlev, swstlev_from_table, qhtab, runoff
contains

   !> Calculate surface water level from surface water storage using table lookup
   !!
   !! This function performs linear interpolation in the storage-level table (`sttab`)
   !! to determine the water level corresponding to a given storage amount.
   !!
   !! Called from: surfacewater.f90
   !!### Algorithm
   !!
   !! 1. Check if storage is within valid table range
   !! 2. Search through table to find bounding entries
   !! 3. Linear interpolation: \( wl = wl_i + \frac{st - st_i}{st_{i+1} - st_i} \cdot (wl_{i+1} - wl_i) \)
   !!
   !!@note
   !! Table `sttab` is ordered with highest storage at index 1, decreasing to index 22.
   !!@endnote
   !!
   !!@warning
   !! Function terminates with fatal error if storage is outside table bounds.
   !!@endwarning
   function wlevst(state, swstor) result(wlevst_r)
      implicit none

      ! Arguments
      type(swap_state_t), intent(in) :: state
      real(real64),       intent(in) :: swstor
      real(real64) :: wlevst_r

      ! Local variables
      integer :: i
      real(real64) :: dswst
      character(len=80) :: messag

      associate(sttab => state%surfacewater%sttab)

      if (swstor < sttab(22,2)) then
         messag = 'Surface water storage below bottom of table'
         call fatalerr_collected('Wlevst', messag)
      end if

      if (swstor > sttab(1,2)) then
         messag = 'Surface water storage above top of table'
         call fatalerr_collected('Wlevst', messag)
      end if

      i = 0
      do
         i = i + 1
         if (swstor >= sttab(i+1,2) .and. swstor <= sttab(i,2)) exit
      end do

      dswst = (swstor - sttab(i+1,2)) / (sttab(i,2) - sttab(i+1,2))
      wlevst_r = sttab(i+1,1) + dswst * (sttab(i,1) - sttab(i+1,1))

      end associate

   end function wlevst


   !> Calculate surface water storage from surface water level — thin wrapper
   !! over [[swstlev_from_table]] that accepts a full `swap_state_t` and pulls
   !! `state%surfacewater%sttab` from it. Canonical implementation lives in
   !! `surfacewater_state_mod` (see `swstlev_from_table` for the table-lookup
   !! algorithm, bounds-check semantics, and the inverse-of-[[wlevst]] note).
   !!
   !! Called from: surfacewater.f90 (compute paths that already have full state).
   function swstlev(state, wlev) result(swstlev_r)
      implicit none
      type(swap_state_t), intent(in) :: state
      real(real64),       intent(in) :: wlev
      real(real64) :: swstlev_r

      swstlev_r = swstlev_from_table(state%surfacewater%sttab, wlev)
   end function swstlev


   !> Calculate surface water discharge from water level using table
   !!
   !! This function performs linear interpolation in the water level-discharge table
   !! (`hqhtab`, `qqhtab`) for the current management period to determine discharge
   !! corresponding to a given water level.
   !!
   !!### Algorithm
   !!
   !! 1. Find table entries bracketing the current water level
   !! 2. Linear interpolation: \( Q = Q_i + \frac{h - h_i}{h_{i-1} - h_i} \cdot (Q_{i-1} - Q_i) \)
   !!
   !!@note
   !! Uses management period index `imper` to select appropriate rating curve.
   !! Different periods can have different level-discharge relationships.
   !!@endnote
   ! SS-SWST Phase 2 Task 11: imper_in passed explicitly (imper removed from globals).
   function qhtab(wlev, imper_in)
      implicit none

      ! Arguments
      real(real64), intent(in) :: wlev
      integer,      intent(in) :: imper_in
      real(real64) :: qhtab

      ! Local variables
      integer :: itab
      real(real64) :: dwl

      itab = 2
      do while (wlev < hqhtab(imper_in,itab))
         itab = itab + 1
      end do

      dwl = (wlev - hqhtab(imper_in,itab)) / (hqhtab(imper_in,itab-1) - hqhtab(imper_in,itab))
      qhtab = qqhtab(imper_in,itab) + dwl * (qqhtab(imper_in,itab-1) - qqhtab(imper_in,itab))

   end function qhtab


   !> Calculate surface runoff from ponded water
   !!
   !! This function calculates runoff flux when ponding exceeds the maximum allowed
   !! ponding depth (`pondmx`). Three calculation modes are supported based on
   !! drainage configuration.
   !!
   !!### Calculation Methods
   !!
   !!#### Mode 1: No surface drainage system (swdra ≠ 2)
   !!
   !! - If resistance negligible (`rsro < 0.001`): Instantaneous drainage
   !!   \[ Q_{ro} = pond - pond_{max} \]
   !!
   !! - If resistance specified: Power law drainage
   !!   \[ Q_{ro} = \frac{\Delta t}{R_{sro}} \cdot (pond - pond_{max})^{E_{sro}} \]
   !!
   !!#### Mode 2: With surface drainage system (swdra = 2)
   !!
   !! **Excess ponding (pond > max(pondmx, wls)):**
   !!   \[ Q_{ro} = \frac{\Delta t}{R_{sro}} \cdot (pond - \max(pond_{max}, wl_s))^{E_{sro}} \]
   !!
   !! **Below surface level (pond < wls):** Inundation from surface water
   !!   \[ Q_{ro} = -\min(inun_{max}, wl_s - \max(pond, pond_{max})) \]
   !!   where \( inun_{max} = swst - swstlev(pond) \)
   !!
   !!@note
   !! Positive runoff indicates drainage from soil surface to surface water system.
   !! Negative runoff indicates inundation from surface water onto soil surface.
   !!@endnote
   function runoff(state)
      implicit none

      ! Arguments
      type(swap_state_t), intent(in) :: state

      ! Function result
      real(real64) :: runoff

      ! Local variables
      real(real64) :: inun_max

      ! SS-TC TC-12: dt read via state%timecontrol tc_* alias.
      associate(tc_dt   => state%timecontrol%dt,      &  ! TC-12
                sw_wls  => state%surfacewater%wls,   &
                sw_swst => state%surfacewater%swst,  &
                ! SS-SWC Phase 2 S-2.8: pond read from state%soilwater
                pond    => state%soilwater%pond)

      runoff = 0.0_real64

      if (pond - pondmx > 0.0_real64 .and. swdra /= 2) then
         if (rsro < 1.0d-3) then
            runoff = pond - pondmx
         else
            runoff = tc_dt / rsro * (pond - pondmx)**rsroexp  ! TC-12
         end if

      else if (swdra == 2) then
         if (pond > pondmx .and. pond > sw_wls) then
            runoff = tc_dt / rsro * (pond - max(pondmx, sw_wls))**rsroexp  ! TC-12
         else if (pond < sw_wls) then
            inun_max = sw_swst - swstlev(state, pond)
            runoff = -min(inun_max, sw_wls - max(pond, pondmx))
         end if
      end if

      end associate

   end function runoff

end module surfacewater_utils
