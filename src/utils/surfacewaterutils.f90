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
   use variables, only: sttab, imper, hqhtab, qqhtab, swdra, pond, pondmx, rsro, rsroexp, dt
   use swap_state_mod, only: swap_state_t

   implicit none

   private
   public :: wlevst, swstlev, qhtab, runoff
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
   function wlevst(swstor)
      implicit none
      
      ! Arguments
      real(real64), intent(in) :: swstor
      real(real64) :: wlevst
      
      ! Local variables
      integer :: i
      real(real64) :: dswst
      character(len=80) :: messag

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
      wlevst = sttab(i+1,1) + dswst * (sttab(i,1) - sttab(i+1,1))
      
   end function wlevst


   !> Calculate surface water storage from surface water level
   !!
   !! This function performs linear interpolation in the level-storage table (`sttab`)
   !! to determine the storage amount corresponding to a given water level.
   !!
      !! Called from: surfacewater.f90, readswap.f90
   !!### Algorithm
   !!
   !! 1. Check if level is within valid table range
   !! 2. Search through table to find bounding entries
   !! 3. Linear interpolation: \( st = st_i + \frac{wl - wl_i}{wl_{i+1} - wl_i} \cdot (st_{i+1} - st_i) \)
   !!
   !!@note
   !! This is the inverse operation of [[wlevst]].
   !!@endnote
   !!
   !!@warning
   !! Function terminates with fatal error if water level is outside table bounds.
   !!@endwarning
   function swstlev(wlev)
      implicit none
      
      ! Arguments
      real(real64), intent(in) :: wlev
      real(real64) :: swstlev
      
      ! Local variables
      integer :: i
      real(real64) :: dwl
      character(len=200) :: messag

      if (wlev < sttab(22,1)) then
         messag = 'Surface water storage below bottom of table'
         call fatalerr_collected('swstlev', messag)
      end if
      
      if (wlev > sttab(1,1)) then
         messag = 'Surface water storage above top of table'
         call fatalerr_collected('swstlev', messag)
      end if

      i = 0
      do
         i = i + 1
         if (wlev >= sttab(i+1,1) .and. wlev <= sttab(i,1)) exit
      end do
      
      dwl = (wlev - sttab(i+1,1)) / (sttab(i,1) - sttab(i+1,1))
      swstlev = sttab(i+1,2) + dwl * (sttab(i,2) - sttab(i+1,2))
      
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
   function qhtab(wlev)
      implicit none
      
      ! Arguments
      real(real64), intent(in) :: wlev
      real(real64) :: qhtab
      
      ! Local variables
      integer :: itab
      real(real64) :: dwl

      itab = 2
      do while (wlev < hqhtab(imper,itab))
         itab = itab + 1
      end do
      
      dwl = (wlev - hqhtab(imper,itab)) / (hqhtab(imper,itab-1) - hqhtab(imper,itab))
      qhtab = qqhtab(imper,itab) + dwl * (qqhtab(imper,itab-1) - qqhtab(imper,itab))
      
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

      associate(sw_wls => state%surfacewater%wls, sw_swst => state%surfacewater%swst)

      runoff = 0.0_real64

      if (pond - pondmx > 0.0_real64 .and. swdra /= 2) then
         if (rsro < 1.0d-3) then
            runoff = pond - pondmx
         else
            runoff = dt / rsro * (pond - pondmx)**rsroexp
         end if

      else if (swdra == 2) then
         if (pond > pondmx .and. pond > sw_wls) then
            runoff = dt / rsro * (pond - max(pondmx, sw_wls))**rsroexp
         else if (pond < sw_wls) then
            inun_max = sw_swst - swstlev(pond)
            runoff = -min(inun_max, sw_wls - max(pond, pondmx))
         end if
      end if

      end associate

   end function runoff

end module surfacewater_utils
