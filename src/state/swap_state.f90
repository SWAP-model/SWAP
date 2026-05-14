!> @file swap_state.f90
!! Top-level aggregator of all per-subsystem typed state records.
!! Phase 1 of the state-migration arc starts with only the
!! surface-water subsystem; subsequent migration arcs (water-balance,
!! crop, drainage, atmosphere, …) add their own fields here.
!!
!! Threaded through subroutine signatures from `swap_main` down,
!! replacing the implicit shared state held in `variables.f90`.

module swap_state_mod
   use iso_c_binding,           only: c_double
   use surfacewater_state_mod,  only: surfacewater_state_t
   use drainage_state_mod,      only: drainage_state_t
   use solute_state_mod,        only: solute_state_t
   use heat_state_mod,          only: heat_state_t
   use soilwater_state_mod,     only: soilwater_state_t
   use atmosphere_state_mod,    only: atmosphere_state_t
   use tillage_state_mod,       only: tillage_state_t
   use timecontrol_state_mod,   only: timecontrol_state_t
   use mesh_state_mod,          only: mesh_state_t
   use crop_state_mod,          only: crop_state_t
   use nutrients_state_mod,     only: nutrients_state_t
   implicit none
   private
   public :: swap_state_t

   type :: swap_state_t
      type(surfacewater_state_t) :: surfacewater
      type(drainage_state_t)     :: drainage
      type(solute_state_t)       :: solute
      type(heat_state_t)         :: heat
      type(soilwater_state_t)    :: soilwater
      type(atmosphere_state_t)   :: atmosphere
      type(tillage_state_t)      :: tillage
      type(timecontrol_state_t)  :: timecontrol
      type(mesh_state_t)         :: mesh
      type(crop_state_t)         :: crop
      type(nutrients_state_t)    :: nutrients

      ! [SS-BMI2] water balance output stream (was: written directly to .inc from outinc)
      real(c_double),    allocatable :: water_balance_row(:)
      character(len=32), allocatable :: water_balance_columns(:)
      integer                        :: water_balance_n_cols = 0

      ! [SS-BMI2] crop output stream (cropoutput/OutCropFixed/OutWofost/OutGrass).
      ! Lives at top-level swap_state_t; no crop_state_t exists today.
      ! N is dynamic: varies by croptype (fixed/wofost/grass) and user config.
      real(c_double),    allocatable :: crop_output_row(:)
      character(len=32), allocatable :: crop_output_columns(:)
      integer                        :: crop_output_n_cols = 0
   end type swap_state_t

end module swap_state_mod
