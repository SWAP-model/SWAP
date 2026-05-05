!> [output.csv] section config: CSV output switches and column lists.
module output_csv_config_mod
   use error_mod, only: error_collection_t
   use validation_mod, only: check_int_enum
   implicit none
   private

   public :: output_csv_config_t

   type :: output_csv_config_t
      !> Enable daily CSV output (legacy SWCSV). Default 1 (on).
      integer :: enabled     = 1
      !> Enable depth-profile CSV output (legacy SWCSV_TZ). Default 0 (off).
      integer :: enabled_tz  = 0
      !> Comma-separated column list for daily CSV (legacy INLIST_CSV).
      !! Default: water-balance summary columns matching hupselbrook baseline.
      character(len=:), allocatable :: inlist
      !> Comma-separated column list for depth-profile CSV (legacy INLIST_CSV_TZ).
      !! Default: volumetric water content, pressure head, solute concentration.
      character(len=:), allocatable :: inlist_tz
   contains
      procedure :: validate => output_csv_config_validate
      procedure :: finalize => output_csv_config_finalize
   end type output_csv_config_t

contains

   subroutine output_csv_config_validate(self, errors)
      class(output_csv_config_t), intent(in)    :: self
      type(error_collection_t),   intent(inout) :: errors

      call check_int_enum(self%enabled,    [0, 1], "output.csv.enabled",    errors)
      call check_int_enum(self%enabled_tz, [0, 1], "output.csv.enabled_tz", errors)
   end subroutine output_csv_config_validate

   subroutine output_csv_config_finalize(self, errors)
      class(output_csv_config_t), intent(inout) :: self
      type(error_collection_t),   intent(inout) :: errors
      ! Apply defaults for unset allocatable fields.
      if (.not. allocated(self%inlist)) then
         self%inlist = 'rain,irrig,interc,runoff,drainage,dstor,epot,eact,tpot,tact,qbottom,gwl'
      end if
      if (.not. allocated(self%inlist_tz)) then
         self%inlist_tz = 'wc,h,conc'
      end if
   end subroutine output_csv_config_finalize

end module output_csv_config_mod
