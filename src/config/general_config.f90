!> [general] section config: project identification, paths, screen/error switches.
module general_config_mod
   use error_mod, only: error_collection_t, ERR_VALIDATION_OUT_OF_RANGE
   use validation_mod, only: check_int_range, check_int_enum, check_not_empty
   implicit none
   private

   public :: general_config_t

   type :: general_config_t
      character(len=:), allocatable :: project
      character(len=:), allocatable :: pathwork
      character(len=:), allocatable :: pathatm
      character(len=:), allocatable :: pathcrop
      character(len=:), allocatable :: pathdrain
      integer :: swscre  = 0   !! 0=no display, 1=wb, 2=daynum
      integer :: swerror = 0   !! 0=no, 1=yes

      !> Optional comma-separated list of variables for the SPECIFIC CSV
      !! output (legacy `INLIST_CSV` in .swp Part 4). When present, the
      !! adapter forwards it to `variables%InList_csv`. When absent, the
      !! adapter falls back to a hard-coded water-balance default. Each
      !! regression case authors its own list so the fixture aggregator
      !! receives the columns it expects (e.g. GRASSDM/MOWDM for case 2).
      character(len=:), allocatable :: inlist_csv
      !> Output-file basename (legacy OUTFIL in .swp Part 1).
      !! Default 'result' matches the legacy hardcode in all regression cases.
      character(len=:), allocatable :: outfil
   contains
      procedure :: validate => general_config_validate
      procedure :: finalize => general_config_finalize
   end type general_config_t

contains

   subroutine general_config_validate(self, errors)
      class(general_config_t),  intent(in)    :: self
      type(error_collection_t), intent(inout) :: errors

      if (allocated(self%project)) then
         call check_not_empty(self%project, "general.project", errors)
      else
         call errors%append(ERR_VALIDATION_OUT_OF_RANGE, &
                            "project is required", "general.project")
      end if
      if (allocated(self%pathwork)) then
         call check_not_empty(self%pathwork, "general.pathwork", errors)
      end if

      call check_int_enum(self%swscre,  [0, 1, 2], "general.swscre",  errors)
      call check_int_enum(self%swerror, [0, 1],    "general.swerror", errors)
   end subroutine general_config_validate

   subroutine general_config_finalize(self, errors)
      class(general_config_t),  intent(inout) :: self
      type(error_collection_t), intent(inout) :: errors
      ! No derivations needed for [general].
   end subroutine general_config_finalize

end module general_config_mod
