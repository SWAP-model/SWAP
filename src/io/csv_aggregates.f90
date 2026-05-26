!> Water-balance aggregates for CSV output. Extracted verbatim from the
!! former inline arithmetic in swap_csv_output.f90's fill_values so it can be
!! unit-tested.
module csv_aggregates_mod
   use iso_fortran_env, only: real64
   implicit none
   private
   public :: water_balance_dev
contains
   pure subroutine water_balance_dev(volact, pond, ssnow, vol_old, pond_old, &
         snow_old, igrai, isnrai, igsnow, igird, irunon, iqssdi, iintc, iruno, &
         irunocn, iqrot, ievap, isubl, iqdra, iqbot, dstor, baldev)
      real(real64), intent(in)  :: volact, pond, ssnow, vol_old, pond_old, snow_old
      real(real64), intent(in)  :: igrai, isnrai, igsnow, igird, irunon, iqssdi
      real(real64), intent(in)  :: iintc, iruno, irunocn, iqrot, ievap, isubl, iqdra, iqbot
      real(real64), intent(out) :: dstor, baldev
      dstor  = (volact + pond + ssnow) - (vol_old + pond_old + snow_old)
      baldev = (igrai + isnrai + igsnow + igird + irunon + iqssdi) - dstor &
             - (iintc + iruno + irunocn + iqrot + ievap + isubl + iqdra + (-1.0_real64*iqbot))
   end subroutine water_balance_dev
end module csv_aggregates_mod
