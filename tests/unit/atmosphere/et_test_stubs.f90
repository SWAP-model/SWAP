! Minimal link-only stubs for legacy subroutines called by PenMon (the
! public wrapper in et_mod) but not by PenMon_calc.  These are never
! reached by the characterization test, which calls PenMon_calc directly.
subroutine astro(iday, lat, avrad, dayl, daylp, sinld, cosld, &
                 difpp, atmtr, dsinbe)
    implicit none
    integer, intent(in)  :: iday
    real(8), intent(in)  :: lat, avrad
    real(8), intent(out) :: dayl, daylp, sinld, cosld
    real(8), intent(out) :: difpp, atmtr, dsinbe
    dayl   = 0.0d0; daylp  = 0.0d0; sinld  = 0.0d0; cosld  = 0.0d0
    difpp  = 0.0d0; atmtr  = 0.0d0; dsinbe = 0.0d0
end subroutine astro

subroutine warn(modul, messag, logf, swscre)
    implicit none
    character(len=*), intent(in) :: modul, messag
    integer, intent(in) :: logf, swscre
end subroutine warn
