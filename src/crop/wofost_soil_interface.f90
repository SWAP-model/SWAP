! File VersionID:
!   $Id: wofost_soil_interface.f90 307 2016-09-12 17:41:52Z kroes006 $
! ----------------------------------------------------------------------
!> WOFOST-SWAP soil interface variables.
!!
!! This module hosts the legacy coupling variables exchanged between crop
!! growth and soil nutrient routines.
!!
!! @note
!! Legacy behavior is preserved. The previous explicit `save` statement was
!! removed because module variables are implicitly static in Fortran.
      module  Wofost_Soil_Interface

      implicit none

      real(8):: idwrt, idwlv, idwst, idwso
      real(8):: iNLOSSL,  iNLOSSR, iNLOSSS, iNLOSSO
      real(8):: NdemandSoil  ! total N demand of leaves, stems and roots (kg/ha N)
      real(8):: NsupplySoil  ! total mineral N from soil and fertiliser (kg/ha N)
      real(8):: Ndemand      ! total N demand of leaves, stems and roots (kg/m2 N)
      real(8):: Nsupply      ! total mineral N from soil and fertiliser (kg/m2 N)
      real(8):: LaiCritNupt  ! Critical LAI value for N-uptake driven by transpiration or dry matter increase (kg/m2 N)

      end module  Wofost_Soil_Interface
