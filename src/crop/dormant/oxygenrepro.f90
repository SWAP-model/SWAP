!> @file src/crop/dormant/oxygenrepro.f90
!! @brief DORMANT — Bartholomeus reproduction-function oxygen stress
!!        (`swoxygen=2/swoxygentype=2` path).
!!
!! ## Status: DORMANT
!!
!! Moved here from `src/crop/oxygenstress.f90` on 2026-05-25. The
!! Bartholomeus reproduction-function oxygen-stress compute path
!! (selected when `swoxygen == 2` and `swoxygentype == 2`) has no live
!! dispatch in the TOML pipeline: rootextraction.f90 already
!! `fatalerr_collected`s when that branch is selected (see the stub
!! around line 152 in src/crop/rootextraction.f90 — Task 3 follow-up).
!! None of the regression cases exercise it.
!!
!! Excluded from `meson.build` and `tests/unit/meson.build` — this file
!! is NOT compiled.
!!
!! ## What it contains
!!
!! - `subroutine oxygen_dat(SwTopSub, NrStaring, OxygenSlope,
!!     OxygenIntercept)` — populates 6-element OxygenSlope /
!!     OxygenIntercept arrays from one of 4 (18×6) hard-coded
!!     coefficient tables (xtop/ytop for topsoil, xsub/ysub for
!!     subsoil, indexed by Staring-series number 1..18).
!! - `subroutine OxygenReproFunction(OxygenSlope, OxygenIntercept,
!!     theta, thetas, tsoil, node, z, dz, rwu_factor, state)` — the
!!     compute body: per-node oxygen-stress rwu_factor as a polynomial
!!     in soil temperature, depth, and (per-node weighted-mean) gas-
!!     filled porosity. The polynomial coefficients (intercept and
!!     slope) come from oxygen_dat above.
!!
!! Live dependencies (re-needed at reactivation):
!! - `state%mesh%zbotcp(node)` (depth normalization in
!!   `OxygenReproFunction`).
!!
!! ## Reactivation checklist
!!
!! 1. Restore the dispatch site in `RootExtraction(state)`: replace the
!!    `fatalerr_collected('RootExtraction', 'swoxygen=2/swoxygentype=2
!!    (OxygenReproFunction) is dormant — ...')` stub (around
!!    rootextraction.f90:155) with the two calls that the original
!!    legacy code performed:
!!       call oxygen_dat(cfg%swtopsub, cfg%nrstaring,
!!                       state%legacy%OxygenSlope,
!!                       state%legacy%OxygenIntercept)
!!       call OxygenReproFunction(state%legacy%OxygenSlope,
!!                                state%legacy%OxygenIntercept,
!!                                state%soilwater%theta,
!!                                state%soilwater%thetas (vg_params),
!!                                state%heat%tsoil, node,
!!                                state%mesh%z, state%mesh%dz,
!!                                alpwet, state)
!! 2. Wire two new fields onto a `crop_oxygen_repro_config_t` (or just
!!    extend `crop_oxygen_state_t`):
!!       state%crop%oxygen%OxygenSlope(6)
!!       state%crop%oxygen%OxygenIntercept(6)
!!    These replace the legacy `OxygenSlope(6)`/`OxygenIntercept(6)`
!!    arrays that were declared in variables.f90 (already retired with
!!    the swoxygen=2 path).
!! 3. Wire `swtopsub` and `nrstaring` from the per-rotation config
!!    (cropfixed_config_t / cropwofost_config_t already have the
!!    slots in their swoxygen=2 sub-record but no live writer
!!    because the validator stub-errors the swoxygen=2 path today).
!! 4. Add this file to `meson.build` and `tests/unit/meson.build`
!!    under the Crop section.
!! 5. Add `use oxygenrepro_dormant_mod, only: oxygen_dat,
!!    OxygenReproFunction` to `rootextraction.f90` and any direct
!!    callers.
!! 6. Run `pixi run check-fast` — the regression suite does NOT cover
!!    this branch; add a dedicated TOML fixture before relying on it.
!!
!! Original SVN revision:
!!   $Id: oxygenstress.f90 378 2018-05-08 13:50:52Z heine003 $

module oxygenrepro_dormant_mod
   use swap_state_mod, only: swap_state_t
   implicit none
   private
   public :: oxygen_dat, OxygenReproFunction

contains

! ----------------------------------------------------------------------
      subroutine oxygen_dat (SwTopSub,NrStaring,OxygenSlope,            &
     &                       OxygenIntercept)
! ----------------------------------------------------------------------
!     date               : december 2009
!     purpose            : set parameter values of metafunction for oxygenstress
! ----------------------------------------------------------------------

! --- global
      integer SwTopSub,NrStaring
      real(8) OxygenSlope(6),OxygenIntercept(6)

! --- local
      integer i
      real(8) xtop(18,6),xsub(18,6),ytop(18,6),ysub(18,6)

      data (xtop(1,i),i=1,6)                                            &
     & /5.07d-03,2.40d+02,-4.39d+00,-6.31d+02,1.67d+00,9.08d+02/
      data (xtop(2,i),i=1,6)                                            &
     & /1.20d-02,3.26d+02,-8.91d+00,-9.29d+02,2.49d+00,1.64d+03/
      data (xtop(3,i),i=1,6)                                            &
     & /1.21d-02,3.64d+02,-9.09d+00,-9.90d+02,2.60d+00,1.70d+03/
      data (xtop(4,i),i=1,6)                                            &
     & /1.67d-02,4.42d+02,-1.21d+01,-1.26d+03,3.34d+00,2.20d+03/
      data (xtop(5,i),i=1,6)                                            &
     & /4.17d-03,1.93d+02,-3.72d+00,-5.28d+02,1.44d+00,7.77d+02/
      data (xtop(6,i),i=1,6)                                            &
     & /2.11d-02,5.02d+02,-1.51d+01,-1.40d+03,3.69d+00,2.71d+03/
      data (xtop(7,i),i=1,6)                                            &
     & /2.86d-02,5.57d+02,-1.97d+01,-1.67d+03,4.50d+00,3.41d+03/
      data (xtop(8,i),i=1,6)                                            &
     & /1.75d-02,5.55d+02,-1.32d+01,-1.34d+03,3.31d+00,2.48d+03/
      data (xtop(9,i),i=1,6)                                            &
     & /2.34d-02,6.00d+02,-1.70d+01,-1.40d+03,3.42d+00,3.09d+03/
      data (xtop(10,i),i=1,6)                                           &
     & /3.12d-02,6.62d+02,-2.20d+01,-1.53d+03,3.65d+00,3.91d+03/
      data (xtop(11,i),i=1,6)                                           &
     & /2.58d-02,6.42d+02,-1.83d+01,-1.61d+03,4.00d+00,3.26d+03/
      data (xtop(12,i),i=1,6)                                           &
     & /2.50d-02,6.53d+02,-1.79d+01,-1.45d+03,3.40d+00,3.21d+03/
      data (xtop(13,i),i=1,6)                                           &
     & /2.53d-02,5.97d+02,-1.79d+01,-1.72d+03,4.58d+00,3.18d+03/
      data (xtop(14,i),i=1,6)                                           &
     & /2.82d-02,7.10d+02,-2.04d+01,-1.54d+03,3.56d+00,3.71d+03/
      data (xtop(15,i),i=1,6)                                           &
     & /2.08d-02,4.72d+02,-1.46d+01,-1.37d+03,3.65d+00,2.57d+03/
      data (xtop(16,i),i=1,6)                                           &
     & /1.99d-02,4.80d+02,-1.39d+01,-1.34d+03,3.47d+00,2.45d+03/
      data (xtop(17,i),i=1,6)                                           &
     & /2.27d-02,6.31d+02,-1.62d+01,-1.58d+03,3.91d+00,2.91d+03/
      data (xtop(18,i),i=1,6)                                           &
     & /2.23d-02,6.50d+02,-1.60d+01,-1.74d+03,4.45d+00,2.89d+03/

      data (xsub(1,i),i=1,6)                                            &
     & /7.21d-04,1.76d+02,-1.59d+00,-4.33d+02,1.16d+00,4.52d+02/
      data (xsub(2,i),i=1,6)                                            &
     & /3.75d-03,2.25d+02,-3.61d+00,-5.96d+02,1.59d+00,7.91d+02/
      data (xsub(3,i),i=1,6)                                            &
     & /7.06d-03,2.86d+02,-5.91d+00,-7.52d+02,2.00d+00,1.19d+03/
      data (xsub(4,i),i=1,6)                                            &
     & /1.29d-02,3.74d+02,-9.74d+00,-1.03d+03,2.72d+00,1.82d+03/
      data (xsub(5,i),i=1,6)                                            &
     & /2.92d-03,1.86d+02,-3.06d+00,-5.07d+02,1.39d+00,6.85d+02/
      data (xsub(6,i),i=1,6)                                            &
     & /3.00d-02,5.41d+02,-2.06d+01,-1.69d+03,4.61d+00,3.55d+03/
      data (xsub(7,i),i=1,6)                                            &
     & /2.76d-02,6.63d+02,-1.95d+01,-1.67d+03,4.15d+00,3.48d+03/
      data (xsub(8,i),i=1,6)                                            &
     & /1.85d-02,4.88d+02,-1.34d+01,-1.34d+03,3.50d+00,2.43d+03/
      data (xsub(9,i),i=1,6)                                            &
     & /2.08d-02,5.46d+02,-1.50d+01,-1.50d+03,3.92d+00,2.72d+03/
      data (xsub(10,i),i=1,6)                                           &
     & /1.85d-02,5.94d+02,-1.39d+01,-1.45d+03,3.60d+00,2.60d+03/
      data (xsub(11,i),i=1,6)                                           &
     & /2.29d-02,6.33d+02,-1.67d+01,-1.65d+03,4.21d+00,3.05d+03/
      data (xsub(12,i),i=1,6)                                           &
     & /2.60d-02,5.89d+02,-1.83d+01,-1.33d+03,3.12d+00,3.26d+03/
      data (xsub(13,i),i=1,6)                                           &
     & /3.36d-02,6.27d+02,-2.29d+01,-1.40d+03,3.26d+00,3.95d+03/
      data (xsub(14,i),i=1,6)                                           &
     & /3.84d-02,6.74d+02,-2.71d+01,-1.40d+03,3.21d+00,4.83d+03/
      data (xsub(15,i),i=1,6)                                           &
     & /2.25d-02,6.16d+02,-1.66d+01,-1.37d+03,3.22d+00,3.05d+03/
      data (xsub(16,i),i=1,6)                                           &
     & /1.88d-02,4.75d+02,-1.32d+01,-1.29d+03,3.31d+00,2.33d+03/
      data (xsub(17,i),i=1,6)                                           &
     & /2.19d-02,5.40d+02,-1.53d+01,-1.50d+03,3.87d+00,2.70d+03/
      data (xsub(18,i),i=1,6)                                           &
     & /1.98d-02,5.00d+02,-1.41d+01,-1.40d+03,3.67d+00,2.52d+03/

      data (ytop(1,i),i=1,6)                                            &
     & /1.89d-04,-3.91d-01,-8.65d-02,2.11d+01,-8.61d-02,7.12d+00/
      data (ytop(2,i),i=1,6)                                            &
     & /5.02d-05,-7.56d-02,-1.01d-02,1.80d+01,-7.27d-02,-2.81d+00/
      data (ytop(3,i),i=1,6)                                            &
     & /5.36d-06,3.54d-01,1.20d-02,1.73d+01,-7.09d-02,-5.51d+00/
      data (ytop(4,i),i=1,6)                                            &
     & /-2.67d-05,4.87d-02,3.03d-02,1.76d+01,-7.02d-02,-7.90d+00/
      data (ytop(5,i),i=1,6)                                            &
     & /2.60d-04,-1.74d+00,-1.16d-01,2.28d+01,-8.99d-02,9.67d+00/
      data (ytop(6,i),i=1,6)                                            &
     & /-1.23d-04,-6.11d-02,8.42d-02,1.66d+01,-6.63d-02,-1.51d+01/
      data (ytop(7,i),i=1,6)                                            &
     & /-1.02d-04,8.54d-03,7.12d-02,1.76d+01,-7.06d-02,-1.33d+01/
      data (ytop(8,i),i=1,6)                                            &
     & /-8.09d-05,2.23d-01,5.65d-02,1.62d+01,-6.61d-02,-1.08d+01/
      data (ytop(9,i),i=1,6)                                            &
     & /-9.68d-05,1.56d-01,6.56d-02,1.67d+01,-6.76d-02,-1.21d+01/
      data (ytop(10,i),i=1,6)                                           &
     & /-1.37d-04,-4.78d-02,8.94d-02,1.81d+01,-7.09d-02,-1.55d+01/
      data (ytop(11,i),i=1,6)                                           &
     & /-1.71d-04,-2.25d-01,1.08d-01,1.68d+01,-6.46d-02,-1.78d+01/
      data (ytop(12,i),i=1,6)                                           &
     & /-1.36d-04,-5.72d-01,8.92d-02,1.02d+01,-3.91d-02,-1.54d+01/
      data (ytop(13,i),i=1,6)                                           &
     & /-1.19d-04,-8.68d-03,8.28d-02,1.77d+01,-6.97d-02,-1.54d+01/
      data (ytop(14,i),i=1,6)                                           &
     & /-9.91d-05,-9.32d-01,7.13d-02,1.35d+01,-5.10d-02,-1.35d+01/
      data (ytop(15,i),i=1,6)                                           &
     & /-7.90d-05,2.75d-02,5.85d-02,1.55d+01,-6.20d-02,-1.16d+01/
      data (ytop(16,i),i=1,6)                                           &
     & /-7.91d-05,2.32d-01,5.71d-02,1.10d+01,-4.40d-02,-1.12d+01/
      data (ytop(17,i),i=1,6)                                           &
     & /-2.21d-04,-5.23d-01,1.39d-01,1.44d+01,-5.42d-02,-2.27d+01/
      data (ytop(18,i),i=1,6)                                           &
     & /-1.16d-04,-3.60d-01,7.85d-02,1.05d+01,-3.96d-02,-1.41d+01/

      data (ysub(1,i),i=1,6)                                            &
     & /4.17d-04,-1.41d+00,-2.06d-01,2.42d+01,-9.84d-02,2.20d+01/
      data (ysub(2,i),i=1,6)                                            &
     & /2.56d-04,-4.40d-01,-1.22d-01,2.19d+01,-8.92d-02,1.16d+01/
      data (ysub(3,i),i=1,6)                                            &
     & /1.93d-04,-2.84d-01,-8.88d-02,1.96d+01,-8.08d-02,7.68d+00/
      data (ysub(4,i),i=1,6)                                            &
     & /5.69d-05,-7.44d-02,-1.48d-02,1.80d+01,-7.36d-02,-2.01d+00/
      data (ysub(5,i),i=1,6)                                            &
     & /5.21d-04,-1.81d+00,-2.61d-01,1.96d+01,-7.85d-02,3.01d+01/
      data (ysub(6,i),i=1,6)                                            &
     & /-1.10d-04,6.88d-02,7.90d-02,1.77d+01,-7.11d-02,-1.48d+01/
      data (ysub(7,i),i=1,6)                                            &
     & /-1.76d-04,-3.81d-01,1.11d-01,1.79d+01,-6.90d-02,-1.85d+01/
      data (ysub(8,i),i=1,6)                                            &
     & /-2.38d-05,5.72d-01,2.47d-02,1.64d+01,-6.82d-02,-6.50d+00/
      data (ysub(9,i),i=1,6)                                            &
     & /-4.68d-05,4.19d-01,3.84d-02,1.74d+01,-7.16d-02,-8.59d+00/
      data (ysub(10,i),i=1,6)                                           &
     & /-1.10d-04,1.03d-01,7.32d-02,1.69d+01,-6.77d-02,-1.32d+01/
      data (ysub(11,i),i=1,6)                                           &
     & /-1.48d-04,-3.36d-01,9.58d-02,1.78d+01,-6.94d-02,-1.64d+01/
      data (ysub(12,i),i=1,6)                                           &
     & /-1.58d-04,1.36d-01,9.92d-02,1.56d+01,-6.18d-02,-1.65d+01/
      data (ysub(13,i),i=1,6)                                           &
     & /-2.69d-04,-6.12d-01,1.66d-01,1.27d+01,-4.80d-02,-2.65d+01/
      data (ysub(14,i),i=1,6)                                           &
     & /-6.34d-05,-4.27d-01,5.05d-02,1.73d+01,-6.82d-02,-1.06d+01/
      data (ysub(15,i),i=1,6)                                           &
     & /3.08d-05,5.25d-02,-6.04d-03,8.44d+00,-3.54d-02,-2.00d+00/
      data (ysub(16,i),i=1,6)                                           &
     & /-6.20d-05,3.40d-01,4.76d-02,9.54d+00,-3.85d-02,-9.98d+00/
      data (ysub(17,i),i=1,6)                                           &
     & /-2.90d-05,3.99d-01,2.75d-02,8.07d+00,-3.30d-02,-6.73d+00/
      data (ysub(18,i),i=1,6)                                           &
     & /-6.76d-05,4.41d-01,5.01d-02,1.50d+01,-6.10d-02,-1.01d+01/

      if (SwTopSub .eq. 1) then
        do i = 1,6
          OxygenSlope(i)     = xtop(NrStaring,i)
          OxygenIntercept(i) = ytop(NrStaring,i)
        enddo
      else
        do i = 1,6
          OxygenSlope(i)     = xsub(NrStaring,i)
          OxygenIntercept(i) = ysub(NrStaring,i)
        enddo
      endif

      return
      end subroutine oxygen_dat

! ----------------------------------------------------------------------
      subroutine OxygenReproFunction (OxygenSlope,OxygenIntercept,      &
     &   theta,thetas,tsoil,node,z,dz,rwu_factor,state)
! ----------------------------------------------------------------------
!     date               : January 2010
!     purpose            : Calculate oxygen stress according to reproduction function
! [GR-BH C7] state added for zbotcp->state%mesh%zbotcp migration.
! ----------------------------------------------------------------------
      use swap_array_dimensions, only: macp
      implicit none

! --- global
      integer node,i
      real(8) OxygenSlope(6),OxygenIntercept(6),theta(macp),thetas(macp)
      real(8) tsoil(macp),z(macp),dz(macp)
      type(swap_state_t), intent(in) :: state  ! [GR-BH C7]

! --- local
      real(8) intercept,slope,sum_porosity
      real(8) gas_filled_porosity
      real(8) soil_temp,depth_ss,mean_gas_filled_porosity
      real(8) rwu_factor

      gas_filled_porosity = thetas(node) - theta(node)
      soil_temp = tsoil(node) + 273.d0
      depth_ss = -z(node) * 0.01d0

      if (gas_filled_porosity .lt. 1.d-10) then
         rwu_factor = 0.d0
         return
      endif

! --- mean gas filled porosity
      sum_porosity = 0.0d0
      do i = 1,node
        sum_porosity = sum_porosity +                                   &
     &                 (thetas(i) - theta(i)) * dz(i)
      enddo
      mean_gas_filled_porosity = sum_porosity /(-state%mesh%zbotcp(node))  ! [GR-BH C7]

      intercept = OxygenIntercept(1)*soil_temp**2 +                     &
     &            OxygenIntercept(2)*depth_ss**2 +                      &
     &            OxygenIntercept(3)*soil_temp +                        &
     &            OxygenIntercept(4)*depth_ss +                         &
     &            OxygenIntercept(5)*soil_temp*depth_ss +               &
     &            OxygenIntercept(6)

      slope = OxygenSlope(1)*soil_temp**2 +                             &
     &        OxygenSlope(2)*depth_ss**2 +                              &
     &        OxygenSlope(3)*soil_temp +                                &
     &        OxygenSlope(4)*depth_ss +                                 &
     &        OxygenSlope(5)*soil_temp*depth_ss +                       &
     &        OxygenSlope(6)

! --- Calculate the sink term (Root Water Uptake) variable due to oxygen stress.
      rwu_factor = intercept + slope*mean_gas_filled_porosity
      if (rwu_factor .gt. 1.d0) then
         rwu_factor = 1.d0
      endif
      if (rwu_factor .lt. 0.d0) then
         rwu_factor = 0.d0
      endif
      return

      end subroutine OxygenReproFunction

end module oxygenrepro_dormant_mod
