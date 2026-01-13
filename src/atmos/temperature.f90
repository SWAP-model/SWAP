! File VersionID:
!   $Id: temperature.f90 362 2018-01-08 13:08:33Z kroes006 $
! ----------------------------------------------------------------------
      subroutine temperature(task)
! ----------------------------------------------------------------------
!     date               : november 2004
!     purpose            : calculate soil temperatures
! ----------------------------------------------------------------------
      use variables
      implicit none 

! --- global
      integer task
! --- local
      integer i,lay, ierror
      real(8) tmpold(macp),tab(mabbc*2),dummy,afgen,gmineral
      real(8) thoma(macp),thomb(macp),thomc(macp),thomf(macp)
      real(8) theave(macp),heacnd(macp)
      real(8) heaconbot,qhbot
      real(8) apar, dzsnw, heaconsnw, Rosnw

      character(len=200) messag

      save
! ----------------------------------------------------------------------

      select case (task)
      case (1)

! === initialization ===================================================

! --- determine initial temperature profile

      if (swcalt.eq.1) then
! ---   analytical solution ---
        do i = 1,numnod
          tsoil(i) = tmean+tampli*(dsin(0.0172d0*(daynr-timref+91.0d0)+ &
     &              z(i)/ddamp)) / dexp(-z(i)/ddamp)
        enddo
      else
! ---   numerical solution, use specified soil temperatures ---
        if (swinco.ne.3) then
          do i = 1, nheat
            tab(i*2) = tsoil(i)
            tab(i*2-1) = dabs(zh(i))
          end do
          do i = 1, numnod
            tsoil(i) = afgen(tab,macp*2,dabs(z(i)))
          end do
        end if
      endif

      if (swcalt.eq.2) then
! ---   initialize dry bulk density and volume fractions sand, clay and organic matter
        do i = 1, numnod
          lay = layer(i)
          dummy = orgmat(lay)/(1.0d0 - orgmat(lay))
          gmineral = (1.0d0 - thetas(i)) / (0.370d0 + 0.714d0*dummy)
          fquartz(i) = (psand(lay) + psilt(lay))*gmineral/2.7d0
          fclay(i) = pclay(lay)*gmineral/2.7d0
          forg(i) = dummy*gmineral/1.4d0
        end do

      endif

      return

      case (2)

! === soil temperature rate and state variables ========================

      if (swcalt .eq. 2) then
!   - numerical solution

! ---   use specified soil surface temperatures as top boundary condition
        if (swtopbhea .eq. 2) then
! ---   use specified soil surface temperatures as top boundary condition
           TeTop = afgen (temtoptab,2*mabbc,t1900+dt)
        elseif (dabs(ssnow).gt.1.0d-10) then
! --- air temperature can not be used with a snow layer, calculate 
! --- temperature on soil- snow interface
           Rosnw = 170.0d0
           heaconsnw = 2.86d-6 * 864.0d0 * Rosnw**2.d0
           dzsnw = ssnow / 0.170d0
           if (heacon(1).lt.1.d-10) heacon(1) = 100.0d0
           apar = (0.5d0*heaconsnw*dz(1)) / (heacon(1)*dzsnw)
           if (flmetdetail) then
             TeTop = (Tsoil(1) + apar*atav(wrecord)) / (1.d0+apar)
           else
             TeTop = (Tsoil(1) + apar*Tav) / (1.d0+apar)
           endif
        else
           if (flmetdetail) then
             TeTop = atav(wrecord)
           else
             TeTop = Tav
           endif
        endif
!
        if (SwBotbHea.eq.1) then
! ---   no heat flow through bottom of profile assumed
           TeBot = Tsoil(Numnod)
        elseif (SwBotbHea.eq.2) then
! ---   bottom temperature is prescribed
           TeBot = afgen (tembtab,2*mabbc,t1900+dt)
        endif

! --- save old temperature profile
        do i = 1,numnod
          tmpold(i) = tsoil(i)
        enddo

! --- compute heat conductivity and capacity ---------------------------

        do i = 1,numnod
          theave(i) = 0.5d0 * (theta(i) + thetm1(i))
        enddo

! --- calculate nodal heat capacity and thermal conductivity 
        call devries(theave,heacap,heacnd)
        heacon(1) = heacnd(1)
        do i = 2,numnod
          heacon(i) = 0.5d0 * (heacnd(i) + heacnd(i-1))
        enddo

! ---   calculate new temperature profile --------------------------------

! ---   calculation of coefficients for node = 1
        i = 1
! ---   temperature fixed at soil surface
        thoma(i) = - dt * heacon(i) / (dz(i) * disnod(i))
        thomc(i) = - dt * heacon(i+1) / (dz(i) * disnod(i+1))
        thomb(i) = heacap(i) - thoma(i) - thomc(i)
        thomf(i) = heacap(i) * tmpold(i) - thoma(i) * TeTop

! ---   calculation of coefficients for 2 < node < numnod
        do i = 2,numnod-1 
          thoma(i) = - dt * heacon(i) / (dz(i) * disnod(i))
          thomc(i) = - dt * heacon(i+1) / (dz(i) * disnod(i+1))
          thomb(i) = heacap(i) - thoma(i) - thomc(i)
          thomf(i) = heacap(i) * tmpold(i)
        enddo

! ---   calculation of coefficients for node = numnod
        i = numnod
        if (SwBotbHea.eq.1) then
! ---   no heat flow through bottom of profile assumed
           qhbot = 0.0d0
           thoma(i) = - dt * heacon(i) / (dz(i) * disnod(i))
           thomb(i) = heacap(i) - thoma(i)
           thomf(i) = heacap(i) * tmpold(i) - (qhbot * dt)/dz(i)
        elseif (SwBotbHea.eq.2) then
! ---   bottom temperature is prescribed
           heaconBot = heacnd(i)
           thoma(i)  = - dt * heacon(i) / (dz(i) * disnod(i))
           thomc(i)  = - dt * heaconBot / (dz(i) * 0.5d0 * dz(i))
           thomb(i)  = heacap(i) - thoma(i) - thomc(i)
           thomf(i)  = heacap(i) * tmpold(i) - thomc(i) * TeBot
        endif

! ---   solve for vector tsoil a tridiagonal linear set
        call tridag (numnod, thoma, thomb, thomc, thomf, tsoil,ierror)
        if(ierror.ne.0)then
           messag =                                                     &
     &       'During a call from Temperature an error occured in TriDag'
           call fatalerr ('Temperature',messag)
        end if
      else

! ---   analytical solution temperature profile ---
        do i = 1,numnod
          tsoil(i) = tmean+tampli*(dsin(0.0172d0*(daynr-timref+91.0d0)+ &
     &              z(i)/ddamp)) / dexp(-z(i)/ddamp)
        enddo

      endif

      case default
         call fatalerr ('Temperature', 'Illegal value for TASK')
      end select

      return
      end

!***********************************************************************
      Subroutine Devries (theta,HeaCap,HeaCon)
!***********************************************************************
!* Purpose:    Calculate soil heat capacity and conductivity for each  *
!*             compartment by full de Vries model                      *     
!* References:                                                         *
!* Description:                                                        *
!* de Vries model for soil heat capcity and thermal conductivity.      * 
!* Heat capacity is calculated as average of heat capacities for each  *
!* soil component. Thermal conductivity is calculated as weighted      *
!* average of conductivities for each component. If theta > 0.05 liquid*
!* water is assumed to be the main transport medium in calculating the *
!* weights. If theta < 0.02 air is assumed to be the main transport    *
!* medium (there is also an empirical adjustment to the conductivity). *
!* For 0.02 < theta < 0.05 conductivity is interpolated.               *
!* See: Heat and water transfer at the bare soil surface, H.F.M Ten    *
!* Berge (pp 48-54 and Appendix 2)                                     *
!***********************************************************************
!* Input:                                                              *
!* NumNod - number of compartments (-)                                 *
!* theta/THETAS - volumetric soil moisture/ saturated vol. s. moist (-)*
!* Fquartz, Fclay and Forg - volume fractions of sand, clay and org.ma.*
!* Output:                                                             *
!* HeaCap - heat capacity (J/m3/K)                                     *
!* HeaCon - thermal conductivity (W/m/K)                               *
!***********************************************************************
      use variables, only: NumNod,THETAS,FQUARTZ,FCLAY,FORG
      Implicit None
      Include 'arrays.fi'
  
!     (i) Global declarations                                          
!     (i.i) Input
      real(8) theta(macp)  ! NOTE: this not the same as theta in VARIABLES, since on input it is some average water content
!     (i.ii)
      real(8) HeaCap(MACP), HeaCon(MACP)
!     (ii) Local declarations                                          
      Integer Node
      real(8) kaw
      real(8), parameter :: kaa = 1.0d0
      real(8), parameter :: kww = 1.0d0
      real(8) fAir(MACP)
      real(8) HeaConDry,HeaConWet
!     (iii) Parameter declarations                                     
!     (iii.i) Physical constants                                       
!     Specific heats (J/kg/K)
      real(8), parameter :: cQuartz =  800.0d0
      real(8), parameter :: cClay   =  900.0d0
      real(8), parameter :: cWat    = 4180.0d0
      real(8), parameter :: cAir    = 1010.0d0
      real(8), parameter :: cOrg    = 1920.0d0
!     Density (kg/m3)
      real(8), parameter :: dQuartz = 2660.0d0
      real(8), parameter :: dClay   = 2650.0d0
      real(8), parameter :: dWat    = 1000.0d0
      real(8), parameter :: dAir    =    1.2d0
      real(8), parameter :: dOrg    = 1300.0d0
!     Thermal conductivities (W/m/K)
      real(8), parameter :: kQuartz = 8.8d0
      real(8), parameter :: kClay   = 2.92d0
      real(8), parameter :: kWat    = 0.57d0
      real(8), parameter :: kAir    = 0.025d0
      real(8), parameter :: kOrg    = 0.25d0
!     
      real(8) GAir,GAirdry
      real(8), parameter :: GQuartz = 0.14d0
      real(8), parameter :: GClay   = 0.125d0
      real(8), parameter :: GWat    = 0.14d0
      real(8), parameter :: GOrg    = 0.5d0
!     (iii.ii) theta 0.02, 0.05 and 1.00
      real(8), parameter :: thetaDry = 0.02d0
      real(8), parameter :: thetaWet = 0.05d0
! ----------------------------------------------------------------------
!     (0) Weights for each component in conductivity calculations      
!     (Calculate these but define as parameters in later version)
      real(8), parameter :: kqw = 0.66d0 / (1.0d0 + ((kQuartz/kWat)-1.0d0) * GQuartz) + 0.33d0 /  &
     &                            (1.0d0 + ((kQuartz/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GQuartz))
      real(8), parameter :: kcw = 0.66d0 / (1.0d0 + ((kClay/kWat) - 1.0d0) * GClay) + 0.33d0 /    &
     &                            (1.0d0 + ((kClay/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GClay))
      real(8), parameter :: kow = 0.66d0 / (1.0d0 + ((kOrg/kWat) - 1.0d0) * GOrg) + 0.33d0 /      &
     &                            (1.0d0 + ((kOrg/kWat) - 1.0d0) * (1.0d0 - 2.0d0 * GOrg))
      real(8), parameter :: kwa = 0.66d0 / (1.0d0 + ((kWat/kAir) - 1.0d0) * GWat) + 0.33d0 /      &
     &                            (1.0d0 + ((kWat/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GWat))
      real(8), parameter :: kqa = 0.66d0 / (1.0d0 + ((kQuartz/kAir)-1.0d0) * GQuartz) + 0.33d0 /  &
     &                            (1.0d0 + ((kQuartz/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GQuartz))
      real(8), parameter :: kca = 0.66d0 / (1.0d0 + ((kClay/kAir) - 1.0d0) * GClay) + 0.33d0 /    &
     &                            (1.0d0 + ((kClay/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GClay))
      real(8), parameter :: koa = 0.66d0 / (1.0d0 + ((kOrg/kAir) - 1.0d0) * GOrg) + 0.33d0 /      &
     &                            (1.0d0 + ((kOrg/kAir) - 1.0d0) * (1.0d0 - 2.0d0 * GOrg))

!     additional constants
      real(8), parameter :: kAirDIVkWat = kAir/kWat
      real(8), parameter :: cdQuartz    = dQuartz*cQuartz
      real(8), parameter :: cdClay      = dClay*cClay
      real(8), parameter :: cdWat       = dWat*cWat
      real(8), parameter :: cdAir       = dAir*cAir
      real(8), parameter :: cdOrg       = dOrg*cOrg
      real(8), parameter :: kqaXkQuartz = kqa*kQuartz
      real(8), parameter :: kcaXkClay   = kca*kClay
      real(8), parameter :: kaaXkAir    = kaa*kAir
      real(8), parameter :: koaXkOrg    = koa*kOrg
      real(8), parameter :: kwaXkWat    = kwa*kWat
      real(8), parameter :: kqwXkQuartz = kqw*kQuartz
      real(8), parameter :: kcwXkClay   = kcw*kClay
      real(8), parameter :: kowXkOrg    = kow*kOrg
      real(8), parameter :: kwwXkWat    = kww*kWat
      
      Do Node = 1,NumNod

!        (1) Air fraction and related parameters
         fAir(Node) = THETAS(Node) - theta(Node)

!        Determine shape factor of air         
         if (theta(node) .gt. thetadry) then
           GAir = 0.333d0 - fair(node)/thetas(node)*0.298d0
         else
           GAirdry = 0.333d0 - fair(node)/thetas(node)*0.298d0
           GAir = 0.013d0 + theta(node)/thetaDry*(GAirdry - 0.013d0)
         endif

!        Determine weighting factor air - water
         kaw = 0.66d0 / (1.0d0 + ((kAirDIVkWat) - 1.0d0) * GAir) + 0.33d0/&
     &       (1.0d0 + ((kAirDIVkWat) - 1.0d0) * (1.0d0 - 2.0d0 * GAir)) 

!        (2) Heat capacity (W/m3/K) is average of heat capacities for  
!        all components (multiplied by density for correct units)      
         HeaCap(Node) = fQuartz(Node)*cdQuartz + fClay(Node)*cdClay +   &
     &                  theta(Node)*cdWat + fAir(Node)*cdAir + fOrg(Node)*cdOrg

!        (3) Thermal conductivity (W/m/K) is weighted average of       
!        conductivities of all components                              
!        (3.1) Dry conditions (include empirical correction) (eq. 3.44)
     
         If (theta(Node).LE.thetaDry) Then
            HeaCon(Node) = 1.25d0 *                                     &
                           (fQuartz(Node)*kqaXkQuartz +                 &
     &                      fClay(Node)*kcaXkClay +                     &
     &                      fAir(Node)*kaaXkAir +                       &
     &                      fOrg(Node)*koaXkOrg +                       &
     &                      theta(Node)*kwaXkWat) /                     &
     &   (kqa * fQuartz(Node) + kca * fClay(Node) + kaa * fAir(Node) +  &
     &    koa * fOrg(Node) + kwa * theta(Node))

!        (3.2) Wet conditions  (eq. 3.43)                              
         Else If (theta(Node).GE.thetaWet) Then
            HeaCon(Node) = (fQuartz(Node)*kqwXkQuartz +                 &
     &                      fClay(Node)*kcwXkClay +                     &
     &                      fAir(Node)*kaw*kAir +                       &
     &                      fOrg(Node)*kowXkOrg +                       &
     &                      theta(Node)*kwwXkWat) /                     &
     &   (kqw * fQuartz(Node) + kcw * fClay(Node) + kaw * fAir(Node) +  &
     &    kow * fOrg(Node) + kww * theta(Node))

!        (3.3) dry < theta < wet (interpolate)
         Else
!           (3.3.1) Conductivity for theta = 0.02                      
            HeaConDry = 1.25d0 *                                        &
                           (fQuartz(Node)*kqaXkQuartz +                 &
     &                      fClay(Node)*kcaXkClay +                     &
     &                      fAir(Node)*kaaXkAir +                       &
     &                      fOrg(Node)*koaXkOrg +                       &
     &                      thetaDry*kwaXkWat) /                        &
     &   (kqa * fQuartz(Node) + kca * fClay(Node) + kaa * fAir(Node) +  &
     &    koa * fOrg(Node) + kwa * thetaDry)
!           (3.3.1) Conductivity for theta = 0.05                      
            HeaConWet = (fQuartz(Node)*kqwXkQuartz +                    &
     &                      fClay(Node)*kcwXkClay +                     &
     &                      fAir(Node)*kaw*kAir +                       &
     &                      fOrg(Node)*kowXkOrg +                       &
     &                      thetaWet*kwwXkWat) /                        &
     &   (kqw * fQuartz(Node) + kcw * fClay(Node) + kaw * fAir(Node) +  &
     &    kow * fOrg(Node) + kww * thetaWet)
!         (3.3.3) Interpolate                                          
          HeaCon(Node) = HeaConDry + (theta(Node)-thetaDry) *           &
     &                   (HeaConWet - HeaConDry)                        &
     &                 / (thetaWet - thetaDry)      
         End If

! ---    conversion of capacity from J/m3/K to J/cm3/K
         HEACAP(NODE) = HEACAP(NODE)*1.0d-6

! ---    conversion of conductivity from W/m/K to J/cm/K/d
         HEACON(NODE) = HEACON(NODE)*864.0d0

      End Do

      Return
      End 

