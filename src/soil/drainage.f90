! File VersionID:
!   $Id: drainage.f90 368 2018-01-11 15:44:15Z heine003 $
! ----------------------------------------------------------------------
      subroutine bocodrb (dh)

! ----------------------------------------------------------------------
!     update             : Dec 2014
!     date               : July 2002
!     purpose            : calculates lateral drainage fluxes
! ----------------------------------------------------------------------  
      use variables, only: dramet,gwl,zbotdr,basegw,l,qdrain,ipos,khtop,khbot,kvtop,kvbot,entres,wetper,zintf,geofac,swdtyp,      &
                           owltab,t1900,swallo,drares,infres,qdrtab,nrlevs,swnrsrf,cofintfl,expintfl,dt,shape,FlMacropore,NumLevRapDra,ZDraBas,swliminf,nowltab
      implicit none
      !include 'arrays.fi'

! --- global
      real(8) dh

! ----------------------------------------------------------------------
! --- local
      real(8) afgen
      integer i,lev

      real(8)  zimp,dbot,pi,totres,x,fx,eqd,rver,rhor,rrad
      real(8)  gwldra   !, temptab(2*maowl)

      logical fldry

      parameter (pi=3.14159d0)

      character(len=200) messag
      include  'params.fi'
! ----------------------------------------------------------------------

      gwldra = gwl

! --- drainage flux calculated according to hooghoudt or ernst
      if (dramet.eq.2) then
        if (shape .gt. small) dh = (gwldra-zbotdr(1)) / shape

! --- contributing layer below drains limited to 1/4 l
        zimp = max (basegw,zbotdr(1)-0.25*l(1))
        dbot = (zbotdr(1)-zimp)  
        if (dbot.lt.0.0d0) then
          messag = 'At the drainage section, the level of the'          &
     &       //' impervious layer is higher than the level of the'      &
     &       //' drain bottom. Adapt drain input!'
          call fatalerr ('Bocodrb',messag)
        endif
              
! --- no infiltration allowed
        if (dh.lt.1.0d-10) then
          qdrain(1) = 0.0d0
          return
        endif

! --- case 1: homogeneous, on top of impervious layer
        if (ipos.eq.1) then

! --- calculation of drainage resistance and drainage flux
          totres = l(1)*l(1)/(4*khtop*abs(dh)) + entres
          qdrain(1) = dh/totres 

! --- case 2,3: in homogeneous profile or at interface of 2 layers
        elseif (ipos.eq.2.or.ipos.eq.3) then

! --- calculation of equivalent depth
          x = 2*pi*dbot/l(1)
          if (x.gt.0.5d0) then
            fx = 0.0d0
            do 10 i = 1,5,2
              fx = fx + (4*exp(-2*i*x))/(i*(1.0d0-exp(-2*i*x)))
   10       continue
            eqd = pi*l(1)/8 / (log(l(1)/wetper(1))+fx)
          else
            if (x.lt.1.0d-6) then
              eqd = dbot
            else
              fx = pi**2/(4*x) + log(x/(2*pi))
              eqd = pi*l(1)/8 / (log(l(1)/wetper(1))+fx)
            endif
          endif 
          if (eqd.gt.dbot) eqd = dbot

! --- calculation of drainage resistance & drainage flux
          if (ipos.eq.2) then
            totres = l(1)*l(1)/(8*khtop*eqd+4*khtop*abs(dh)) + entres 
          elseif (ipos.eq.3) then
            totres = l(1)*l(1)/(8*khbot*eqd+4*khtop*abs(dh)) + entres
          endif
          qdrain(1) = dh/totres 

! --- case 4: drain in bottom layer
        elseif (ipos.eq.4) then
          if (zbotdr(1).gt.zintf) then
            messag = 'At the drainage section, the level of the'        &
     &         //' impervious layer is higher than the level of the'    &
     &         //' drain bottom. Adapt drain input!'
            call fatalerr ('bocodrb',messag)
          endif 
            rver = max (gwldra-zintf,0.0d0) / kvtop +                   &
     &                (min (zintf,gwldra) -zbotdr(1))/kvbot
          rhor = l(1)*l(1)/(8*khbot*dbot) 
          rrad = l(1)/(pi*dsqrt(khbot*kvbot)) * log(dbot/wetper(1))
          totres = rver+rhor+rrad+entres
          qdrain(1) = dh/totres 

! --- case 5 : drain in top layer
        elseif (ipos.eq.5) then
          if (zbotdr(1).lt.zintf) then
            messag = 'At the drainage section, the level of the'        &
     &         //' impervious layer is higher than the level of the'    &
     &         //' drain bottom. Adapt drain input!'
            call fatalerr ('bocodrb',messag)
          endif
          rver = (gwldra-zbotdr(1))/kvtop
          rhor = l(1)*l(1) / (8*khtop*(zbotdr(1)-zintf) +               &
     &                             8*khbot*(zintf-zimp))
          rrad = l(1)/(pi*dsqrt(khtop*kvtop))*log((geofac*               &
     &                (zbotdr(1)-zintf))/wetper(1))
          totres = rver+rhor+rrad+entres
          qdrain(1) = dh/totres
        endif

! --- drainage flux calc. using given drainage/infiltration resistance
      elseif (dramet.eq.3) then

        do lev = 1,nrlevs
          fldry = .false.

! ---     first copy to 1-dimensional table temptab
!          do i = 1,2*nowltab(lev)
!          do i = 1,2*maowl 
!            temptab(i) = owltab(lev,i)
!          end do
!          x = afgen(temptab,2*nowltab(lev),t1900+dt-1.d0)
!          x = afgen(temptab,2*maowl,t1900+dt-1.d0)
          x = afgen(owltab(lev,1:2*nowltab(lev)),2*nowltab(lev),t1900+dt-1.d0)
!          x = afgen(owltab(lev,1:2*maowl),2*maowl,t1900+dt-1.d0)
          if ( (x-zbotdr(lev)) .lt. 1.0d-3) then
            fldry = .true.
          endif
          dh = gwldra- x
          if (fldry) dh = gwldra-zbotdr(lev)
! ---     drainage basis for rapid drainage through macropores
          if (FlMacropore .and. lev.eq.NumLevRapdra .and. swdtyp(lev).eq.2) then
            ZDraBas= dmax1(x,zbotdr(NumLevRapDra))
          endif

! ---     drainage              
          if (dh.ge.0.0d0) then

! ---       interflow flux calculated by a power function
            if ((lev.eq.nrlevs).and.(swnrsrf.eq.1)) then
              qdrain(lev) = cofintfl*dh**expintfl              
            else
              qdrain(lev) = dh/drares(lev)
              if (swallo(lev).eq.2 ) qdrain(lev) = 0.0d0
            endif

! ---     infiltration
          else
! Pvw_begin , implemented by KRO_20170907
!           Limit the infiltration head (-dh, dh<0. here) to the waterdepth in the channel
            if(swdtyp(lev).eq.2) then
              if(swliminf.eq.1) then
                dh = max(dh,(zbotdr(lev)-x))
              endif
            endif
! Pvw_end   , implemented by KRO_20170907
            qdrain(lev) = dh/infres(lev)
            if (swallo(lev).eq.3.or.fldry ) qdrain(lev) = 0.0d0
          endif
        end do


! --- drainage flux from table with gwlevel - flux data pairs
      elseif (dramet.eq.1) then
        qdrain(1) = afgen (qdrtab,50,abs(gwldra))
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine drainage 
! ----------------------------------------------------------------------
!     Date               : Aug 2004   
!     Purpose            : calculate drainage rate and state variables 
! ----------------------------------------------------------------------

      use variables
      implicit none

!     local
      integer node,level
      real(8) zCum,difzTopDisLay(madr),ratio,ratiodz,sumqdr(madr),dh
      real(8) afgen   !, temptab(2*maowl)
      integer nodeTopDisLay(madr)
      CHARACTER(len=33) messag

!   - In case of macropores: initialise drainage basis for rapid drainage through macropores
      if (flInitDraBas) then
         if (NumLevRapDra.gt.nrlevs) then
            messag = ' NUMLEVRAPDRA greater then NRLEVS'
            call fatalerr('MacroRead',messag)
         endif
!
         if (dramet.lt.3) then
            ZDraBas = zbotdr(1)
         else
            if (swdtyp(NumLevRapDra).eq.1) then
               ZDraBas = zbotdr(NumLevRapDra)
            else
               !do i = 1,2*maowl 
               !   temptab(i) = owltab(NumLevRapDra,i)
               !end do
               !ZDraBas = afgen (temptab,2*maowl,t1900)
               ZDraBas = afgen(owltab(NumLevRapDra,1:2*nowltab(NumLevRapDra)),2*nowltab(NumLevRapDra),t1900)
            endif
         endif
!
         flInitDraBas = .false.
!
         Return
!
      endif

! --- reset intermediate soil water fluxes
      if (flzerointr) then
        do node = 1,numnod
          do level = 1,nrlevs
            inqdra(level,node)     = 0.0d0
            inqdra_in(level,node)  = 0.0d0
            inqdra_out(level,node) = 0.0d0
          enddo
        enddo
        iqdra = 0.0d0
      endif

! --- reset cumulative soil water fluxes
      if (flzerocumu) then
        cqdra = 0.0d0
        do level = 1,nrlevs
          cqdrain(level) = 0.0d0
          cqdrainin(level) = 0.0d0
          cqdrainout(level) = 0.0d0
        enddo
      endif

! --- reset to zero if groundwater level under soil profile and return
      if (gwl.gt.998.0d0) then
        do level = 1,nrlevs
          qdrain(level) = 0.0d0
        end do
        return
      endif

! --- calculate total drainage rate and state variables
      call bocodrb (dh)

! --- partition drainage flux over compartments
      if (swdivd.eq.1) then
         call divdra (numnod,nrlevs,dz,ksatfit,ksatexm,fluseksatexm,    &
     &      layer,cofani,gwl,l,qdrain,qdra,Swdivdinf,Swnrsrf,           &
     &      SwTopnrsrf,Zbotdr,dt,FacDpthInf,owltab,t1900)                  !  Divdra, infiltration
!       redistribute qdrain with new top boundary for discharge layers
        if(swdislay.eq.2) then
            do level=1,nrlevs
               if(swtopdislay(level).eq.1)  then
                  zTopDisLay(level) = fTopDisLay(level) * gwl  +        &
     &                       (1.0d0-fTopDisLay(level)) * (gwl-dh)
               end if
            end do
        end if
        if(swdislay.eq.1 .or. swdislay.eq.2) then
            do level=1,nrlevs
               if(swtopdislay(level).eq.1)  then
!                 find node nr of new top of discharge layer
                  nodeTopDisLay(level) = 1
                  zCum               = - dz(1)
                  do while (zTopDisLay(level) .lt. zCum)
                     nodeTopDisLay(level) = nodeTopDisLay(level) + 1
                     zCum              = zCum - dz(nodeTopDisLay(level))
                  enddo
!                 saturated part (difzTopDisLay(lev)) of compartment containing waterlevel
                  difzTopDisLay(level) = zTopDisLay(level) - zCum
                  ratiodz =                                             &
     &                     difzTopDisLay(level)/dz(nodeTopDisLay(level))
                  sumqdr(level) =                                       &
     &                        ratiodz * qdra(level,nodeTopDisLay(level))
                  do node = nodeTopDisLay(level)+1,numnod
                     sumqdr(level) =  sumqdr(level) + qdra(level,node)
                  end do
                  if( dabs(sumqdr(level)) .lt. 1.0d-8)then
                     ratio = 1.0d0
                  else
                     ratio = qdrain(level)/sumqdr(level)
                  end if
!                 redistribute drainwater fluxes
                  do node = 1,nodeTopDisLay(level)-1
                     qdra(level,node) = 0.0d0
                  end do           
                  qdra(level,nodeTopDisLay(level)) =                    &
     &                qdra(level,nodeTopDisLay(level)) * ratio * ratiodz
                  do node = nodeTopDisLay(level)+1,numnod
                     qdra(level,node) = qdra(level,node)* ratio
                  end do           
               endif
            end do           
         endif
      else
! --- drainage flux through lowest compartment 
        do level = 1,nrlevs
          do node = 1,numnod-1
            qdra(level,node) = 0.0d0
          end do
          qdra(level,numnod) = qdrain(level)
        end do
      endif

      qdrtot = 0.0d0
      do level=1,nrlevs
          qdrtot = qdrtot + qdrain(level)
      end do

      return
      end 

! ----------------------------------------------------------------------
      subroutine bocodre (dh)
! ----------------------------------------------------------------------
!     Date               : 5/5/2001
!     Purpose            :        
! --- Calculate drainage/infiltration fluxes for all levels: qdrain(level).
! --- Summate fluxes to secondary system: qdrd
! --- Given present storage (swst), qdrd and wscap(imper), check if the
! --- system falls dry. If so reduce drainage fluxes proportionally in 
! --- such a way that total drainage flux equals available amount 
! --- ( = swst + wscap).
!     Subroutines called : -                                           
!     Functions called   : -                                     
!     File usage         : -  Error handling
! ----------------------------------------------------------------------  
      use variables, only: swsec,swsrf,nrlevs,nrpri,gwl,zbotdr,taludr,widthr,pond,pondmx,swdtyp,dt,wls,wlp,drainl,l,rdrain,rinfi,      &
                           rentry,rexit,gwlinf,wetper,qdrain,qdrd,impend,nmper,wscap,swst,swnrsrf,rsurfdeep,rsurfshallow,cofintfl,          &
                           expintfl,t1900,FlMacropore,NumLevRapdra,ZDraBas
      IMPLICIT NONE

! --- global
      real(8) dh

! --- local
      integer level,imper
      real(8) qdrdm,qdratio,swdepth,swexbrd,dvmax,swstmax,wl,rd,re
      character(len=200) messag

      save
! ----------------------------------------------------------------------

! --- summate fluxes for use by swballev and swlevbal
      qdrd = 0.0d0

      do 500 level = 1,nrlevs

! --- surface water level
        if (swsrf.ge.2) then
          if (level.gt.nrpri) then
            wl = wls
          else
            wl = wlp
          endif
        endif

! --- drainage fluxes are set to zero if both groundwater level and surface 
!     water level are above ponding sill (so the nonzero drainage flux 
!     is only computed if either the gwl or the wl is below pondmx)
        if (wl .lt. pondmx .or. gwl .lt. pondmx) then

! --- channel is active medium if either groundwater or surface water
!     level is above channel bottom
          if (gwl.gt.(zbotdr(level)+0.001d0) .or.                       &
     &        wl.gt.(zbotdr(level)+0.001d0)) then
            if (wl .le. (zbotdr(level)+0.001d0).or.swsrf.eq.1) then

! --- only groundw. level above channel bottom; bottom is dr. base
              drainl(level) = zbotdr(level)

! --- wetted perimeter only computed for open channels
!     (for drains it is input) 
              if (swdtyp(level) .eq. 0) then
                wetper(level) = widthr(level)
              endif
            else

! --- surface water level above channel bottom
              drainl(level) = wl
              if (swdtyp(level) .eq. 0) then
                swdepth = wl-zbotdr(level)
                swexbrd = (wl-zbotdr(level))/taludr(level)
                wetper(level) = widthr(level) +                         &
     &            2*dsqrt(swdepth**2 + swexbrd**2) 
              endif
            endif

! --- drainage flux (cm/d)
            ! calculate head difference
            dh = gwl - drainl(level)
            if (gwl.gt.-0.1d0) dh = dh+pond
            if (dh.lt.0.0d0 .and.gwl.lt.gwlinf(level)) then
              dh = gwlinf(level)-drainl(level)
            endif
            ! interflow flux calculated by a power function, 
            if ((level.eq.nrlevs).and.(swnrsrf.eq.2)) then               
              qdrain(level) = cofintfl * dh ** expintfl                 
            else 
              if (dh.gt.0.0d0) then
                rd = rdrain(level)
                re = rentry(level)
! ---           surface drainage (vaccuum cleaner)
                if ((level.eq.nrlevs).and.(swnrsrf.eq.1)) then
                  rd = rsurfdeep - dh
                  rd = max(rd,rsurfshallow)
                endif
              else
                rd = rinfi(level)
                re = rexit(level)
              endif 
              if (swdtyp(level).eq.0) then
                qdrain(level) = dh / ((rd + re*l(level)/wetper(level)))
              else
                qdrain(level) = dh/rd
              endif 
            endif
          else
            if (swdtyp(level) .eq. 0) then
              wetper(level) = 0.0d0
            endif
            dh = 0.0d0
            qdrain(level) = 0.0d0

!   - for determining drainage basis for rapid drainage through macropores
            drainl(level) = zbotdr(level)

          endif
!
        else
          qdrain(level) = 0.0d0

!   - for determining drainage basis for rapid drainage through macropores
          drainl(level) = wl

        endif
!
        if (swsrf.ge.2.and.level.gt.nrpri) then

! --- qdrd is total flux to or from secondary system
          qdrd = qdrd + qdrain(level)
        endif

 500  continue

!   - drainage basis for rapid drainage through macropores
      if (FlMacropore) then
        if(swdtyp(NumLevRapDra).ne.1) then
          ZDraBas = drainl(NumLevRapDra)
        endif
      endif

! ----------------------------------------------------------------------
! --- check for system falling dry (only for swsec = 2):
      if (swsec .eq. 1) return
      if (swsrf.eq.1) then
        do 10 level = 1,nrlevs
          if (qdrain(level).lt.0.0d0) then
            qdrain(level) = 0.0d0
          endif
 10     continue
      elseif (swsrf.ge.2) then

! --- determine which management period the model is in:
        imper = 0
 800    imper = imper + 1

! ---   Error handling
        if (imper .gt. nmper) then
          messag = 'sw-management periods(IMPER), more than defined'
          call fatalerr ('Bocodre',messag)
        endif
        if (t1900-1.d0+0.1d-10 .gt. impend(imper)) goto 800

! ---   determine whether the system will become empty
        dvmax = (qdrd + wscap(imper)) * dt
        swstmax = swst + dvmax

        if (swstmax .lt. 0.0d0) then
! ---     storage decreases to below zero, then the surface water system 
!         falls dry; make the total infiltration exactly equal to the
!         available amount:
          qdrdm = - (swst + wscap(imper)*dt) / dt
          qdratio = qdrdm / qdrd

! ---     Error handling
          if (qdratio .gt. 1.0d0 .or. qdratio .lt. 0.0d0) then
            messag = 'sw-management error with storage (qdratio)'
            call fatalerr ('Bocodre',messag)
          endif

          do 820 level=1+NRPRI,nrlevs
            qdrain(level) = qdrain(level)*qdratio
 820       continue
          qdrd = qdrdm
        endif
      endif

      return
      end

