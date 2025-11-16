! File VersionID:
!   $Id: calcgwl.f90 341 2017-09-29 18:12:25Z kroes006 $
! ----------------------------------------------------------------------
      subroutine calcgwl ()
! ----------------------------------------------------------------------
!     date               : july 2002, updated april 2008
!     purpose            : search for the watertable and perched
!                          watertable (if existing).
!
!     update             : an unsaturated zone embedded in a saturated soil
!                          column should contain at least a total of 'CritAir' cm 
!                          of air to be recognized as really unsaturated
! ----------------------------------------------------------------------
      use variables, only: disnod,logf,swscre,swbotb,flmacropore,numnod,gwlinp,h,z,pond,t1900,  &
                           gwl,nodgwl,bpegwl,npegwl,pegwl,nodgwlflcpzo,gwlflcpzo,CritUndSatVol

      implicit none
! --- global

! --- local
      integer   i, node, nodhlp, nodheq1
      real(8)   level
      logical   flsat,flunsat
      character(len=200) messag
      character(len=19) datexti

      save
! ----------------------------------------------------------------------

! --- set initial values
      gwl       = 999.0d0
      pegwl     = 999.0d0
      flsat     = .false.
      nodgwl    = numnod+1
      nodhlp    = numnod
      nodheq1   = numnod

! --- search for groundwater table
      if (h(numnod).ge.0.0d0) flsat  = .true.

      node = numnod
      nodgwlflcpzo = numnod + 1
      gwlflcpzo    = gwl
      do while (flsat .and. node.gt.1)
         node = node - 1 
         if(swbotb.eq.1)then
            if (h(node) .lt. 0.0d0) then
               gwl = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
               flsat   =.false.
               nodgwl  = node
            endif
         else 
!
            if (h(node) .lt. 1.0d0 .and. nodheq1.eq.numnod) nodheq1 = node
!
            if (h(node) .lt. 0.0d0) then
               if (.not.flmacropore) then
                  flsat  = .false.
                  nodgwl = node
                  gwl    = level (1,node,nodheq1)
               elseif (flmacropore) then
                  if (gwl.gt.990.0d0) then
                     nodgwl = node
                     gwl    = level (2,node,nodheq1)
                  endif
                  call watertable (node,nodgwlflcpzo,nodhlp,nodheq1,0.0d0,flsat,gwlflcpzo)
               endif
            endif
         endif
      end do

!   - whole profile saturated, then add ponding layer to groundwater level
      if (flsat)then
         if(h(1) .gt. 0.0d0)then
            if (pond .lt. 1.d-8) then
               gwl = min(z(1)+h(1),pond)
            else
               gwl = pond
            endif
         else
            gwl = 0.0d0
         end if
         nodgwl = 1
         if (flmacropore) then
            nodgwlflcpzo = 1
            gwlflcpzo    = gwl
         endif         
      endif         
 
! --- search for perched groundwater table

!   - first, search for first saturated compartment (i) above groundwater level
      i = nodhlp
      flunsat = .true.
      do while (flunsat .and. i.ge.1)
         if (h(i).ge.0.0d0) flunsat = .false.
         i = i - 1 
      enddo
!
!   - if saturated compartment above gwl exists, then find perched groundwater table
      if (i.ne.0) then
         flsat  = .true.
         bpegwl = i
         node   = bpegwl
         nodheq1 = bpegwl

         do while (flsat .and. node.gt.1)
            node = node - 1 
!
            if (h(node) .lt. 1.0d0 .and. nodheq1.eq.bpegwl) nodheq1 = node
!
            if (h(node) .lt. 0.0d0) then
               if (.not.flmacropore) then
                  flsat = .false.
                  npegwl = node
                  pegwl  = level (1,node,nodheq1)
               elseif (flmacropore) then
                  call watertable (node,npegwl,nodhlp,nodheq1,CritUndSatVol,flsat,pegwl)
               endif
            endif
         end do
!
!   - whole profile saturated, then add ponding layer to perched groundwater level
         if (flsat)then
            if(h(1) .gt. 0.0d0)then
               if (pond .lt. 1.d-8) then
                  pegwl = min(z(1)+h(1),pond)
               else
                  pegwl = pond
               endif
            else
               pegwl = 0.0d0
            end if
            npegwl = 1
         endif 
      else
         bpegwl = -1
         npegwl = -1
      endif

! --- fatal error if gwl below profile and flux has to be calculated
      if ((swbotb.eq.3.or.swbotb.eq.4).and.gwl.gt.998.0d0) then
          messag = 'The groundwater level descends below the lower'     &
     &     //' boundary. This conflicts with bottom boundary'           &
     &     //' condition 3 and 4. Extend soil profile!'
         call fatalerr ('calcgwl',messag)
      endif

! --- warning error if there is inconsistency between defined gwl and soil physics
      if (swbotb.eq.1 .and. (gwlinp .ge.z(1) .or. gwl.gt.998.0d0)) then
!         determine date and date-time
         call dtdpst('year-month-day,hour:minute:seconds',t1900,datexti)
         write(messag,'(6a)')                                           &
     &         'No groundwater level because unsaturation at bottom ',  &
     &         'compartment ( ', datexti,  ' ). ',                      &
     &         'This is caused by inconsistency between ',              &
     &         'given gwl and soil physical parameters '
         call warn ('Calcgwl',messag,logf,swscre)
      endif

      return
      end


!-----------------------------------------------------------------------
      real(8) function level (swoptlev,node,nodheq1)
! ----------------------------------------------------------------------
!     Date               : april 2008
!     Purpose            : calc. water level from pressure head
! ----------------------------------------------------------------------
      use variables, only: numnod, disnod, dz, h, z, zbotcp
      implicit none
      
! --- global 
      integer node, nodheq1, swoptlev

! --- local
      integer i
      real(8) levm1, levp1
! ----------------------------------------------------------------------

      if (swoptlev.eq.1) then
! - groundwater level equals elevation head where h = 0
      if (h(node+1).ge.0.0d0)then
         level = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
      else   
         level = zbotcp(node) - h(node)
         level = min(z(node),max(zbotcp(node),level))
      end if
!
      elseif (swoptlev.eq.2) then
! - groundwater level equals average of elevation heads of h = -1 and h = +1
!   - elevation head of h = +1
         i = nodheq1
         if (nodheq1.eq.numnod) then
            levp1 = z(i) - 0.5d0 * dz(i)
         else
            levp1 = z(i) - (z(i) - z(i+1)) * (1.d0-h(i)) / (h(i+1)-h(i))
         endif
!   - elevation head of h = -1
         i = node
         do while (h(i).gt.-1.d0 .and. i.gt.1)
            i = i - 1
         enddo
         if (i.eq.1 .and. h(1).gt.-1.d0 .and. node.gt.2) then
!  - no compartment with pressure head < -1 cm in top of profile:
!    use elevation head of h = 0 as estimation for groundwater level
            levm1 = z(node+1) + h(node+1) / (h(node+1)-h(node)) * disnod(node+1)
            levp1 = levm1
         else
            levm1 = z(i+1) + (z(i) - z(i+1)) * (1.d0+h(i+1)) / (h(i+1)-h(i))
         endif
!   - groundwater level = average of levp1 and levm1
         level = (levp1 + levm1) / 2.d0
      endif

      return
      end

! ----------------------------------------------------------------------
      subroutine watertable (node,nodlev,nodhlp,nodheq1,CritUndSatVol,flsat,waterlevel)
! ----------------------------------------------------------------------
!     date               : april 2008
!     purpose            : search for watertable and perched
!                          watertable (if existing).
! ----------------------------------------------------------------------
      use variables, only: numnod,dz,h,Theta,ThetaS,z
      implicit none

! --- global                                                          in
      integer node,nodhlp,nodheq1,nodlev
      logical flsat
!                                                                    out
      real(8) waterlevel
! ----------------------------------------------------------------------
! --- local
      integer i
      real(8) CritUndSatVol, level, TotUndSatVol
      logical flsat2
! ----------------------------------------------------------------------
      TotUndSatVol = 0.0d0
      flsat2 = .false. 
      i = node
      do while (TotUndSatVol.lt.CritUndSatVol .and. .not.flsat2 .and. i.ge.1)
         TotUndSatVol = TotUndSatVol + (ThetaS(i) - Theta(i)) * dz(i)
         if (h(i).gt.-1.d-7) flsat2 = .true.
         i = i - 1
      enddo
!
      if (i.eq.0 .or. TotUndSatVol.gt.CritUndSatVol-1.d-8) then
         flsat = .false.
!!!!!!!! NODGWL is NOT node with GWL, but DEEPEST UNSATURATED NODE !!!!!!!!!!!!!
         nodlev = node
         nodhlp = i
      elseif(flsat2) then
         node   = i + 1
      endif
!
      if (.not.flsat) then        
! find groundwater level containing node   
         if (CritUndSatVol.gt.0.d0) then  
            waterlevel = level(1,node,nodheq1)
         else
            waterlevel = level(2,node,nodheq1)
         endif
         i = max(node-2,1)
         do while(z(i)-0.5d0*dz(i).gt.waterlevel .and. i.gt.2 .and. i.lt.numnod)
            i = i + 1
         enddo
         nodlev = min(max(i,1),numnod)
      endif

      return
      end
