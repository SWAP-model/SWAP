! File VersionID:
!   $Id: frozencond.f90 368 2018-01-11 15:44:15Z heine003 $
! ----------------------------------------------------------------------
      subroutine FrozenCond(heat, soil)
! ----------------------------------------------------------------------
!     date               : sept 2005
!     purpose            : if Soil temperatures are simulated, determine 
!                          the reduction factors and frozen depth for 
!                          frozen conditions
! ----------------------------------------------------------------------
! global
      use swap_state_mod, only: heat_state_t, soil_state_t
      implicit none

      type(heat_state_t), intent(inout) :: heat
      type(soil_state_t), intent(in)    :: soil

! local
      integer node
      logical flthaw

! ----------------------------------------------------------------------

!    reduction factor
      do node=1,soil%numnod
         heat%rfcp(node) = 1.0d0
         if (heat%swfrost.eq.1)then
           if(heat%tsoil(node).ge.heat%tfroststa)then
              heat%rfcp(node) = 1.0d0
           else if(heat%tsoil(node).le.heat%tfrostend) then
              heat%rfcp(node) = 0.0d0
           else if(heat%tsoil(node).lt.heat%tfroststa .and.             &
     &             heat%tsoil(node).gt.heat%tfrostend) then 
              heat%rfcp(node) = (heat%tsoil(node)-heat%tfrostend)/      &
     &                          (heat%tfroststa-heat%tfrostend)
           endif
         endif
      end do

 
!     frozen soil : frozen depth z and frozen node nr
      flthaw              = .true.
      heat%nodfrostbot    = -1
      heat%zfrostbot      = 0.0d0
!      nodfrosttop  = -1
      heat%zfrosttop      = 0.0d0

      node = soil%numnod
      do while (flthaw .and. node.gt.1)
         node = node - 1 
         if(heat%tsoil(node) .le. heat%tfrostend+1.0d-6)then
            heat%zfrostbot = soil%z(node+1) + soil%disnod(node+1) *     &
     &           (heat%tfrostend-heat%tsoil(node+1)) /                  &
     &           (heat%tsoil(node)-heat%tsoil(node+1))
            flthaw             =.false.
            heat%nodfrostbot   = node
         endif
      end do

      if(.not.flthaw)then
         flthaw  = .true.
         node = 0
         do while (flthaw .and. node.lt.heat%nodfrostbot)
            node = node + 1 
            if(heat%tsoil(node) .le. heat%tfrostend+1.0d-6)then
               if(node.eq.1) then
                  if(heat%tetop.le.heat%tfrostend) then
                     heat%zfrosttop = 0.0d0
                  else
                     heat%zfrosttop = soil%z(node) -                    &
     &             (soil%z(node) - 0.0d0) *                             &
     &             (heat%tsoil(node)-heat%tfrostend) /                  &
     &             (heat%tsoil(node)-heat%tetop)
                  endif
               else
                  heat%zfrosttop = soil%z(node) + soil%disnod(node) *   &
     &             (heat%tsoil(node)-heat%tfrostend) /                  &
     &             (heat%tsoil(node)-heat%tsoil(node-1))
               endif
               heat%zfrosttop = min(0.0d0,heat%zfrosttop)
               flthaw      =.false.
            endif
         end do
      end if

      return
      end
! ----------------------------------------------------------------------
      subroutine FrozenBounds
! ----------------------------------------------------------------------
!     date               : 20070206
!     purpose            : reduce or stop boundary (drainage and bottom) 
!                          fluxes under frost conditions
! ----------------------------------------------------------------------
!     Swap modules for data communication
      use variables
      use distribute_drainage, only: DIVDRA
      implicit none

!     global - in
!      integer numnod,nodfrostbot,nodfrosttop,nrlevs,swnrsrf,layer(macp)
!      real(8) rfcp(macp),zfrostbot,zfrosttop,ksatfit(maho),thetas(macp)
!      real(8) dz(macp),cofani(maho),qbot,Zbotdr(Madr),L(Madr),theta(macp)
!      real(8) qbot_nonfrozen
!     global - out
!     global - in/out
!      real(8) qdra(Madr,macp),qdrain(Madr),qdrtot

!     local
      integer node,level,layercp(macp),leveldeepest
      real(8) volair,ksatcp(macp),cofanicp(macp),qdratot
      real(8) zdeepest,ztop
      logical frozencomp
      
      ! Hydraulic conductivity for complete frozen soils (constant)
      real(8), parameter :: hconode_vsmall = 1.0d-10

! ----------------------------------------------------------------------

!     initialize qbot
      qbot = qbot_nonfrozen


! --  verify available air volume

      node = numnod
      volair = 0.0d0
      frozencomp = .true.
      do while (frozencomp)
         volair = volair + (thetas(node)-theta(node))*dz(node)
         node = node - 1
         if(node.eq.0)then
            frozencomp = .false.
         else
            if(rfcp(node) .le. 0.01d0)then
               frozencomp = .false.
            end if
         end if 
      end do

! --  consider reduction when volair is very low
!     reduction of drainage only when systems are present
      if(swdra.eq.0) then
         if(nodfrostbot.gt.1 .and. volair.lt.0.01d0)then
            qbot = 0.0d0
         endif
      else
         if(nodfrostbot.gt.1 .and. volair.lt.0.01d0)then

            leveldeepest = 0
            zdeepest     = 0.0d0
            do level=1,nrlevs
               if(zbotdr(level).lt.zdeepest) then
                  leveldeepest = level
                  zdeepest     = zbotdr(level)
               endif
            enddo

            do node=1,numnod
               if(fluseksatexm(node))then
                 ksatcp(node)  = ksatexm(layer(node))*rfcp(node) +      &
     &                      (1.0d0-rfcp(node))*hconode_vsmall
               else
                 ksatcp(node)  = ksatfit(layer(node))*rfcp(node) +      &
     &                      (1.0d0-rfcp(node))*hconode_vsmall
               endif

               cofanicp(node) = cofani(layer(node))
               layercp(node) = node
               do level=1,nrlevs
                  if(zfrostbot.lt.zbotdr(level)) then
                     qdra(level,node) = 0.0d0
                     qdrain(level) = 0.0d0
                  endif
               enddo
            enddo

            qdratot = 0.0d0
            do level = 1,nrlevs
               qdratot = qdratot + qdrain(level)
            end do

            if(abs(qdratot).lt.1.0d-6) then
               if(zfrostbot.lt.zbotdr(leveldeepest)) then
                  qbot = 0.0d0
               else
                  qdrain(leveldeepest) = qbot 
               endif
            else
               do level = 1,nrlevs
                  qdrain(level) = qdrain(level) * (1.0d0 + qbot/qdratot)
               end do
            end if

            if (swdivd.eq.1) then
               ztop = min(gwl,zfrostbot)
               call divdra (numnod,nrlevs,dz,ksatcp,ksatcp,fluseksatexm,&
     &                      layercp,cofanicp,ztop,L,qdrain,qdra,        &
     &                      Swdivdinf,Swnrsrf,SwTopnrsrf,Zbotdr,        &
     &                      dt,FacDpthInf,owltab,t1900)
            endif
         else

            do level = 1,nrlevs
               qdrain(level) = 0.0d0
               do node = 1,numnod
                  qdra(level,node) = qdra(level,node)*rfcp(node)
                  qdrain(level) = qdrain(level) + qdra(level,node)
               end do
            end do

         endif

         qdrtot = 0.0d0
         do level=1,nrlevs
             qdrtot = qdrtot + qdrain(level)
         end do

      endif


!     write(981,'(i5,3('','',f16.10) , 99('','',f16.6:))')               &
!    &      daycum, tcum, t1900, dt, zfrosttop, zfrostbot,gwl,qbot,      &
!    &      (qdrain(level),level=1,Madr)


      return
      end
