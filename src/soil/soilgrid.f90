!> Module for grid and discretization operations
!!
!! This module provides subroutines for calculating grid parameters and
!! converting vertical discretization in soil profiles.
module soilgrid_mod
    implicit none
    private
    public :: calcgrid, convertdiscrvert
  public :: calcgrid_state, convertdiscrvert_state
contains

      !> Calculate grid parameters for soil compartments
      !!
      !! Computes grid parameters including nodal positions, distances between nodes,
      !! layer assignments, and interpolation coefficients for the soil profile.
      !!
      !! @note
      !! Date: Aug 2004
      !! Purpose: calculate grid parameters
      !! @endnote
      subroutine calcgrid 

      use variables
      implicit none

      integer i,j,lay,node,layold
      character(len=200) messag
      character(len=11) tmp

      ! check correct input of number and height of soil compartments
      do i = 1,nsublay
        if (abs(ncomp(i)*hcomp(i) - hsublay(i)) .gt. 1.d-5) then
          ! error in input data
          write(tmp,'(i11)') i
          tmp = adjustl(tmp)
          messag ='At the soil water section, part 4, at layer '//      &
     &   trim(tmp)//' the height of this soil layer hsublay corresponds'&
     &    //' not to the product of height and number of compartments'
          call fatalerr ('calcgrid',messag)
        endif
      end do

      ! position of nodal points and distances between them; also layer of each node
      node = 0
      do i = 1,nsublay
        do j = 1,ncomp(i)
          node = node + 1
          dz(node) = hcomp(i)          
          if (node .eq. 1) then
            z(node) = - 0.5 * dz(node)
            disnod(node) = - z(node)
            layer(node) = isoillay(i)
          else
            z(node) = z(node-1) - 0.5*(dz(node-1)+dz(node))
            disnod(node) = z(node-1)-z(node)
            layer(node) = isoillay(i)
          endif
        end do
      end do
      numnod = node
      disnod(numnod+1) = 0.5*dz(numnod)

      ! store top and bottom depths of each compartment (2019-05-17, MH)
      do i = 1, numnod
         if (i == 1) then
            ztopcp(i) = 0.0d0
            zbotcp(i) = -dz(i)
         else
            ztopcp(i) = zbotcp(i-1)
            zbotcp(i) = zbotcp(i-1) - dz(i)
         end if
      end do      
      
      ! determine bottom compartment of each soil layer
      layold = 1
      do node = 1, numnod
        if (layer(node) .gt. layold) then
          botcom(layold) = node - 1
          layold = layold + 1
        endif
      end do
      numlay = layold
      botcom(numlay) = numnod

      ! linear interpolation values between nodes
      inpolb(1) = 0.5*dz(1)/disnod(2)
      do node = 2,numnod-1
        inpola(node) = 0.5*dz(node)/disnod(node)
        inpolb(node) = 0.5*dz(node)/disnod(node+1)
      end do
      inpola(numnod) = 0.5*dz(numnod)/disnod(numnod)

      ! find first Node of the Layer
      do lay = 1,numlay
        Node = 1
        do while(Layer(Node).ne.lay)
           Node = Node + 1
        enddo
        nod1lay(lay) = node
      enddo


      return
      end subroutine calcgrid

    
      !> Convert vertical discretization
      !!
      !! Converts soil profile discretization from one grid to another, redistributing
      !! state variables and fluxes according to weighted averages or integration.
      !!
      !! @param[in] part Part indicator (1=initial, 2=dynamic)
      !! @param[in] swop Macropore option switch
      !! @param[out] botcomnew Bottom compartment indices for new discretization
      !! @param[out] hnew Pressure head for new discretization
      !! @param[out] thetanew Volumetric water content for new discretization
      !! @param[out] inqnew Water flux for new discretization
      !! @param[out] inqrotnew Root water uptake for new discretization
      !! @param[out] inqdranew Drainage flux for new discretization
      !! @param[out] ithetabegnew Initial water content for new discretization
      !! @param[in] tsoil Soil temperature (old discretization)
      !! @param[out] tsoilnew Soil temperature for new discretization
      !! @param[out] dipocpnew Macropore diameter for new discretization
      !! @param[out] iavfrmpwlwtdm1new Average macropore fraction (domain 1) for new discretization
      !! @param[out] iavfrmpwlwtdm2new Average macropore fraction (domain 2) for new discretization
      !! @param[out] iqexcmtxdm1cpnew Matrix exchange flux (domain 1) for new discretization
      !! @param[out] iqexcmtxdm2cpnew Matrix exchange flux (domain 2) for new discretization
      !! @param[out] iqoutdrrapcpnew Rapid drainage outflow for new discretization
      !! @param[out] vlmpstdm1new Macropore storage volume (domain 1) for new discretization
      !! @param[out] vlmpstdm2new Macropore storage volume (domain 2) for new discretization
      !!
      !! @note
      !! Routine: ConvertDiscrVert()
      !! Purpose: converts vertical descritization 
      !! Usage: Call ConvertDiscrVert(......)
      !! 
      !! Arguments:
      !! I/O Data_type Name         Meaning
      !!  I     I4     SwDiscrVert  
      !!  I     I4     part         
      !!  I     I4     Swop  
      !!  I     I4     logf         Unit number
      !!  I     I4     nrlevs       
      !!  I     I4     numlay       
      !!  I     I4     botcom(maho) 
      !!  I     I4     numnod       
      !!  I     R8     dz(macp)     
      !!  I     R8     theta(macp)  
      !!  I     R8     inq(macp)    
      !!  I     R8     inqrot(macp) 
      !!  I     R8     inqdra(macp) 
      !!  I     R8     Tsoil(0:macp)  
      !!  I     I4     numnodNew       
      !!  I     R8     dzNew(macp) 
      !!  I     R8     DiPoCp(macp)
      !!  I     R8     IAvFrMpWlWtDm1(macp)
      !!  I     R8     IAvFrMpWlWtDm2(macp)
      !!  I     R8     IQExcMtxDm1Cp(macp)
      !!  I     R8     IQExcMtxDm2Cp(macp)
      !!  I     R8     IQOutDrRapCp(macp)
      !!  I     R8     VlMpStDm1
      !!  I     R8     VlMpStDm2
      !!  O     I4     botcomNew(maho) 
      !!  O     R8     thetaNew(macp)  
      !!  O     R8     inqNew(macp)    
      !!  O     R8     inqrotNew(macp) 
      !!  O     R8     inqdraNew(macp) 
      !!  O     R8     TsoilNew(0:macp) 
      !!  O     R8     DiPoCpNew(macp)
      !!  O     R8     IAvFrMpWlWtDm1New(macp)
      !!  O     R8     IAvFrMpWlWtDm2New(macp)
      !!  O     R8     IQExcMtxDm1CpNew(macp)
      !!  O     R8     IQExcMtxDm2CpNew(macp)
      !!  O     R8     IQOutDrRapCpNew(macp)
      !!  O     R8     VlMpStDm1New(macp)
      !!  O     R8     VlMpStDm2New(macp)
      !! @endnote
      subroutine ConvertDiscrVert(part,swop,botcomnew,hnew,thetanew,inqnew,            &
                                  inqrotnew,inqdranew,ithetabegnew,tsoil,tsoilnew,     &
                                  dipocpnew,iavfrmpwlwtdm1new,iavfrmpwlwtdm2new,       &
                                  iqexcmtxdm1cpnew,iqexcmtxdm2cpnew,iqoutdrrapcpnew,   &
                                  vlmpstdm1new,vlmpstdm2new)

      !---- Declarations
      use variables, only: SwDiscrvert,nrlevs,numlay,botcom,numnod,dz,h,theta,inq,inqrot,inqdra,                           &
                           IThetaBeg,cofgen,numnodNew,dzNew,DiPoCp,FrArMtrx,IAvFrMpWlWtDm1,IAvFrMpWlWtDm2,           &
                           IQExcMtxDm1Cp,IQExcMtxDm2Cp,IQOutDrRapCp,VlMpStDm1,VlMpStDm2
      use soilhydraulics_utils, only: prhead
      use swap_array_dimensions, only: macp, maho, madr
      IMPLICIT NONE

      ! global variables
      integer   part,Swop
      integer   botcomNew(maho)
      real(8)   hNew(macp),thetaNew(macp),inqrotNew(macp),inqNew(macp+1),inqdraNew(Madr,macp),Tsoil(0:macp),TsoilNew(0:macp)
      real(8)   IThetaBegNew(MaCp)
      real(8)   DiPoCpNew(macp), IAvFrMpWlWtDm1New(macp) 
      real(8)   IAvFrMpWlWtDm2New(macp), IQExcMtxDm1CpNew(macp)
      real(8)   IQExcMtxDm2CpNew(macp), IQOutDrRapCpNew(macp)
      real(8)   VlMpStDm1New(macp), VlMpStDm2New(macp)

      ! local variables
      integer   lay,node,nodeN,nodeNew(macp,2),i,level
      real(8)   disnodNew(macp+1),total,zNew(macp)
      real(8)   CofgenNew(21,macp)
      character(len=80) Message
      character(len=*), parameter :: ModuleName = 'ConvertDiscrVert'

      ! SAVE statement removed - local variables are temporary work arrays
      ! ModuleName converted to parameter for proper initialization

      ! error in call of part
      if (part.lt.1 .or. part.gt.2) then
        write(message,*) 'fatal error in variabel PART '
        call fatalerr(ModuleName,message)
      endif

      ! always convert to maintain identic writing

      if(SwDiscrVert.eq.0) then
        ! Simply copy
        numnodNew = numnod
        do lay = 1,numlay
          botcomNew(lay) = botcom(lay)
        enddo
        do node = 1,numnodNew
          dzNew(node) = dz(node)
          hNew(node) = h(node)
          thetaNew(node) = theta(node)
          IThetaBegNew(node) = IThetaBeg(node)
          inqNew(node) = inq(node)
          inqrotNew(node) = inqrot(node)
          do level=1,nrlevs
            inqdraNew(level,node) = inqdra(level,node)
          enddo

          if (swop.eq.2) then
             VlMpStDm1New(node)     = VlMpStDm1(node)
             VlMpStDm2New(node)     = VlMpStDm2(node)
             DiPoCpNew(node)        = DiPoCp(node)
             IQExcMtxDm1CpNew(node) = IQExcMtxDm1Cp(node)
             IQOutDrRapCpNew(node)  = IQOutDrRapCp(node)
             IAvFrMpWlWtDm1New(node)= IAvFrMpWlWtDm1(node)
             IQExcMtxDm2CpNew(node) = IQExcMtxDm2Cp(node)
             IAvFrMpWlWtDm2New(node)= IAvFrMpWlWtDm2(node)

             ThetaNew(node)         = FrArMtrx(node)*Theta(node)
             IThetaBegNew(node)     = FrArMtrx(node)*IThetaBeg(node)
          endif
        enddo
        do node = 0,numnodNew
           TsoilNew(node) = Tsoil(node)
        enddo
        inqNew(numnodNew+1) = inq(numnod+1)

      else if(SwDiscrVert.eq.1) then
        ! initial part (part 1)
        if (part.eq.1) then

          ! calculate zNew = depth at bottom of new compartment
          Total = 0.0d0
          Do node = 1,numnodNew
            Total = Total + dzNew(node)
            zNew(node) = Total
          enddo

          ! Fill array NodeNew(Numnod,1) = top old compartment and NodeNew(Numnod,2) = bottom old comp. in new comp. NodeN
          NodeN = 1
          Node = 1
          Total = 0.0d0
          do while (node.le.numnod) 
            Total = Total + dz(node) 
            If (ABS(zNew(NodeN)-Total).LT.1.0D-6) then
              NodeNew(NodeN,2) = Node
              NodeN = NodeN+1
            endif
            node = node+1
          enddo
          NodeNew(1,1) = 1
          Do node = 2,numnodNew
            NodeNew(Node,1) = NodeNew(Node-1,2)+1
          Enddo

          ! botcomNew
          do lay = 1, numlay
            do node = 1,numnodNew
              if (NodeNew(Node,2).eq.botcom(lay)) then
                botcomNew(lay) = node
              else
                continue
              endif
            Enddo
          Enddo

          ! for hNew: position of nodal points and distances between them
          zNew(1) = - 0.5d0*dzNew(1)
          disnodNew(1) = - zNew(1)
          do node = 2,numnodNew
            zNew(node) = zNew(node-1)-0.5d0*(dzNew(node-1)+dzNew(node))
            disnodNew(node) = zNew(node-1)-zNew(node)
          enddo

        endif

        ! In both part 1 and 2: convert theta to thetaNew based on weighed average
        Do node = 1,NumNodNew
          Total = 0.0d0
          Do i = NodeNew(node,1),NodeNew(node,2)   
            Total = Total+dz(i)
          Enddo
          thetaNew(node) = 0.0d0
          IThetaBegNew(node) = 0.d0
          Do i = NodeNew(node,1),NodeNew(node,2)
            if (Swop.ne.2) then
               thetaNew(node) = thetaNew(node)+theta(i)*dz(i)/Total
               IThetaBegNew(node) = IThetaBegNew(node)+IThetaBeg(i)*    &
     &                              dz(i)/Total
            else
               thetaNew(node) = thetaNew(node)+FrArMtrx(i)*theta(i)*    &
     &                          dz(i)/Total
               IThetaBegNew(node) = IThetaBegNew(node)+FrArMtrx(i)*     &
     &                              IThetaBeg(i)*dz(i)/Total
            endif
          Enddo
        Enddo

        ! convert Tsoil to TSoilNew based on weighed average
        Do node = 1,NumNodNew
          Total = 0.0d0
          Do i = NodeNew(node,1),NodeNew(node,2)   
            Total = Total+dz(i)
          Enddo
          TsoilNew(node) = 0.0d0
          Do i = NodeNew(node,1),NodeNew(node,2)
            TsoilNew(node) = TsoilNew(node)+Tsoil(i)*dz(i)/Total
          Enddo               
        Enddo

        if (swop.eq.2) then
          ! convert VlMpStDm to VlMpStDmNew based on integration
           Do node = 1,NumNodNew
             VlMpStDm1New(node) = 0.0d0
             VlMpStDm2New(node) = 0.0d0
             Do i = NodeNew(node,1),NodeNew(node,2)
               VlMpStDm1New(node) = VlMpStDm1New(node) + VlMpStDm1(i)
               VlMpStDm2New(node) = VlMpStDm2New(node) + VlMpStDm2(i)
             Enddo               
           Enddo

           ! convert DiPoCp to DiPoCpNew based on weighed average
           Do node = 1,NumNodNew
             Total = 0.0d0
             Do i = NodeNew(node,1),NodeNew(node,2)   
               Total = Total+dz(i)
             Enddo
             DiPoCpNew(node) = 0.0d0
             Do i = NodeNew(node,1),NodeNew(node,2)
               DiPoCpNew(node) = DiPoCpNew(node)+DiPoCp(i)*dz(i)/Total
             Enddo               
           Enddo
        endif

        ! only in dynamic part (part 2)
        if (part.eq.2) then

          ! Determine hNew
          do node = 1,numnodNew

             ! Bug fixed: Rob H 28 okt 2011 - this output for tabulated soil physics is suppressed
             do i = 1, 9
                cofgenNew(i,node) = cofgen(i,NodeNew(node,2))  ! CofgenNew Newnode = Cofgen bottom old node
             enddo
             
            hNew(node) = prhead (node,disnodNew(node),thetaNew(node),cofgenNew,hNew)
          enddo

          ! convert inqrot to inqrotNew based on integration
          Do node = 1,NumNodNew
            inqrotNew(node) = 0.0d0
            Do i = NodeNew(node,1),NodeNew(node,2)
              inqrotNew(node) = inqrotNew(node)+inqrot(i)
            Enddo               
          Enddo

          ! convert inq to inqNew
          Do node = 1,NumNodNew
            inqNew(node) = inq(NodeNew(node,1))
          Enddo
          inqNew(NumNodNew+1) = inq(numnod+1)

          ! convert inqdra to inqdraNew based on integration
          Do level = 1,nrlevs
            Do node = 1,NumNodNew
              inqdraNew(level,node) = 0.0d0
              Do i = NodeNew(node,1),NodeNew(node,2)
                inqdraNew(level,node) = inqdraNew(level,node) +         &
     &                                 inqdra(level,i)
              Enddo               
            Enddo
          Enddo
        endif

        ! convert macropore fluxes to New fluxes based on integration, and
        ! IAvFrMpWlWtDm1 to IAvFrMpWlWtDm1New, based on weighed average
        if (swop.eq.2) then
          Do node = 1,NumNodNew
            Total = 0.0d0
            Do i = NodeNew(node,1),NodeNew(node,2)   
              Total = Total+dz(i)
            Enddo
            IQExcMtxDm1CpNew(node) = 0.0d0
            IQExcMtxDm2CpNew(node) = 0.0d0
            IQOutDrRapCpNew(node)  = 0.0d0
            IAvFrMpWlWtDm1New(node) = 0.0d0
            IAvFrMpWlWtDm2New(node) = 0.0d0
            Do i = NodeNew(node,1),NodeNew(node,2)
              IQExcMtxDm1CpNew(node) =                                  &
     &           IQExcMtxDm1CpNew(node) + IQExcMtxDm1Cp(i)
              IQExcMtxDm2CpNew(node) =                                  &
     &           IQExcMtxDm2CpNew(node) + IQExcMtxDm2Cp(i)
              IQOutDrRapCpNew(node) =                                   &
     &           IQOutDrRapCpNew(node)  + IQOutDrRapCp(i) 
              IAvFrMpWlWtDm1New(node) =                                 &
     &           IAvFrMpWlWtDm1New(node) + IAvFrMpWlWtDm1(i)*dz(i)/Total
              IAvFrMpWlWtDm2New(node) =                                 &
     &           IAvFrMpWlWtDm2New(node) + IAvFrMpWlWtDm2(i)*dz(i)/Total
            Enddo               
          Enddo
        endif

      else
        write(message,*) 'fatal error in variable SwDiscrVert'
        call fatalerr(ModuleName,message)
      endif

      return
      end subroutine ConvertDiscrVert

      !> State-aware wrapper for `calcgrid`
      !!
      !! Executes the legacy grid setup and synchronizes grid dimensions and,
      !! when allocated, grid arrays to explicit state.
      !!
      !! @param[inout] state SWAP model state container
      subroutine calcgrid_state(state)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilgrid_outputs_from_variables
      use variables, only: numnod, numlay
      implicit none

      type(swap_state_t), intent(inout) :: state

      call calcgrid()

      state%numnod = numnod
      state%numlay = numlay

      if (allocated(state%soil%z)) then
         call soilgrid_outputs_from_variables(state%soil, numnod, numlay)
      end if
      end subroutine calcgrid_state

      !> State-aware wrapper for `ConvertDiscrVert`
      !!
      !! Executes legacy discretization conversion and then synchronizes current
      !! soil balance outputs that may be used downstream.
      !!
      !! @param[inout] state SWAP model state container
      !! @param[in]    part Part indicator (1=initial, 2=dynamic)
      !! @param[in]    swop Macropore option switch
      !! @param[out]   botcomnew Bottom compartment indices for new discretization
      !! @param[out]   hnew Pressure head for new discretization
      !! @param[out]   thetanew Volumetric water content for new discretization
      !! @param[out]   inqnew Water flux for new discretization
      !! @param[out]   inqrotnew Root water uptake for new discretization
      !! @param[out]   inqdranew Drainage flux for new discretization
      !! @param[out]   ithetabegnew Initial water content for new discretization
      !! @param[in]    tsoil Soil temperature (old discretization)
      !! @param[out]   tsoilnew Soil temperature for new discretization
      !! @param[out]   dipocpnew Macropore diameter for new discretization
      !! @param[out]   iavfrmpwlwtdm1new Average macropore fraction (domain 1)
      !! @param[out]   iavfrmpwlwtdm2new Average macropore fraction (domain 2)
      !! @param[out]   iqexcmtxdm1cpnew Matrix exchange flux (domain 1)
      !! @param[out]   iqexcmtxdm2cpnew Matrix exchange flux (domain 2)
      !! @param[out]   iqoutdrrapcpnew Rapid drainage outflow
      !! @param[out]   vlmpstdm1new Macropore storage volume (domain 1)
      !! @param[out]   vlmpstdm2new Macropore storage volume (domain 2)
      subroutine convertdiscrvert_state(state, part, swop, botcomnew, hnew, thetanew, inqnew,     &
                    inqrotnew, inqdranew, ithetabegnew, tsoil, tsoilnew,      &
                    dipocpnew, iavfrmpwlwtdm1new, iavfrmpwlwtdm2new,           &
                    iqexcmtxdm1cpnew, iqexcmtxdm2cpnew, iqoutdrrapcpnew,       &
                    vlmpstdm1new, vlmpstdm2new)
      use swap_state_mod, only: swap_state_t
      use swap_state_sync, only: soilwaterbalance_outputs_from_variables
      use swap_array_dimensions, only: macp, maho, madr
      implicit none

      type(swap_state_t), intent(inout) :: state
      integer, intent(in) :: part, swop
      integer, intent(out) :: botcomnew(maho)
      real(8), intent(out) :: hnew(macp), thetanew(macp), inqrotnew(macp), inqnew(macp+1)
      real(8), intent(out) :: inqdranew(madr,macp), ithetabegnew(macp)
      real(8), intent(in) :: tsoil(0:macp)
      real(8), intent(out) :: tsoilnew(0:macp)
      real(8), intent(out) :: dipocpnew(macp), iavfrmpwlwtdm1new(macp), iavfrmpwlwtdm2new(macp)
      real(8), intent(out) :: iqexcmtxdm1cpnew(macp), iqexcmtxdm2cpnew(macp), iqoutdrrapcpnew(macp)
      real(8), intent(out) :: vlmpstdm1new(macp), vlmpstdm2new(macp)

      call ConvertDiscrVert(part, swop, botcomnew, hnew, thetanew, inqnew, inqrotnew, inqdranew,   &
                ithetabegnew, tsoil, tsoilnew, dipocpnew, iavfrmpwlwtdm1new,            &
                iavfrmpwlwtdm2new, iqexcmtxdm1cpnew, iqexcmtxdm2cpnew, iqoutdrrapcpnew, &
                vlmpstdm1new, vlmpstdm2new)

      call soilwaterbalance_outputs_from_variables(state%soil, state%numnod)
      end subroutine convertdiscrvert_state

end module soilgrid_mod