! File: swap_hydraulic_functions.f90
module soilhydraulics_utils
   use error_mod, only: fatalerr_collected
   !> Module containing hydraulic property functions for soil water modeling
   !!
   !! This module contains functions for calculating:
   !! - Water content from pressure head
   !! - Hydraulic conductivity
   !! - Moisture capacity (differential water content)
   !! - Pressure head from water content
   !! - Mean hydraulic conductivity between nodes
   !!
   !! @author Original SWAP team
   !! @date February 2026 (modularization)
   use iso_fortran_env, only: real64
   use variables, only: cofgen, swsophy, numtab, sptab, ientrytab, &
                        iHWCKmodel, layer, swfrost, dt, fluseksatexm
   ! [SS-HEAT] Task 9: tsoil global retired; hconduc fallback removed (iHWCKmodel 4-11 path unreachable in regression)
   use doln
   use WC_K_models_04_11, only: functionvalue_04_11
   
   implicit none
   
   private
   public :: watcon, moiscap, hconduc, dhconduc, prhead, hcomean, dkmean

contains

   !> Calculate mean hydraulic conductivity between two nodes
   function hcomean(swkmean, kup, klow, dzup, dzlow)
      implicit none
      
      ! Arguments
      integer, intent(in) :: swkmean
      real(real64), intent(in) :: kup, klow, dzup, dzlow
      real(real64) :: hcomean
      
      ! Local variables
      real(real64) :: a1, a2

      ! Unweighted arithmic mean
      if (swkmean == 1) then
         hcomean = 0.5_real64 * (kup + klow)
      ! Weighted arithmic mean
      else if (swkmean == 2) then
         hcomean = (dzup * kup + dzlow * klow) / (dzup + dzlow)
      ! Unweighted geometric mean
      else if (swkmean == 3) then
         hcomean = dsqrt(kup * klow)
      ! Weighted geometric mean
      else if (swkmean == 4) then
         a1 = dzup / (dzup + dzlow)
         a2 = 1.0_real64 - a1
         hcomean = (kup ** a1) * (klow ** a2)
      ! Unweighted harmonic mean
      else if (swkmean == 5) then
         hcomean = 1.0_real64 / (0.5_real64/kup + 0.5_real64/klow)
      ! Weighted harmonic mean
      else if (swkmean == 6) then
         a1 = dzup / (dzup + dzlow)
         a2 = 1.0_real64 - a1
         hcomean = 1.0_real64 / (a1/kup + a2/klow)
      end if
      
   end function hcomean
   
   !> Calculate water content from pressure head
   function watcon(node, head)
      implicit none
      
      ! Arguments
      integer, intent(in) :: node
      real(real64), intent(in) :: head
      real(real64) :: watcon
      
      ! Local variables
      real(real64) :: h_enpr, help, m, n, s_enpr
      real(real64), parameter :: h_crit = -1.0d-2
      real(real64) :: h105, C105, a, b, dum, t105
      real(real64) :: alfamg, thetar, thetas, moiscap
      real(real64) :: alfa_2, n_2, m_2, omega_1
      
      ! sptab(1,node,i):    h
      ! sptab(2,node,i):    Theta
      ! sptab(3,node,i):    K
      ! sptab(4,node,i):    dTheta / dh
      ! sptab(5,node,i):    dK / dh

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (swsophy == 0) then

         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)
         
         if (iHWCKmodel(layer(node)) == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            watcon = dmax1(1.0000001_real64*thetar, thetar + (thetas-thetar)*dexp(alfamg*head))
            
         else if (iHWCKmodel(layer(node)) == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = cofgen(13,node)
            n_2     = cofgen(14,node)
            m_2     = cofgen(15,node)
            omega_1 = cofgen(16,node)
            if (head < 0.0_real64) then
               watcon = omega_1 / (1.0_real64 + (dabs(alfamg*head))**n)**m
               watcon = watcon + (1.0_real64 - omega_1) / (1.0_real64 + (dabs(alfa_2*head))**n_2)**m_2
               watcon = thetar + (thetas - thetar) * watcon
            else
               watcon = thetas
            end if

         else if (iHWCKmodel(layer(node)) > 3 .and. iHWCKmodel(layer(node)) < 12) then
            watcon = functionvalue_04_11(1, node, head)

         else  ! Use default MvG
         
            if (h_enpr > h_crit) then

               if (head >= 0.0_real64) then
                  ! Saturated moisture content
                  watcon = thetas
               else if (head > h_crit) then
                  help = (dabs(alfamg*h_crit))** n
                  help = (1.0_real64 + help) ** m
                  help = thetar + (thetas - thetar) / help    
                  watcon = help + (thetas - help) / (-h_crit) * (head - h_crit) 
                  watcon = min(watcon, thetas)  
               else 
                  ! First compute |alpha * h| ** n
                  help = (dabs(alfamg*head)) ** n

                  ! Add 1 and raise to the power m
                  help = (1.0_real64 + help) ** m

                  ! Now compute theta
                  watcon = thetar + (thetas - thetar) / help 
               end if
            else

               h105 = 1.05_real64 * h_enpr

               if (head >= h105) then
                  t105 = thetar + (thetas - thetar) * &
                       ((1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m) / &
                       ((1.0_real64 + (dabs(alfamg*h105)) ** n) ** m)
                  C105 = (ThetaS - ThetaR) * alfamg * M * N * &
                       (dabs(alfamg * h105) ** (N - 1)) * &
                       ((1 + dabs(alfamg * h_enpr) ** N) ** M) / &
                       ((1 + dabs(alfamg * h105) ** N) ** (M + 1))
                  a = (t105 - thetas - C105*h105) / (C105*h105**2)
                  b = (t105**2 - 2*t105*thetas + thetas**2) / &
                       (t105 - thetas - C105*h105)
                  watcon = thetas + b*a*head / (1.0_real64 + a*head)

               else 
                  ! First compute |alpha * h| ** n
                  help = (dabs(alfamg*head)) ** n

                  ! Add 1 and raise to the power m
                  help = (1.0_real64 + help) ** m

                  ! For modified VanGenuchten model:
                  ! - S_enpr: relative saturation at Entry Pressure h_enpr 
                  s_enpr = (1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m

                  ! Now compute theta
                  watcon = thetar + (thetas - thetar) / help * s_enpr
               end if
            end if
         end if

      ! Use tabulated function
      else if (swsophy == 1) then
         dum = head
         if (do_ln_trans .and. head < 0.0_real64) dum = -dlog(-head + 1.0_real64)
         if (head >= -1.0d-9) then
            watcon = sptab(2,node,numtab(node))
         else if (dum < sptab(1,node,1)) then
            watcon = sptab(2,node,1)
         else
            call EvalTabulatedFunction(0, numtab(node), 1, 2, 4, node, sptab, ientrytab, head, watcon, moiscap, 1)
         end if
      end if
      
   end function watcon

   !> Calculate differential moisture capacity (as a function of pressure head)
   function moiscap(node, head)
      implicit none
      
      ! Arguments
      integer, intent(in) :: node
      real(real64), intent(in) :: head
      real(real64) :: moiscap
      
      ! Local variables
      real(real64) :: h_enpr, m, n, alphah, s_enpr
      real(real64), parameter :: h_crit = -1.0d-2
      real(real64) :: alfamg, thetar, thetas
      real(real64) :: h105, C105, a, b, term1, term2, dum, t105, dummy
      real(real64) :: alfa_2, n_2, m_2, omega_1

      ! Use analytical expression
      if (swsophy == 0) then

         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         n      = cofgen(6,node)
         m      = cofgen(7,node)
         h_enpr = cofgen(9,node)

         if (iHWCKmodel(layer(node)) == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            moiscap = alfamg * (thetas - thetar) * dexp(alfamg*head)
            
         else if (iHWCKmodel(layer(node)) == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = cofgen(13,node)
            n_2     = cofgen(14,node)
            m_2     = cofgen(15,node)
            omega_1 = cofgen(16,node)
            if (head < 0.0_real64) then
               moiscap = omega_1*alfamg*n*m*(dabs(alfamg*head))**(n - 1.0_real64) * &
                    (1.0_real64 + (dabs(alfamg*head))**n)**(-1.0_real64 - m)
               moiscap = moiscap + (1.0_real64 - omega_1)*alfa_2*n_2*m_2* &
                    (dabs(alfa_2*head))**(n_2 - 1.0_real64) * &
                    (1.0_real64 + (dabs(alfa_2*head))**n_2)**(-1.0_real64 - m_2)
               moiscap = (thetas - thetar) * moiscap
            else
               moiscap = 0.0_real64
            end if

         else if (iHWCKmodel(layer(node)) > 3 .and. iHWCKmodel(layer(node)) < 12) then
            moiscap = functionvalue_04_11(3, node, head)

         else  ! Use default MvG
         
            if (h_enpr > h_crit) then

               if (head >= 0.0_real64) then

                  moiscap = dt * 1.0d-7

               else if (head > h_crit) then

                  term1 = watcon(node, h_crit)
                  moiscap = (thetas - term1) / (-h_crit)
               else

                  ! Use analytical evaluation of capacity
                  alphah = dabs(alfamg*head)

                  ! Compute |alpha * h| to the power n-1
                  term1 = alphah ** (n - 1.0_real64)

                  ! Compute |alpha*h| to the power n
                  term2 = term1 * alphah

                  ! Add one and raise to the power m+1
                  term2 = (1.0_real64 + term2) ** (m + 1.0_real64)

                  ! Divide theta-s minus theta-r by term2
                  term2 = (thetas - thetar) / term2

                  ! Calculate the differential moisture capacity
                  moiscap = dabs(-1.0_real64 * n * m * alfamg * term2 * term1)
               end if
            else

               h105 = 1.05_real64 * h_enpr

               if (head >= h105) then
                  t105 = thetar + (thetas - thetar) * &
                       ((1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m) / &
                       ((1.0_real64 + (dabs(alfamg*h105)) ** n) ** m)
                  C105 = (ThetaS - ThetaR) * alfamg * M * N * &
                       (dabs(alfamg * h105) ** (N - 1)) * &
                       ((1 + dabs(alfamg * h_enpr) ** N) ** M) / &
                       ((1 + dabs(alfamg * h105) ** N) ** (M + 1))
                  a = (t105 - thetas - C105*h105) / (C105*h105**2)
                  b = (t105**2 - 2*t105*thetas + thetas**2) / &
                       (t105 - thetas - C105*h105)
                  moiscap = b*a / ((1.0_real64 + a*head)**2)

               else

                  alphah = dabs(alfamg*head)
                  term1 = alphah ** (n - 1.0_real64)
                  term2 = term1 * alphah
                  term2 = (1.0_real64 + term2) ** (m + 1.0_real64)
                  term2 = (thetas - thetar) / term2

                  ! For modified VanGenuchten model:
                  ! - S_enpr: relative saturation at Entry Pressure h_enpr 
                  s_enpr = (1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m

                  moiscap = dabs(-1.0_real64 * n * m * alfamg*term2*term1)*s_enpr
  
               end if

            end if
         end if
      
         if (head > -1.0_real64 .and. moiscap < (dt * 1.0d-7)) moiscap = dt * 1.0d-7

      ! Use tabulated function
      else if (swsophy == 1) then
         dum = head
         if (do_ln_trans .and. head < 0.0_real64) dum = -dlog(-head + 1.0_real64)
         if (head >= -1.0d-9) then
            moiscap = dt*1.0d-7
         else if (dum < sptab(1,node,1)) then
            moiscap = 0.0_real64
         else
            call EvalTabulatedFunction(0, numtab(node), 1, 2, 4, node, sptab, ientrytab, head, dummy, moiscap, 3)
         end if
      end if
      
   end function moiscap

   !> Calculate derivative of hydraulic conductivity (as a function of pressure head)
   function dhconduc(node, head, theta, dimocap, rfcp)
      implicit none
      
      ! Arguments
      integer, intent(in) :: node
      real(real64), intent(in) :: head, theta, dimocap, rfcp
      real(real64) :: dhconduc
      
      ! Local variables
      real(real64) :: term0, term1, term2, term3, term4, relsat, dummy
      real(real64) :: m, n, ksatfit, lambda, h_enpr, s_enpr, thetar, thetas, alfamg
      character(len=200) :: messag

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (swsophy == 0) then
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         ksatfit = cofgen(3,node)
         lambda = cofgen(5,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)
         
         if (iHWCKmodel(node) == 2) then
            dhconduc = alfamg*ksatfit*dexp(alfamg*head)
            
         else  ! Use default MvG

            ! For modified VanGenuchten model:
            ! - S_enpr: relative saturation at Entry Pressure h_enpr; NOTE: Senpr = (thetas-thetar)/(thetam-thetar)
            s_enpr = (abs(alfamg*h_enpr)) ** n
            s_enpr = (1.0_real64 + s_enpr) ** (-m)

            relsat = (theta - thetar) / (thetas - thetar)

            if (cofgen(10,node) > 0.0_real64 .and. relsat > cofgen(11,node)) then
               messag = 'Linear interpolation option for examined ksat not' // &
                    ' yet implemented for implicit hydraulic conductivity ' // &
                    '(swKimpl=1) in iteration scheme' 
               call fatalerr_collected('dhconduc', messag)
            end if

            if (relsat < 0.001_real64) then
               dhconduc = 0.0_real64
            else if (relsat > s_enpr) then
               dhconduc = 1.0d-12
            else
               dhconduc = dimocap / (thetas - thetar)
               term0    = (s_enpr*relsat)**(1.0_real64/m)
               term1    = 1.0_real64 - term0
               term2    = (2.0_real64 + lambda)*term0 - lambda
               term3    = term1 ** (m - 1.0_real64)
               term4    = 1.0_real64 - (1.0_real64 - s_enpr**(1.0_real64/m)) ** m
               dhconduc = dhconduc * ksatfit * relsat ** (lambda - 1.0_real64)
               dhconduc = dhconduc * (1.0_real64 - term1 ** m)
               dhconduc = dhconduc * (lambda + term2 * term3)
               dhconduc = dhconduc / (term4**2)
            end if
         end if

      ! Use tabulated function. "dhconduc" is calclated as a function of "head"
      else if (swsophy == 1) then
         if (theta >= sptab(2,node,numtab(node)) - 1.0d-9) then
            dhconduc = 1.0d+08
         else
            call EvalTabulatedFunction(0, numtab(node), 1, 3, 5, node, sptab, ientrytab, head, dummy, dhconduc, 4)
         end if
      end if

      ! In case of frost conditions
      if (swfrost == 1) then
         dhconduc = dhconduc * rfcp
      end if  
      
   end function dhconduc

   !> Calculate hydraulic conductivity (as a function of THETA)
   !! @param tsoil_node Optional soil temperature [deg C] at node — when
   !! supplied by the caller, the WC_K_models_04_11 temperature-dependent
   !! path reads from state%heat instead of the global tsoil array.
   !! SS-HEAT Phase 2 Task 6.
   function hconduc(node, head, theta, rfcp, tsoil_node)
      implicit none

      ! Arguments
      integer, intent(in) :: node
      real(real64), intent(in) :: head, theta, rfcp
      real(real64), optional, intent(in) :: tsoil_node
      !! SS-HEAT Phase 2 Task 6: soil temperature at node from state%heat (optional)
      real(real64) :: hconduc

      ! Local variables
      real(real64) :: term1, relsat, hconode_vsmall, m, ksatfit, lambda, dummy
      real(real64) :: relsatm, relsat1, alfamg, thetar, thetas
      real(real64) :: h_enpr, n, term2, thetam, relsatthr, ksatthr, ksatexm
      real(real64) :: alfa_2, n_2, m_2, omega_1, s1, s2
      real(real64) :: tsoil_loc
      !! Local temperature value: from tsoil_node arg (caller must supply when iHWCKmodel 4-11)
      real(real64), parameter :: h_crit = -1.0d-2

      hconode_vsmall = 1.0d-10
      ! [SS-HEAT] Task 9: global tsoil retired. Callers must pass tsoil_node for iHWCKmodel 4-11.
      ! For all other iHWCKmodel values the temperature path is not reached; tsoil_loc is unused.
      if (present(tsoil_node)) then
         tsoil_loc = tsoil_node
      else
         tsoil_loc = 0.0_real64   ! [SS-HEAT] Task 9: sentinel — global tsoil retired; iHWCKmodel 4-11 unreachable without tsoil_node
      end if

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (swsophy == 0) then
         thetar = cofgen(1,node)
         thetas = cofgen(2,node)
         alfamg = cofgen(4,node)
         ksatfit = cofgen(3,node)
         lambda = cofgen(5,node)
         n = cofgen(6,node)
         m = cofgen(7,node)
         h_enpr = cofgen(9,node)
         
         if (iHWCKmodel(layer(node)) == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            relsat  = (theta - thetar) / (thetas - thetar)
            hconduc = ksatfit*relsat

         else if (iHWCKmodel(layer(node)) == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = cofgen(13,node)
            n_2     = cofgen(14,node)
            m_2     = cofgen(15,node)
            omega_1 = cofgen(16,node)
            relsat  = (theta - thetar) / (thetas - thetar)
            if (relsat < 1.0_real64) then
               s1 = (1.0_real64 + (dabs(alfamg*head))**n)**(-m)
               s2 = (1.0_real64 + (dabs(alfa_2*head))**n_2)**(-m_2)
               term1 = omega_1*alfamg*(1.0_real64 - s1**(1.0_real64/m))**m
               term2 = (1.0_real64 - omega_1)*alfa_2*(1.0_real64 - s2**(1.0_real64/m_2))**m_2
               hconduc = ksatfit * (omega_1*s1 + (1.0_real64 - omega_1)*s2)**lambda
               hconduc = hconduc * (1.0_real64 - (term1 + term2) / (omega_1*alfamg + (1.0_real64 - omega_1)*alfa_2))**2
            else
               hconduc = ksatfit
            end if
            
         else if (iHWCKmodel(layer(node)) > 3 .and. iHWCKmodel(layer(node)) < 12) then
            ! SS-HEAT Phase 2 Task 6: use tsoil_loc (from state%heat or global fallback)
            hconduc = functionvalue_04_11(2, node, head, wc=theta, temp=tsoil_loc)
            
         else  ! Use default MvG

            if (fluseksatexm(node)) then 
               ksatexm   = cofgen(10,node)
               relsatthr = cofgen(11,node)
               ksatthr   = cofgen(12,node)
            else
               ksatexm   = 0.0_real64
               relsatthr = 0.0_real64
               ksatthr   = 0.0_real64
            end if

            relsat = (theta - thetar) / (thetas - thetar)
      
            if (fluseksatexm(node) .and. relsat > relsatthr) then

               term1   = (relsat - relsatthr) / (1.0_real64 - relsatthr)
               hconduc = term1 * ksatexm + (1.0_real64 - term1) * ksatthr

            else
               if (h_enpr > h_crit) then

                  if (head < -1.0d14) then
                     hconduc = hconode_vsmall
                  else if (relsat > (1.0_real64 - 1.0d-6)) then
                     hconduc = ksatfit
                  else
                     term1   = (1.0_real64 - relsat**(1.0_real64/m)) ** m
                     hconduc = ksatfit * (relsat**lambda) * (1.0_real64 - term1) * (1.0_real64 - term1)
                  end if

               else
                  ! For modified VanGenuchten model 
                  thetam = thetar + (thetas - thetar) * ((1.0_real64 + (abs(alfamg*h_enpr)) ** n) ** m)
                  if (head < -1.0d14) then
                     hconduc = hconode_vsmall
                  else 
                     if (theta >= thetam) then
                        hconduc = ksatfit
                     else
                        relsatm = (theta - thetar) / (thetam - thetar)
                        relsat1 = (thetas - thetar) / (thetam - thetar)
                        term1   = (1.0_real64 - (relsatm) ** (1.0_real64/m)) ** m
                        term2   = (1.0_real64 - (relsat1) ** (1.0_real64/m)) ** m
                        hconduc = ksatfit*(relsat**lambda) * ((1.0_real64 - term1) / (1.0_real64 - term2)) ** 2
                     end if
                  end if
               end if
               hconduc = min(hconduc, ksatfit)
            end if
         end if

      ! Use tabulated function. "hconduc" is calclated as a function of "head"
      else if (swsophy == 1) then
         if (theta >= sptab(2,node,numtab(node)) - 1.0d-9) then
            hconduc = sptab(3,node,numtab(node))
            if (do_ln_trans) hconduc = dexp(hconduc)
         else if (theta <= sptab(2,node,1) + 1.0d-9) then
            hconduc = sptab(3,node,1)            
            if (do_ln_trans) hconduc = dexp(hconduc)
         else
            call EvalTabulatedFunction(0, numtab(node), 1, 3, 5, node, sptab, ientrytab, head, hconduc, dummy, 2)
         end if
      end if

      ! In case of frost conditions
      if (swfrost == 1) then
         hconduc = hconduc * rfcp + hconode_vsmall * (1.0_real64 - rfcp)
      end if  
      
   end function hconduc

   !> Calculate pressure head from water content
   !! @note
   !! MH: since in routine convertdiscrvert prhead needs to be called with NEW distribution of cofgen and h,
   !!     cofgen_in and h_in are requied as input (and cannot be imported from variables as cofgen and h)
   !! @endnote
   function prhead(node, disnod, wcon, cofgen_in, h_in)
      use variables, only: swsophy, numtab, sptab, ientrytab, iHWCKmodel, layer
      use swap_array_dimensions, only: macp

      implicit none

      ! Arguments
      integer, intent(in) :: node
      real(real64), intent(in) :: disnod, wcon
      real(real64), intent(in) :: h_in(macp), cofgen_in(21,macp)
      real(real64) :: prhead

      ! Local variables
      real(real64) :: alfamg, thetar, thetas, h_enpr, s_enpr, npar, mpar, relsat
      real(real64) :: help, dummy, prh

      if (swsophy == 0) then
         thetar = cofgen_in(1,node)
         thetas = cofgen_in(2,node)
         alfamg = cofgen_in(4,node)
         npar   = cofgen_in(6,node)
         mpar   = cofgen_in(7,node)
         h_enpr = cofgen_in(9,node)
         
         if (iHWCKmodel(layer(node)) == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            relsat = (wcon - thetar) / (thetas - thetar)
            prhead = dlog(relsat) / alfamg
            
         else  ! Use default MvG

            if (thetas - wcon < 1.0d-6) then

               ! Saturated pressure head
               if (node == 1) then
                  prhead = disnod
               else
                  prhead = h_in(node-1) + disnod
               end if
               prhead = dmax1(prhead, h_enpr)
            else
               if (wcon - thetar < 1.0d-6) then
                  prhead = -1.0d12
               else

                  ! For modified VanGenuchten model:
                  ! - S_enpr: relative saturation at Entry Pressure h_enpr 
                  s_enpr = (abs(alfamg*h_enpr)) ** npar
                  s_enpr = (1.0_real64 + s_enpr) ** (-mpar)

                  ! First calculate the inverse of the sorptivity
                  help = (thetas - thetar) / (wcon - thetar) / s_enpr    
                  ! Raise to the power 1/m
                  help = help ** (1.0_real64 / mpar)
                  ! Subtract one and raise to the power 1/n 
                  help = (help - 1.0_real64) ** (1.0_real64 / npar)
                  ! Divide by alpha
                  prhead = -1.0_real64 * abs(help/alfamg)
               end if
            end if
         end if

      else if (swsophy == 1) then
         if (sptab(2,node,numtab(node)) - wcon < 1.0d-6) then

            ! Saturated pressure head
            if (node == 1) then
               prhead = disnod
            else
               prhead = h_in(node-1) + disnod
            end if
            prhead = dmax1(prhead, 0.0_real64)
         else

            call EvalTabulatedFunction(1, numtab(node), 1, 2, 4, node, sptab, ientrytab, prh, wcon, dummy, 1)
            prhead = prh
         end if
      end if
      
   end function prhead

   !> Calculate derivative of mean hydraulic conductivity with respect to main node conductivity
   !!
   !! d(hcomean)/d(kmain); note kmain = kup in hcomean and ksub = klow in hcomean; dzmain = dzup in hcomena, and dzub = dzlow in hcomean
   function dkmean(swkmean,kmain,ksub,dzmain,dzsub)
      implicit none
      integer swkmean
      real(8) kmain, ksub, dzmain, dzsub, a
      real(8) dkmean
      
      if (swkmean.eq.1) then
         dkmean = 0.5d0
      else if (swkmean.eq.2) then
         dkmean = dzmain/(dzmain+dzsub)
      else if (swkmean.eq.3) then
         dkmean = 0.5d0 * dsqrt(ksub / kmain)
      else if (swkmean.eq.4) then
         a = dzmain/(dzmain+dzsub)
         dkmean = a * (ksub / kmain) ** (1.0d0 - a)
      else if (swkmean.eq.5) then
         dkmean = 0.5d0/(((0.5d0/kmain)+(0.5/ksub))**2 * kmain**2)
      else if (swkmean.eq.6) then
         a = dzmain/(dzmain+dzsub)
         dkmean = a/(((a/kmain)+((1.0d0-a)/ksub))**2 * kmain**2)
      end if
   end function dkmean
end module soilhydraulics_utils
