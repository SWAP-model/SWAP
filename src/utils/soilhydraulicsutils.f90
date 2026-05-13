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
   ! [SS-SWC S-2.12B] cofgen/fluseksatexm retired from variables — bound via state%soilwater
   ! SS-TC TC-12: dt retired from only-list; bound via bind_tc_target (module-level pointer).
   use variables, only: swsophy, numtab, sptab, ientrytab, &
                        iHWCKmodel, layer, swfrost
   ! [SS-HEAT] Task 9: tsoil global retired; hconduc fallback removed (iHWCKmodel 4-11 path unreachable in regression)
   use doln
   use WC_K_models_04_11, only: functionvalue_04_11

   implicit none

   ! [SS-SWC S-2.12B] module-level pointers; bound once by bind_state_targets(state)
   real(real64), pointer :: cofgen(:,:) => null()
   logical,      pointer :: fluseksatexm(:) => null()
   ! SS-TC TC-12: module-level pointer for dt; bound once by bind_tc_target
   real(real64), pointer :: tc_dt_ptr => null()

   private
   public :: watcon, moiscap, hconduc, dhconduc, prhead, hcomean, dkmean
   public :: bind_state_targets
   public :: bind_tc_target  ! SS-TC TC-12

contains

   !> [SS-SWC S-2.12B] Bind module-level pointers to state%soilwater.
   !! Must be called once after state%soilwater is allocated (e.g. in SoilHydraulics(1)
   !! before any reader uses cofgen/fluseksatexm here).
   subroutine bind_state_targets(sw_cofgen_in, sw_fluseksatexm_in)
      real(real64), target, intent(in) :: sw_cofgen_in(:,:)
      logical,      target, intent(in) :: sw_fluseksatexm_in(:)
      cofgen       => sw_cofgen_in
      fluseksatexm => sw_fluseksatexm_in
   end subroutine bind_state_targets

   !> [SS-TC TC-12] Bind module-level tc_dt_ptr to state%timecontrol%dt.
   !! Must be called once after state is allocated (e.g. in swap.f90 init block
   !! alongside bind_state_targets). Used by moiscap() which has no state arg.
   subroutine bind_tc_target(tc_dt_in)
      real(real64), target, intent(in) :: tc_dt_in
      tc_dt_ptr => tc_dt_in
   end subroutine bind_tc_target


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
   !! [SS-GR-UTILS Task 5] New signature: explicit vg record + model + soilwater
   function watcon(head, vg, model, node, soilwater) result(theta)
      use soilwater_state_mod,  only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t
      implicit none

      real(real64),                  intent(in) :: head
      type(vanGenuchten_params_t),   intent(in) :: vg
      integer,                       intent(in) :: model       ! iHWCKmodel value for this node's layer
      integer,                       intent(in) :: node        ! for tabulated branch + functionvalue_04_11
      type(soilwater_state_t),       intent(in) :: soilwater   ! for tabulated branch (swsophy/sptab/numtab)
      real(real64) :: theta

      ! Local variables
      real(real64) :: h_enpr, help, m, n, s_enpr
      real(real64), parameter :: h_crit = -1.0d-2
      real(real64) :: h105, C105, a, b, dum, t105
      real(real64) :: alfamg, thetar, thetas, moiscap_loc
      real(real64) :: alfa_2, n_2, m_2, omega_1

      ! sptab(1,node,i):    h
      ! sptab(2,node,i):    Theta
      ! sptab(3,node,i):    K
      ! sptab(4,node,i):    dTheta / dh
      ! sptab(5,node,i):    dK / dh

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (soilwater%swsophy == 0) then

         thetar = vg%thetar
         thetas = vg%thetas
         alfamg = vg%alpha
         n = vg%npar
         m = vg%mpar
         h_enpr = vg%h_enpr

         if (model == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            theta = dmax1(1.0000001_real64*thetar, thetar + (thetas-thetar)*dexp(alfamg*head))

         else if (model == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = vg%alpha_2
            n_2     = vg%npar_2
            m_2     = vg%mpar_2
            omega_1 = vg%omega_1
            if (head < 0.0_real64) then
               theta = omega_1 / (1.0_real64 + (dabs(alfamg*head))**n)**m
               theta = theta + (1.0_real64 - omega_1) / (1.0_real64 + (dabs(alfa_2*head))**n_2)**m_2
               theta = thetar + (thetas - thetar) * theta
            else
               theta = thetas
            end if

         else if (model > 3 .and. model < 12) then
            theta = functionvalue_04_11(1, head, vg, model, &
                        soilwater%BiModal(soilwater%layer(node)), &
                        soilwater%NoVap(soilwater%layer(node)))

         else  ! Use default MvG

            if (h_enpr > h_crit) then

               if (head >= 0.0_real64) then
                  ! Saturated moisture content
                  theta = thetas
               else if (head > h_crit) then
                  help = (dabs(alfamg*h_crit))** n
                  help = (1.0_real64 + help) ** m
                  help = thetar + (thetas - thetar) / help
                  theta = help + (thetas - help) / (-h_crit) * (head - h_crit)
                  theta = min(theta, thetas)
               else
                  ! First compute |alpha * h| ** n
                  help = (dabs(alfamg*head)) ** n

                  ! Add 1 and raise to the power m
                  help = (1.0_real64 + help) ** m

                  ! Now compute theta
                  theta = thetar + (thetas - thetar) / help
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
                  theta = thetas + b*a*head / (1.0_real64 + a*head)

               else
                  ! First compute |alpha * h| ** n
                  help = (dabs(alfamg*head)) ** n

                  ! Add 1 and raise to the power m
                  help = (1.0_real64 + help) ** m

                  ! For modified VanGenuchten model:
                  ! - S_enpr: relative saturation at Entry Pressure h_enpr
                  s_enpr = (1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m

                  ! Now compute theta
                  theta = thetar + (thetas - thetar) / help * s_enpr
               end if
            end if
         end if

      ! Use tabulated function
      else if (soilwater%swsophy == 1) then
         dum = head
         if (do_ln_trans .and. head < 0.0_real64) dum = -dlog(-head + 1.0_real64)
         if (head >= -1.0d-9) then
            theta = soilwater%sptab(2,node,soilwater%numtab(node))
         else if (dum < soilwater%sptab(1,node,1)) then
            theta = soilwater%sptab(2,node,1)
         else
            call EvalTabulatedFunction(0, soilwater%numtab(node), 1, 2, 4, node, &
                                       soilwater%sptab, soilwater%ientrytab, head, theta, moiscap_loc, 1)
         end if
      end if

   end function watcon

   !> Calculate differential moisture capacity (as a function of pressure head)
   !! [SS-GR-UTILS Task 7] New signature: explicit vg record + model + dt + soilwater
   function moiscap(h, vg, model, dt, node, soilwater) result(c)
      use soilwater_state_mod,  only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t
      implicit none

      real(real64),                  intent(in) :: h
      type(vanGenuchten_params_t),   intent(in) :: vg
      integer,                       intent(in) :: model
      real(real64),                  intent(in) :: dt
      integer,                       intent(in) :: node
      type(soilwater_state_t),       intent(in) :: soilwater
      real(real64) :: c

      ! Local variables
      real(real64) :: h_enpr, m, n, alphah, s_enpr
      real(real64), parameter :: h_crit = -1.0d-2
      real(real64) :: alfamg, thetar, thetas
      real(real64) :: h105, C105, a, b, term1, term2, dum, t105, dummy
      real(real64) :: alfa_2, n_2, m_2, omega_1

      ! Use analytical expression
      if (soilwater%swsophy == 0) then

         thetar = vg%thetar
         thetas = vg%thetas
         alfamg = vg%alpha
         n      = vg%npar
         m      = vg%mpar
         h_enpr = vg%h_enpr

         if (model == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            c = alfamg * (thetas - thetar) * dexp(alfamg*h)

         else if (model == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = vg%alpha_2
            n_2     = vg%npar_2
            m_2     = vg%mpar_2
            omega_1 = vg%omega_1
            if (h < 0.0_real64) then
               c = omega_1*alfamg*n*m*(dabs(alfamg*h))**(n - 1.0_real64) * &
                    (1.0_real64 + (dabs(alfamg*h))**n)**(-1.0_real64 - m)
               c = c + (1.0_real64 - omega_1)*alfa_2*n_2*m_2* &
                    (dabs(alfa_2*h))**(n_2 - 1.0_real64) * &
                    (1.0_real64 + (dabs(alfa_2*h))**n_2)**(-1.0_real64 - m_2)
               c = (thetas - thetar) * c
            else
               c = 0.0_real64
            end if

         else if (model > 3 .and. model < 12) then
            c = functionvalue_04_11(3, h, vg, model, &
                    soilwater%BiModal(soilwater%layer(node)), &
                    soilwater%NoVap(soilwater%layer(node)))

         else  ! Use default MvG

            if (h_enpr > h_crit) then

               if (h >= 0.0_real64) then

                  c = dt * 1.0d-7  ! TC-12

               else if (h > h_crit) then

                  ! [SS-GR-UTILS Task 5] inlined: watcon at h_crit (swsophy==0, default MvG branch)
                  term1 = (dabs(alfamg*h_crit)) ** n
                  term1 = thetar + (thetas - thetar) / ((1.0_real64 + term1) ** m)
                  c = (thetas - term1) / (-h_crit)
               else

                  ! Use analytical evaluation of capacity
                  alphah = dabs(alfamg*h)

                  ! Compute |alpha * h| to the power n-1
                  term1 = alphah ** (n - 1.0_real64)

                  ! Compute |alpha*h| to the power n
                  term2 = term1 * alphah

                  ! Add one and raise to the power m+1
                  term2 = (1.0_real64 + term2) ** (m + 1.0_real64)

                  ! Divide theta-s minus theta-r by term2
                  term2 = (thetas - thetar) / term2

                  ! Calculate the differential moisture capacity
                  c = dabs(-1.0_real64 * n * m * alfamg * term2 * term1)
               end if
            else

               h105 = 1.05_real64 * h_enpr

               if (h >= h105) then
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
                  c = b*a / ((1.0_real64 + a*h)**2)

               else

                  alphah = dabs(alfamg*h)
                  term1 = alphah ** (n - 1.0_real64)
                  term2 = term1 * alphah
                  term2 = (1.0_real64 + term2) ** (m + 1.0_real64)
                  term2 = (thetas - thetar) / term2

                  ! For modified VanGenuchten model:
                  ! - S_enpr: relative saturation at Entry Pressure h_enpr
                  s_enpr = (1.0_real64 + (dabs(alfamg*h_enpr)) ** n) ** m

                  c = dabs(-1.0_real64 * n * m * alfamg*term2*term1)*s_enpr

               end if

            end if
         end if

         if (h > -1.0_real64 .and. c < (dt * 1.0d-7)) c = dt * 1.0d-7  ! TC-12

      ! Use tabulated function
      else if (soilwater%swsophy == 1) then
         dum = h
         if (do_ln_trans .and. h < 0.0_real64) dum = -dlog(-h + 1.0_real64)
         if (h >= -1.0d-9) then
            c = dt*1.0d-7  ! TC-12
         else if (dum < soilwater%sptab(1,node,1)) then
            c = 0.0_real64
         else
            call EvalTabulatedFunction(0, soilwater%numtab(node), 1, 2, 4, node, &
                                       soilwater%sptab, soilwater%ientrytab, h, dummy, c, 3)
         end if
      end if

   end function moiscap

   !> Calculate derivative of hydraulic conductivity (as a function of pressure head)
   !! [SS-GR-UTILS Task 6] New signature: explicit vg record + model + soilwater
   function dhconduc(h, theta, dimoca, rfcp, vg, model, node, soilwater) result(dkdh)
      use soilwater_state_mod,  only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t
      implicit none

      real(real64),                  intent(in) :: h, theta, dimoca, rfcp
      type(vanGenuchten_params_t),   intent(in) :: vg
      integer,                       intent(in) :: model
      integer,                       intent(in) :: node
      type(soilwater_state_t),       intent(in) :: soilwater
      real(real64) :: dkdh

      ! Local variables
      real(real64) :: term0, term1, term2, term3, term4, relsat, dummy
      real(real64) :: m, n, ksatfit, lambda, h_enpr, s_enpr, thetar, thetas, alfamg
      character(len=200) :: messag

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (soilwater%swsophy == 0) then
         thetar  = vg%thetar
         thetas  = vg%thetas
         alfamg  = vg%alpha
         ksatfit = vg%ksat
         lambda  = vg%lpar
         n       = vg%npar
         m       = vg%mpar
         h_enpr  = vg%h_enpr

         if (model == 2) then
            dkdh = alfamg*ksatfit*dexp(alfamg*h)

         else  ! Use default MvG

            ! For modified VanGenuchten model:
            ! - S_enpr: relative saturation at Entry Pressure h_enpr; NOTE: Senpr = (thetas-thetar)/(thetam-thetar)
            s_enpr = (abs(alfamg*h_enpr)) ** n
            s_enpr = (1.0_real64 + s_enpr) ** (-m)

            relsat = (theta - thetar) / (thetas - thetar)

            if (vg%ksatexm > 0.0_real64 .and. relsat > vg%relsatthr) then
               messag = 'Linear interpolation option for examined ksat not' // &
                    ' yet implemented for implicit hydraulic conductivity ' // &
                    '(swKimpl=1) in iteration scheme'
               call fatalerr_collected('dhconduc', messag)
            end if

            if (relsat < 0.001_real64) then
               dkdh = 0.0_real64
            else if (relsat > s_enpr) then
               dkdh = 1.0d-12
            else
               dkdh  = dimoca / (thetas - thetar)
               term0 = (s_enpr*relsat)**(1.0_real64/m)
               term1 = 1.0_real64 - term0
               term2 = (2.0_real64 + lambda)*term0 - lambda
               term3 = term1 ** (m - 1.0_real64)
               term4 = 1.0_real64 - (1.0_real64 - s_enpr**(1.0_real64/m)) ** m
               dkdh  = dkdh * ksatfit * relsat ** (lambda - 1.0_real64)
               dkdh  = dkdh * (1.0_real64 - term1 ** m)
               dkdh  = dkdh * (lambda + term2 * term3)
               dkdh  = dkdh / (term4**2)
            end if
         end if

      ! Use tabulated function. "dhconduc" is calclated as a function of "head"
      else if (soilwater%swsophy == 1) then
         if (theta >= soilwater%sptab(2,node,soilwater%numtab(node)) - 1.0d-9) then
            dkdh = 1.0d+08
         else
            call EvalTabulatedFunction(0, soilwater%numtab(node), 1, 3, 5, node, &
                                       soilwater%sptab, soilwater%ientrytab, h, dummy, dkdh, 4)
         end if
      end if

      ! In case of frost conditions
      if (swfrost == 1) then
         dkdh = dkdh * rfcp
      end if

   end function dhconduc

   !> Calculate hydraulic conductivity (as a function of THETA)
   !! [SS-GR-UTILS Task 6] New signature: explicit vg record + model + use_ksatexm + soilwater
   function hconduc(h, theta, rfcp, tsoil_node, vg, model, use_ksatexm, node, soilwater) result(k)
      use soilwater_state_mod,  only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t
      implicit none

      real(real64),                  intent(in) :: h, theta, rfcp, tsoil_node
      type(vanGenuchten_params_t),   intent(in) :: vg
      integer,                       intent(in) :: model       ! iHWCKmodel value
      logical,                       intent(in) :: use_ksatexm ! fluseksatexm flag
      integer,                       intent(in) :: node        ! for swsophy=1 / functionvalue_04_11
      type(soilwater_state_t),       intent(in) :: soilwater   ! for swsophy=1 (sptab/numtab)
      real(real64) :: k

      ! Local variables
      real(real64) :: term1, relsat, hconode_vsmall, m, ksatfit, lambda, dummy
      real(real64) :: relsatm, relsat1, alfamg, thetar, thetas
      real(real64) :: h_enpr, n, term2, thetam, relsatthr, ksatthr, ksatexm
      real(real64) :: alfa_2, n_2, m_2, omega_1, s1, s2
      real(real64), parameter :: h_crit = -1.0d-2

      hconode_vsmall = 1.0d-10

      ! Use analytical expression. "hconduc" is calclated as a function of "watcon"
      if (soilwater%swsophy == 0) then
         thetar  = vg%thetar
         thetas  = vg%thetas
         alfamg  = vg%alpha
         ksatfit = vg%ksat
         lambda  = vg%lpar
         n       = vg%npar
         m       = vg%mpar
         h_enpr  = vg%h_enpr

         if (model == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            relsat = (theta - thetar) / (thetas - thetar)
            k      = ksatfit*relsat

         else if (model == 3) then
            ! Bi-modal MvG relationships; basic form without air-entry h_enpr or h_crit
            alfa_2  = vg%alpha_2
            n_2     = vg%npar_2
            m_2     = vg%mpar_2
            omega_1 = vg%omega_1
            relsat  = (theta - thetar) / (thetas - thetar)
            if (relsat < 1.0_real64) then
               s1    = (1.0_real64 + (dabs(alfamg*h))**n)**(-m)
               s2    = (1.0_real64 + (dabs(alfa_2*h))**n_2)**(-m_2)
               term1 = omega_1*alfamg*(1.0_real64 - s1**(1.0_real64/m))**m
               term2 = (1.0_real64 - omega_1)*alfa_2*(1.0_real64 - s2**(1.0_real64/m_2))**m_2
               k     = ksatfit * (omega_1*s1 + (1.0_real64 - omega_1)*s2)**lambda
               k     = k * (1.0_real64 - (term1 + term2) / (omega_1*alfamg + (1.0_real64 - omega_1)*alfa_2))**2
            else
               k = ksatfit
            end if

         else if (model > 3 .and. model < 12) then
            ! SS-HEAT Phase 2 Task 6: use tsoil_node (from state%heat)
            k = functionvalue_04_11(2, h, vg, model, &
                    soilwater%BiModal(soilwater%layer(node)), &
                    soilwater%NoVap(soilwater%layer(node)), &
                    wc=theta, temp=tsoil_node)

         else  ! Use default MvG

            if (use_ksatexm) then
               ksatexm   = vg%ksatexm
               relsatthr = vg%relsatthr
               ksatthr   = vg%ksatthr
            else
               ksatexm   = 0.0_real64
               relsatthr = 0.0_real64
               ksatthr   = 0.0_real64
            end if

            relsat = (theta - thetar) / (thetas - thetar)

            if (use_ksatexm .and. relsat > relsatthr) then

               term1 = (relsat - relsatthr) / (1.0_real64 - relsatthr)
               k     = term1 * ksatexm + (1.0_real64 - term1) * ksatthr

            else
               if (h_enpr > h_crit) then

                  if (h < -1.0d14) then
                     k = hconode_vsmall
                  else if (relsat > (1.0_real64 - 1.0d-6)) then
                     k = ksatfit
                  else
                     term1 = (1.0_real64 - relsat**(1.0_real64/m)) ** m
                     k     = ksatfit * (relsat**lambda) * (1.0_real64 - term1) * (1.0_real64 - term1)
                  end if

               else
                  ! For modified VanGenuchten model
                  thetam = thetar + (thetas - thetar) * ((1.0_real64 + (abs(alfamg*h_enpr)) ** n) ** m)
                  if (h < -1.0d14) then
                     k = hconode_vsmall
                  else
                     if (theta >= thetam) then
                        k = ksatfit
                     else
                        relsatm = (theta - thetar) / (thetam - thetar)
                        relsat1 = (thetas - thetar) / (thetam - thetar)
                        term1   = (1.0_real64 - (relsatm) ** (1.0_real64/m)) ** m
                        term2   = (1.0_real64 - (relsat1) ** (1.0_real64/m)) ** m
                        k       = ksatfit*(relsat**lambda) * ((1.0_real64 - term1) / (1.0_real64 - term2)) ** 2
                     end if
                  end if
               end if
               k = min(k, ksatfit)
            end if
         end if

      ! Use tabulated function. "hconduc" is calclated as a function of "head"
      else if (soilwater%swsophy == 1) then
         if (theta >= soilwater%sptab(2,node,soilwater%numtab(node)) - 1.0d-9) then
            k = soilwater%sptab(3,node,soilwater%numtab(node))
            if (do_ln_trans) k = dexp(k)
         else if (theta <= soilwater%sptab(2,node,1) + 1.0d-9) then
            k = soilwater%sptab(3,node,1)
            if (do_ln_trans) k = dexp(k)
         else
            call EvalTabulatedFunction(0, soilwater%numtab(node), 1, 3, 5, node, &
                                       soilwater%sptab, soilwater%ientrytab, h, k, dummy, 2)
         end if
      end if

      ! In case of frost conditions
      if (swfrost == 1) then
         k = k * rfcp + hconode_vsmall * (1.0_real64 - rfcp)
      end if

   end function hconduc

   !> Calculate pressure head from water content
   !! @note
   !! New typed signature (SS-GR-UTILS Task 8). Optional vg_in allows
   !! callers to pass a tentative / locally-modified vg struct (e.g.
   !! soilgrid grid-redistribution, hysteresis parameter update).
   !! When vg_in is absent the node's vg_params entry from state is used.
   !! @endnote
   function prhead(disnod, wcon, h_in, model, node, soilwater, vg_in) result(h_out)
      use soilwater_state_mod,  only: soilwater_state_t
      use hydraulic_params_mod, only: vanGenuchten_params_t

      implicit none

      ! Arguments
      real(real64),                intent(in)           :: disnod
      real(real64),                intent(in)           :: wcon
      real(real64),                intent(in)           :: h_in(:)
      integer,                     intent(in)           :: model
      integer,                     intent(in)           :: node
      type(soilwater_state_t),     intent(in)           :: soilwater
      type(vanGenuchten_params_t), intent(in), optional :: vg_in
      real(real64) :: h_out

      ! Local variables
      type(vanGenuchten_params_t) :: vg
      real(real64) :: alfamg, thetar, thetas, h_enpr, s_enpr, npar, mpar, relsat
      real(real64) :: help, dummy, prh

      if (present(vg_in)) then
         vg = vg_in
      else
         vg = soilwater%vg_params(node)
      end if

      if (soilwater%swsophy == 0) then
         thetar = vg%thetar
         thetas = vg%thetas
         alfamg = vg%alpha
         npar   = vg%npar
         mpar   = vg%mpar
         h_enpr = vg%h_enpr

         if (model == 2) then
            ! Exponential relationships; special for testing against analytical solutions
            relsat = (wcon - thetar) / (thetas - thetar)
            h_out = dlog(relsat) / alfamg

         else  ! Use default MvG

            if (thetas - wcon < 1.0d-6) then

               ! Saturated pressure head
               if (node == 1) then
                  h_out = disnod
               else
                  h_out = h_in(node-1) + disnod
               end if
               h_out = dmax1(h_out, h_enpr)
            else
               if (wcon - thetar < 1.0d-6) then
                  h_out = -1.0d12
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
                  h_out = -1.0_real64 * abs(help/alfamg)
               end if
            end if
         end if

      else if (soilwater%swsophy == 1) then
         if (soilwater%sptab(2,node,soilwater%numtab(node)) - wcon < 1.0d-6) then

            ! Saturated pressure head
            if (node == 1) then
               h_out = disnod
            else
               h_out = h_in(node-1) + disnod
            end if
            h_out = dmax1(h_out, 0.0_real64)
         else

            call EvalTabulatedFunction(1, soilwater%numtab(node), 1, 2, 4, node, &
                                       soilwater%sptab, soilwater%ientrytab, prh, wcon, dummy, 1)
            h_out = prh
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
