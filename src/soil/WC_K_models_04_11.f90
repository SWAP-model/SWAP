module WC_K_models_04_11
   use error_mod, only: fatalerr_collected
   use iso_fortran_env, only: real64
   use hydraulic_params_mod, only: vanGenuchten_params_t

implicit none

!***************************************************************************************************************************
! Soil hydraulic property library (Mualem-van Genuchten + PDI model variants, models 4-11).
!
! Given a pressure head h (cm), functionvalue_04_11 returns one of:
!   iType=1 : water content            theta(h)
!   iType=2 : hydraulic conductivity   K(h)
!   iType=3 : differential moisture capacity  C(h) = d(theta)/dh
!
! for the following parameterisations (selected by `model`):
!   4  MvG uni-modal              5  MvG uni-modal, saturation-corrected (_s)
!   6  MvG bi-modal               7  MvG bi-modal, saturation-corrected
!   8  PDI uni-modal              9  PDI uni-modal, saturation-corrected
!   10 PDI bi-modal               11 PDI bi-modal, saturation-corrected
! The PDI variants add an adsorptive water term (Sad) plus film and vapour
! conductivity contributions on top of the capillary (MvG) term.
!
! This routine is called per node per Newton iteration via watcon/hconduc/moiscap.
!
! REENTRANT / THREAD-SAFE: all curve parameters travel through `vg`
! (vanGenuchten_params_t, the dispatcher's input) and the model flags
! is_bimodal/no_vap are passed explicitly; intermediate "help" quantities are
! function-local. The module holds NO mutable state.
!***************************************************************************************************************************

! functionvalue_04_11 is the sole public symbol
private
public    :: functionvalue_04_11

contains

function functionvalue_04_11(iType, h, vg, model, is_bimodal, no_vap, wc, temp) result(val)
   implicit none

   integer,                     intent(in)           :: iType
   real(real64),                intent(in)           :: h
   type(vanGenuchten_params_t), intent(in)           :: vg
   integer,                     intent(in)           :: model       ! iHWCKmodel(layer(iNode))
   logical,                     intent(in)           :: is_bimodal  ! BiModal(layer(iNode))
   logical,                     intent(in)           :: no_vap      ! NoVap(layer(iNode))
   real(real64),                intent(in), optional :: wc, temp
   real(real64) :: val

! local
real(8), parameter :: dummy = 0.0d0

select case (iType)
case (1)
   select case (model)
      case (4);  val = WC_MvG (h, vg)
      case (5);  val = WC_MvG_s (h, vg)
      case (6);  val = WC_MvG_2 (h, vg)
      case (7);  val = WC_MvG_2_s (h, vg)
      case (8);  val = WC_PDI (h, vg, is_bimodal)
      case (9);  val = WC_PDI_s (h, vg, is_bimodal)
      case (10); val = WC_PDI_2 (h, vg, is_bimodal)
      case (11); val = WC_PDI_2_s (h, vg, is_bimodal)
   end select

case (2)
   select case (model)
      case (4);  val = K_MvG (h, vg)
      case (5);  val = K_MvG_s (h, vg)
      case (6);  val = K_MvG_2 (h, vg)
      case (7);  val = K_MvG_2_s (h, vg)
      case (8)
         if (no_vap) then
            val = K_PDI (h,dummy,dummy, vg, is_bimodal, no_vap)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr_collected ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            val = K_PDI (h,wc,temp, vg, is_bimodal, no_vap)
         end if

      case (9)
         if (no_vap) then
            val = K_PDI_s (h,dummy,dummy, vg, is_bimodal, no_vap)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr_collected ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            val = K_PDI_s (h,wc,temp, vg, is_bimodal, no_vap)
         end if

      case (10)
         if (no_vap) then
            val = K_PDI_2 (h,dummy,dummy, vg, is_bimodal, no_vap)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr_collected ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            val = K_PDI_2 (h,wc,temp, vg, is_bimodal, no_vap)
         end if

      case (11)
         if (no_vap) then
            val = K_PDI_2_s (h,dummy,dummy, vg, is_bimodal, no_vap)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr_collected ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            val = K_PDI_2_s (h,wc,temp, vg, is_bimodal, no_vap)
         end if

   end select

case (3)
   select case (model)
      case (4);  val = C_MvG (h, vg)
      case (5);  val = C_MvG_s (h, vg)
      case (6);  val = C_MvG_2 (h, vg)
      case (7);  val = C_MvG_2_s (h, vg)
      case (8);  val = C_PDI (h, vg, is_bimodal)
      case (9);  val = C_PDI_s (h, vg, is_bimodal)
      case (10); val = C_PDI_2 (h, vg, is_bimodal)
      case (11); val = C_PDI_2_s (h, vg, is_bimodal)
   end select

case default
   call fatalerr_collected ('functionvalue_04_11', 'Illegal iType; allowed values [1,2,3]')

end select
return
end function functionvalue_04_11

!************** P R I V A T E **********************************************************************************************
! for each of WC, K and C calculations 8 functions are given: MvG-1, MvG-1-s, MvG-2, MvG-2-s, PDI-1, PDI-1-s, PDI-2, PDI-2-s
! There are 7 help functions (Gamm1, Gamm2, b, Sad, Kvap_func, C1, C2)
!***************************************************************************************************************************

function Gamma1 (h, vg)
real(8) :: h, Gamma1
type(vanGenuchten_params_t), intent(in) :: vg
Gamma1 = (1.0d0 + (vg%alpha*h)**(vg%npar))**(-vg%mpar)
end function Gamma1

function Gamma2 (h, vg)
real(8) :: h, Gamma2
type(vanGenuchten_params_t), intent(in) :: vg
Gamma2 = (1.0d0 + (vg%alpha_2*h)**(vg%npar_2))**(-vg%mpar_2)
end function Gamma2

function b (vg, is_bimodal)
real(8) :: b
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in) :: is_bimodal
real(8) :: nn
if (.not. is_bimodal) then
   b  = 0.1d0 + 0.2d0/vg%npar**2 * (1.0d0 - dexp(-((vg%thetar/(vg%thetas-vg%thetar))**2)))
else
   nn = vg%npar
   if (vg%alpha_2 > vg%alpha) nn = vg%npar_2
   b  = 0.1d0 + 0.2d0/nn**2 * (1.0d0 - dexp(-((vg%thetar/(vg%thetas-vg%thetar))**2)))
end if
end function b

function Sad (h, vg, is_bimodal)
real(8) :: Sad, h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in) :: is_bimodal
real(8) :: x, xa, x0, bb
xa = dlog10(vg%ha)
x0 = dlog10(vg%h0)
x  = dlog10(h)
bb = b(vg, is_bimodal)
Sad = 1.0d0 + (x - xa + bb*dlog(1.0d0 + dexp((xa-x)/bb))) / (xa - x0)
end function Sad

function dSad_dh (h, vg, is_bimodal)
real(8) :: dSad_dh, h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in) :: is_bimodal
real(8) :: x, xa, x0, bb
xa = dlog10(vg%ha)
x0 = dlog10(vg%h0)
x  = dlog10(h)
bb = b(vg, is_bimodal)
dSad_dh = -1.0d0/(h*dlog(10.0d0)*(xa-x0) * (1.0d0 + dexp((xa-x)/bb)))
end function dSad_dh

function Kvap_func (WC, h, Temp, vg)
! Temp in degree Celsius
real(8), intent(in)    :: WC, h, Temp
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: Kvap_func, ksi, D, Hr
real(8)                :: fKvap, Da, MgRT, Rho_sv
real(8), parameter     :: p = 7.0d0/3.0d0
! M : molecular weight of water;  kg/mol
! g : gravitational acceleration; m/s2
! R : universal gas constant;     J/mol/K; J = kg.m2/s2
real(8), parameter     :: MgR = 0.018015d0*9.81d0/8.314d0  ! (kg/mol * m/s2) / (kg.m2/s2)
real(8), parameter     :: Rho_w = 1000.0d0                 ! density of water; kg/m3

MgRT      = MgR/(Temp+273.15d0)
Da        = 2.14d-5*((Temp+273.15d0)/273.15d0)**2                                ! diffusivity of water vapor in air; m2/s
Rho_sv    = 1.0d-3*dexp(31.3716d0 - 6014.79d0/Temp - 7.92495d-3*Temp)/Temp       ! saturated vapor density; kg/m3
fKvap     = Rho_sv/Rho_w * MgRT
ksi       = (vg%thetas-WC)**p/vg%thetas**2
D         = ksi*(vg%thetas-WC)*Da
Hr        = dexp(h/100.0d0*MgRT)     ! h must be in m, thus h (cm) is idvided by 100
Kvap_func = fKvap*D*Hr
end function Kvap_func

function C1 (h, vg)
real(8) :: h, C1
type(vanGenuchten_params_t), intent(in) :: vg
C1 = vg%alpha*vg%npar*vg%mpar*(vg%alpha*dabs(h))**(vg%npar-1.0d0)*(1.0d0+(vg%alpha*dabs(h))**vg%npar)**(-vg%mpar-1.0d0)
end function C1

function C2 (h, vg)
real(8) :: h, C2
type(vanGenuchten_params_t), intent(in) :: vg
C2 = vg%alpha_2*vg%npar_2*vg%mpar_2*(vg%alpha_2*dabs(h))**(vg%npar_2-1.0d0)*(1.0d0+(vg%alpha_2*dabs(h))**vg%npar_2)**(-vg%mpar_2-1.0d0)
end function C2

function WC_MvG (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: WC_MvG
real(8)                :: Scap

if (h >= 0.0d0) then
   WC_MvG = vg%thetas
else
   Scap = Gamma1 (dabs(h), vg)
   WC_MvG = vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_MvG

function WC_MvG_s (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: WC_MvG_s
real(8)                :: Gam01, Gamh1, Scap

if (h >= 0.0d0) then
   WC_MvG_s = vg%thetas
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Gamh1 = Gamma1 (dabs(h), vg)
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   WC_MvG_s = vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_MvG_s

function WC_MvG_2 (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: WC_MvG_2
real(8)                :: Gamh1, Gamh2, Scap

if (h >= 0.0d0) then
   WC_MvG_2 = vg%thetas
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   Scap = vg%omega_1 * Gamh1 + vg%omega_2 * Gamh2
   WC_MvG_2 = vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_MvG_2

function WC_MvG_2_s (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: WC_MvG_2_s
real(8)                :: Gam01, Gamh1, Gam02, Gamh2, Scap

if (h >= 0.0d0) then
   WC_MvG_2_s = vg%thetas
else
   Gam01 = vg%omega_1 * Gamma1 (dabs(vg%h0), vg)
   Gamh1 = vg%omega_1 * Gamma1 (dabs(h), vg)
   Gam02 = vg%omega_2 * Gamma2 (dabs(vg%h0), vg)
   Gamh2 = vg%omega_2 * Gamma2 (dabs(h), vg)
   Scap = (Gamh1 + Gamh2 - Gam01 - Gam02) / (1.0d0 - Gam01 - Gam02)
   WC_MvG_2_s = vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_MvG_2_s

function K_MvG (h, vg)
implicit none
real(8), intent(in)          :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                      :: K_MvG
real(8)                      :: Scap

if (h >= 0.0d0) then
   K_MvG = vg%ksat
else
   Scap = Gamma1 (dabs(h), vg)
   K_MvG = vg%ksat * Scap**vg%lpar * (1.0d0 - (1.0d0 - Scap**(1.0d0/vg%mpar))**vg%mpar)**2
end if

end function K_MvG

function K_MvG_s (h, vg)
implicit none
real(8), intent(in)          :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                      :: K_MvG_s
real(8)                      :: Gam01, Gamh1, Scap

if (h >= 0.0d0) then
   K_MvG_s = vg%ksat
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Gamh1 = Gamma1 (dabs(h), vg)
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   K_MvG_s = vg%ksat*Scap**vg%lpar * (1.0d0 - ((1.0d0-Gamh1**(1.0d0/vg%mpar))/(1.0d0-Gam01**(1.0d0/vg%mpar)))**vg%mpar)**2
end if

end function K_MvG_s

function K_MvG_2 (h, vg)
implicit none
real(8), intent(in)          :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                      :: K_MvG_2
real(8)                      :: Gamh1, Gamh2, t1, t2, t3

if (h >= 0.0d0) then
   K_MvG_2 = vg%ksat
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   t1 = (vg%omega_1*Gamh1 + vg%omega_2*Gamh2)**vg%lpar
   t2 = vg%omega_1*vg%alpha*(1.0d0-Gamh1**(1.0d0/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gamh2**(1.0d0/vg%mpar_2))**vg%mpar_2
   t3 = vg%omega_1*vg%alpha + vg%omega_2*vg%alpha_2
   K_MvG_2 = vg%ksat*t1*(1.0d0-t2/t3)**2
end if

end function K_MvG_2

function K_MvG_2_s (h, vg)
implicit none
real(8), intent(in)          :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                      :: K_MvG_2_s
real(8)                      :: Gam01, Gamh1, Gam02, Gamh2, Scap1, Scap2, t1, t2, t3

if (h >= 0.0d0) then
   K_MvG_2_s = vg%ksat
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Gamh1 = Gamma1 (dabs(h), vg)
   Gam02 = Gamma2 (dabs(vg%h0), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   Scap1 = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Scap2 = (Gamh2 - Gam02) / (1.0d0 - Gam02)
   t1 = (vg%omega_1*Scap1 + vg%omega_2*Scap2)**vg%lpar
   t2 = vg%omega_1*vg%alpha*(1.0d0-Gamh1**(1.0d0/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gamh2**(1.0d0/vg%mpar_2))**vg%mpar_2
   t3 = vg%omega_1*vg%alpha*(1.0d0-Gam01**(1.0d0/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gam02**(1.0d0/vg%mpar_2))**vg%mpar_2
   K_MvG_2_s = vg%ksat*t1*(1.0d0-t2/t3)**2
end if

end function K_MvG_2_s

function WC_PDI (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: WC_PDI
real(8)                :: Scap

if (h >= 0.0d0) then
   WC_PDI = vg%thetas
else
   Scap = Gamma1 (dabs(h), vg)
   WC_PDI = Sad(dabs(h), vg, is_bimodal)*vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_PDI

function WC_PDI_s (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: WC_PDI_s
real(8)                :: Gamh1, Gam01, Scap

if (h >= 0.0d0) then
   WC_PDI_s = vg%thetas
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   WC_PDI_s = Sad(dabs(h), vg, is_bimodal)*vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_PDI_s

function WC_PDI_2 (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: WC_PDI_2
real(8)                :: Gamh1, Gamh2, Scap

if (h >= 0.0d0) then
   WC_PDI_2 = vg%thetas
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   Scap = vg%omega_1 * Gamh1 + vg%omega_2 * Gamh2
   WC_PDI_2 = Sad(dabs(h), vg, is_bimodal)*vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_PDI_2

function WC_PDI_2_s (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: WC_PDI_2_s
real(8)                :: Gam01, Gamh1, Gam02, Gamh2, Scap

if (h >= 0.0d0) then
   WC_PDI_2_s = vg%thetas
else
   Gam01 = vg%omega_1*Gamma1 (dabs(vg%h0), vg)
   Gamh1 = vg%omega_1*Gamma1 (dabs(h), vg)
   Gam02 = vg%omega_2*Gamma2 (dabs(vg%h0), vg)
   Gamh2 = vg%omega_2*Gamma2 (dabs(h), vg)
   Scap = (Gamh1 + Gamh2 - (Gam01 + Gam02)) / (1.0d0 - (Gam01 + Gam02))
   WC_PDI_2_s = Sad(dabs(h), vg, is_bimodal)*vg%thetar + Scap*(vg%thetas-vg%thetar)
end if

end function WC_PDI_2_s

function K_PDI (h, WC, Temp, vg, is_bimodal, no_vap)
implicit none
real(8), intent(in)          :: h, WC, Temp
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)          :: is_bimodal, no_vap
real(8)                      :: K_PDI
real(8)                      :: Scap, Kcap, Kfilm, Kvap
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI = vg%ksat
else
   Scap = Gamma1 (dabs(h), vg)
   Kcap = Scap**vg%lpar*(1.0d0 - (1.0d0 - Scap**(1.0d0/vg%mpar))**vg%mpar)**2
   Kfilm = (vg%h0/vg%ha)**(vg%apar*(1.0d0-Sad(dabs(h), vg, is_bimodal)))
   if (no_vap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp, vg) * Conv
   end if
   K_PDI = vg%ksat*((1.0d0-vg%omega_k)*Kcap + vg%omega_k*Kfilm) + Kvap
end if

end function K_PDI

function K_PDI_s (h, WC, Temp, vg, is_bimodal, no_vap)
implicit none
real(8), intent(in)          :: h, WC, Temp
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)          :: is_bimodal, no_vap
real(8)                      :: K_PDI_s
real(8)                      :: Gamh1, Gam01, Scap, Kcap, Kfilm, Kvap
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_s = vg%ksat
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Kcap = Scap**vg%lpar*(1.0d0 - ((1.0d0-Gamh1**(1.0d0/vg%mpar))/(1.0d0-Gam01**(1.0d0/vg%mpar)))**vg%mpar)**2
   Kfilm = (vg%h0/vg%ha)**(vg%apar*(1.0d0-Sad(dabs(h), vg, is_bimodal)))
   if (no_vap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp, vg) * Conv
   end if
   K_PDI_s = vg%ksat*((1.0d0-vg%omega_k)*Kcap + vg%omega_k*Kfilm) + Kvap
end if

end function K_PDI_s

function K_PDI_2 (h, WC, Temp, vg, is_bimodal, no_vap)
implicit none
real(8), intent(in)          :: h, WC, Temp
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)          :: is_bimodal, no_vap
real(8)                      :: K_PDI_2
real(8)                      :: Gamh1, Gamh2, t1, t2, t3, Kcap, Kfilm, Kvap
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_2 = vg%ksat
else
   Gamh1 = Gamma1 (dabs(h), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   t1 = (vg%omega_1*Gamh1 + vg%omega_2*Gamh2)**vg%lpar
   t2 = vg%omega_1*vg%alpha*(1.0d0-Gamh1**(1/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gamh2**(1/vg%mpar_2))**vg%mpar_2
   t3 = vg%omega_1*vg%alpha + vg%omega_2*vg%alpha_2
   Kcap = t1*(1.0d0-t2/t3)**2
   Kfilm = (vg%h0/vg%ha)**(vg%apar*(1.0d0-Sad(dabs(h), vg, is_bimodal)))
   if (no_vap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp, vg) * Conv
   end if
   K_PDI_2 = vg%ksat*((1.0d0-vg%omega_k)*Kcap + vg%omega_k*Kfilm) + Kvap
end if

end function K_PDI_2

function K_PDI_2_s (h, WC, Temp, vg, is_bimodal, no_vap)
implicit none
real(8), intent(in)          :: h, WC, Temp
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)          :: is_bimodal, no_vap
real(8)                      :: K_PDI_2_s
real(8)                      :: Gam01, Gamh1, Gam02, Gamh2, Scap1, Scap2, t1, t2, t3, Kcap, Kfilm, Kvap
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_2_s = vg%ksat
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Gamh1 = Gamma1 (dabs(h), vg)
   Gam02 = Gamma2 (dabs(vg%h0), vg)
   Gamh2 = Gamma2 (dabs(h), vg)
   Scap1 = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Scap2 = (Gamh2 - Gam02) / (1.0d0 - Gam02)
   t1 = (vg%omega_1*Scap1 + vg%omega_2*Scap2)**vg%lpar
   t2 = vg%omega_1*vg%alpha*(1.0d0-Gamh1**(1.0d0/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gamh2**(1.0d0/vg%mpar_2))**vg%mpar_2
   t3 = vg%omega_1*vg%alpha*(1.0d0-Gam01**(1.0d0/vg%mpar))**vg%mpar + vg%omega_2*vg%alpha_2*(1.0d0-Gam02**(1.0d0/vg%mpar_2))**vg%mpar_2
   Kcap = t1*(1.0d0-t2/t3)**2
   Kfilm = (vg%h0/vg%ha)**(vg%apar*(1.0d0-Sad(dabs(h), vg, is_bimodal)))
   if (no_vap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp, vg) * Conv
   end if
   K_PDI_2_s = vg%ksat*((1.0d0-vg%omega_k)*Kcap + vg%omega_k*Kfilm) + Kvap
end if

end function K_PDI_2_s

function C_MvG (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: C_MvG
if (h >= 0.0d0) then
   C_MvG = 0.0d0
else
   C_MvG = (vg%thetas-vg%thetar) * C1(h, vg)
end if
end function C_MvG

function C_MvG_s (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: C_MvG_s
real(8)                :: Gam01
if (h >= 0.0d0) then
   C_MvG_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   C_MvG_s = (vg%thetas-vg%thetar)/(1.0d0-Gam01) * C1(h, vg)
end if
end function C_MvG_s

function C_MvG_2 (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: C_MvG_2
if (h >= 0.0d0) then
   C_MvG_2 = 0.0d0
else
   C_MvG_2 = (vg%thetas-vg%thetar)*(vg%omega_1*C1(h, vg) + vg%omega_2*C2(h, vg))
end if
end function C_MvG_2

function C_MvG_2_s (h, vg)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
real(8)                :: C_MvG_2_s
real(8)                :: Gam01, Gam02
if (h >= 0.0d0) then
   C_MvG_2_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   Gam02 = Gamma2 (dabs(vg%h0), vg)
   C_MvG_2_s = (vg%thetas-vg%thetar)*(vg%omega_1*C1(h, vg)/(1.0d0-Gam01) + vg%omega_2*C2(h, vg)/(1.0d0-Gam02))
end if
end function C_MvG_2_s

function C_PDI (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: C_PDI
real(8)                :: SSad
if (h >= 0.0d0) then
   C_PDI = 0.0d0
else
   SSad = dSad_dh (dabs(h), vg, is_bimodal)
   C_PDI = (vg%thetas-vg%thetar)*C1(h, vg) + vg%thetar*SSad
end if
end function C_PDI

function C_PDI_s (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: C_PDI_s
real(8)                :: Gam01, SSad
if (h >= 0.0d0) then
   C_PDI_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(vg%h0), vg)
   SSad  = dSad_dh (dabs(h), vg, is_bimodal)
   C_PDI_s = (vg%thetas-vg%thetar)/(1.0d0-Gam01)*C1(h, vg) + vg%thetar*SSad
end if
end function C_PDI_s

function C_PDI_2 (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: C_PDI_2
real(8)                :: SSad
if (h >= 0.0d0) then
   C_PDI_2 = 0.0d0
else
   SSad = dSad_dh (dabs(h), vg, is_bimodal)
   C_PDI_2 = (vg%thetas-vg%thetar)*(vg%omega_1*C1(h, vg) + vg%omega_2*C2(h, vg)) + vg%thetar*SSad
end if
end function C_PDI_2

function C_PDI_2_s (h, vg, is_bimodal)
real(8), intent(in)    :: h
type(vanGenuchten_params_t), intent(in) :: vg
logical, intent(in)    :: is_bimodal
real(8)                :: C_PDI_2_s, Gam0
real(8)                :: Gam01, Gam02, SSad
if (h >= 0.0d0) then
   C_PDI_2_s = 0.0d0
else
   Gam01 = vg%omega_1*Gamma1 (dabs(vg%h0), vg)
   Gam02 = vg%omega_2*Gamma2 (dabs(vg%h0), vg)
   Gam0  = Gam01 + Gam02
   SSad  = dSad_dh (dabs(h), vg, is_bimodal)
   C_PDI_2_s = (vg%thetas-vg%thetar)*(vg%omega_1*C1(h, vg)/(1.0d0-Gam0) + vg%omega_2*C2(h, vg)/(1.0d0-Gam0)) + vg%thetar*SSad
end if
end function C_PDI_2_s

end module WC_K_models_04_11
