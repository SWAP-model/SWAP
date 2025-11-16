module WC_K_models_04_11

use variables, only: cofgen, iHWCKmodel, BiModal, NoVap, layer
implicit none

! curve parameters
real(8)  :: WCr, WCs, Alpha1, Npar1, Mpar1, Alpha2, Npar2, Mpar2, Lpar, Ksat
real(8)  :: ha, h0, Apar
real(8)  :: Omega1, Omega2, OmegaK

! local help variables
real(8) :: nn, x, xa, x0, bb, t1, t2, t3, SSad
real(8) :: Scap, Scap1, Scap2
real(8) :: Gamh1, Gamh2, Gam01, Gam02
real(8) :: Kcap, Kfilm, Kvap

logical :: L_BiModal, L_NoVap

! all functions are private, except functionvalue_04_11
private
public    :: functionvalue_04_11

contains

function functionvalue_04_11 (iType, iNode, h, wc, temp)
integer, intent(in)             :: iType, iNode
real(8), intent(in)             :: h
real(8), intent(in), optional   :: wc, temp
real(8)                         :: functionvalue_04_11
! local
integer                         :: iModel
real(8), parameter              :: dummy = 0.0d0

! local information
iModel    = iHWCKmodel(layer(iNode))
L_BiModal = BiModal(layer(iNode))
L_NoVap   = NoVap(layer(iNode))

! for all cases
WCr    = cofgen(1,iNode)
WCs    = cofgen(2,iNode)
Ksat   = cofgen(3,iNode)
Alpha1 = cofgen(4,iNode)
Lpar   = cofgen(5,iNode)
Npar1  = cofgen(6,iNode)
Mpar1  = cofgen(7,iNode)

select case (iModel)

   case (4)
      ! no additional settings needed
      continue

   case (5)
      h0 = cofgen(18,iNode)

   case (6)
      Alpha2 = cofgen(13,iNode)
      Npar2  = cofgen(14,iNode)
      Mpar2  = cofgen(15,iNode)
      Omega1 = cofgen(16,iNode)
      Omega2 = cofgen(17,iNode)

   case (7)
      Alpha2 = cofgen(13,iNode)
      Npar2  = cofgen(14,iNode)
      Mpar2  = cofgen(15,iNode)
      Omega1 = cofgen(16,iNode)
      Omega2 = cofgen(17,iNode)
      h0     = cofgen(18,iNode)

   case (8,9)
      h0      = cofgen(18,iNode)
      ha      = cofgen(19,iNode)
      Apar    = cofgen(20,iNode)
      OmegaK  = cofgen(21,iNode)

   case (10,11)
      Alpha2 = cofgen(13,iNode)
      Npar2  = cofgen(14,iNode)
      Mpar2  = cofgen(15,iNode)
      Omega1 = cofgen(16,iNode)
      Omega2 = cofgen(17,iNode)
      h0     = cofgen(18,iNode)
      ha     = cofgen(19,iNode)
      Apar   = cofgen(20,iNode)
      OmegaK = cofgen(21,iNode)

end select

select case (iType)
case (1)
   select case (iModel)
      case (4);  functionvalue_04_11 = WC_MvG (h)
      case (5);  functionvalue_04_11 = WC_MvG_s (h)
      case (6);  functionvalue_04_11 = WC_MvG_2 (h)
      case (7);  functionvalue_04_11 = WC_MvG_2_s (h)
      case (8);  functionvalue_04_11 = WC_PDI (h)
      case (9);  functionvalue_04_11 = WC_PDI_s (h)
      case (10); functionvalue_04_11 = WC_PDI_2 (h)
      case (11); functionvalue_04_11 = WC_PDI_2_s (h)
   end select

case (2)
   select case (iModel)
      case (4);  functionvalue_04_11 = K_MvG (h)
      case (5);  functionvalue_04_11 = K_MvG_s (h)
      case (6);  functionvalue_04_11 = K_MvG_2 (h)
      case (7);  functionvalue_04_11 = K_MvG_2_s (h)
      case (8)
         if (L_NoVap) then
            functionvalue_04_11 = K_PDI (h,dummy,dummy)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            functionvalue_04_11 = K_PDI (h,wc,temp)
         end if

      case (9)
         if (L_NoVap) then
            functionvalue_04_11 = K_PDI_s (h,dummy,dummy)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            functionvalue_04_11 = K_PDI_s (h,wc,temp)
         end if
      
      case (10)
         if (L_NoVap) then
            functionvalue_04_11 = K_PDI_2 (h,dummy,dummy)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            functionvalue_04_11 = K_PDI_2 (h,wc,temp)
         end if
      
      case (11)
         if (L_NoVap) then
            functionvalue_04_11 = K_PDI_2_s (h,dummy,dummy)
         else
            ! both wc and temp must be present as input
            if (.not.present(wc) .or. .not.present(temp)) call fatalerr ('functionvalue_04_11','For Kvap both WC and TEMP must be given as arguments')
            functionvalue_04_11 = K_PDI_2_s (h,wc,temp)
         end if
            
   end select

case (3)
   select case (iModel)
      case (4);  functionvalue_04_11 = C_MvG (h)
      case (5);  functionvalue_04_11 = C_MvG_s (h)
      case (6);  functionvalue_04_11 = C_MvG_2 (h)
      case (7);  functionvalue_04_11 = C_MvG_2_s (h)
      case (8);  functionvalue_04_11 = C_PDI (h)
      case (9);  functionvalue_04_11 = C_PDI_s (h)
      case (10); functionvalue_04_11 = C_PDI_2 (h)
      case (11); functionvalue_04_11 = C_PDI_2_s (h)
   end select

case default
   call fatalerr ('functionvalue_04_11', 'Illegal iType; allowed values [1,2,3]')

end select
return
end function functionvalue_04_11

!************** P R I V A T E **********************************************************************************************
! for each of WC, K and C calculations 8 functions are given: MvG-1, MvG-1-s, MvG-2, MvG-2-s, PDI-1, PDI-1-s, PDI-2, PDI-2-s
! There are 7 help functions (Gamm1, Gamm2, b, Sad, Kvap_func, C1, C2)
!***************************************************************************************************************************

function Gamma1 (h)
real(8) :: h, Gamma1
Gamma1 = (1.0d0 + (Alpha1*h)**(Npar1))**(-Mpar1)
end function Gamma1

function Gamma2 (h)
real(8) :: h, Gamma2
Gamma2 = (1.0d0 + (Alpha2*h)**(Npar2))**(-Mpar2)
end function Gamma2

function b ()
real(8) :: b
if (.not. L_BiModal) then
   b  = 0.1d0 + 0.2d0/Npar1**2 * (1.0d0 - dexp(-((WCr/(WCs-WCr))**2)))
else
   nn = Npar1
   if (Alpha2 > Alpha1) nn = Npar2
   b  = 0.1d0 + 0.2d0/nn**2 * (1.0d0 - dexp(-((WCr/(WCs-WCr))**2)))
end if
end function b

function Sad (h)
real(8) :: Sad, h
xa = dlog10(ha)
x0 = dlog10(h0)
x  = dlog10(h)
bb = b()
Sad = 1.0d0 + (x - xa + bb*dlog(1.0d0 + dexp((xa-x)/bb))) / (xa - x0)
end function Sad

function dSad_dh (h)
real(8) :: dSad_dh, h
xa = dlog10(ha)
x0 = dlog10(h0)
x  = dlog10(h)
bb = b()
dSad_dh = -1.0d0/(h*dlog(10.0d0)*(xa-x0) * (1.0d0 + dexp((xa-x)/bb)))
end function dSad_dh

function Kvap_func (WC, h, Temp)
! Temp in degree Celsius
real(8), intent(in)    :: WC, h, Temp
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
ksi       = (WCs-WC)**p/WCs**2
D         = ksi*(WCs-WC)*Da
Hr        = dexp(h/100.0d0*MgRT)     ! h must be in m, thus h (cm) is idvided by 100
Kvap_func = fKvap*D*Hr
end function Kvap_func

function C1 (h)
real(8) :: h, C1
C1 = Alpha1*Npar1*Mpar1*(Alpha1*dabs(h))**(Npar1-1.0d0)*(1.0d0+(Alpha1*dabs(h))**Npar1)**(-Mpar1-1.0d0)
end function C1

function C2 (h)
real(8) :: h, C2
C2 = Alpha2*Npar2*Mpar2*(Alpha2*dabs(h))**(Npar2-1.0d0)*(1.0d0+(Alpha2*dabs(h))**Npar2)**(-Mpar2-1.0d0)
end function C2
   
function WC_MvG (h)
real(8), intent(in)    :: h
real(8)                :: WC_MvG

if (h >= 0.0d0) then
   WC_MvG = WCs
else
   Scap = Gamma1 (dabs(h))
   WC_MvG = WCr + Scap*(WCs-WCr)
end if

end function WC_MvG

function WC_MvG_s (h)
real(8), intent(in)    :: h
real(8)                :: WC_MvG_s

if (h >= 0.0d0) then
   WC_MvG_s = WCs
else
   Gam01 = Gamma1 (dabs(h0))
   Gamh1 = Gamma1 (dabs(h))
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   WC_MvG_s = WCr + Scap*(WCs-WCr)
end if

end function WC_MvG_s

function WC_MvG_2 (h)
real(8), intent(in)    :: h
real(8)                :: WC_MvG_2

if (h >= 0.0d0) then
   WC_MvG_2 = WCs
else
   Gamh1 = Gamma1 (dabs(h))
   Gamh2 = Gamma2 (dabs(h))
   Scap = Omega1 * Gamh1 + Omega2 * Gamh2
   WC_MvG_2 = WCr + Scap*(WCs-WCr)
end if

end function WC_MvG_2

function WC_MvG_2_s (h)
real(8), intent(in)    :: h
real(8)                :: WC_MvG_2_s

if (h >= 0.0d0) then
   WC_MvG_2_s = WCs
else
   Gam01 = Omega1 * Gamma1 (dabs(h0))
   Gamh1 = Omega1 * Gamma1 (dabs(h))
   Gam02 = Omega2 * Gamma2 (dabs(h0))
   Gamh2 = Omega2 * Gamma2 (dabs(h))
   Scap = (Gamh1 + Gamh2 - Gam01 - Gam02) / (1.0d0 - Gam01 - Gam02)
   WC_MvG_2_s = WCr + Scap*(WCs-WCr)
end if

end function WC_MvG_2_s

function K_MvG (h)
implicit none
real(8), intent(in)          :: h
real(8)                      :: K_MvG

if (h >= 0.0d0) then
   K_MvG = Ksat
else
   Scap = Gamma1 (dabs(h))
   K_MvG = Ksat * Scap**Lpar * (1.0d0 - (1.0d0 - Scap**(1.0d0/Mpar1))**Mpar1)**2
end if

end function K_MvG

function K_MvG_s (h)
implicit none
real(8), intent(in)          :: h
real(8)                      :: K_MvG_s

if (h >= 0.0d0) then
   K_MvG_s = Ksat
else
   Gam01 = Gamma1 (dabs(h0))
   Gamh1 = Gamma1 (dabs(h))
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   K_MvG_s = Ksat*Scap**Lpar * (1.0d0 - ((1.0d0-Gamh1**(1.0d0/Mpar1))/(1.0d0-Gam01**(1.0d0/Mpar1)))**Mpar1)**2
end if

end function K_MvG_s

function K_MvG_2 (h)
implicit none
real(8), intent(in)          :: h
real(8)                      :: K_MvG_2

if (h >= 0.0d0) then
   K_MvG_2 = Ksat
else
   Gamh1 = Gamma1 (dabs(h))
   Gamh2 = Gamma2 (dabs(h))
   t1 = (Omega1*Gamh1 + Omega2*Gamh2)**Lpar
   t2 = Omega1*Alpha1*(1.0d0-Gamh1**(1.0d0/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gamh2**(1.0d0/Mpar2))**Mpar2
   t3 = Omega1*Alpha1 + Omega2*Alpha2
   K_MvG_2 = Ksat*t1*(1.0d0-t2/t3)**2
end if

end function K_MvG_2

function K_MvG_2_s (h)
implicit none
real(8), intent(in)          :: h
real(8)                      :: K_MvG_2_s

if (h >= 0.0d0) then
   K_MvG_2_s = Ksat
else
   Gam01 = Gamma1 (dabs(h0))
   Gamh1 = Gamma1 (dabs(h))
   Gam02 = Gamma2 (dabs(h0))
   Gamh2 = Gamma2 (dabs(h))
   Scap1 = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Scap2 = (Gamh2 - Gam02) / (1.0d0 - Gam02)
   t1 = (Omega1*Scap1 + Omega2*Scap2)**Lpar
   t2 = Omega1*Alpha1*(1.0d0-Gamh1**(1.0d0/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gamh2**(1.0d0/Mpar2))**Mpar2
   t3 = Omega1*Alpha1*(1.0d0-Gam01**(1.0d0/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gam02**(1.0d0/Mpar2))**Mpar2
   K_MvG_2_s = Ksat*t1*(1.0d0-t2/t3)**2
end if

end function K_MvG_2_s

function WC_PDI (h)
real(8), intent(in)    :: h
real(8)                :: WC_PDI

if (h >= 0.0d0) then
   WC_PDI = WCs
else
   Scap = Gamma1 (dabs(h))
   WC_PDI = Sad(dabs(h))*WCr + Scap*(WCs-WCr)
end if

end function WC_PDI

function WC_PDI_s (h)
real(8), intent(in)    :: h
real(8)                :: WC_PDI_s

if (h >= 0.0d0) then
   WC_PDI_s = WCs
else
   Gamh1 = Gamma1 (dabs(h))
   Gam01 = Gamma1 (dabs(h0))
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   WC_PDI_s = Sad(dabs(h))*WCr + Scap*(WCs-WCr)
end if

end function WC_PDI_s

function WC_PDI_2 (h)
real(8), intent(in)    :: h
real(8)                :: WC_PDI_2

if (h >= 0.0d0) then
   WC_PDI_2 = WCs
else
   Gamh1 = Gamma1 (dabs(h))
   Gamh2 = Gamma2 (dabs(h))
   Scap = Omega1 * Gamh1 + Omega2 * Gamh2
   WC_PDI_2 = Sad(dabs(h))*WCr + Scap*(WCs-WCr)
end if

end function WC_PDI_2

function WC_PDI_2_s (h)
real(8), intent(in)    :: h
real(8)                :: WC_PDI_2_s

if (h >= 0.0d0) then
   WC_PDI_2_s = WCs
else
   Gam01 = Omega1*Gamma1 (dabs(h0))
   Gamh1 = Omega1*Gamma1 (dabs(h))
   Gam02 = Omega2*Gamma2 (dabs(h0))
   Gamh2 = Omega2*Gamma2 (dabs(h))
   Scap = (Gamh1 + Gamh2 - (Gam01 + Gam02)) / (1.0d0 - (Gam01 + Gam02))
   WC_PDI_2_s = Sad(dabs(h))*WCr + Scap*(WCs-WCr)
end if

end function WC_PDI_2_s

function K_PDI (h, WC, Temp)
implicit none
real(8), intent(in)          :: h, WC, Temp
real(8)                      :: K_PDI
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI = Ksat
else
   Scap = Gamma1 (dabs(h))
   Kcap = Scap**Lpar*(1.0d0 - (1.0d0 - Scap**(1.0d0/Mpar1))**Mpar1)**2
   Kfilm = (h0/ha)**(Apar*(1.0d0-Sad(dabs(h))))
   if (L_NoVap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp) * Conv
   end if
   K_PDI = Ksat*((1.0d0-OmegaK)*Kcap + OmegaK*Kfilm) + Kvap
end if

end function K_PDI

function K_PDI_s (h, WC, Temp)
implicit none
real(8), intent(in)          :: h, WC, Temp
real(8)                      :: K_PDI_s
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_s = Ksat
else
   Gamh1 = Gamma1 (dabs(h))
   Gam01 = Gamma1 (dabs(h0))
   Scap = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Kcap = Scap**Lpar*(1.0d0 - ((1.0d0-Gamh1**(1.0d0/Mpar1))/(1.0d0-Gam01**(1.0d0/Mpar1)))**Mpar1)**2
   Kfilm = (h0/ha)**(Apar*(1.0d0-Sad(dabs(h))))
   if (L_NoVap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp) * Conv
   end if
   K_PDI_s = Ksat*((1.0d0-OmegaK)*Kcap + OmegaK*Kfilm) + Kvap
end if

end function K_PDI_s

function K_PDI_2 (h, WC, Temp)
implicit none
real(8), intent(in)          :: h, WC, Temp
real(8)                      :: K_PDI_2
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_2 = Ksat
else
   Gamh1 = Gamma1 (dabs(h))
   Gamh2 = Gamma2 (dabs(h))
   t1 = (Omega1*Gamh1 + Omega2*Gamh2)**Lpar
   t2 = Omega1*Alpha1*(1.0d0-Gamh1**(1/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gamh2**(1/Mpar2))**Mpar2
   t3 = Omega1*Alpha1 + Omega2*Alpha2
   Kcap = t1*(1.0d0-t2/t3)**2
   Kfilm = (h0/ha)**(Apar*(1.0d0-Sad(dabs(h))))
   if (L_NoVap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp) * Conv
   end if
   K_PDI_2 = Ksat*((1.0d0-OmegaK)*Kcap + OmegaK*Kfilm) + Kvap
end if

end function K_PDI_2

function K_PDI_2_s (h, WC, Temp)
implicit none
real(8), intent(in)          :: h, WC, Temp
real(8)                      :: K_PDI_2_s
real(8), parameter           :: Conv = 100.0d0 * 86400.0d0    ! to convert m/s to cm/d; 100 cm in 1 m, 86400 sec in 1 day

if (h >= 0.0d0) then
   K_PDI_2_s = Ksat
else
   Gam01 = Gamma1 (dabs(h0))
   Gamh1 = Gamma1 (dabs(h))
   Gam02 = Gamma2 (dabs(h0))
   Gamh2 = Gamma2 (dabs(h))
   Scap1 = (Gamh1 - Gam01) / (1.0d0 - Gam01)
   Scap2 = (Gamh2 - Gam02) / (1.0d0 - Gam02)
   t1 = (Omega1*Scap1 + Omega2*Scap2)**Lpar
   t2 = Omega1*Alpha1*(1.0d0-Gamh1**(1.0d0/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gamh2**(1.0d0/Mpar2))**Mpar2
   t3 = Omega1*Alpha1*(1.0d0-Gam01**(1.0d0/Mpar1))**Mpar1 + Omega2*Alpha2*(1.0d0-Gam02**(1.0d0/Mpar2))**Mpar2
   Kcap = t1*(1.0d0-t2/t3)**2
   Kfilm = (h0/ha)**(Apar*(1.0d0-Sad(dabs(h))))
   if (L_NoVap) then
      Kvap = 0.0d0
   else
      Kvap = Kvap_func (WC, dabs(h), Temp) * Conv
   end if
   K_PDI_2_s = Ksat*((1.0d0-OmegaK)*Kcap + OmegaK*Kfilm) + Kvap
end if

end function K_PDI_2_s

function C_MvG (h)
real(8), intent(in)    :: h
real(8)                :: C_MvG
if (h >= 0.0d0) then
   C_MvG = 0.0d0
else
   C_MvG = (WCs-WCr) * C1(h)
end if
end function C_MvG

function C_MvG_s (h)
real(8), intent(in)    :: h
real(8)                :: C_MvG_s
if (h >= 0.0d0) then
   C_MvG_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(h0))
   C_MvG_s = (WCs-WCr)/(1.0d0-Gam01) * C1(h)
end if
end function C_MvG_s

function C_MvG_2 (h)
real(8), intent(in)    :: h
real(8)                :: C_MvG_2
if (h >= 0.0d0) then
   C_MvG_2 = 0.0d0
else
   C_MvG_2 = (WCs-WCr)*(Omega1*C1(h) + Omega2*C2(h))
end if
end function C_MvG_2

function C_MvG_2_s (h)
real(8), intent(in)    :: h
real(8)                :: C_MvG_2_s
if (h >= 0.0d0) then
   C_MvG_2_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(h0))
   Gam02 = Gamma2 (dabs(h0))
   C_MvG_2_s = (WCs-WCr)*(Omega1*C1(h)/(1.0d0-Gam01) + Omega2*C2(h)/(1.0d0-Gam02))
end if
end function C_MvG_2_s

function C_PDI (h)
real(8), intent(in)    :: h
real(8)                :: C_PDI
if (h >= 0.0d0) then
   C_PDI = 0.0d0
else
   SSad = dSad_dh (dabs(h))
   C_PDI = (WCs-WCr)*C1(h) + WCr*SSad
end if
end function C_PDI

function C_PDI_s (h)
real(8), intent(in)    :: h
real(8)                :: C_PDI_s
if (h >= 0.0d0) then
   C_PDI_s = 0.0d0
else
   Gam01 = Gamma1 (dabs(h0))
   SSad  = dSad_dh (dabs(h))
   C_PDI_s = (WCs-WCr)/(1.0d0-Gam01)*C1(h) + WCr*SSad
end if
end function C_PDI_s

function C_PDI_2 (h)
real(8), intent(in)    :: h
real(8)                :: C_PDI_2
if (h >= 0.0d0) then
   C_PDI_2 = 0.0d0
else
   SSad = dSad_dh (dabs(h))
   C_PDI_2 = (WCs-WCr)*(Omega1*C1(h) + Omega2*C2(h)) + WCr*SSad
end if
end function C_PDI_2

function C_PDI_2_s (h)
real(8), intent(in)    :: h
real(8)                :: C_PDI_2_s, Gam0
if (h >= 0.0d0) then
   C_PDI_2_s = 0.0d0
else
   Gam01 = Omega1*Gamma1 (dabs(h0))
   Gam02 = Omega2*Gamma2 (dabs(h0))
   Gam0  = Gam01 + Gam02
   SSad  = dSad_dh (dabs(h))
   C_PDI_2_s = (WCs-WCr)*(Omega1*C1(h)/(1.0d0-Gam0) + Omega2*C2(h)/(1.0d0-Gam0)) + WCr*SSad
end if
end function C_PDI_2_s

end module WC_K_models_04_11
