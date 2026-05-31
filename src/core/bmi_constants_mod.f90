!> BMI string-length constants exported as C globals (xmipy reads them via
!! c_int.in_dll). Must be bind(C) module variables to have external linkage.
module bmi_constants_mod
   use iso_c_binding, only: c_int
   implicit none
   integer(c_int), bind(C, name='BMI_LENCOMPONENTNAME') :: BMI_LENCOMPONENTNAME = 256
   integer(c_int), bind(C, name='BMI_LENVERSION') :: BMI_LENVERSION = 256
   integer(c_int), bind(C, name='BMI_LENVARADDRESS') :: BMI_LENVARADDRESS = 256
   integer(c_int), bind(C, name='BMI_LENVARTYPE') :: BMI_LENVARTYPE = 256
   integer(c_int), bind(C, name='BMI_LENGRIDTYPE') :: BMI_LENGRIDTYPE = 256
   integer(c_int), bind(C, name='BMI_LENERRMESSAGE') :: BMI_LENERRMESSAGE = 1024
end module bmi_constants_mod
