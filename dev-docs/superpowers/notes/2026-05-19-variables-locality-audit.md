# variables.f90 audit

Total active declarations: 606

- **config-side**: 221
- **runtime-no-home**: 174
- **orphan-read-only**: 141
- **state-homed-runtime**: 65
- **state-homed-unknown**: 5

## Locality candidates
Symbols touched by exactly 1 subroutine in 1 file: 136

| Symbol | Type | Category | Sub | File |
|---|---|---|---|---|
| `Agedrain` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `Ageirr` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `Agepond` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `Agepondm1` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `Agepre` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `co2ppm` | real(8) | orphan-read-only | src/crop/cropgrowth_helpers.f90:facco2 | src/crop/cropgrowth_helpers.f90 |
| `co2year` | integer | orphan-read-only | src/crop/cropgrowth_helpers.f90:facco2 | src/crop/cropgrowth_helpers.f90 |
| `cQMpLatSs` | real(8) | runtime-no-home | src/soil/soilhydraulics.f90:soilwater | src/soil/soilhydraulics.f90 |
| `CritDevMasBal` | real(8) | orphan-read-only | src/soil/waterbalance.f90:checkmassbal | src/soil/waterbalance.f90 |
| `dayfix` | integer | runtime-no-home | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `daygrowth` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `daygrowthpot` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `days_counter_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `days_interval_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `DaysGrazingtab` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `dev_cmb` | integer | orphan-read-only | src/soil/waterbalance.f90:checkmassbal | src/soil/waterbalance.f90 |
| `DiPoCp` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `ditab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `drbl` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `drblpot` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `drfil` | character(len=16) | config-side | src/io/toml/read_drainage_toml.f90:read_drainage_inner | src/io/toml/read_drainage_toml.f90 |
| `dropr` | real(8) | orphan-read-only | src/drainage/surfacewater.f90:wlevbal | src/drainage/surfacewater.f90 |
| `dwatab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `fidtab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `flanthesis` | logical | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `flCropOpenFile` | logical | runtime-no-home | src/crop/cropgrowth_helpers.f90:cropoutput | src/crop/cropgrowth_helpers.f90 |
| `flCropReadFile` | logical | state-homed-runtime | src/crop/cropgrowth.f90:cropgrowth | src/crop/cropgrowth.f90 |
| `FlDecMpRat` | logical | orphan-read-only | src/core/timecontrol_mod.f90:timecontrol_reduce_dt | src/core/timecontrol_mod.f90 |
| `fldumpconvcrit` | logical | orphan-read-only | src/soil/soilhydraulics.f90:headcalc | src/soil/soilhydraulics.f90 |
| `flGrazing` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flGrazingpot` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flHarvest` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flHarvestpot` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flhrvendact` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flhrvendpot` | logical | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `flrunon` | logical | config-side | src/boundary/boundtop.f90:boundtop | src/boundary/boundtop.f90 |
| `flwarn_hc` | logical | runtime-no-home | src/soil/soilhydraulics.f90:headcalc | src/soil/soilhydraulics.f90 |
| `gasst` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `gasstpot` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `grzdm` | real(8) | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `gwlcrit` | real(8) | orphan-read-only | src/drainage/surfacewater.f90:wlevbal | src/drainage/surfacewater.f90 |
| `haqtab` | real(8) | config-side | src/boundary/boundbottom.f90:boundbottom | src/boundary/boundbottom.f90 |
| `hbotab` | real(8) | config-side | src/boundary/boundbottom.f90:boundbottom | src/boundary/boundbottom.f90 |
| `hcritab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `hqhtab` | real(8) | orphan-read-only | src/utils/surfacewaterutils.f90:qhtab | src/utils/surfacewaterutils.f90 |
| `IAvFrMpWlWtDm1` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `IAvFrMpWlWtDm2` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `icAgetopdwn` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `icAgetopupw` | real(8) | runtime-no-home | src/solute/agetracer.f90:agetracer | src/solute/agetracer.f90 |
| `icn_atm` | integer | state-homed-runtime | src/atmosphere/meteoday.f90:cnmethod | src/atmosphere/meteoday.f90 |
| `idaysgraz` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `idaysgrazpot` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `ientrytablay` | integer | orphan-read-only | src/soil/soilhydraulics.f90:soilwater | src/soil/soilhydraulics.f90 |
| `iharvest` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `iQExcMtxDm1Cp` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `iQExcMtxDm2Cp` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `iqinfmax` | real(8) | orphan-read-only | src/io/swap_csv_output.f90:set_values | src/io/swap_csv_output.f90 |
| `iQOutDrRapCp` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `irconc` | real(8) | config-side | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `irdate` | real(8) | config-side | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `irdepth` | real(8) | config-side | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `irtype` | integer | config-side | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `iseqgm` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `iseqgmpot` | integer | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `iwarn_hc` | integer | runtime-no-home | src/soil/soilhydraulics.f90:headcalc | src/soil/soilhydraulics.f90 |
| `LossGrazingtab` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `lossgrztab` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `lossmowtab` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `mrest` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `mrestpot` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `mrftb` | real(8) | orphan-read-only | src/crop/cropfixed_runtime.f90:cropfixed | src/crop/cropfixed_runtime.f90 |
| `nhead` | integer | config-side | src/soil/soilhydraulics.f90:soilwater | src/soil/soilhydraulics.f90 |
| `nirri_ssdi_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `nod10_cn` | integer | state-homed-runtime | src/atmosphere/meteoday.f90:cnmethod | src/atmosphere/meteoday.f90 |
| `nod_ssdi_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `nod_ssdi_sensor_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `nphase` | integer | orphan-read-only | src/drainage/surfacewater.f90:wlevbal | src/drainage/surfacewater.f90 |
| `nsublay` | integer | config-side | src/soil/soilgrid.f90:calcgrid | src/soil/soilgrid.f90 |
| `o2_capac_term` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_d_soil_term1` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_d_soil_term2` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_gfp100` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_ini_stress` | logical | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_mplus1` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `o2_nmin1` | real(8) | runtime-no-home | src/crop/oxygenstress.f90:oxygenstress | src/crop/oxygenstress.f90 |
| `outdat` | real(8) | orphan-read-only | src/core/timecontrol_mod.f90:timecontrol_advance | src/core/timecontrol_mod.f90 |
| `outdatint` | real(8) | config-side | src/core/timecontrol_mod.f90:timecontrol_advance | src/core/timecontrol_mod.f90 |
| `pgrzdm` | real(8) | runtime-no-home | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `pondmxtab` | real(8) | orphan-read-only | src/boundary/boundtop.f90:pondrunoff | src/boundary/boundtop.f90 |
| `qdraincomp` | real(8) | runtime-no-home | src/soil/waterbalance.f90:integral | src/soil/waterbalance.f90 |
| `qdrtab` | real(8) | orphan-read-only | src/drainage/drainage.f90:bocodrb | src/drainage/drainage.f90 |
| `QExcMpMtx` | real(8) | orphan-read-only | src/soil/waterbalance.f90:fluxes | src/soil/waterbalance.f90 |
| `qimmob` | real(8) | orphan-read-only | src/soil/waterbalance.f90:fluxes | src/soil/waterbalance.f90 |
| `QMaPo` | real(8) | orphan-read-only | src/soil/waterbalance.f90:fluxes | src/soil/waterbalance.f90 |
| `raintab` | real(8) | config-side | src/atmosphere/meteodt.f90:processrainevents | src/atmosphere/meteodt.f90 |
| `rawtab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `rot` | integer | orphan-read-only | src/io/swapoutput.f90:outrot | src/io/swapoutput.f90 |
| `runonarr` | real(8) | orphan-read-only | src/boundary/boundtop.f90:boundtop | src/boundary/boundtop.f90 |
| `siccaptb` | real(8) | orphan-read-only | src/atmosphere/meteoday.f90:processmeteoday | src/atmosphere/meteoday.f90 |
| `snw` | integer | orphan-read-only | src/io/swapoutput.f90:snowoutput | src/io/swapoutput.f90 |
| `sptablay` | real(8) | orphan-read-only | src/soil/soilhydraulics.f90:soilwater | src/soil/soilhydraulics.f90 |
| `ssdi_amount_f_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_amount_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_appl_rate_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_date_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_rate_f_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_sched_type_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_schedule_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_threshold_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `ssdi_threshold_z_irr` | real(8) | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `sw_interval_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `swcaprise` | logical | config-side | src/soil/soilhydraulics.f90:headcalc | src/soil/soilhydraulics.f90 |
| `swpondmx` | integer | orphan-read-only | src/boundary/boundtop.f90:pondrunoff | src/boundary/boundtop.f90 |
| `swssdi_irr` | integer | runtime-no-home | src/crop/irrigation.f90:ssdi_irrigation | src/crop/irrigation.f90 |
| `swtopsub` | integer | orphan-read-only | src/crop/oxygenstress.f90:oxygen_dat | src/crop/oxygenstress.f90 |
| `tadw` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `tadwpot` | real(8) | runtime-no-home | src/crop/cropwofost_runtime.f90:wofost | src/crop/cropwofost_runtime.f90 |
| `tau` | real(8) | orphan-read-only | src/soil/soilhydraulics.f90:hysteresis | src/soil/soilhydraulics.f90 |
| `tawtab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `tcritab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `tendirrig` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `tmnr` | real(8) | runtime-no-home | src/crop/cropgrowth.f90:cropgrowth | src/crop/cropgrowth.f90 |
| `treltab` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `tstairrig` | real(8) | orphan-read-only | src/crop/irrigation.f90:irrigation | src/crop/irrigation.f90 |
| `tsunrise_atm` | real(8) | runtime-no-home | src/atmosphere/meteodt.f90:etsine | src/atmosphere/meteodt.f90 |
| `tsunset_atm` | real(8) | runtime-no-home | src/atmosphere/meteodt.f90:etsine | src/atmosphere/meteodt.f90 |
| `tz_z1_z2` | real(8) | orphan-read-only | src/io/swap_csv_output.f90:csv_out_tz | src/io/swap_csv_output.f90 |
| `UptGrazingtab` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `VlMpStDm1` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `VlMpStDm2` | real(8) | orphan-read-only | src/soil/soilgrid.f90:convertdiscrvert | src/soil/soilgrid.f90 |
| `wlsman` | real(8) | orphan-read-only | src/drainage/surfacewater.f90:wlevbal | src/drainage/surfacewater.f90 |
| `wlstab` | real(8) | orphan-read-only | src/drainage/surfacewater.f90:wballev | src/drainage/surfacewater.f90 |
| `wrtb` | real(8) | orphan-read-only | src/crop/cropfixed_runtime.f90:cropfixed | src/crop/cropfixed_runtime.f90 |
| `zgrz` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |
| `zi` | real(8) | config-side | src/soil/soilhydraulics.f90:soilwater | src/soil/soilhydraulics.f90 |
| `zmow` | real(8) | orphan-read-only | src/crop/cropgrass_runtime.f90:grass | src/crop/cropgrass_runtime.f90 |

## Per-symbol table

| Symbol | Type | Category | State home | NReadFiles | NWriteFiles | NSubs |
|---|---|---|---|---|---|---|
| `adcrh` | real(8) | config-side |  | 7 | 3 | 9 |
| `adcrl` | real(8) | config-side |  | 7 | 3 | 9 |
| `aeratecrit` | real(8) | config-side |  | 7 | 3 | 9 |
| `air_filled_root_por` | real(8) | config-side |  | 4 | 1 | 4 |
| `alt` | real(8) | config-side |  | 4 | 0 | 4 |
| `altw` | real(8) | config-side |  | 4 | 0 | 4 |
| `amaxtb` | real(8) | config-side |  | 6 | 2 | 6 |
| `angstroma` | real(8) | config-side |  | 3 | 0 | 3 |
| `angstromb` | real(8) | config-side |  | 3 | 0 | 3 |
| `aqamp` | real(8) | config-side |  | 2 | 0 | 2 |
| `aqave` | real(8) | config-side |  | 2 | 0 | 2 |
| `aqper` | real(8) | config-side |  | 2 | 0 | 2 |
| `aqtmax` | real(8) | config-side |  | 2 | 0 | 2 |
| `basegw` | real(8) | config-side |  | 3 | 0 | 3 |
| `bexp` | real(8) | config-side |  | 2 | 0 | 2 |
| `cfevappond` | real(8) | config-side |  | 3 | 0 | 2 |
| `cofintfl` | real(8) | config-side |  | 2 | 0 | 3 |
| `cofqha` | real(8) | config-side |  | 2 | 0 | 2 |
| `cofqhb` | real(8) | config-side |  | 2 | 0 | 2 |
| `cofqhc` | real(8) | config-side |  | 2 | 0 | 2 |
| `cofred` | real(8) | config-side |  | 1 | 0 | 3 |
| `cpre` | real(8) | config-side |  | 3 | 0 | 3 |
| `cref` | real(8) | config-side |  | 2 | 0 | 2 |
| `CritDevPondDt` | real(8) | config-side |  | 2 | 0 | 2 |
| `cropfil` | character(len=40) | config-side |  | 3 | 0 | 2 |
| `croptype` | integer | config-side |  | 6 | 0 | 7 |
| `cseeptab` | real(8) | config-side |  | 2 | 0 | 2 |
| `cvl` | real(8) | config-side |  | 7 | 2 | 7 |
| `cvo` | real(8) | config-side |  | 4 | 1 | 4 |
| `cvr` | real(8) | config-side |  | 7 | 2 | 7 |
| `cvs` | real(8) | config-side |  | 7 | 2 | 7 |
| `daquif` | real(8) | config-side |  | 2 | 0 | 2 |
| `dcritrtz` | real(8) | config-side |  | 6 | 1 | 8 |
| `ddamp` | real(8) | config-side |  | 2 | 0 | 2 |
| `ddif` | real(8) | config-side |  | 3 | 0 | 3 |
| `decpot` | real(8) | config-side |  | 2 | 0 | 2 |
| `decsat` | real(8) | config-side |  | 2 | 0 | 2 |
| `DelayRegrowthTab` | real(8) | config-side |  | 2 | 1 | 2 |
| `dlc` | real(8) | config-side |  | 3 | 1 | 3 |
| `dlo` | real(8) | config-side |  | 3 | 1 | 3 |
| `dmmowtb` | real(8) | config-side |  | 3 | 1 | 3 |
| `dramet` | integer | config-side |  | 4 | 0 | 6 |
| `drfil` | character(len=16) | config-side |  | 1 | 0 | 1 |
| `dry_mat_cont_roots` | real(8) | config-side |  | 4 | 1 | 4 |
| `dtsmtb` | real(8) | config-side |  | 3 | 1 | 3 |
| `dvsend` | real(8) | config-side |  | 6 | 2 | 6 |
| `dvsnlt` | real(8) | config-side |  | 3 | 1 | 3 |
| `eff` | real(8) | config-side |  | 6 | 2 | 8 |
| `entres` | real(8) | config-side |  | 3 | 0 | 3 |
| `fdepth` | real(8) | config-side |  | 2 | 0 | 2 |
| `flCO2` | logical | config-side |  | 4 | 1 | 3 |
| `flrunon` | logical | config-side |  | 1 | 0 | 1 |
| `fltb` | real(8) | config-side |  | 7 | 2 | 7 |
| `fotb` | real(8) | config-side |  | 4 | 1 | 4 |
| `fraharlosorm_lv` | real(8) | config-side |  | 3 | 1 | 3 |
| `frexp` | real(8) | config-side |  | 2 | 0 | 2 |
| `frtb` | real(8) | config-side |  | 7 | 2 | 7 |
| `fstb` | real(8) | config-side |  | 7 | 2 | 7 |
| `ftopdislay` | real(8) | config-side |  | 3 | 0 | 3 |
| `gampar` | real(8) | config-side |  | 2 | 0 | 2 |
| `geofac` | real(8) | config-side |  | 2 | 0 | 3 |
| `gwlconv` | real(8) | config-side |  | 3 | 0 | 3 |
| `gwli` | real(8) | config-side |  | 3 | 0 | 3 |
| `gwltab` | real(8) | config-side |  | 2 | 0 | 2 |
| `haqtab` | real(8) | config-side |  | 1 | 0 | 1 |
| `hbotab` | real(8) | config-side |  | 1 | 0 | 1 |
| `hcomp` | real(8) | config-side |  | 2 | 0 | 1 |
| `hdrain` | real(8) | config-side |  | 2 | 0 | 2 |
| `hlim1` | real(8) | config-side |  | 7 | 3 | 9 |
| `hlim2l` | real(8) | config-side |  | 7 | 3 | 9 |
| `hlim2u` | real(8) | config-side |  | 7 | 3 | 9 |
| `hlim3h` | real(8) | config-side |  | 7 | 3 | 9 |
| `hlim3l` | real(8) | config-side |  | 7 | 3 | 9 |
| `hlim4` | real(8) | config-side |  | 7 | 3 | 9 |
| `hplate` | real(8) | config-side |  | 2 | 0 | 2 |
| `hsublay` | real(8) | config-side |  | 2 | 0 | 2 |
| `idev` | integer | config-side |  | 4 | 1 | 4 |
| `idsl` | integer | config-side |  | 3 | 1 | 3 |
| `ilnmxl` | integer | config-side |  | 2 | 1 | 1 |
| `impend` | real(8) | config-side |  | 3 | 0 | 3 |
| `InList_csv_tz` | character(len=1024) | config-side |  | 1 | 0 | 2 |
| `intwl` | integer | config-side |  | 2 | 0 | 2 |
| `irconc` | real(8) | config-side |  | 1 | 0 | 1 |
| `irdate` | real(8) | config-side |  | 1 | 0 | 1 |
| `irdepth` | real(8) | config-side |  | 1 | 0 | 1 |
| `irtype` | integer | config-side |  | 1 | 0 | 1 |
| `isoillay` | integer | config-side |  | 2 | 0 | 2 |
| `kf` | real(8) | config-side |  | 2 | 0 | 2 |
| `kfsat` | real(8) | config-side |  | 2 | 0 | 2 |
| `khbot` | real(8) | config-side |  | 2 | 0 | 2 |
| `khtop` | real(8) | config-side |  | 2 | 0 | 2 |
| `kvbot` | real(8) | config-side |  | 2 | 0 | 2 |
| `kvtop` | real(8) | config-side |  | 2 | 0 | 2 |
| `lat` | real(8) | config-side |  | 7 | 0 | 9 |
| `ldis` | real(8) | config-side |  | 3 | 0 | 3 |
| `lrnr` | real(8) | config-side |  | 5 | 1 | 6 |
| `MaxBackTr` | integer | config-side |  | 2 | 0 | 2 |
| `metfil` | character(len=200) | config-side |  | 2 | 0 | 2 |
| `ncomp` | integer | config-side |  | 2 | 0 | 2 |
| `nconc` | integer | config-side |  | 2 | 0 | 2 |
| `nhead` | integer | config-side |  | 1 | 0 | 1 |
| `nlai` | real(8) | config-side |  | 3 | 1 | 3 |
| `nlue` | real(8) | config-side |  | 5 | 1 | 4 |
| `nmetdetail` | integer | config-side |  | 7 | 0 | 8 |
| `nmxlv` | real(8) | config-side |  | 5 | 1 | 6 |
| `nowltab` | integer | config-side |  | 3 | 0 | 4 |
| `nrstaring` | integer | config-side |  | 2 | 0 | 2 |
| `nsla` | real(8) | config-side |  | 3 | 1 | 3 |
| `nsublay` | integer | config-side |  | 1 | 0 | 1 |
| `NumLevRapDra` | integer | config-side |  | 2 | 0 | 4 |
| `osswlm` | real(8) | config-side |  | 2 | 0 | 2 |
| `outdatint` | real(8) | config-side |  | 1 | 0 | 1 |
| `outfil` | character(len=16) | config-side |  | 7 | 0 | 15 |
| `pathatm` | character(len=80) | config-side |  | 5 | 0 | 4 |
| `pathcrop` | character(len=80) | config-side |  | 3 | 0 | 3 |
| `pathdrain` | character(len=80) | config-side |  | 2 | 0 | 2 |
| `pathwork` | character(len=80) | config-side |  | 7 | 0 | 16 |
| `perdl` | real(8) | config-side |  | 6 | 2 | 7 |
| `poros` | real(8) | config-side |  | 2 | 0 | 2 |
| `project` | character(len=80) | config-side |  | 10 | 0 | 21 |
| `q10` | real(8) | config-side |  | 8 | 2 | 9 |
| `q10_microbial` | real(8) | config-side |  | 4 | 1 | 8 |
| `qbotab` | real(8) | config-side |  | 2 | 0 | 2 |
| `raintab` | real(8) | config-side |  | 1 | 0 | 1 |
| `rdctb` | real(8) | config-side |  | 8 | 3 | 10 |
| `rdmax` | real(8) | config-side |  | 4 | 0 | 4 |
| `rdrrtb` | real(8) | config-side |  | 6 | 2 | 6 |
| `rdrstb` | real(8) | config-side |  | 6 | 2 | 6 |
| `rdtb` | real(8) | config-side |  | 9 | 1 | 9 |
| `rfsetb` | real(8) | config-side |  | 7 | 2 | 8 |
| `rgrlai` | real(8) | config-side |  | 6 | 2 | 7 |
| `rimlay` | real(8) | config-side |  | 3 | 0 | 3 |
| `rlwtb` | real(8) | config-side |  | 6 | 2 | 6 |
| `rml` | real(8) | config-side |  | 7 | 2 | 7 |
| `rmo` | real(8) | config-side |  | 4 | 1 | 4 |
| `rmr` | real(8) | config-side |  | 7 | 2 | 8 |
| `rms` | real(8) | config-side |  | 7 | 2 | 7 |
| `root_radiusO2` | real(8) | config-side |  | 4 | 1 | 4 |
| `rsigni` | real(8) | config-side |  | 3 | 0 | 3 |
| `rsoil` | real(8) | config-side |  | 4 | 0 | 4 |
| `rsurfshallow` | real(8) | config-side |  | 2 | 0 | 2 |
| `rsw` | real(8) | config-side |  | 8 | 3 | 8 |
| `rtheta` | real(8) | config-side |  | 2 | 0 | 2 |
| `salthead` | real(8) | config-side |  | 5 | 1 | 5 |
| `saltmax` | real(8) | config-side |  | 5 | 2 | 7 |
| `saltslope` | real(8) | config-side |  | 5 | 2 | 7 |
| `schedule` | integer | config-side |  | 12 | 3 | 9 |
| `shape` | real(8) | config-side |  | 12 | 0 | 6 |
| `sinamp` | real(8) | config-side |  | 2 | 0 | 2 |
| `sinave` | real(8) | config-side |  | 2 | 0 | 2 |
| `sinmax` | real(8) | config-side |  | 2 | 0 | 2 |
| `slatb` | real(8) | config-side |  | 6 | 2 | 6 |
| `snowcoef` | real(8) | config-side |  | 2 | 0 | 2 |
| `spa` | real(8) | config-side |  | 3 | 1 | 3 |
| `span` | real(8) | config-side |  | 6 | 2 | 8 |
| `spec_weight_root_tissue` | real(8) | config-side |  | 4 | 1 | 4 |
| `specific_resp_humus` | real(8) | config-side |  | 4 | 1 | 8 |
| `srl` | real(8) | config-side |  | 4 | 1 | 4 |
| `ssa` | real(8) | config-side |  | 6 | 2 | 6 |
| `sw2` | integer | config-side |  | 2 | 0 | 2 |
| `sw3` | integer | config-side |  | 2 | 0 | 2 |
| `sw4` | integer | config-side |  | 3 | 0 | 3 |
| `swbotb3Impl` | integer | config-side |  | 2 | 0 | 2 |
| `swbotbc` | integer | config-side |  | 2 | 0 | 2 |
| `swbotbhea` | integer | config-side |  | 2 | 0 | 2 |
| `swbr` | integer | config-side |  | 2 | 0 | 2 |
| `swcalt` | integer | config-side |  | 3 | 0 | 3 |
| `swcaprise` | logical | config-side |  | 1 | 0 | 1 |
| `swcofqhc` | integer | config-side |  | 2 | 0 | 2 |
| `swcompensate` | integer | config-side |  | 7 | 3 | 9 |
| `swdiscrvert` | integer | config-side |  | 4 | 0 | 3 |
| `swdislay` | integer | config-side |  | 4 | 0 | 4 |
| `swdivide` | integer | config-side |  | 4 | 0 | 4 |
| `swdmi2rd` | integer | config-side |  | 9 | 3 | 9 |
| `swetr` | integer | config-side |  | 4 | 0 | 4 |
| `swetsine` | integer | config-side |  | 4 | 0 | 4 |
| `swgc` | integer | config-side |  | 6 | 1 | 6 |
| `swharv` | integer | config-side |  | 8 | 2 | 8 |
| `swhea` | integer | config-side |  | 3 | 0 | 2 |
| `swhyst` | integer | config-side |  | 4 | 0 | 4 |
| `swinter` | integer | config-side |  | 12 | 3 | 12 |
| `swirfix` | integer | config-side |  | 3 | 0 | 3 |
| `swkmean` | integer | config-side |  | 4 | 0 | 6 |
| `swliminf` | integer | config-side |  | 2 | 0 | 2 |
| `swoxygen` | integer | config-side |  | 8 | 3 | 10 |
| `swoxygentype` | integer | config-side |  | 3 | 1 | 5 |
| `swqhbot` | integer | config-side |  | 2 | 0 | 2 |
| `swqhr` | integer | config-side |  | 2 | 0 | 2 |
| `swrd` | integer | config-side |  | 9 | 3 | 9 |
| `swrdc` | integer | config-side |  | 9 | 3 | 8 |
| `swredu` | integer | config-side |  | 4 | 0 | 3 |
| `swrootradius` | integer | config-side |  | 4 | 1 | 4 |
| `swsalinity` | integer | config-side |  | 7 | 3 | 10 |
| `swsolu` | integer | config-side |  | 4 | 0 | 4 |
| `swstressor` | integer | config-side |  | 7 | 3 | 9 |
| `swsublim` | integer | config-side |  | 2 | 0 | 2 |
| `swtopbhea` | integer | config-side |  | 2 | 0 | 2 |
| `swtopdislay` | integer | config-side |  | 3 | 0 | 3 |
| `swWrtNonox` | integer | config-side |  | 7 | 3 | 9 |
| `taccur` | real(8) | config-side |  | 2 | 0 | 4 |
| `tampli` | real(8) | config-side |  | 2 | 0 | 2 |
| `tbase` | real(8) | config-side |  | 9 | 3 | 10 |
| `tdwi` | real(8) | config-side |  | 6 | 2 | 6 |
| `tembtab` | real(8) | config-side |  | 2 | 0 | 2 |
| `temtoptab` | real(8) | config-side |  | 2 | 0 | 2 |
| `tfrostend` | real(8) | config-side |  | 3 | 0 | 3 |
| `tfroststa` | real(8) | config-side |  | 3 | 0 | 3 |
| `till_swtill` | integer | config-side |  | 1 | 0 | 0 |
| `timref` | real(8) | config-side |  | 2 | 0 | 2 |
| `tmean` | real(8) | config-side |  | 2 | 0 | 2 |
| `tmnftb` | real(8) | config-side |  | 5 | 2 | 5 |
| `tmpftb` | real(8) | config-side |  | 5 | 2 | 5 |
| `tscf` | real(8) | config-side |  | 2 | 0 | 2 |
| `tsumam` | real(8) | config-side |  | 6 | 2 | 6 |
| `tsumea` | real(8) | config-side |  | 6 | 2 | 6 |
| `var_a` | real(8) | config-side |  | 4 | 1 | 4 |
| `wldip` | real(8) | config-side |  | 2 | 0 | 2 |
| `wrtmax` | real(8) | config-side |  | 6 | 2 | 6 |
| `zc` | real(8) | config-side |  | 2 | 0 | 2 |
| `zi` | real(8) | config-side |  | 1 | 0 | 1 |
| `zintf` | real(8) | config-side |  | 2 | 0 | 2 |
| `ad` | integer | orphan-read-only |  | 1 | 0 | 3 |
| `am` | integer | orphan-read-only |  | 1 | 0 | 3 |
| `ArMpTp` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `BiModal` | logical | orphan-read-only |  | 3 | 0 | 4 |
| `cirrs` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `cirrthres` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `co2amaxtb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `co2efftb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `co2ppm` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `co2tratb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `co2year` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `CritDevMasBal` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `CriterHr` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `CritUndSatVol` | real(8) | orphan-read-only |  | 1 | 0 | 2 |
| `crp` | integer | orphan-read-only |  | 18 | 0 | 9 |
| `DaysGrazingtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `dcrit` | real(8) | orphan-read-only |  | 2 | 0 | 1 |
| `dev_cmb` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `dewrest` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `DiPoCp` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `ditab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `dmgrztb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `dra` | integer | orphan-read-only |  | 1 | 0 | 0 |
| `dropr` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `dwatab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `fbltb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `fidtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `FlDecMpRat` | logical | orphan-read-only |  | 1 | 0 | 1 |
| `fldumpconvcrit` | logical | orphan-read-only |  | 1 | 0 | 1 |
| `flSwapShared` | logical | orphan-read-only |  | 1 | 0 | 3 |
| `gctb` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `gwlcrit` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `hcritab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `hPrep` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `hqhtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `hSow` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `IAvFrMpWlWtDm1` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `IAvFrMpWlWtDm2` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `ientrytablay` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `inc` | integer | orphan-read-only |  | 1 | 0 | 2 |
| `iQExcMtxDm1Cp` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `iQExcMtxDm2Cp` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `iqinfmax` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `iQMpOutDrRap` | real(8) | orphan-read-only |  | 2 | 0 | 4 |
| `iQOutDrRapCp` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `isuas` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `Itnumb` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `Kroot` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `kstem` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `logf` | integer | orphan-read-only |  | 14 | 0 | 15 |
| `LossGrazingtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `lossgrztab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `lossmowtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `lsda` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `MaxPrepDelay` | integer | orphan-read-only |  | 3 | 0 | 3 |
| `MaxSowDelay` | integer | orphan-read-only |  | 3 | 0 | 3 |
| `mrftb` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `NoVap` | logical | orphan-read-only |  | 3 | 0 | 4 |
| `nphase` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `numtablay` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `o2_bunsencoeff` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_c_macro` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_c_min_micro` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_ctopnode` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_d_o2inwater` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_d_root` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_d_soil` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_depth` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_gas_filled_porosity` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_perc_org_mat` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_r_microbial_z0` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_root_radius` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_sat_water_cont` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_shape_factor_microbialr` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_soil_density` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_soil_temp` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_w_root` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_w_root_z0` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `o2_waterfilm_thickness` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `outdat` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `perirrsurp` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `pld` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `pondmxtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `qdrtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `QExcMpMtx` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `qimmob` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `QMaPo` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `QRapDra` | real(8) | orphan-read-only |  | 2 | 0 | 3 |
| `raithreshold` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `rawtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `remoc` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `rootcoefa` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `rooteff` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `rootradius` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `rot` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `runonarr` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `Rxylem` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `siccaplai` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `siccaptb` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `snw` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `sptablay` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `StepHr` | real(8) | orphan-read-only |  | 1 | 0 | 3 |
| `SwBotb3ResVert` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `swcirrthres` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `swpondmx` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `swtopsub` | integer | orphan-read-only |  | 1 | 0 | 1 |
| `swuseCN` | integer | orphan-read-only |  | 2 | 0 | 2 |
| `tau` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `tawtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `tcritab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `tem` | integer | orphan-read-only |  | 1 | 0 | 2 |
| `TempSow` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `tendirrig` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `till_i_n_model` | integer | orphan-read-only |  | 1 | 0 | 0 |
| `till_iRedist` | integer | orphan-read-only |  | 1 | 0 | 0 |
| `till_Max_Z_tillage` | real(8) | orphan-read-only |  | 1 | 0 | 0 |
| `till_Ntill` | integer | orphan-read-only |  | 1 | 0 | 0 |
| `till_Ntypes` | integer | orphan-read-only |  | 1 | 0 | 0 |
| `treltab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `tstairrig` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `tsumdepth` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `tsumtemp` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `tsumtime` | integer | orphan-read-only |  | 3 | 0 | 3 |
| `tz_z1_z2` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `UptGrazingtab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `vernbase` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `verndvs` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `vernrtb` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `vernsat` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `VlMpStDm1` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `VlMpStDm2` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `wiltpoint` | real(8) | orphan-read-only |  | 4 | 0 | 7 |
| `wlsman` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `wlstab` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `wrtb` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `Z_Tp` | real(8) | orphan-read-only |  | 2 | 0 | 2 |
| `zgrz` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `zmow` | real(8) | orphan-read-only |  | 1 | 0 | 1 |
| `zPrep` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `zSow` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `zTempSow` | real(8) | orphan-read-only |  | 3 | 0 | 3 |
| `Agedrain` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `Ageirr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `Agepond` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `Agepondm1` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `Agepre` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `agerm` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `alphacrit` | real(8) | runtime-no-home |  | 6 | 3 | 8 |
| `atmtr` | real(8) | runtime-no-home |  | 4 | 1 | 5 |
| `bdens` | real(8) | runtime-no-home |  | 5 | 1 | 9 |
| `bgerm` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `botcom` | integer | runtime-no-home |  | 3 | 1 | 6 |
| `c_mroot` | real(8) | runtime-no-home |  | 1 | 1 | 4 |
| `c_top` | real(8) | runtime-no-home |  | 2 | 1 | 3 |
| `cgerm` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `cirr` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `cml` | real(8) | runtime-no-home |  | 6 | 2 | 9 |
| `cmsy` | real(8) | runtime-no-home |  | 2 | 1 | 4 |
| `cQMpLatSs` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `CritDevh1Cp` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `CritDevh2Cp` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dayfix` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `daygrowth` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `daygrowthpot` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `daylp` | real(8) | runtime-no-home |  | 5 | 2 | 7 |
| `days_counter_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `days_interval_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `dethum` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `detrad` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `detrain` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `detrecord` | integer | runtime-no-home |  | 2 | 1 | 2 |
| `dettav` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dettime` | real(8) | runtime-no-home |  | 2 | 1 | 3 |
| `detwind` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dhPrep` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dhSow` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `difpp` | real(8) | runtime-no-home |  | 4 | 1 | 6 |
| `drares` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `drbl` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `drblpot` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `dsinbe` | real(8) | runtime-no-home |  | 4 | 1 | 6 |
| `dt_SSDI_event` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `dtempSow` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dtEventRain` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `dznew` | real(8) | runtime-no-home |  | 4 | 1 | 4 |
| `f_senes` | real(8) | runtime-no-home |  | 1 | 1 | 4 |
| `fbl` | real(8) | runtime-no-home |  | 1 | 1 | 3 |
| `finterception` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `flanthesis` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flCropOpenFile` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flGrazing` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flGrazingpot` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flHarvest` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flHarvestpot` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flhrvendact` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `flhrvendpot` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `FlHydrLift` | logical | runtime-no-home |  | 5 | 3 | 6 |
| `flwarn_hc` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `fstr` | real(8) | runtime-no-home |  | 3 | 2 | 3 |
| `gasst` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `gasstpot` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `gc` | real(8) | runtime-no-home |  | 3 | 2 | 2 |
| `glaiex` | real(8) | runtime-no-home |  | 2 | 2 | 3 |
| `glaiexpot` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `grzdm` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `gwlinf` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `gwrt` | real(8) | runtime-no-home |  | 3 | 2 | 3 |
| `hdrygerm` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `hwetgerm` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `icAgetopdwn` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `icAgetopupw` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `idaysgraz` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `idaysgrazpot` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `idregr` | integer | runtime-no-home |  | 2 | 1 | 2 |
| `idregrpot` | integer | runtime-no-home |  | 2 | 1 | 1 |
| `ientrytab` | integer | runtime-no-home |  | 3 | 1 | 7 |
| `iharvest` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `iHWCKmodel` | integer | runtime-no-home |  | 15 | 2 | 18 |
| `ilvold` | integer | runtime-no-home |  | 2 | 2 | 6 |
| `ilvoldpot` | integer | runtime-no-home |  | 2 | 2 | 2 |
| `infres` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `InList_csv` | character(len=1024) | runtime-no-home |  | 1 | 1 | 2 |
| `inpola` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `inpolb` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `ipos` | integer | runtime-no-home |  | 3 | 1 | 5 |
| `irectotal` | integer | runtime-no-home |  | 2 | 2 | 2 |
| `irrigevent` | integer | runtime-no-home |  | 1 | 1 | 2 |
| `iseqgm` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `iseqgmpot` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `issnowbeg` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `iwarn_hc` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `laiem` | real(8) | runtime-no-home |  | 6 | 4 | 6 |
| `laiexp` | real(8) | runtime-no-home |  | 2 | 2 | 3 |
| `laiexppot` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `laimax` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `lv` | real(8) | runtime-no-home |  | 2 | 2 | 6 |
| `lvage` | real(8) | runtime-no-home |  | 2 | 2 | 5 |
| `lvagepot` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `lvpot` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `max_resp_factor` | real(8) | runtime-no-home |  | 4 | 2 | 6 |
| `mowdm` | real(8) | runtime-no-home |  | 3 | 1 | 5 |
| `mrest` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `mrestpot` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `nird` | real(8) | runtime-no-home |  | 8 | 2 | 10 |
| `nirri` | integer | runtime-no-home |  | 2 | 2 | 3 |
| `nirri_ssdi_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `nmper` | integer | runtime-no-home |  | 3 | 1 | 5 |
| `nod1lay` | integer | runtime-no-home |  | 3 | 1 | 3 |
| `nod_ssdi_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `nod_ssdi_sensor_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `numbit` | integer | runtime-no-home |  | 2 | 1 | 2 |
| `numlay` | integer | runtime-no-home |  | 5 | 2 | 6 |
| `numnodnew` | integer | runtime-no-home |  | 4 | 1 | 4 |
| `numtab` | integer | runtime-no-home |  | 2 | 1 | 6 |
| `o2_capac_term` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `o2_d_soil_term1` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `o2_d_soil_term2` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `o2_gfp100` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `o2_ini_stress` | logical | runtime-no-home |  | 1 | 1 | 1 |
| `o2_mplus1` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `o2_nmin1` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `owltab` | real(8) | runtime-no-home |  | 5 | 2 | 7 |
| `OxygenIntercept` | real(8) | runtime-no-home |  | 2 | 1 | 5 |
| `OxygenSlope` | real(8) | runtime-no-home |  | 2 | 1 | 5 |
| `paramvg` | real(8) | runtime-no-home |  | 3 | 1 | 5 |
| `pgrzdm` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `pmowdm` | real(8) | runtime-no-home |  | 3 | 1 | 5 |
| `q10_root` | real(8) | runtime-no-home |  | 1 | 1 | 4 |
| `qdraincomp` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `qssdi` | real(8) | runtime-no-home |  | 4 | 1 | 5 |
| `qssdisum` | real(8) | runtime-no-home |  | 2 | 1 | 3 |
| `rad` | real(8) | runtime-no-home |  | 8 | 3 | 7 |
| `reltr` | real(8) | runtime-no-home |  | 3 | 3 | 4 |
| `rid` | real(8) | runtime-no-home |  | 2 | 1 | 3 |
| `shape_factor_rootr` | real(8) | runtime-no-home |  | 1 | 1 | 4 |
| `sla` | real(8) | runtime-no-home |  | 2 | 2 | 4 |
| `slapot` | real(8) | runtime-no-home |  | 2 | 2 | 2 |
| `sptab` | real(8) | runtime-no-home |  | 4 | 1 | 8 |
| `ssdi_amount_f_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_amount_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_appl_rate_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_date_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_rate_f_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_sched_type_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_schedule_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_threshold_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `ssdi_threshold_z_irr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `sw_interval_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `swallo` | integer | runtime-no-home |  | 2 | 1 | 2 |
| `swdrought` | integer | runtime-no-home |  | 12 | 4 | 13 |
| `swinco` | integer | runtime-no-home |  | 15 | 3 | 14 |
| `swkimpl` | integer | runtime-no-home |  | 3 | 1 | 3 |
| `swmetdetail` | integer | runtime-no-home |  | 7 | 1 | 8 |
| `swrain` | integer | runtime-no-home |  | 6 | 1 | 8 |
| `swsnow` | integer | runtime-no-home |  | 7 | 1 | 7 |
| `swssdi_irr` | integer | runtime-no-home |  | 1 | 1 | 1 |
| `swtsum` | integer | runtime-no-home |  | 4 | 2 | 4 |
| `tadw` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `tadwpot` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `TBASEM` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `TEFFMX` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `tmn` | real(8) | runtime-no-home |  | 4 | 1 | 4 |
| `tmnr` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `tmx` | real(8) | runtime-no-home |  | 4 | 1 | 3 |
| `tsoil` | real(8) | runtime-no-home |  | 21 | 1 | 26 |
| `tsumemeopt` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `tsumgerm` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `tsunrise_atm` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `tsunset_atm` | real(8) | runtime-no-home |  | 1 | 1 | 1 |
| `twilt` | real(8) | runtime-no-home |  | 5 | 3 | 6 |
| `w_root_ss` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `widthr` | real(8) | runtime-no-home |  | 2 | 1 | 2 |
| `wrtmin` | real(8) | runtime-no-home |  | 3 | 2 | 3 |
| `zgerm` | real(8) | runtime-no-home |  | 3 | 1 | 3 |
| `zh` | real(8) | runtime-no-home |  | 2 | 1 | 0 |
| `aetr` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 3 |
| `ahum` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 2 | 4 |
| `arad` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 2 | 4 |
| `arai` | real(8) | state-homed-runtime | atmosphere_state.f90 | 4 | 2 | 5 |
| `atav` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 3 |
| `atmin7` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 2 |
| `atmn` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 4 |
| `atmx` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 4 |
| `awin` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 3 |
| `cdrain` | real(8) | state-homed-runtime | solute_state.f90 | 3 | 1 | 3 |
| `cfeic` | real(8) | state-homed-runtime | crop_fixed_state.f90 | 5 | 4 | 5 |
| `cfeictb` | real(8) | state-homed-runtime | crop_fixed_state.f90 | 6 | 2 | 5 |
| `cftb` | real(8) | state-homed-runtime | crop_fixed_state.f90 | 10 | 4 | 10 |
| `chtb` | real(8) | state-homed-runtime | crop_fixed_state.f90 | 10 | 4 | 10 |
| `CNdry` | real(8) | state-homed-runtime | atmosphere_state.f90 | 2 | 1 | 1 |
| `cropend` | real(8) | state-homed-runtime | crop_common_state.f90 | 3 | 1 | 2 |
| `cropstart` | real(8) | state-homed-runtime | crop_common_state.f90 | 7 | 1 | 6 |
| `cumdens` | real(8) | state-homed-runtime | crop_common_state.f90 | 9 | 5 | 9 |
| `cuptgraz` | real(8) | state-homed-runtime | crop_common_state.f90 | 4 | 2 | 4 |
| `cuptgrazpot` | real(8) | state-homed-runtime | crop_common_state.f90 | 4 | 2 | 4 |
| `dateharvest` | real(8) | state-homed-runtime | crop_grass_state.f90 | 3 | 2 | 4 |
| `daycrop` | integer | state-homed-runtime | crop_common_state.f90 | 8 | 3 | 12 |
| `daynrfirst` | integer | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 4 |
| `daynrlast` | integer | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 4 |
| `dtsolu` | real(8) | state-homed-runtime | solute_state.f90 | 2 | 2 | 2 |
| `epot` | real(8) | state-homed-runtime | atmosphere_state.f90 | 6 | 1 | 7 |
| `flCropCalendar` | logical | state-homed-runtime | crop_common_state.f90 | 6 | 3 | 8 |
| `flCropEmergence` | logical | state-homed-runtime | crop_state.f90 | 5 | 2 | 4 |
| `flCropGerm` | logical | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 2 |
| `flCropHarvest` | logical | state-homed-runtime | crop_common_state.f90 | 6 | 2 | 6 |
| `flCropNut` | logical | state-homed-runtime | crop_common_state.f90 | 6 | 2 | 8 |
| `flCropOutput` | logical | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 3 |
| `flCropPrep` | logical | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 2 |
| `flCropReadFile` | logical | state-homed-runtime | crop_common_state.f90 | 1 | 1 | 1 |
| `flCropSow` | logical | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 2 |
| `flHarvestDay` | logical | state-homed-runtime | crop_common_state.f90 | 3 | 2 | 4 |
| `grain` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 2 |
| `h_enpr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 | 4 | 3 | 9 |
| `icn_atm` | integer | state-homed-runtime | atmosphere_state.f90 | 1 | 1 | 1 |
| `icrop` | integer | state-homed-runtime | crop_common_state.f90 | 14 | 3 | 14 |
| `isqbot` | real(8) | state-homed-runtime | solute_state.f90 | 2 | 2 | 2 |
| `isqtop` | real(8) | state-homed-runtime | solute_state.f90 | 2 | 2 | 2 |
| `isua` | integer | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 4 |
| `ksatthr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 | 3 | 2 | 2 |
| `nod10_cn` | integer | state-homed-runtime | atmosphere_state.f90 | 1 | 1 | 1 |
| `noddrz` | integer | state-homed-runtime | crop_common_state.f90 | 7 | 2 | 9 |
| `nofd` | integer | state-homed-runtime | atmosphere_state.f90 | 4 | 3 | 3 |
| `nrain` | real(8) | state-homed-runtime | atmosphere_state.f90 | 3 | 1 | 2 |
| `PrepDelay` | integer | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 2 |
| `relsatthr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 | 3 | 2 | 3 |
| `rottot` | real(8) | state-homed-runtime | solute_state.f90 | 3 | 2 | 3 |
| `samini` | real(8) | state-homed-runtime | solute_state.f90 | 2 | 2 | 2 |
| `seqgrazmow` | integer | state-homed-runtime | crop_grass_state.f90 | 4 | 2 | 4 |
| `seqgrazmowpot` | integer | state-homed-runtime | crop_grass_state.f90 | 2 | 2 | 2 |
| `SowDelay` | integer | state-homed-runtime | crop_common_state.f90 | 2 | 2 | 2 |
| `sqdra` | real(8) | state-homed-runtime | solute_state.f90 | 3 | 2 | 3 |
| `swbulb` | integer | state-homed-runtime | crop_wofost_state.f90 | 6 | 3 | 6 |
| `swcf` | integer | state-homed-runtime | crop_state.f90 | 12 | 4 | 12 |
| `swcfbs` | integer | state-homed-runtime | crop_state.f90 | 3 | 1 | 3 |
| `swcrp` | integer | state-homed-runtime | crop_common_state.f90 | 2 | 1 | 4 |
| `swfrost` | integer | state-homed-runtime | soilwater_state.f90 | 6 | 1 | 10 |
| `swpotrelmf` | integer | state-homed-runtime | crop_grass_state.f90 | 6 | 3 | 6 |
| `tav` | real(8) | state-homed-runtime | atmosphere_state.f90 | 12 | 1 | 16 |
| `tpot` | real(8) | state-homed-runtime | atmosphere_state.f90 | 6 | 1 | 7 |
| `wet` | real(8) | state-homed-runtime | atmosphere_state.f90 | 10 | 1 | 5 |
| `avevaptb` | real(8) | state-homed-unknown | atmosphere_state.f90 | 2 | 0 | 1 |
| `avprectb` | real(8) | state-homed-unknown | atmosphere_state.f90 | 2 | 0 | 1 |
| `pfreetb` | real(8) | state-homed-unknown | atmosphere_state.f90 | 2 | 0 | 1 |
| `pstemtb` | real(8) | state-homed-unknown | atmosphere_state.f90 | 2 | 0 | 1 |
| `scanopytb` | real(8) | state-homed-unknown | atmosphere_state.f90 | 2 | 0 | 1 |
