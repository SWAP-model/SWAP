# variables.f90 audit

Total active declarations: 616

- **config-side**: 222
- **runtime-no-home**: 173
- **orphan-read-only**: 141
- **state-homed-runtime**: 74
- **state-homed-config-only**: 4
- **dead-no-readers**: 2

## Per-symbol table

| Symbol | Type | Category | State home |
|---|---|---|---|
| `adcrh` | real(8) | config-side |  |
| `adcrl` | real(8) | config-side |  |
| `aeratecrit` | real(8) | config-side |  |
| `air_filled_root_por` | real(8) | config-side |  |
| `alt` | real(8) | config-side |  |
| `altw` | real(8) | config-side |  |
| `amaxtb` | real(8) | config-side |  |
| `angstroma` | real(8) | config-side |  |
| `angstromb` | real(8) | config-side |  |
| `aqamp` | real(8) | config-side |  |
| `aqave` | real(8) | config-side |  |
| `aqper` | real(8) | config-side |  |
| `aqtmax` | real(8) | config-side |  |
| `basegw` | real(8) | config-side |  |
| `bexp` | real(8) | config-side |  |
| `cfevappond` | real(8) | config-side |  |
| `cofintfl` | real(8) | config-side |  |
| `cofqha` | real(8) | config-side |  |
| `cofqhb` | real(8) | config-side |  |
| `cofqhc` | real(8) | config-side |  |
| `cofred` | real(8) | config-side |  |
| `cpre` | real(8) | config-side |  |
| `cref` | real(8) | config-side |  |
| `CritDevPondDt` | real(8) | config-side |  |
| `cropfil` | character(len=40) | config-side |  |
| `croptype` | integer | config-side |  |
| `cseeptab` | real(8) | config-side |  |
| `cvl` | real(8) | config-side |  |
| `cvo` | real(8) | config-side |  |
| `cvr` | real(8) | config-side |  |
| `cvs` | real(8) | config-side |  |
| `daquif` | real(8) | config-side |  |
| `dcritrtz` | real(8) | config-side |  |
| `ddamp` | real(8) | config-side |  |
| `ddif` | real(8) | config-side |  |
| `decpot` | real(8) | config-side |  |
| `decsat` | real(8) | config-side |  |
| `DelayRegrowthTab` | real(8) | config-side |  |
| `dlc` | real(8) | config-side |  |
| `dlo` | real(8) | config-side |  |
| `dmmowtb` | real(8) | config-side |  |
| `dramet` | integer | config-side |  |
| `drfil` | character(len=16) | config-side |  |
| `dry_mat_cont_roots` | real(8) | config-side |  |
| `dtsmtb` | real(8) | config-side |  |
| `dvsend` | real(8) | config-side |  |
| `dvsnlt` | real(8) | config-side |  |
| `eff` | real(8) | config-side |  |
| `entres` | real(8) | config-side |  |
| `fdepth` | real(8) | config-side |  |
| `flCO2` | logical | config-side |  |
| `flrunon` | logical | config-side |  |
| `fltb` | real(8) | config-side |  |
| `fotb` | real(8) | config-side |  |
| `fraharlosorm_lv` | real(8) | config-side |  |
| `frexp` | real(8) | config-side |  |
| `frtb` | real(8) | config-side |  |
| `fstb` | real(8) | config-side |  |
| `ftopdislay` | real(8) | config-side |  |
| `gampar` | real(8) | config-side |  |
| `geofac` | real(8) | config-side |  |
| `gwlconv` | real(8) | config-side |  |
| `gwli` | real(8) | config-side |  |
| `gwltab` | real(8) | config-side |  |
| `haqtab` | real(8) | config-side |  |
| `hbotab` | real(8) | config-side |  |
| `hcomp` | real(8) | config-side |  |
| `hdrain` | real(8) | config-side |  |
| `hlim1` | real(8) | config-side |  |
| `hlim2l` | real(8) | config-side |  |
| `hlim2u` | real(8) | config-side |  |
| `hlim3h` | real(8) | config-side |  |
| `hlim3l` | real(8) | config-side |  |
| `hlim4` | real(8) | config-side |  |
| `hplate` | real(8) | config-side |  |
| `hsublay` | real(8) | config-side |  |
| `idev` | integer | config-side |  |
| `idsl` | integer | config-side |  |
| `ilnmxl` | integer | config-side |  |
| `impend` | real(8) | config-side |  |
| `InList_csv_tz` | character(len=1024) | config-side |  |
| `intwl` | integer | config-side |  |
| `irconc` | real(8) | config-side |  |
| `irdate` | real(8) | config-side |  |
| `irdepth` | real(8) | config-side |  |
| `irtype` | integer | config-side |  |
| `isoillay` | integer | config-side |  |
| `kf` | real(8) | config-side |  |
| `kfsat` | real(8) | config-side |  |
| `khbot` | real(8) | config-side |  |
| `khtop` | real(8) | config-side |  |
| `kvbot` | real(8) | config-side |  |
| `kvtop` | real(8) | config-side |  |
| `lat` | real(8) | config-side |  |
| `ldis` | real(8) | config-side |  |
| `lrnr` | real(8) | config-side |  |
| `MaxBackTr` | integer | config-side |  |
| `metfil` | character(len=200) | config-side |  |
| `ncomp` | integer | config-side |  |
| `nconc` | integer | config-side |  |
| `nhead` | integer | config-side |  |
| `nlai` | real(8) | config-side |  |
| `nlue` | real(8) | config-side |  |
| `nmetdetail` | integer | config-side |  |
| `nmxlv` | real(8) | config-side |  |
| `nowltab` | integer | config-side |  |
| `nrstaring` | integer | config-side |  |
| `nsla` | real(8) | config-side |  |
| `nsublay` | integer | config-side |  |
| `NumLevRapDra` | integer | config-side |  |
| `osswlm` | real(8) | config-side |  |
| `outdatint` | real(8) | config-side |  |
| `outfil` | character(len=16) | config-side |  |
| `pathatm` | character(len=80) | config-side |  |
| `pathcrop` | character(len=80) | config-side |  |
| `pathdrain` | character(len=80) | config-side |  |
| `pathwork` | character(len=80) | config-side |  |
| `perdl` | real(8) | config-side |  |
| `poros` | real(8) | config-side |  |
| `project` | character(len=80) | config-side |  |
| `q10` | real(8) | config-side |  |
| `q10_microbial` | real(8) | config-side |  |
| `qbotab` | real(8) | config-side |  |
| `raintab` | real(8) | config-side |  |
| `rdctb` | real(8) | config-side |  |
| `rdmax` | real(8) | config-side |  |
| `rdrrtb` | real(8) | config-side |  |
| `rdrstb` | real(8) | config-side |  |
| `rdtb` | real(8) | config-side |  |
| `rfsetb` | real(8) | config-side |  |
| `rgrlai` | real(8) | config-side |  |
| `rimlay` | real(8) | config-side |  |
| `rlwtb` | real(8) | config-side |  |
| `rml` | real(8) | config-side |  |
| `rmo` | real(8) | config-side |  |
| `rmr` | real(8) | config-side |  |
| `rms` | real(8) | config-side |  |
| `root_radiusO2` | real(8) | config-side |  |
| `rsigni` | real(8) | config-side |  |
| `rsoil` | real(8) | config-side |  |
| `rsurfshallow` | real(8) | config-side |  |
| `rsw` | real(8) | config-side |  |
| `rtheta` | real(8) | config-side |  |
| `salthead` | real(8) | config-side |  |
| `saltmax` | real(8) | config-side |  |
| `saltslope` | real(8) | config-side |  |
| `schedule` | integer | config-side |  |
| `shape` | real(8) | config-side |  |
| `sinamp` | real(8) | config-side |  |
| `sinave` | real(8) | config-side |  |
| `sinmax` | real(8) | config-side |  |
| `slatb` | real(8) | config-side |  |
| `snowcoef` | real(8) | config-side |  |
| `spa` | real(8) | config-side |  |
| `span` | real(8) | config-side |  |
| `spec_weight_root_tissue` | real(8) | config-side |  |
| `specific_resp_humus` | real(8) | config-side |  |
| `srl` | real(8) | config-side |  |
| `ssa` | real(8) | config-side |  |
| `sw2` | integer | config-side |  |
| `sw3` | integer | config-side |  |
| `sw4` | integer | config-side |  |
| `swbotb3Impl` | integer | config-side |  |
| `swbotbc` | integer | config-side |  |
| `swbotbhea` | integer | config-side |  |
| `swbr` | integer | config-side |  |
| `swcalt` | integer | config-side |  |
| `swcaprise` | logical | config-side |  |
| `swcofqhc` | integer | config-side |  |
| `swcompensate` | integer | config-side |  |
| `swdiscrvert` | integer | config-side |  |
| `swdislay` | integer | config-side |  |
| `swdivide` | integer | config-side |  |
| `swdmi2rd` | integer | config-side |  |
| `swetr` | integer | config-side |  |
| `swetsine` | integer | config-side |  |
| `swgc` | integer | config-side |  |
| `swharv` | integer | config-side |  |
| `swhea` | integer | config-side |  |
| `swhyst` | integer | config-side |  |
| `swinter` | integer | config-side |  |
| `swirfix` | integer | config-side |  |
| `swkmean` | integer | config-side |  |
| `swliminf` | integer | config-side |  |
| `swoxygen` | integer | config-side |  |
| `swoxygentype` | integer | config-side |  |
| `swqhbot` | integer | config-side |  |
| `swqhr` | integer | config-side |  |
| `swrd` | integer | config-side |  |
| `swrdc` | integer | config-side |  |
| `swredu` | integer | config-side |  |
| `swrootradius` | integer | config-side |  |
| `swsalinity` | integer | config-side |  |
| `swsolu` | integer | config-side |  |
| `swsrf` | integer | config-side |  |
| `swstressor` | integer | config-side |  |
| `swsublim` | integer | config-side |  |
| `swtopbhea` | integer | config-side |  |
| `swtopdislay` | integer | config-side |  |
| `swWrtNonox` | integer | config-side |  |
| `taccur` | real(8) | config-side |  |
| `tampli` | real(8) | config-side |  |
| `tbase` | real(8) | config-side |  |
| `tdwi` | real(8) | config-side |  |
| `tembtab` | real(8) | config-side |  |
| `temtoptab` | real(8) | config-side |  |
| `tfrostend` | real(8) | config-side |  |
| `tfroststa` | real(8) | config-side |  |
| `till_swtill` | integer | config-side |  |
| `timref` | real(8) | config-side |  |
| `tmean` | real(8) | config-side |  |
| `tmnftb` | real(8) | config-side |  |
| `tmpftb` | real(8) | config-side |  |
| `tscf` | real(8) | config-side |  |
| `tsumam` | real(8) | config-side |  |
| `tsumea` | real(8) | config-side |  |
| `var_a` | real(8) | config-side |  |
| `wldip` | real(8) | config-side |  |
| `wrtmax` | real(8) | config-side |  |
| `zc` | real(8) | config-side |  |
| `zi` | real(8) | config-side |  |
| `zintf` | real(8) | config-side |  |
| `cQMpOutDrRap` | real(8) | dead-no-readers |  |
| `inifil` | character(len=200) | dead-no-readers |  |
| `ad` | integer | orphan-read-only |  |
| `am` | integer | orphan-read-only |  |
| `ArMpTp` | real(8) | orphan-read-only |  |
| `BiModal` | logical | orphan-read-only |  |
| `cirrs` | real(8) | orphan-read-only |  |
| `cirrthres` | real(8) | orphan-read-only |  |
| `co2amaxtb` | real(8) | orphan-read-only |  |
| `co2efftb` | real(8) | orphan-read-only |  |
| `co2ppm` | real(8) | orphan-read-only |  |
| `co2tratb` | real(8) | orphan-read-only |  |
| `co2year` | integer | orphan-read-only |  |
| `CritDevMasBal` | real(8) | orphan-read-only |  |
| `CriterHr` | real(8) | orphan-read-only |  |
| `CritUndSatVol` | real(8) | orphan-read-only |  |
| `crp` | integer | orphan-read-only |  |
| `DaysGrazingtab` | real(8) | orphan-read-only |  |
| `dcrit` | real(8) | orphan-read-only |  |
| `dev_cmb` | integer | orphan-read-only |  |
| `dewrest` | real(8) | orphan-read-only |  |
| `DiPoCp` | real(8) | orphan-read-only |  |
| `ditab` | real(8) | orphan-read-only |  |
| `dmgrztb` | real(8) | orphan-read-only |  |
| `dra` | integer | orphan-read-only |  |
| `dropr` | real(8) | orphan-read-only |  |
| `dwatab` | real(8) | orphan-read-only |  |
| `fbltb` | real(8) | orphan-read-only |  |
| `fidtab` | real(8) | orphan-read-only |  |
| `FlDecMpRat` | logical | orphan-read-only |  |
| `fldumpconvcrit` | logical | orphan-read-only |  |
| `flSwapShared` | logical | orphan-read-only |  |
| `gctb` | real(8) | orphan-read-only |  |
| `gwlcrit` | real(8) | orphan-read-only |  |
| `hcritab` | real(8) | orphan-read-only |  |
| `hPrep` | real(8) | orphan-read-only |  |
| `hqhtab` | real(8) | orphan-read-only |  |
| `hSow` | real(8) | orphan-read-only |  |
| `IAvFrMpWlWtDm1` | real(8) | orphan-read-only |  |
| `IAvFrMpWlWtDm2` | real(8) | orphan-read-only |  |
| `ientrytablay` | integer | orphan-read-only |  |
| `inc` | integer | orphan-read-only |  |
| `iQExcMtxDm1Cp` | real(8) | orphan-read-only |  |
| `iQExcMtxDm2Cp` | real(8) | orphan-read-only |  |
| `iqinfmax` | real(8) | orphan-read-only |  |
| `iQMpOutDrRap` | real(8) | orphan-read-only |  |
| `iQOutDrRapCp` | real(8) | orphan-read-only |  |
| `isuas` | integer | orphan-read-only |  |
| `Itnumb` | integer | orphan-read-only |  |
| `Kroot` | real(8) | orphan-read-only |  |
| `kstem` | real(8) | orphan-read-only |  |
| `logf` | integer | orphan-read-only |  |
| `LossGrazingtab` | real(8) | orphan-read-only |  |
| `lossgrztab` | real(8) | orphan-read-only |  |
| `lossmowtab` | real(8) | orphan-read-only |  |
| `lsda` | real(8) | orphan-read-only |  |
| `MaxPrepDelay` | integer | orphan-read-only |  |
| `MaxSowDelay` | integer | orphan-read-only |  |
| `mrftb` | real(8) | orphan-read-only |  |
| `NoVap` | logical | orphan-read-only |  |
| `nphase` | integer | orphan-read-only |  |
| `numtablay` | integer | orphan-read-only |  |
| `o2_bunsencoeff` | real(8) | orphan-read-only |  |
| `o2_c_macro` | real(8) | orphan-read-only |  |
| `o2_c_min_micro` | real(8) | orphan-read-only |  |
| `o2_ctopnode` | real(8) | orphan-read-only |  |
| `o2_d_o2inwater` | real(8) | orphan-read-only |  |
| `o2_d_root` | real(8) | orphan-read-only |  |
| `o2_d_soil` | real(8) | orphan-read-only |  |
| `o2_depth` | real(8) | orphan-read-only |  |
| `o2_gas_filled_porosity` | real(8) | orphan-read-only |  |
| `o2_perc_org_mat` | real(8) | orphan-read-only |  |
| `o2_r_microbial_z0` | real(8) | orphan-read-only |  |
| `o2_root_radius` | real(8) | orphan-read-only |  |
| `o2_sat_water_cont` | real(8) | orphan-read-only |  |
| `o2_shape_factor_microbialr` | real(8) | orphan-read-only |  |
| `o2_soil_density` | real(8) | orphan-read-only |  |
| `o2_soil_temp` | real(8) | orphan-read-only |  |
| `o2_w_root` | real(8) | orphan-read-only |  |
| `o2_w_root_z0` | real(8) | orphan-read-only |  |
| `o2_waterfilm_thickness` | real(8) | orphan-read-only |  |
| `outdat` | real(8) | orphan-read-only |  |
| `perirrsurp` | real(8) | orphan-read-only |  |
| `pld` | real(8) | orphan-read-only |  |
| `pondmxtab` | real(8) | orphan-read-only |  |
| `qdrtab` | real(8) | orphan-read-only |  |
| `QExcMpMtx` | real(8) | orphan-read-only |  |
| `qimmob` | real(8) | orphan-read-only |  |
| `QMaPo` | real(8) | orphan-read-only |  |
| `QRapDra` | real(8) | orphan-read-only |  |
| `raithreshold` | real(8) | orphan-read-only |  |
| `rawtab` | real(8) | orphan-read-only |  |
| `remoc` | real(8) | orphan-read-only |  |
| `rootcoefa` | real(8) | orphan-read-only |  |
| `rooteff` | real(8) | orphan-read-only |  |
| `rootradius` | real(8) | orphan-read-only |  |
| `rot` | integer | orphan-read-only |  |
| `runonarr` | real(8) | orphan-read-only |  |
| `Rxylem` | real(8) | orphan-read-only |  |
| `siccaplai` | real(8) | orphan-read-only |  |
| `siccaptb` | real(8) | orphan-read-only |  |
| `snw` | integer | orphan-read-only |  |
| `sptablay` | real(8) | orphan-read-only |  |
| `StepHr` | real(8) | orphan-read-only |  |
| `SwBotb3ResVert` | integer | orphan-read-only |  |
| `swcirrthres` | integer | orphan-read-only |  |
| `swpondmx` | integer | orphan-read-only |  |
| `swtopsub` | integer | orphan-read-only |  |
| `swuseCN` | integer | orphan-read-only |  |
| `tau` | real(8) | orphan-read-only |  |
| `tawtab` | real(8) | orphan-read-only |  |
| `tcritab` | real(8) | orphan-read-only |  |
| `tem` | integer | orphan-read-only |  |
| `TempSow` | real(8) | orphan-read-only |  |
| `tendirrig` | real(8) | orphan-read-only |  |
| `till_i_n_model` | integer | orphan-read-only |  |
| `till_iRedist` | integer | orphan-read-only |  |
| `till_Max_Z_tillage` | real(8) | orphan-read-only |  |
| `till_Ntill` | integer | orphan-read-only |  |
| `till_Ntypes` | integer | orphan-read-only |  |
| `treltab` | real(8) | orphan-read-only |  |
| `tstairrig` | real(8) | orphan-read-only |  |
| `tsumdepth` | real(8) | orphan-read-only |  |
| `tsumtemp` | real(8) | orphan-read-only |  |
| `tsumtime` | integer | orphan-read-only |  |
| `tz_z1_z2` | real(8) | orphan-read-only |  |
| `UptGrazingtab` | real(8) | orphan-read-only |  |
| `vernbase` | real(8) | orphan-read-only |  |
| `verndvs` | real(8) | orphan-read-only |  |
| `vernrtb` | real(8) | orphan-read-only |  |
| `vernsat` | real(8) | orphan-read-only |  |
| `VlMpStDm1` | real(8) | orphan-read-only |  |
| `VlMpStDm2` | real(8) | orphan-read-only |  |
| `wiltpoint` | real(8) | orphan-read-only |  |
| `wlsman` | real(8) | orphan-read-only |  |
| `wlstab` | real(8) | orphan-read-only |  |
| `wrtb` | real(8) | orphan-read-only |  |
| `Z_Tp` | real(8) | orphan-read-only |  |
| `zgrz` | real(8) | orphan-read-only |  |
| `zmow` | real(8) | orphan-read-only |  |
| `zPrep` | real(8) | orphan-read-only |  |
| `zSow` | real(8) | orphan-read-only |  |
| `zTempSow` | real(8) | orphan-read-only |  |
| `Agedrain` | real(8) | runtime-no-home |  |
| `Ageirr` | real(8) | runtime-no-home |  |
| `Agepond` | real(8) | runtime-no-home |  |
| `Agepondm1` | real(8) | runtime-no-home |  |
| `Agepre` | real(8) | runtime-no-home |  |
| `agerm` | real(8) | runtime-no-home |  |
| `alphacrit` | real(8) | runtime-no-home |  |
| `atmtr` | real(8) | runtime-no-home |  |
| `bdens` | real(8) | runtime-no-home |  |
| `bgerm` | real(8) | runtime-no-home |  |
| `botcom` | integer | runtime-no-home |  |
| `c_mroot` | real(8) | runtime-no-home |  |
| `c_top` | real(8) | runtime-no-home |  |
| `cgerm` | real(8) | runtime-no-home |  |
| `cirr` | real(8) | runtime-no-home |  |
| `cml` | real(8) | runtime-no-home |  |
| `cmsy` | real(8) | runtime-no-home |  |
| `cQMpLatSs` | real(8) | runtime-no-home |  |
| `CritDevh1Cp` | real(8) | runtime-no-home |  |
| `CritDevh2Cp` | real(8) | runtime-no-home |  |
| `dayfix` | integer | runtime-no-home |  |
| `daygrowth` | integer | runtime-no-home |  |
| `daygrowthpot` | integer | runtime-no-home |  |
| `daylp` | real(8) | runtime-no-home |  |
| `days_counter_irr` | integer | runtime-no-home |  |
| `days_interval_irr` | integer | runtime-no-home |  |
| `dethum` | real(8) | runtime-no-home |  |
| `detrad` | real(8) | runtime-no-home |  |
| `detrain` | real(8) | runtime-no-home |  |
| `detrecord` | integer | runtime-no-home |  |
| `dettav` | real(8) | runtime-no-home |  |
| `dettime` | real(8) | runtime-no-home |  |
| `detwind` | real(8) | runtime-no-home |  |
| `dhPrep` | real(8) | runtime-no-home |  |
| `dhSow` | real(8) | runtime-no-home |  |
| `difpp` | real(8) | runtime-no-home |  |
| `drares` | real(8) | runtime-no-home |  |
| `drbl` | real(8) | runtime-no-home |  |
| `drblpot` | real(8) | runtime-no-home |  |
| `dsinbe` | real(8) | runtime-no-home |  |
| `dt_SSDI_event` | real(8) | runtime-no-home |  |
| `dtempSow` | real(8) | runtime-no-home |  |
| `dtEventRain` | real(8) | runtime-no-home |  |
| `dznew` | real(8) | runtime-no-home |  |
| `f_senes` | real(8) | runtime-no-home |  |
| `fbl` | real(8) | runtime-no-home |  |
| `finterception` | real(8) | runtime-no-home |  |
| `flanthesis` | logical | runtime-no-home |  |
| `flCropOpenFile` | logical | runtime-no-home |  |
| `flGrazing` | logical | runtime-no-home |  |
| `flGrazingpot` | logical | runtime-no-home |  |
| `flHarvest` | logical | runtime-no-home |  |
| `flHarvestpot` | logical | runtime-no-home |  |
| `flhrvendact` | logical | runtime-no-home |  |
| `flhrvendpot` | logical | runtime-no-home |  |
| `FlHydrLift` | logical | runtime-no-home |  |
| `flwarn_hc` | logical | runtime-no-home |  |
| `fstr` | real(8) | runtime-no-home |  |
| `gasst` | real(8) | runtime-no-home |  |
| `gasstpot` | real(8) | runtime-no-home |  |
| `gc` | real(8) | runtime-no-home |  |
| `glaiex` | real(8) | runtime-no-home |  |
| `glaiexpot` | real(8) | runtime-no-home |  |
| `grzdm` | real(8) | runtime-no-home |  |
| `gwlinf` | real(8) | runtime-no-home |  |
| `gwrt` | real(8) | runtime-no-home |  |
| `hdrygerm` | real(8) | runtime-no-home |  |
| `hwetgerm` | real(8) | runtime-no-home |  |
| `icAgetopdwn` | real(8) | runtime-no-home |  |
| `icAgetopupw` | real(8) | runtime-no-home |  |
| `idaysgraz` | integer | runtime-no-home |  |
| `idaysgrazpot` | integer | runtime-no-home |  |
| `idregr` | integer | runtime-no-home |  |
| `idregrpot` | integer | runtime-no-home |  |
| `ientrytab` | integer | runtime-no-home |  |
| `iharvest` | integer | runtime-no-home |  |
| `iHWCKmodel` | integer | runtime-no-home |  |
| `ilvold` | integer | runtime-no-home |  |
| `ilvoldpot` | integer | runtime-no-home |  |
| `infres` | real(8) | runtime-no-home |  |
| `InList_csv` | character(len=1024) | runtime-no-home |  |
| `inpola` | real(8) | runtime-no-home |  |
| `inpolb` | real(8) | runtime-no-home |  |
| `ipos` | integer | runtime-no-home |  |
| `irectotal` | integer | runtime-no-home |  |
| `irrigevent` | integer | runtime-no-home |  |
| `iseqgm` | integer | runtime-no-home |  |
| `iseqgmpot` | integer | runtime-no-home |  |
| `issnowbeg` | real(8) | runtime-no-home |  |
| `iwarn_hc` | integer | runtime-no-home |  |
| `laiem` | real(8) | runtime-no-home |  |
| `laiexp` | real(8) | runtime-no-home |  |
| `laiexppot` | real(8) | runtime-no-home |  |
| `laimax` | real(8) | runtime-no-home |  |
| `lv` | real(8) | runtime-no-home |  |
| `lvage` | real(8) | runtime-no-home |  |
| `lvagepot` | real(8) | runtime-no-home |  |
| `lvpot` | real(8) | runtime-no-home |  |
| `max_resp_factor` | real(8) | runtime-no-home |  |
| `mowdm` | real(8) | runtime-no-home |  |
| `mrest` | real(8) | runtime-no-home |  |
| `mrestpot` | real(8) | runtime-no-home |  |
| `nird` | real(8) | runtime-no-home |  |
| `nirri` | integer | runtime-no-home |  |
| `nirri_ssdi_irr` | integer | runtime-no-home |  |
| `nmper` | integer | runtime-no-home |  |
| `nod1lay` | integer | runtime-no-home |  |
| `nod_ssdi_irr` | integer | runtime-no-home |  |
| `nod_ssdi_sensor_irr` | integer | runtime-no-home |  |
| `numbit` | integer | runtime-no-home |  |
| `numlay` | integer | runtime-no-home |  |
| `numnodnew` | integer | runtime-no-home |  |
| `numtab` | integer | runtime-no-home |  |
| `o2_capac_term` | real(8) | runtime-no-home |  |
| `o2_d_soil_term1` | real(8) | runtime-no-home |  |
| `o2_d_soil_term2` | real(8) | runtime-no-home |  |
| `o2_gfp100` | real(8) | runtime-no-home |  |
| `o2_ini_stress` | logical | runtime-no-home |  |
| `o2_mplus1` | real(8) | runtime-no-home |  |
| `o2_nmin1` | real(8) | runtime-no-home |  |
| `owltab` | real(8) | runtime-no-home |  |
| `OxygenIntercept` | real(8) | runtime-no-home |  |
| `OxygenSlope` | real(8) | runtime-no-home |  |
| `paramvg` | real(8) | runtime-no-home |  |
| `pgrzdm` | real(8) | runtime-no-home |  |
| `pmowdm` | real(8) | runtime-no-home |  |
| `q10_root` | real(8) | runtime-no-home |  |
| `qdraincomp` | real(8) | runtime-no-home |  |
| `qssdi` | real(8) | runtime-no-home |  |
| `qssdisum` | real(8) | runtime-no-home |  |
| `rad` | real(8) | runtime-no-home |  |
| `reltr` | real(8) | runtime-no-home |  |
| `rid` | real(8) | runtime-no-home |  |
| `shape_factor_rootr` | real(8) | runtime-no-home |  |
| `sla` | real(8) | runtime-no-home |  |
| `slapot` | real(8) | runtime-no-home |  |
| `sptab` | real(8) | runtime-no-home |  |
| `ssdi_amount_f_irr` | real(8) | runtime-no-home |  |
| `ssdi_amount_irr` | real(8) | runtime-no-home |  |
| `ssdi_appl_rate_irr` | real(8) | runtime-no-home |  |
| `ssdi_date_irr` | real(8) | runtime-no-home |  |
| `ssdi_rate_f_irr` | real(8) | runtime-no-home |  |
| `ssdi_sched_type_irr` | integer | runtime-no-home |  |
| `ssdi_schedule_irr` | integer | runtime-no-home |  |
| `ssdi_threshold_irr` | real(8) | runtime-no-home |  |
| `ssdi_threshold_z_irr` | real(8) | runtime-no-home |  |
| `sw_interval_irr` | integer | runtime-no-home |  |
| `swdrought` | integer | runtime-no-home |  |
| `swinco` | integer | runtime-no-home |  |
| `swkimpl` | integer | runtime-no-home |  |
| `swmetdetail` | integer | runtime-no-home |  |
| `swrain` | integer | runtime-no-home |  |
| `swsnow` | integer | runtime-no-home |  |
| `swssdi_irr` | integer | runtime-no-home |  |
| `swtsum` | integer | runtime-no-home |  |
| `tadw` | real(8) | runtime-no-home |  |
| `tadwpot` | real(8) | runtime-no-home |  |
| `TBASEM` | real(8) | runtime-no-home |  |
| `TEFFMX` | real(8) | runtime-no-home |  |
| `tmn` | real(8) | runtime-no-home |  |
| `tmnr` | real(8) | runtime-no-home |  |
| `tmx` | real(8) | runtime-no-home |  |
| `tsoil` | real(8) | runtime-no-home |  |
| `tsumemeopt` | real(8) | runtime-no-home |  |
| `tsumgerm` | real(8) | runtime-no-home |  |
| `tsunrise_atm` | real(8) | runtime-no-home |  |
| `tsunset_atm` | real(8) | runtime-no-home |  |
| `twilt` | real(8) | runtime-no-home |  |
| `w_root_ss` | real(8) | runtime-no-home |  |
| `widthr` | real(8) | runtime-no-home |  |
| `wrtmin` | real(8) | runtime-no-home |  |
| `zgerm` | real(8) | runtime-no-home |  |
| `zh` | real(8) | runtime-no-home |  |
| `rsro` | real(8) | state-homed-config-only | surfacewater_state.f90 |
| `rsroexp` | real(8) | state-homed-config-only | surfacewater_state.f90 |
| `swdra` | integer | state-homed-config-only | surfacewater_state.f90 |
| `swsophy` | integer | state-homed-config-only | soilwater_state.f90 |
| `aetr` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `ahum` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `arad` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `arai` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `atav` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `atmin7` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `atmn` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `atmx` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `avevaptb` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `avprectb` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `awin` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `cdrain` | real(8) | state-homed-runtime | solute_state.f90 |
| `cfeic` | real(8) | state-homed-runtime | crop_fixed_state.f90 |
| `cfeictb` | real(8) | state-homed-runtime | crop_fixed_state.f90 |
| `cftb` | real(8) | state-homed-runtime | crop_fixed_state.f90 |
| `chtb` | real(8) | state-homed-runtime | crop_fixed_state.f90 |
| `CNdry` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `CNrefTAB` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `cropend` | real(8) | state-homed-runtime | crop_common_state.f90 |
| `cropstart` | real(8) | state-homed-runtime | crop_common_state.f90 |
| `cumdens` | real(8) | state-homed-runtime | crop_common_state.f90 |
| `cuptgraz` | real(8) | state-homed-runtime | crop_common_state.f90 |
| `cuptgrazpot` | real(8) | state-homed-runtime | crop_common_state.f90 |
| `dateharvest` | real(8) | state-homed-runtime | crop_grass_state.f90 |
| `daycrop` | integer | state-homed-runtime | crop_common_state.f90 |
| `daynrfirst` | integer | state-homed-runtime | atmosphere_state.f90 |
| `daynrlast` | integer | state-homed-runtime | atmosphere_state.f90 |
| `dtsolu` | real(8) | state-homed-runtime | solute_state.f90 |
| `epot` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `flCropCalendar` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropEmergence` | logical | state-homed-runtime | crop_state.f90 |
| `flCropGerm` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropHarvest` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropNut` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropOutput` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropPrep` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropReadFile` | logical | state-homed-runtime | crop_common_state.f90 |
| `flCropSow` | logical | state-homed-runtime | crop_common_state.f90 |
| `flHarvestDay` | logical | state-homed-runtime | crop_common_state.f90 |
| `grain` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `h_enpr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 |
| `icn_atm` | integer | state-homed-runtime | atmosphere_state.f90 |
| `iCNtab` | integer | state-homed-runtime | atmosphere_state.f90 |
| `icrop` | integer | state-homed-runtime | crop_common_state.f90 |
| `isqbot` | real(8) | state-homed-runtime | solute_state.f90 |
| `isqtop` | real(8) | state-homed-runtime | solute_state.f90 |
| `isua` | integer | state-homed-runtime | atmosphere_state.f90 |
| `ksatthr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 |
| `nod10_cn` | integer | state-homed-runtime | atmosphere_state.f90 |
| `noddrz` | integer | state-homed-runtime | crop_common_state.f90 |
| `nofd` | integer | state-homed-runtime | atmosphere_state.f90 |
| `nrain` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `pfreetb` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `pondmx` | real(8) | state-homed-runtime | surfacewater_state.f90 |
| `PrepDelay` | integer | state-homed-runtime | crop_common_state.f90 |
| `pstemtb` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `relsatthr` | real(8) | state-homed-runtime | hydraulic_params_mod.f90 |
| `rottot` | real(8) | state-homed-runtime | solute_state.f90 |
| `samini` | real(8) | state-homed-runtime | solute_state.f90 |
| `scanopytb` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `seqgrazmow` | integer | state-homed-runtime | crop_grass_state.f90 |
| `seqgrazmowpot` | integer | state-homed-runtime | crop_grass_state.f90 |
| `SowDelay` | integer | state-homed-runtime | crop_common_state.f90 |
| `sqdra` | real(8) | state-homed-runtime | solute_state.f90 |
| `swbulb` | integer | state-homed-runtime | crop_wofost_state.f90 |
| `swcf` | integer | state-homed-runtime | crop_state.f90 |
| `swcfbs` | integer | state-homed-runtime | crop_state.f90 |
| `swcrp` | integer | state-homed-runtime | crop_common_state.f90 |
| `swfrost` | integer | state-homed-runtime | soilwater_state.f90 |
| `swpotrelmf` | integer | state-homed-runtime | crop_grass_state.f90 |
| `tav` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `tpot` | real(8) | state-homed-runtime | atmosphere_state.f90 |
| `wc_cor` | integer | state-homed-runtime | atmosphere_state.f90 |
| `wet` | real(8) | state-homed-runtime | atmosphere_state.f90 |
