# Crop Subsystem Symbol Ledger — Pre-Flight Discovery (2026-05-25)

> **Status**: pre-flight discovery completed  
> **Commit**: [awaiting ledger commit]  
> **Purpose**: inventory of all symbols imported via `use variables` across 12 crop files for retirement sweep orchestration.

## Section 1: Per-File Summary

| File | use-blocks | Total imports | Estimated class distribution |
|---|---|---|---|
| cropfixed_init            | 1 |   1 | B(1) |
| cropfixed_runtime         | 1 |  11 | B(1), C(7), D(3) |
| cropgrass_init            | 1 |   6 | B(2), D(4) |
| cropgrass_runtime         | 1 |  14 | A(1), C(8), D(5) |
| cropgrowth                | 1 |  58 | A(17), C(21), D(20) |
| cropgrowth_helpers        | 1 |  46 | A(6), C(8), D(32) |
| cropwofost_init           | 1 |  33 | A(2), C(23), D(8) |
| cropwofost_runtime        | 1 |  62 | A(6), B(1), C(43), D(12) |
| irrigation                | 1 |  46 | A(3), B(2), C(24), E(17) |
| oxygenstress              | 1 |  58 | A(2), B(7), C(49) |
| rootextraction            | 1 |  15 | A(1), D(14) |
| tillage                   | 1 |  39 | B(4), C(25), D(10) |

**Totals**: 389 imports across 12 files

## Section 2: Full Symbol Ledger

| Symbol | Imported by | Count | Class | Target home | Notes |
|---|---|---|---|---|---|
| Date_tillage              | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| I_tillage                 | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| MaxPrepDelay              | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| MaxSowDelay               | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| Max_Z_tillage             | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| Ntill                     | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| Ntypes                    | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| ParamVG                   | tillage                                  | 1 | B | state%cfg% (direct read)            | cross=N |
| PrepDelay                 | cropgrowth_helpers, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| SowDelay                  | cropgrowth_helpers, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| SwDiscrvert               | tillage                                  | 1 | B | state%cfg% (direct read)            | cross=N |
| TAB_K_R_cons              | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| TAB_N_match               | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| TAB_Rho_cons              | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| TAB_Rho_match             | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| TAB_Rho_tillage           | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| TBASEM                    | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| TEFFMX                    | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| TempSow                   | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| Type_Tillage              | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| Z_tillage                 | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| agerm                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| anlv                      | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| anst                      | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| bgerm                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| bunsencoeff               | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen | cross=N |
| c_macro                   | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen | cross=N |
| c_min_micro               | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen | cross=N |
| c_mroot                   | oxygenstress, cropfixed_init/runtime     | 3 | C | state%crop%oxygen (seeded from legacy at OxygenStress entry; cross-file writer retires after Task 6) | cross=Y |
| c_top                     | oxygenstress + swap_csv_output           | 2 | E | retained — CSV output buffer (out of scope) | cross=Y |
| cgerm                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| cirrs                     | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| cirrthres                 | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| co2amaxtb                 | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| co2efftb                  | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| co2ppm                    | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| co2tratb                  | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| co2year                   | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| criterhr                  | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| cropend                   | cropwofost_runtime, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| cropfil                   | cropgrowth_helpers, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| cropstart                 | cropgrowth                               | 1 | A | state%crop% (exists)                | cross=N |
| crp                       | cropgrowth_helpers                       | 1 | C | state%crop% (new field)             | cross=N |
| ctopnode                  | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| d_o2inwater               | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| d_root                    | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| d_soil                    | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| daycrop                   | oxygenstress, cropgrass_runtime, cropwofost_init, +2 | 5 | A | state%crop% (exists)                | cross=Y |
| dayfix                    | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| days_counter_irr          | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| days_interval_irr         | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| depth                     | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| dhPrep                    | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| dhSow                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| ditab                     | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| dlc                       | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| dlo                       | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| drbl                      | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| drblpot                   | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| dt_SSDI_event             | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| dtempSow                  | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| dummy_tsoil_alg_          | cropgrowth_helpers                       | 1 | C | state%crop% (new field)             | cross=N |
| dummy_tsoil_cg_           | cropgrowth                               | 1 | C | state%crop% (new field)             | cross=N |
| dummy_tsoil_gr_           | cropgrass_runtime                        | 1 | C | state%crop% (new field)             | cross=N |
| dummy_tsoil_sumttd_       | cropgrowth_helpers                       | 1 | C | state%crop% (new field)             | cross=N |
| dvsnlt                    | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| dvsnt                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| dwatab                    | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| f_senes                   | oxygenstress                             | 1 | B | state%cfg% (direct read)            | cross=N |
| fbl                       | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| fidtab                    | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| flCropCalendar            | irrigation, cropgrowth                   | 2 | A | state%crop% (exists)                | cross=Y |
| flCropEmergence           | cropgrowth                               | 1 | A | state%crop% (exists)                | cross=N |
| flCropGerm                | cropgrowth_helpers, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| flCropHarvest             | irrigation, cropwofost_runtime, cropgrowth | 3 | A | state%crop% (exists)                | cross=Y |
| flCropNut                 | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | A | state%crop% (exists)                | cross=Y |
| flCropOpenFile            | cropgrowth_helpers                       | 1 | A | state%crop% (exists)                | cross=N |
| flCropPrep                | cropgrowth_helpers, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| flCropReadFile            | cropgrowth                               | 1 | A | state%crop% (exists)                | cross=N |
| flCropSow                 | cropgrowth_helpers, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| flHarvestDay              | cropwofost_runtime, cropgrowth           | 2 | A | state%crop% (exists)                | cross=Y |
| flanthesis                | cropwofost_runtime                       | 1 | A | state%crop% (exists)                | cross=N |
| flhydrlift                | rootextraction, cropfixed_runtime, cropgrass_runtime, +1 | 4 | D | state%cfg% (new config)             | cross=Y |
| fntrt                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| fraharlosorm_lv           | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| fraharlosorm_so           | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| fraharlosorm_st           | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| frnx                      | cropwofost_init, cropgrowth              | 2 | C | state%crop% (new field)             | cross=Y |
| fstr                      | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| gas_filled_porosity       | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| gasst                     | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| gasstpot                  | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| gwrt                      | cropgrowth_helpers, cropgrass_runtime, cropwofost_runtime | 3 | D | state%cfg% (new config)             | cross=Y |
| hPrep                     | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| hSow                      | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| hcritab                   | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| hdrygerm                  | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| hprep                     | cropgrowth                               | 1 | C | state%crop% (new field)             | cross=N |
| hwetgerm                  | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| iRedist                   | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| iTT1                      | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| iTT2                      | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| iType_Tillage             | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| i_n_model                 | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| icrop                     | oxygenstress, cropgrowth                 | 2 | A | state%crop% (exists)                | cross=Y |
| idsl                      | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| ilnmxl                    | cropwofost_init                          | 1 | C | state%crop% (new field)             | cross=N |
| irconc                    | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| irdate                    | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| irdepth                   | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| irrigevent                | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| irtype                    | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| isua                      | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| isuas                     | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| kroot                     | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| kstem                     | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| lrnr                      | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | C | state%crop% (new field)             | cross=Y |
| lsnr                      | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | C | state%crop% (new field)             | cross=Y |
| macp                      | cropgrass_runtime, cropwofost_runtime    | 2 | C | state%crop% (new field)             | cross=Y |
| magrs                     | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime | 3 | C | state%crop% (new field)             | cross=Y |
| max_resp_factor           | oxygenstress, cropfixed_init, cropfixed_runtime | 3 | B | state%cfg% (direct read)            | cross=Y |
| mayrs                     | cropgrowth_helpers                       | 1 | D | state%cfg% (new config)             | cross=N |
| mrest                     | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| mrestpot                  | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| mrftb                     | cropfixed_runtime                        | 1 | C | state%crop% (new field)             | cross=N |
| nfixf                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| nirri                     | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| nirri_ssdi_irr            | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| nlai                      | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| nlue                      | cropwofost_init, cropgrowth              | 2 | C | state%crop% (new field)             | cross=Y |
| nmaxlv                    | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| nmaxrt                    | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| nmaxso                    | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| nmaxst                    | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| nmxlv                     | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | C | state%crop% (new field)             | cross=Y |
| nni                       | cropwofost_runtime, cropgrowth           | 2 | C | state%crop% (new field)             | cross=Y |
| nod_ssdi_irr              | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| nod_ssdi_sensor_irr       | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| noddrz                    | RETIRED                                  | 0 | A | state%crop%common%noddrz            | RETIRED 2026-05-25 (Task 3) |
| nofd                      | cropwofost_init                          | 1 | C | state%crop% (new field)             | cross=N |
| npart                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| nsla                      | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| o2_bunsencoeff            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_c_macro                | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_c_min_micro            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_capac_term             | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_ctopnode               | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_d_o2inwater            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_d_root                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_d_soil                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_d_soil_term1           | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_d_soil_term2           | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_depth                  | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_gas_filled_porosity    | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_gfp100                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_ini_stress             | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_mplus1                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_nmin1                  | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_perc_org_mat           | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_r_microbial_z0         | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_root_radius            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_sat_water_cont         | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_shape_factor_microbialr | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_soil_density           | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_soil_temp              | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_w_root                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_w_root_z0              | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| o2_waterfilm_thickness    | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| outfil                    | cropgrowth_helpers, cropwofost_runtime   | 2 | C | state%crop% (new field)             | cross=Y |
| oxygenintercept           | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| oxygenslope               | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| pathcrop                  | cropgrowth                               | 1 | C | state%crop% (new field)             | cross=N |
| pathwork                  | cropgrowth_helpers, cropwofost_runtime   | 2 | C | state%crop% (new field)             | cross=Y |
| perc_org_mat              | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| perirrsurp                | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| pld                       | cropgrowth                               | 1 | C | state%crop% (new field)             | cross=N |
| project                   | cropgrowth_helpers, cropwofost_runtime   | 2 | C | state%crop% (new field)             | cross=Y |
| q10_microbial             | oxygenstress, cropgrass_init             | 2 | B | state%cfg% (direct read)            | cross=Y |
| q10_root                  | oxygenstress                             | 1 | B | state%cfg% (direct read)            | cross=N |
| r_microbial_z0            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| raithreshold              | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| rawtab                    | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| rdmax                     | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime | 3 | C | state%crop% (new field)             | cross=Y |
| rdrns                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| reltr                     | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime | 3 | C | state%crop% (new field)             | cross=Y |
| remoc                     | cropgrowth                               | 1 | C | state%crop% (new field)             | cross=N |
| rid                       | oxygenstress, cropgrass_runtime          | 2 | C | state%crop% (new field)             | cross=Y |
| rnflv                     | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | C | state%crop% (new field)             | cross=Y |
| rnfrt                     | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| rnfst                     | cropwofost_init, cropwofost_runtime, cropgrowth | 3 | C | state%crop% (new field)             | cross=Y |
| root_radius               | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| rootcoefa                 | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| rooteff                   | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| rootradius                | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| rxylem                    | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| sat_water_cont            | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| shape_factor_microbialr   | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| shape_factor_rootr        | oxygenstress                             | 1 | C | state%crop% (new field)             | cross=N |
| siccaplai                 | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime | 3 | C | state%crop% (new field)             | cross=Y |
| soil_density              | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| soil_temp                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| specific_resp_humus       | oxygenstress, cropgrass_init             | 2 | B | state%cfg% (direct read)            | cross=Y |
| ssdi_amount_f_irr         | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_amount_irr           | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_appl_rate_irr        | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_date_irr             | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_rate_f_irr           | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_sched_type_irr       | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_schedule_irr         | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_threshold_irr        | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| ssdi_threshold_z_irr      | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| stephr                    | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| sw_interval_irr           | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| swbulb                    | cropwofost_runtime                       | 1 | B | state%cfg% (direct read)            | cross=N |
| swcirrthres               | irrigation                               | 1 | C | state%crop% (new field)             | cross=N |
| swcrp                     | cropgrowth                               | 1 | A | state%crop% (exists)                | cross=N |
| swfrost                   | rootextraction                           | 1 | D | state%cfg% (new config)             | cross=N |
| swirfix                   | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| swpotrelmf                | cropgrass_init, cropwofost_init          | 2 | D | state%cfg% (new config)             | cross=Y |
| swsolu                    | tillage, irrigation                      | 2 | B | state%cfg% (direct read)            | cross=Y |
| swssdi_irr                | irrigation                               | 1 | B | state%cfg% (direct read)            | cross=N |
| swtill                    | tillage                                  | 1 | B | state%cfg% (direct read)            | cross=N |
| tadw                      | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| tadwpot                   | cropwofost_runtime                       | 1 | C | state%crop% (new field)             | cross=N |
| tawtab                    | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| tcnt                      | cropwofost_init, cropwofost_runtime      | 2 | C | state%crop% (new field)             | cross=Y |
| tcritab                   | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| tendirrig                 | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| till_Date_tillage         | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_I_tillage            | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_Max_Z_tillage        | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_Ntill                | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_Ntypes               | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_TAB_K_R_cons         | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| till_TAB_N_match          | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| till_TAB_Rho_cons         | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| till_TAB_Rho_match        | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| till_TAB_Rho_tillage      | tillage                                  | 1 | D | state%cfg% (new config)             | cross=N |
| till_Type_Tillage         | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_Z_tillage            | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_iRedist              | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_iTT1                 | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_iTT2                 | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_iType_Tillage        | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_i_n_model            | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| till_swtill               | tillage                                  | 1 | C | state%crop% (new field)             | cross=N |
| treltab                   | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| tsoil                     | oxygenstress, cropgrowth_helpers, cropgrass_runtime, +1 | 4 | C | state%crop% (new field)             | cross=Y |
| tstairrig                 | irrigation                               | 1 | E | retire (orphan)                     | cross=N |
| tsumdepth                 | cropgrowth_helpers, cropgrass_init       | 2 | D | state%cfg% (new config)             | cross=Y |
| tsumemeopt                | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| tsumgerm                  | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| tsumtemp                  | cropgrowth_helpers, cropgrass_init       | 2 | D | state%cfg% (new config)             | cross=Y |
| tsumtime                  | cropgrowth_helpers, cropgrass_init       | 2 | D | state%cfg% (new config)             | cross=Y |
| twilt                     | rootextraction, cropfixed_runtime, cropgrass_runtime, +1 | 4 | D | state%cfg% (new config)             | cross=Y |
| vernbase                  | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| verndvs                   | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| vernrtb                   | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| vernsat                   | cropwofost_init, cropwofost_runtime      | 2 | D | state%cfg% (new config)             | cross=Y |
| w_root                    | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| w_root_ss                 | oxygenstress, cropfixed_runtime          | 2 | C | state%crop% (new field)             | cross=Y |
| w_root_z0                 | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| waterfilm_thickness       | oxygenstress                             | 1 | C | RETIRED 2026-05-25 — state%crop%oxygen      | cross=N |
| wiltpoint                 | rootextraction, cropfixed_runtime, cropgrass_runtime, +1 | 4 | D | state%cfg% (new config)             | cross=Y |
| wrtb                      | cropfixed_runtime                        | 1 | C | state%crop% (new field)             | cross=N |
| wrtmin                    | cropgrowth_helpers, cropgrass_runtime, cropwofost_runtime | 3 | D | state%cfg% (new config)             | cross=Y |
| zPrep                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| zSow                      | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| zTempSow                  | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |
| zgerm                     | cropgrowth_helpers, cropgrowth           | 2 | D | state%cfg% (new config)             | cross=Y |

## Section 3: Cross-File Symbols (2+ files)

| Symbol | Last Consumer (proposal) | Current consumers |
|---|---|---|
| MaxPrepDelay              | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| MaxSowDelay               | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| PrepDelay                 | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| SowDelay                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| TBASEM                    | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| TEFFMX                    | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| TempSow                   | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| agerm                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| anlv                      | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| anst                      | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| bgerm                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| cgerm                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| cropend                   | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| cropfil                   | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| daycrop                   | oxygenstress                   | cropgrass_runtime, cropgrowth, cropwofost_init, cropwofost_runtime, oxygenstress |
| dhPrep                    | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| dhSow                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| dlc                       | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| dlo                       | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| dtempSow                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| dvsnlt                    | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| dvsnt                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| flCropCalendar            | irrigation                     | cropgrowth, irrigation |
| flCropGerm                | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| flCropHarvest             | irrigation                     | cropgrowth, cropwofost_runtime, irrigation |
| flCropNut                 | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| flCropPrep                | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| flCropSow                 | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| flHarvestDay              | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| flhydrlift                | rootextraction                 | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime, rootextraction |
| fntrt                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| fraharlosorm_lv           | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| fraharlosorm_so           | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| fraharlosorm_st           | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| frnx                      | cropwofost_init                | cropgrowth, cropwofost_init |
| fstr                      | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| gwrt                      | cropwofost_runtime             | cropgrass_runtime, cropgrowth_helpers, cropwofost_runtime |
| hSow                      | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| hdrygerm                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| hwetgerm                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| icrop                     | oxygenstress                   | cropgrowth, oxygenstress |
| idsl                      | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| lrnr                      | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| lsnr                      | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| macp                      | cropwofost_runtime             | cropgrass_runtime, cropwofost_runtime |
| magrs                     | cropwofost_runtime             | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime |
| max_resp_factor           | oxygenstress                   | cropfixed_init, cropfixed_runtime, oxygenstress |
| nfixf                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| nlai                      | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| nlue                      | cropwofost_init                | cropgrowth, cropwofost_init |
| nmaxlv                    | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| nmaxrt                    | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| nmaxso                    | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| nmaxst                    | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| nmxlv                     | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| nni                       | cropwofost_runtime             | cropgrowth, cropwofost_runtime |
| noddrz                    | RETIRED (Task 3 2026-05-25)    | declaration removed; all reads/writes via state%crop%common%noddrz |
| npart                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| nsla                      | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| outfil                    | cropwofost_runtime             | cropgrowth_helpers, cropwofost_runtime |
| pathwork                  | cropwofost_runtime             | cropgrowth_helpers, cropwofost_runtime |
| project                   | cropwofost_runtime             | cropgrowth_helpers, cropwofost_runtime |
| q10_microbial             | oxygenstress                   | cropgrass_init, oxygenstress |
| rdmax                     | cropwofost_runtime             | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime |
| rdrns                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| reltr                     | cropwofost_runtime             | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime |
| rid                       | oxygenstress                   | cropgrass_runtime, oxygenstress |
| rnflv                     | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| rnfrt                     | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| rnfst                     | cropwofost_runtime             | cropgrowth, cropwofost_init, cropwofost_runtime |
| siccaplai                 | cropwofost_runtime             | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime |
| specific_resp_humus       | oxygenstress                   | cropgrass_init, oxygenstress |
| swpotrelmf                | cropwofost_init                | cropgrass_init, cropwofost_init |
| swsolu                    | tillage                        | irrigation, tillage |
| tcnt                      | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| tsoil                     | oxygenstress                   | cropgrass_runtime, cropgrowth, cropgrowth_helpers, oxygenstress |
| tsumdepth                 | cropgrowth_helpers             | cropgrass_init, cropgrowth_helpers |
| tsumemeopt                | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| tsumgerm                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| tsumtemp                  | cropgrowth_helpers             | cropgrass_init, cropgrowth_helpers |
| tsumtime                  | cropgrowth_helpers             | cropgrass_init, cropgrowth_helpers |
| twilt                     | rootextraction                 | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime, rootextraction |
| vernbase                  | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| verndvs                   | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| vernrtb                   | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| vernsat                   | cropwofost_runtime             | cropwofost_init, cropwofost_runtime |
| w_root_ss                 | oxygenstress                   | cropfixed_runtime, oxygenstress |
| wiltpoint                 | rootextraction                 | cropfixed_runtime, cropgrass_runtime, cropwofost_runtime, rootextraction |
| wrtmin                    | cropwofost_runtime             | cropgrass_runtime, cropgrowth_helpers, cropwofost_runtime |
| zPrep                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| zSow                      | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| zTempSow                  | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |
| zgerm                     | cropgrowth_helpers             | cropgrowth, cropgrowth_helpers |

## Section 4: Summary Statistics

- **Total unique symbols**: 270
- **Cross-file symbols**: 93
- **Single-file symbols**: 177

**Classification breakdown**:

| Class | A (state-rebind) | B (config-direct) | C (state-migrate) | D (config-migrate) | E (retired-zero) | F (dormant) |
|---|---|---|---|---|---|---|
| Count | 19 | 13 | 157 | 64 | 17 | 0 |

- **Cheap retirements** (A+B, state already exists): 32
- **Expensive retirements** (C+D, new state/config needed): 221
- **Orphan/dormant** (E+F, no home): 17

**Next steps**:
1. Verify classification against actual state/config files
2. Confirm `noddrz` last-consumer is rootextraction.f90
3. Confirm `tsoil` is dormant config-staging buffer (retire renaming only)
4. Dispatch subagent sub-arcs bottom-up (tillage → irrigation → ... → cropgrowth)

