# ADR 0046 — `state%cfg` pointer retirement

**Status:** Accepted (2026-05-28)
**Closes:** the broader `state%cfg` reach-through anti-pattern flagged during the 2026-05-28 model-setup-layer analysis. Direct follow-on to ADR 0045.

## Context

`state%cfg => config` was a back-pointer set in `swap_init_body` so any compute path could read configuration via `state%cfg%X%Y` without threading `config` through signatures. It enabled fast strangler-pattern migration but left an architectural debt: `state%cfg` reads were the equivalent of bare global reads at the config layer — implicit dependencies invisible at the call site.

At the start of this arc, 159 live `state%cfg` reads existed across 25 non-dormant files, the heaviest concentrations in `src/soil/soilhydraulics.f90` (40), `src/io/csv_output.f90` (12), `src/state/legacy_state.f90` (16), `src/core/timecontrol_mod.f90` (11), and various crop runtime modules (~30).

## Decision

Retire the `state%cfg` pointer. Migrate every live read using one of three patterns:

1. **Direct config-parameter substitution:** when the consumer already takes `config` as a parameter (or trivially can), replace `state%cfg%X%Y` with `config%X%Y`.
2. **Snapshot to state at init:** when the consumer takes only `state`, snapshot the needed config values onto a sensible state sub-record at init time (e.g. `state%timecontrol%swssdi/swrain/swmetdetail`, `state%atmosphere%lat`, `state%soilwater%swkmean`, `state%crop%rotation_wofost`).
3. **Signature expansion:** when the consumer is init-only and the config sub-tree is small, expand its signature to take the specific sub-record directly.

After all reads migrated, delete the `cfg :: type(swap_config_t), pointer` field from `swap_state_t` and the corresponding assignment in `swap_init_body`.

## Cluster summary

Each cluster shipped as one commit:

1. **swap_mod + atmosphere** (~20 reads): substitution. Config already in scope.
2. **timecontrol_mod** (~11 reads): signature expansion (`timecontrol_init` gains `config` parameter; `target` attribute required to prevent gfortran heap corruption); state%timecontrol snapshots for runtime reads.
3. **drainage cluster** (~9 reads): snapshot to state%drainage / state%surfacewater / state%soilwater%swkmean.
4. **csv_output** (~12 reads): snapshot to state%timecontrol (pathwork, outfil, project, csv_enabled, csv_inlist, csv_tz_z1_z2).
5. **crop subsystem** (~23 reads): snapshot rotation arrays to state%crop; substitution where possible.
6. **soilhydraulics + waterbalance** (~41 reads): snapshot of numerical, bottom-boundary, and soil-hydraulics sub-records to state%soilwater; atmosphere%swsnow for waterbalance.
7. **Mop-up of deferred items** (8 reads): per-rotation crop arrays, populate_hydraulic_params signature expansion, meteodt raintab snapshot.

Final commit (this one): delete the field, the assignment, and the stub references.

## Consequences

+ **No more bare-global-equivalent reads at the config layer.** Every config dependency is visible at the call site (as a parameter) or owned by an init that documents it.
+ **State subsystems own their derived data.** A maintainer reading a compute routine knows the data lives on state, populated by the corresponding init.
+ **Multi-instance safety.** With no global pointer, the model is closer to supporting multiple concurrent `swap_state_t` instances (the broader goal of the rescue effort).
+ **Init signatures grew.** `state%timecontrol%init`, `state%soilwater%init`, `state%crop%init` now take more config sub-trees or the full config object. Trade explicitness vs. brevity — explicitness wins.

- **State records grew.** ~30 new snapshot fields across the state hierarchy.
- **Some init-time work moved.** The `populate_hydraulic_params` init helper signature expanded; the `timecontrol_init` body folded into the type-bound init.

## Verification

- 159 → 0 live `state%cfg` reads in src/.
- check-fast PASS: 4/4 byte-identical.
- check-full PASS: 5/5 non-xfail byte-identical.
- 9 commits on branch `state-cfg-retirement`.

## What's next

The state-init layer is now `state%cfg`-free, every config dependency explicit. The remaining architectural items from the 2026-05-28 model-setup-layer analysis (W5-W12) are smaller cleanup items:
- W5: init-order dependencies encoded only in comments — could be made explicit via a docstring or compile-time check.
- W7/W8: CAPI path resolution (hardcoded `base_dir = './'` and two path-resolution conventions).
- W9: TOML reader code duplication (`read_table_2d` × 5).
- Other minor items.

The major architectural arc that started with the meteo typed-CSV pilot is complete: typed CSV records on state (ADR 0044), orchestrator dissolved (ADR 0045), config-as-pointer retired (this ADR 0046).
