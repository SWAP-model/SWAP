# src/driver

## Responsibility

Top-level lifecycle orchestration and the standalone executable. The apex of
the dependency graph — depends on every compute subsystem.

- **`swap_mod.f90`** (`swap_mod`) — the three named lifecycle procedures
  `swap_init` / `swap_run_step` / `swap_close` that replace the legacy
  `(iCaller, iTask)` dispatch. The time loop lives in the caller; state and
  config are threaded explicitly.
- **`swap_main.f90`** — the `swap` executable entry point.
- **`swap_ensemble_mod.f90`** — N-instance ensemble driver used by the XMI
  coupling facade.

## Dependencies

Depends on `core/`, `config/`, `state/`, `io/`, `timecontrol/`, and every
compute feature folder. Wrapped by `bindings/` for the C-ABI.
