# Refactoring Pitfalls & FAQ

This document captures common mistakes and their solutions encountered during the Phase 2 state integration refactoring.

---

## Pitfall #1: Bidirectional Sync at Entry/Exit Causes Infinite Loops

### Symptom
- Simulation appears to run but time never advances (e.g., stuck at `day=1, year=1997`)
- No output files are produced
- Log shows the same day/year repeating indefinitely

### Root Cause
When integrating a state type (e.g., `time_state_t`) into a subroutine using the sync bridge pattern, we initially added **bidirectional sync**:

```fortran
subroutine TimeControl(task, tstate)
    ! WRONG: Sync state → variables at entry
    if (present(tstate)) call time_state_to_variables(tstate)
    
    ! ... subroutine updates daynr, iyear, etc. ...
    
    ! Sync variables → state at exit
    if (present(tstate)) call time_state_from_variables(tstate)
end subroutine
```

**Problem**: The entry sync (`time_state_to_variables`) overwrites the module variables with old state values before the subroutine can update them. On the next call:
1. `daynr` was updated from 1 → 2 during the previous call
2. Exit sync captured `daynr=2` into state
3. But on the NEXT call, entry sync copies the OLD state value back
4. Since we're syncing BEFORE the update, `daynr` is reset to 1
5. Time never advances

### Solution
The sync bridge pattern should be **ONE-WAY at exit only**:

```fortran
subroutine TimeControl(task, tstate)
    ! NO sync at entry - module variables retain their current values
    
    ! ... subroutine updates daynr, iyear, etc. ...
    
    ! ONLY sync variables → state at exit (captures updated values)
    if (present(tstate)) call time_state_from_variables(tstate)
end subroutine
```

### When to Use Entry Sync
Entry sync (`state_to_variables`) is ONLY needed when:
1. The state is the **source of truth** and the subroutine needs to read from it
2. This is typically in **multi-instance scenarios** where each instance has its own state

During the transition phase (sync bridge), the module variables are still the source of truth, so:
- ✅ Sync variables → state at exit (capture updates)
- ❌ Do NOT sync state → variables at entry (would overwrite updates)

---

## Pitfall #2: State Initialization Order Matters

### Symptom
- Error: `"ERROR in DTDPAR: date before 1-jan-0001"`
- State contains zeros instead of initialized values

### Root Cause
`swap_state_init` was called BEFORE `TimeControl(1)`, so when `time_state_from_variables` was called, the module variables hadn't been initialized yet.

### Solution
Initialize time state AFTER the subroutine that populates the module variables:

```fortran
! 1. TimeControl(1) initializes module variables (daynr, iyear, etc.)
call TimeControl(1)

! 2. Now we can populate state from initialized module variables
call time_state_from_variables(state%time)

! 3. CalcGrid determines numnod/numlay
call CalcGrid

! 4. Initialize array-based states (soil, drainage, etc.)
call swap_state_init(state, numnod, numlay)
```

**Key insight**: `time_state_t` only contains scalar values with defaults (no allocation needed), so we can populate it before `swap_state_init`. The array-based states need `numnod`/`numlay` which aren't known until after `CalcGrid`.

---

## Pitfall #3: Regression Tests May Pass Despite Model Failure

### Symptom
- All regression tests report "PASS"
- But manual inspection shows no output or wrong output

### Root Cause
The CSV comparison script may pass if:
1. The reference file and result file are both empty
2. The comparison tolerance is too loose
3. The script doesn't check if files exist or have content

### Solution
- Always manually verify at least one test case after major changes
- Check that output files have expected content, not just that comparison passes
- Consider adding assertions for minimum file size or row count

---

## General Guidelines for State Integration

### The Correct Pattern

1. **Add optional state argument** to subroutine:
   ```fortran
   subroutine ProcessName(task, pstate)
       type(process_state_t), intent(inout), optional :: pstate
   ```

2. **Sync at exit only** (at ALL return points):
   ```fortran
   ! Before each RETURN statement:
   if (present(pstate)) call process_state_from_variables(pstate)
   return
   ```

3. **Populate state after initialization**:
   ```fortran
   call ProcessName(1)  ! Initialize module variables
   call process_state_from_variables(state%process)  ! Capture into state
   ```

4. **Test after each integration**:
   - Build: `pixi run build-linux`
   - Unit tests: `pixi run test-unit`
   - Regression tests: `pixi run regression_all`
   - Manual verification of at least one test case output

### Order of Integration
1. Time state (scalar, no arrays)
2. Soil state (needs numnod, numlay)
3. Drainage state (needs numnod, nrlevs)
4. Other states...
5. Crop state (most complex, do last)
