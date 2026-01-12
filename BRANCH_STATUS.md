# SWAP Library Branch Status

## Current Branch: `swaplib-simple`
✅ Clean baseline from commit 267202d ("Add BMI interface for SWAP library")

### What's Working:
- Basic BMI interface (`src/swap_bmi.f90`)
- Library compiles successfully (`libswap.so`)
- All 9 BMI functions export correctly:
  - `initialize` / `initialize_memory` / `finalize`
  - `update` / `update_day`
  - `get_current_time` / `get_grid_size`
  - `get_value_ptr_theta` / `get_value_gwl`
- Test suite passes (`test_bmi_lib.py`)

### Build Command:
```bash
pixi run bash -c "source /opt/intel/oneapi/setvars.sh --force && FC=ifx meson compile -C builddir"
```

### Test Command:
```bash
pixi run bash -c "source /opt/intel/oneapi/setvars.sh --force && python test_bmi_lib.py"
```

---

## Reference Branch: `swaplib`
Contains experimental memory-based configuration work (NOT currently usable).

See `MEMORY_WORK_README.txt` on that branch for details.

### Key Stashes:
Run `git stash list` to see stashed work:
- Interface design patterns document
- Memory-based implementation experiments

---

## Git Commands Reference:

### View stashed work:
```bash
git stash list                        # List all stashes
git stash show -p stash@{0}          # View latest stash contents
```

### Switch branches:
```bash
git checkout swaplib-simple          # Use this for new work
git checkout swaplib                 # Reference branch with experiments
```

### Retrieve stashed work (if needed later):
```bash
git stash apply stash@{0}            # Apply without removing from stash
git stash pop stash@{0}              # Apply and remove from stash
```
