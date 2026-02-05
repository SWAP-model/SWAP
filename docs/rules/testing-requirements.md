# Testing Requirements

## Mandatory Testing for Every Change

Every code modification MUST pass ALL of the following:

### 1. Build Verification
```bash
pixi run clean
pixi run build-linux
# etc...
```
Success criteria:

    No compilation errors

    Zero new compiler warnings

    All modules link successfully

### 3. Integration Tests

```bash
pixi run test-linux-hupselbrook
pixi run test-linux-grassgrowth
pixi run test-linux-surfacewater
pixi run test-linux-oxygenstress
pixi run test-linux-macropore
```

Success criteria:

    Simulation completes without errors

### 4. Regression Validation

```bash
pixi run regression_all
```

Success criteria:

    All simulations complete successfully
    Only the oxygenstress can have higher discrepancies which are expected.
    The rest of the tests should end successfully 

### 5. Unit tests (pHUnit)

```bash
pixi run test-unit-all
```
