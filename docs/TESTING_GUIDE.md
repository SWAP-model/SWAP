# Developer Guide: Writing Unit Tests for SWAP

This guide explains how to add new unit tests to the SWAP testing infrastructure using pFUnit.

## Table of Contents

1. [Testing Infrastructure Overview](#testing-infrastructure-overview)
2. [Prerequisites](#prerequisites)
3. [Quick Start](#quick-start)
4. [Writing Test Files](#writing-test-files)
5. [Test Patterns and Best Practices](#test-patterns-and-best-practices)
6. [Adding Tests to Build System](#adding-tests-to-build-system)
7. [Running Tests](#running-tests)
8. [Troubleshooting](#troubleshooting)
9. [Examples](#examples)

## Testing Infrastructure Overview

SWAP uses **pFUnit 4.15** (Fortran unit testing framework) with:
- **Build system**: Meson + Ninja
- **Compiler**: gfortran with `-std=legacy` for compatibility
- **Test location**: `tests/unit/`
- **Test file extension**: `.pf` (pFUnit preprocessed files)

### Current Test Coverage

- **32 tests** across 4 test suites:
  - `test_functions` - SWAP interpolation functions (5 tests)
  - `test_ttutil` - ttutil utility functions (9 tests)
  - `test_van_genuchten` - Soil hydraulic functions (4 tests)
  - `test_datetime` - Date/time utilities (14 tests)

## Prerequisites

### Environment Setup

```bash
cd /home/zawadzkim/Code/swap
# Ensure pFUnit is built with gfortran
# Location: tests/pFUnit/build/install_gfortran/PFUNIT-4.15
```

### Required Knowledge

- Basic Fortran syntax
- Understanding of the function/module you're testing
- Familiarity with unit testing concepts

## Quick Start

**Create a new test suite in 5 minutes:**

```fortran
! File: tests/unit/test_mymodule.pf
module test_mymodule
    use funit
    implicit none
contains

    @test
    subroutine test_my_function()
        real(8), external :: my_function
        real(8) :: result
        
        result = my_function(2.0d0, 3.0d0)
        
        @assertEqual(5.0d0, result, tolerance=1.0d-10)
    end subroutine

end module test_mymodule
```

Then add to build system (see [Adding Tests to Build System](#adding-tests-to-build-system)).

## Writing Test Files

### File Structure

Test files use `.pf` extension (pFUnit preprocessor format):

```fortran
! File: tests/unit/test_<name>.pf

! 1. Module declaration (MUST match filename without extension)
module test_<name>
    ! 2. Import pFUnit
    use funit
    implicit none

contains

    ! 3. Test subroutines with @test directive
    @test
    subroutine test_<feature>()
        ! Test code here
    end subroutine

    ! 4. Helper functions (optional, no @test directive)
    real(8) function helper_function(x)
        real(8), intent(in) :: x
        helper_function = x * 2.0d0
    end function

end module test_<name>
```

### Critical Naming Rule

**Module name MUST match filename** (without `.pf`):
- File: `test_mymodule.pf` → Module: `test_mymodule` ✅
- File: `test_mymodule.pf` → Module: `test_mymodule_suite` ❌

pFUnit preprocessor will fail if names don't match.

### pFUnit Assertions

Available assertion macros:

```fortran
! Equality assertions
@assertEqual(expected, actual)
@assertEqual(expected, actual, tolerance=1.0e-6)

! Boolean assertions
@assertTrue(condition)
@assertFalse(condition)

! Comparison assertions
@assertGreaterThan(value1, value2)
@assertLessThan(value1, value2)

! String assertions
@assertEqual('expected', trim(actual_string))
```

### Important: No Inline Comments on Assertions

**Legacy Fortran compatibility issue** - do NOT add comments on same line as assertions:

```fortran
! ❌ WRONG - Will cause compilation errors
@assertEqual(5.0d0, result)  ! Check result is 5

! ✅ CORRECT - Comment on separate line
! Check result is 5
@assertEqual(5.0d0, result)
```

## Test Patterns and Best Practices

### 1. Test Pure Functions First

**Easiest to test** - no global dependencies:

```fortran
@test
subroutine test_pure_calculation()
    real(8), external :: calculate_something
    real(8) :: input, output
    
    input = 10.0d0
    output = calculate_something(input)
    
    @assertEqual(20.0d0, output, tolerance=1.0d-10)
end subroutine
```

### 2. Test Edge Cases

```fortran
@test
subroutine test_boundary_conditions()
    real(8), external :: safe_divide
    
    ! Test division by very small number
    @assertEqual(1.0d6, safe_divide(1.0d0, 1.0d-6), tolerance=1.0)
    
    ! Test zero input
    @assertEqual(0.0d0, safe_divide(0.0d0, 5.0d0), tolerance=1.0d-10)
end subroutine
```

### 3. Test Monotonicity and Physical Constraints

```fortran
@test
subroutine test_water_content_monotonic()
    real(8) :: theta1, theta2, theta3
    
    ! Water content should decrease as pressure becomes more negative
    theta1 = vg_theta(-10.0d0, theta_r, theta_s, alpha, n, m)
    theta2 = vg_theta(-100.0d0, theta_r, theta_s, alpha, n, m)
    theta3 = vg_theta(-1000.0d0, theta_r, theta_s, alpha, n, m)
    
    @assertTrue(theta1 > theta2)
    @assertTrue(theta2 > theta3)
    @assertTrue(theta3 >= theta_r)  ! Physical constraint
end subroutine
```

### 4. Use Helper Functions for Complex Setup

```fortran
module test_soil_properties
    use funit
    implicit none
contains

    ! Helper function
    subroutine setup_sandy_soil(theta_r, theta_s, alpha, n, m)
        real(8), intent(out) :: theta_r, theta_s, alpha, n, m
        theta_r = 0.045d0
        theta_s = 0.43d0
        alpha = 0.145d0
        n = 2.68d0
        m = 1.0d0 - 1.0d0/n
    end subroutine

    @test
    subroutine test_sandy_soil_retention()
        real(8) :: theta_r, theta_s, alpha, n, m, theta
        
        call setup_sandy_soil(theta_r, theta_s, alpha, n, m)
        theta = vg_theta(-100.0d0, theta_r, theta_s, alpha, n, m)
        
        @assertTrue(theta > 0.1d0)
        @assertTrue(theta < 0.3d0)
    end subroutine

end module test_soil_properties
```

### 5. Group Related Tests

```fortran
@test
subroutine test_leap_year_regular()
    ! Test regular leap years
end subroutine

@test
subroutine test_leap_year_century_not_divisible_400()
    ! Test century exception
end subroutine

@test
subroutine test_leap_year_century_divisible_400()
    ! Test 400-year rule
end subroutine
```

## Adding Tests to Build System

### Step 1: Create Test File

```bash
cd tests/unit
vim test_mymodule.pf  # Create your test file
```

### Step 2: Register Test Suite

Edit `tests/unit/testSuites.inc`:

```fortran
ADD_TEST_SUITE(test_functions_suite)
ADD_TEST_SUITE(test_ttutil_suite)
ADD_TEST_SUITE(test_van_genuchten_suite)
ADD_TEST_SUITE(test_datetime_suite)
ADD_TEST_SUITE(test_mymodule_suite)  ! Add your new suite
```

**Note**: Suite name = `<module_name>_suite`

### Step 3: Add to Meson Build

Edit `tests/unit/meson.build`:

```meson
# Preprocess test files
pp_sources = pfunit_pp.process(
  'test_functions.pf',
  'test_ttutil.pf',
  'test_van_genuchten.pf',
  'test_datetime.pf',
  'test_mymodule.pf'  # Add your test file
)
```

### Step 4: Add Source Dependencies (if needed)

If your tests need SWAP source files:

```meson
# Add SWAP source files needed for testing
swap_pure_functions = files('swap_pure_functions.f90')
my_module_sources = files('../../src/path/to/mymodule.f90')  # Add if needed

# Build test executable
unit_tests = executable(
  'unit-swap-tests',
  sources: [
    swap_pure_functions,
    my_module_sources,  # Include your sources
    pp_sources,
    driver_src,
  ],
  # ... rest of configuration
)
```

**Important**: Only add sources with **no global dependencies** to avoid compilation errors.

## Running Tests

### Build and Run All Tests

```bash
cd /home/zawadzkim/Code/swap
pixi run meson compile -C builddir
pixi run meson test -C builddir
```

### Run Tests Directly (More Verbose)

```bash
cd builddir
./tests/unit/unit-swap-tests          # Run all tests
./tests/unit/unit-swap-tests -v       # Verbose output
```

### Expected Output

```
................................
Time:         0.000 seconds
  
 OK
 (32 tests)
```

### Run Specific Test Suite

```bash
# Not directly supported, but you can use test names in code
# Best practice: Run all tests (they're fast!)
```

## Troubleshooting

### Error: Module name doesn't match filename

```
Exception: pFUnit preprocessor: module name (test_foo_suite) and 
file name (test_foo) do not match (ignoring case).
```

**Solution**: Rename module to match filename:
```fortran
! File: test_foo.pf
module test_foo  ! Not test_foo_suite
```

### Error: Cannot open module file for reading

```
Fatal Error: Cannot open module file 'variables.mod' for reading
```

**Solution**: Your test depends on a module with global state. Either:
1. Extract the function to a standalone file (preferred)
2. Mock the dependencies
3. Add the module to sources (may cause cascading dependencies)

### Error: Invalid character in name

```
Error: Invalid character in name at (1)
/path/to/test.pf:42:10:
   42 |         @assertEqual(5, result)  ! This is wrong
```

**Solution**: Remove inline comment:
```fortran
! This is correct
@assertEqual(5, result)
```

### Tests Fail with Wrong Values

**Debug with verbose output**:
```bash
./builddir/tests/unit/unit-swap-tests -v
```

Check actual vs expected values and verify:
1. Correct function behavior
2. Appropriate tolerance
3. Physical reasonableness of test inputs

### Compilation Fails with Legacy Code

If testing legacy Fortran, you might need to:
1. Extract the function to clean code
2. Add to `swap_pure_functions.f90` as wrapper
3. Update function signatures for modern Fortran

## Examples

### Example 1: Testing Interpolation

```fortran
module test_interpolation
    use funit
    implicit none
contains

    @test
    subroutine test_linear_interpolation()
        real(8), external :: interpolate
        real(8) :: x1, y1, x2, y2, x, y
        
        ! Set up line from (0,0) to (10,100)
        x1 = 0.0d0
        y1 = 0.0d0
        x2 = 10.0d0
        y2 = 100.0d0
        
        ! Test midpoint
        x = 5.0d0
        y = interpolate(x1, y1, x2, y2, x)
        
        @assertEqual(50.0d0, y, tolerance=1.0d-10)
    end subroutine

end module test_interpolation
```

### Example 2: Testing with Multiple Cases

```fortran
module test_soil_types
    use funit
    implicit none
contains

    @test
    subroutine test_sandy_soil()
        ! Test parameters for sandy soil
    end subroutine
    
    @test
    subroutine test_clay_soil()
        ! Test parameters for clay soil
    end subroutine
    
    @test
    subroutine test_loamy_soil()
        ! Test parameters for loamy soil
    end subroutine

end module test_soil_types
```

### Example 3: Testing Logical Functions

```fortran
module test_date_logic
    use funit
    implicit none
contains

    @test
    subroutine test_is_leap_year()
        logical, external :: is_leap_year
        
        ! Regular leap years
        @assertTrue(is_leap_year(2000))
        @assertTrue(is_leap_year(2024))
        
        ! Non-leap years
        @assertFalse(is_leap_year(2001))
        @assertFalse(is_leap_year(1900))
    end subroutine

end module test_date_logic
```

## Advanced Topics

### Testing Functions with Array Arguments

```fortran
@test
subroutine test_array_function()
    real(8), external :: sum_array
    real(8) :: array(5)
    real(8) :: result
    
    array = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0]
    result = sum_array(array, 5)
    
    @assertEqual(15.0d0, result, tolerance=1.0d-10)
end subroutine
```

### Testing with Tolerances

For floating-point comparisons, always use tolerance:

```fortran
! Too strict - may fail due to rounding
@assertEqual(0.3333333333333333d0, result)

! Better - appropriate tolerance
@assertEqual(0.33333333d0, result, tolerance=1.0d-8)

! Physical quantities - match measurement precision
@assertEqual(expected_gwl, actual_gwl, tolerance=0.1d0)  ! cm precision
```

### Incrementally Testing Complex Modules

When testing modules with dependencies:

1. **Start simple** - test pure calculations
2. **Extract functions** - create standalone versions
3. **Add to test suite** - build up gradually
4. **Refactor module** - reduce dependencies over time

## Next Steps

After adding your tests:

1. ✅ Run tests locally: `pixi run meson test -C builddir`
2. ✅ Verify all pass
3. ✅ Commit test files to git
4. 🔄 Consider CI/CD integration
5. 📝 Update [TODOs.md](TODOs.md) with refactoring opportunities

## Resources

- **pFUnit Documentation**: Tests in `tests/unit/` serve as examples
- **Existing Tests**: Study `test_functions.pf`, `test_datetime.pf` for patterns
- **SWAP Source**: Browse `src/` for functions to test
- **TODOs**: See [TODOs.md](TODOs.md) for improvement suggestions

## Questions?

The test infrastructure is designed to be simple and maintainable. If something isn't clear:

1. Look at existing test files as examples
2. Check error messages carefully (they're usually specific)
3. Start with simple pure functions
4. Gradually tackle more complex modules

Happy testing! 🧪
