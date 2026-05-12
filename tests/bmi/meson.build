# ============================================================================
# BMI hello-world integration test
# ============================================================================
# Runs the Python cffi driver against libswap_bmi.so in the
# hupselbrook case directory. Registered in the `bmi` suite so it
# can be run independently of pFUnit and the regression suite.
python3 = import('python').find_installation('python3')

test('bmi-hello-world',
    python3,
    args: [
        meson.current_source_dir() / 'hello_swap.py',
        swap_bmi_lib.full_path(),
    ],
    workdir: meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
    suite: 'bmi',
    depends: [swap_bmi_lib],
    timeout: 60,
)
