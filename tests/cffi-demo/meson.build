# ============================================================================
# SS-BMI2 Phase 2 end-to-end demo
# ============================================================================
# Loads SWAP TOML in Python, passes via cffi, runs headless, reads
# state via capi accessors. Proves the in-memory workflow.

python3 = import('python').find_installation('python3')

test('cffi-demo-hupselbrook',
    python3,
    args: [
        meson.current_source_dir() / 'run_ensemble.py',
        meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
        swap_lib.full_path(),
    ],
    workdir: meson.project_source_root(),
    suite: 'cffi-demo',
    depends: [swap_lib],
    timeout: 120,
)
