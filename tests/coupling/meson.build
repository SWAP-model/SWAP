python3 = import('python').find_installation('python3')

test('swap-xmi-smoke',
    python3,
    args: [
        meson.current_source_dir() / 'test_swap_xmi_smoke.py',
        swap_xmi_lib.full_path(),
        meson.project_source_root() / 'tests/swap-cases/toml/1.hupselbrook',
        meson.current_build_dir(),
    ],
    suite: 'coupling',
    depends: [swap_xmi_lib],
    timeout: 600,
)
