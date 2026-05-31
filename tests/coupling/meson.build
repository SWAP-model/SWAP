python3 = import('python').find_installation('python3')

test('swap-xmi-smoke',
    python3,
    args: [
        meson.current_source_dir() / 'test_swap_xmi_smoke.py',
        swap_xmi_lib.full_path(),
        meson.current_source_dir() / 'hupselbrook_coupled',
        meson.current_build_dir(),
    ],
    suite: 'coupling',
    depends: [swap_xmi_lib],
    timeout: 600,
)

# Coupled SWAP<->MODFLOW6 qualitative smoke (short run: ncol=10, nper=30).
# libmf6.so lives in the pixi test environment; the path is discovered at
# configure time from the project root. If it is absent (e.g. building without
# the test env), the test is skipped rather than failing the suite. The smoke
# definition of done is qualitative -- see tests/coupling/README.md.
fs = import('fs')
libmf6 = meson.project_source_root() / '.pixi/envs/test/lib/libmf6.so'
if fs.exists(libmf6)
    test('coupled-smoke',
        python3,
        args: [
            meson.current_source_dir() / 'test_coupled_smoke.py',
            swap_xmi_lib.full_path(),
            libmf6,
            meson.current_source_dir() / 'hupselbrook_coupled',
        ],
        suite: 'coupling',
        depends: [swap_xmi_lib],
        timeout: 1200,
    )
else
    message('coupled-smoke: libmf6.so not found at ' + libmf6 + ' -- test not registered (run tests/coupling/test_coupled_smoke.py manually).')
endif
