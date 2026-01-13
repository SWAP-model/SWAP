project('swap', 'fortran', version : '4.2.0', default_options : ['warning_level=0'])

fc = meson.get_compiler('fortran')
fs = import('fs')
src_inc = include_directories(
    'src',
    'src/core',
    'src/io',
    'src/soil',
    'src/crop',
    'src/atmos',
    'src/solute'
)

# Intel Fortran compiler settings
if fc.get_id() == 'intel'
    add_project_arguments([
        '-O2',
        '-fpscomp', 'general',
        '-warn', 'declarations',
        '-warn', 'unused',
        '-warn', 'uncalled', 
        '-warn', 'interfaces',
        '-init=zero',
        '-fpe0',
        '-fp-model', 'source'
    ], language : 'fortran')
    
    # NOTE: Windows Intel Fortran (ifx) static linking remains problematic.
    # Despite multiple attempts (/MT, -static-intel, /NODEFAULTLIB, etc.),
    # the executable still requires Intel runtime DLLs (libmmd.dll).
    # This appears to be a known limitation of ifx on Windows.
    # Use MinGW cross-compilation for fully static Windows executables.
    
    # Commented out Windows Intel build - uncomment to re-enable testing
    # if host_machine.system() == 'windows'
    #     add_project_arguments([
    #         '/MT',                           # Static multithreaded C runtime
    #     ], language : 'fortran')
    #     add_project_link_arguments([
    #         '-static-intel'                  # Attempt static Intel libraries
    #     ], language : 'fortran')
    # else
    if host_machine.system() != 'windows'
        # Linux: Statically link Intel runtime libraries
        # System libraries (libc, libm) remain dynamically linked
        add_project_link_arguments([
            '-static-intel',    # Statically link Intel-specific libraries
        ], language : 'fortran')
    endif
    
# MinGW Fortran compiler settings for Windows cross-compilation 
elif fc.get_id() == 'gcc' and host_machine.system() == 'windows'
    add_project_arguments([
        '-O2',
        '-ffree-line-length-none',  # Remove line length limit
        '-Wno-line-truncation',
        '-Wno-compare-reals', 
        '-Wno-missing-include-dirs',  # Ignore missing include directories
        '-Wno-unused-variable',       # Ignore unused variable warnings  
        '-Wno-unused-dummy-argument', # Ignore unused dummy argument warnings
        '-std=legacy',
        '-mconsole',
        '-static-libgfortran',
        '-static-libgcc'
    ], language : 'fortran')
    add_project_link_arguments(['-static', '-mconsole'], language : 'fortran')
endif

# Source files grouped by subfolder (preserve original order where possible)
sources = [
    'src/core/variables.f90',
    'src/core/swap_state.f90',              # NEW: State object module (after variables)
    'src/io/swap_csv_output.f90', 
    'src/crop/sptabulated.f90',
    'src/crop/wofost_soil_declarations.f90',
    'src/crop/wofost_soil_interface.f90',
    'src/soil/WC_K_models_04_11.f90',
    'src/crop/tillage.f90',
    'src/crop/wofostnut.f90',
    'src/soil/tridag.f90',
    'src/core/checkmassbal.f90',
    'src/crop/wofost_soil_watern.f90',
    'src/soil/fluxes.f90',
    'src/soil/calcgwl.f90',
    'src/soil/hysteresis.f90',
    'src/atmos/penmon.f90',
    'src/soil/divdra.f90',
    'src/io/sharedsimulation.f90',
    'src/soil/watstor.f90',
    'src/io/sharedexchange.f90',
    'src/core/convertdiscrvert.f90',
    'src/atmos/meteoday.f90',
    'src/soil/oxygenstress.f90',
    'src/crop/irrigation.f90',
    'src/crop/rootextraction.f90',
    'src/soil/frozencond.f90',
    'src/solute/solute.f90',
    'src/core/integral.f90',
    'src/crop/wofost_soil_cropresidues.f90',
    'src/io/readswap.f90',
    'src/core/timecontrol.f90',
    'src/atmos/temperature.f90',
    'src/crop/management_soil.f90',
    'src/crop/wofost_soil_rateconstants.f90',
    'src/soil/soilwater.f90',
    'src/core/initialize.f90',
    'src/crop/wofost_soil_balancecheck.f90',
    'src/io/readmeteo.f90',
    'src/soil/drainage.f90',
    'src/io/swapoutput.f90',
    'src/soil/macrorate.f90',
    'src/soil/calcgrid.f90',
    'src/crop/wofost_soil_parameters.f90',
    'src/atmos/snow.f90',
    'src/soil/boundtop.f90',
    'src/crop/wofost_soil_amendments.f90',
    'src/core/swap.f90',
    'src/soil/macroporeoutput.f90',
    'src/crop/wofost_soil_orgmatn.f90',
    'src/soil/headcalc.f90',
    'src/soil/boundbottom.f90',
    'src/crop/cropgrowth.f90',
    'src/soil/surfacewater.f90',
    'src/soil/macropore.f90',
    'src/core/functions.f90',
    'src/atmos/meteodt.f90',
    'src/core/swap_main.f90'
]

# Download ttutil library dependency from GitHub releases
ttutil_version = '4.2.7'
ttutil_lib_dir = join_paths(meson.current_build_dir(), 'lib')

# Determine platform-specific library name
if host_machine.system() == 'windows' and fc.get_id() != 'gcc'
    # Native Windows compilation with Intel Fortran
    ttutil_filename = 'libttutil@0@-windows.a'.format(ttutil_version)
elif host_machine.system() == 'windows' and fc.get_id() == 'gcc'
    # MinGW cross-compilation 
    ttutil_filename = 'libttutil@0@-mingw.a'.format(ttutil_version)
else
    # Linux compilation
    ttutil_filename = 'libttutil@0@-linux.a'.format(ttutil_version)
endif

ttutil_lib_path = join_paths(ttutil_lib_dir, 'libttutil.a')
ttutil_url = 'https://github.com/SWAP-model/ttutil/releases/download/v@0@/@1@'.format(ttutil_version, ttutil_filename)

# Download ttutil library if not present
if not fs.exists(ttutil_lib_path)
    # Create directory
    run_command('mkdir', '-p', ttutil_lib_dir, check: true)
    run_command('curl', '-L', ttutil_url, '-o', ttutil_lib_path, check: true)
endif

# Create ttutil dependency
ttutil_dep = declare_dependency(
    link_args: [ttutil_lib_path]
)

# Build fully static executable with downloaded TTUTIL library
if host_machine.system() == 'windows'
    if fc.get_id() == 'gcc'
        # Windows cross-compilation with MinGW
        link_args = ['-static']
        executable_name = 'swap'
    else
        # Windows native compilation with Intel
        # No additional link args needed - static libs already specified above
        link_args = ['/SUBSYSTEM:CONSOLE']
        executable_name = 'swap'
    endif
else
    # Linux native compilation with Intel
    link_args = ['-static']
    executable_name = 'swap'
endif

executable(executable_name, sources,
    dependencies: [ttutil_dep],
    link_args : link_args,
    include_directories: src_inc,
    install : false
)

# Build SWAP BMI shared library for Python
bmi_sources = sources + ['src/core/swap_bmi.f90']

shared_library('swap',
    bmi_sources,
    dependencies: [ttutil_dep],
    include_directories: src_inc,
    install: true,
    install_dir: get_option('libdir')
)