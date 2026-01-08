project('swap', 'fortran', version : '4.2.0', default_options : ['warning_level=0'])

fc = meson.get_compiler('fortran')
fs = import('fs')

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
    
    if host_machine.system() == 'windows'
        # Windows: Force static linking for ifx compiler
        add_project_link_arguments([
            '/Qstatic',           # Force static linking of Intel runtime libraries
            '/SUBSYSTEM:CONSOLE', # Console application
            '/INCREMENTAL:NO',    # No incremental linking
        ], language : 'fortran')
    else
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

# Source files in the exact order from the original Intel compilation script
sources = [
  'src/variables.f90',
  'src/swap_csv_output.f90', 
  'src/sptabulated.f90',
  'src/wofost_soil_declarations.f90',
  'src/wofost_soil_interface.f90',
  'src/WC_K_models_04_11.f90',
  'src/tillage.f90',
  'src/wofostnut.f90',
  'src/tridag.f90',
  'src/checkmassbal.f90',
  'src/wofost_soil_watern.f90',
  'src/fluxes.f90',
  'src/calcgwl.f90',
  'src/hysteresis.f90',
  'src/penmon.f90',
  'src/divdra.f90',
  'src/sharedsimulation.f90',
  'src/watstor.f90',
  'src/sharedexchange.f90',
  'src/convertdiscrvert.f90',
  'src/meteoday.f90',
  'src/oxygenstress.f90',
  'src/irrigation.f90',
  'src/rootextraction.f90',
  'src/frozencond.f90',
  'src/solute.f90',
  'src/integral.f90',
  'src/wofost_soil_cropresidues.f90',
  'src/readswap.f90',
  'src/timecontrol.f90',
  'src/temperature.f90',
  'src/management_soil.f90',
  'src/wofost_soil_rateconstants.f90',
  'src/soilwater.f90',
  'src/initialize.f90',
  'src/wofost_soil_balancecheck.f90',
  'src/readmeteo.f90',
  'src/drainage.f90',
  'src/swapoutput.f90',
  'src/macrorate.f90',
  'src/calcgrid.f90',
  'src/wofost_soil_parameters.f90',
  'src/snow.f90',
  'src/boundtop.f90',
  'src/wofost_soil_amendments.f90',
  'src/swap.f90',
  'src/macroporeoutput.f90',
  'src/wofost_soil_orgmatn.f90',
  'src/headcalc.f90',
  'src/boundbottom.f90',
  'src/cropgrowth.f90',
  'src/surfacewater.f90',
  'src/macropore.f90',
  'src/functions.f90',
  'src/meteodt.f90',
  'src/swap_main.f90'
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
        link_args = ['Qstatic']
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
    install : false
)