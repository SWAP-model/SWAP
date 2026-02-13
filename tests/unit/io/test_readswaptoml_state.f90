!> Standalone unit test for TOML-based SWAP configuration reader.
!!
!! Verifies representative key mappings from a temporary TOML file into
!! `swap_state_t` using `ReadSwapToml_state`.
program test_readswaptoml_state
    use swap_state_mod, only: swap_state_t, swap_state_init, swap_state_finalize
    use readswaptoml_mod, only: ReadSwapToml_state
    implicit none

    call test_readswaptoml_maps_core_fields()
    call test_readswaptoml_reads_drainage_file()
    print *, 'test_readswaptoml_state: PASS'

contains

!> Validate that representative TOML keys are mapped to `swap_state_t`.
    subroutine test_readswaptoml_maps_core_fields()
        type(swap_state_t) :: state
        integer :: iu
        character(len=*), parameter :: cfg = 'unit_readswaptoml_config.toml'

        call swap_state_init(state, numnod=3, numlay=1, nrlevs=1, ncrop=1)

        open(newunit=iu, file=cfg, status='replace', action='write')
        write(iu,'(A)') '[general]'
        write(iu,'(A)') 'project = "unit_project"'
        write(iu,'(A)') 'swscre = 2'
        write(iu,'(A)') ''
        write(iu,'(A)') '[general.paths]'
        write(iu,'(A)') 'work = "./work"'
        write(iu,'(A)') 'atmosphere = "./atm"'
        write(iu,'(A)') 'crop = "./crop"'
        write(iu,'(A)') 'drain = "./drain"'
        write(iu,'(A)') ''
        write(iu,'(A)') '[simulation]'
        write(iu,'(A)') 'start_date = 2025-01-02'
        write(iu,'(A)') 'end_date = 2025-01-05'
        write(iu,'(A)') ''
        write(iu,'(A)') '[output]'
        write(iu,'(A)') 'file_prefix = "unit_out"'
        write(iu,'(A)') 'nprintday = 1'
        write(iu,'(A)') 'swheader = 1'
        write(iu,'(A)') ''
        write(iu,'(A)') '[output.timing]'
        write(iu,'(A)') 'swmonth = 0'
        write(iu,'(A)') 'period = 7'
        write(iu,'(A)') 'swres = 1'
        write(iu,'(A)') 'swodat = 0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[meteorology]'
        write(iu,'(A)') 'file = "meteo.bbc"'
        write(iu,'(A)') 'lat = 52.1'
        write(iu,'(A)') ''
        write(iu,'(A)') '[meteorology.evapotranspiration]'
        write(iu,'(A)') 'swetr = 1'
        write(iu,'(A)') 'alt = 10.0'
        write(iu,'(A)') 'altw = 2.0'
        write(iu,'(A)') 'angstrom_a = 0.25'
        write(iu,'(A)') 'angstrom_b = 0.50'
        write(iu,'(A)') 'swdivide = 0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[meteorology.temporal]'
        write(iu,'(A)') 'swmetdetail = 0'
        write(iu,'(A)') 'nmetdetail = 1'
        write(iu,'(A)') 'swetsine = 0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[meteorology.rainfall]'
        write(iu,'(A)') 'swrain = 1'
        write(iu,'(A)') 'rainfall_file = "rain.bbc"'
        write(iu,'(A)') ''
        write(iu,'(A)') '[crop.rotation]'
        write(iu,'(A)') 'crop_file = ["maizes.toml", "potatod.toml"]'
        write(iu,'(A)') 'start_date = ["2025-04-01", "2025-09-01"]'
        write(iu,'(A)') 'end_date = ["2025-08-15", "2026-02-01"]'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage]'
        write(iu,'(A)') 'drainage_file = "swap.dra.toml"'
        write(iu,'(A)') ''
        write(iu,'(A)') '[boundary.bottom]'
        write(iu,'(A)') 'swbotb = 6'
        write(iu,'(A)') 'swqhbot = 2'
        write(iu,'(A)') 'deepgw = -150.0'
        write(iu,'(A)') 'hbot = -120.0'
        write(iu,'(A)') 'qbot = -0.01'
        write(iu,'(A)') ''
        write(iu,'(A)') '[boundary.top]'
        write(iu,'(A)') 'swpondmx = 1'
        write(iu,'(A)') 'pondmx = 2.5'
        write(iu,'(A)') 'rsro = 0.8'
        write(iu,'(A)') 'runon = 0.12'
        write(iu,'(A)') 'qtop = 0.03'
        write(iu,'(A)') ''
        write(iu,'(A)') '[soil.initial]'
        write(iu,'(A)') 'gwli = -75.0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[soil.surface]'
        write(iu,'(A)') 'pondmx = 1.5'
        write(iu,'(A)') ''
        write(iu,'(A)') '[soil.evaporation]'
        write(iu,'(A)') 'cofredbl = 0.35'
        write(iu,'(A)') ''
        write(iu,'(A)') '[soil.hysteresis]'
        write(iu,'(A)') 'swhyst = 2'
        write(iu,'(A)') 'tau = 0.2'
        write(iu,'(A)') ''
        write(iu,'(A)') '[soil.numerical]'
        write(iu,'(A)') 'gwlconv = 88.0'
        write(iu,'(A)') 'critdevh1cp = 0.02'
        write(iu,'(A)') 'critdevh2cp = 0.12'
        write(iu,'(A)') 'maxit = 45'
        close(iu)

        call ReadSwapToml_state(state, cfg)

        call assert_equal_char('unit_project', trim(state%time%project), 'state%time%project')
        call assert_equal_int(2, state%time%swscre, 'state%time%swscre')
        call assert_equal_char('./work', trim(state%time%pathwork), 'state%time%pathwork')
        call assert_equal_char('./atm', trim(state%atm%pathatm), 'state%atm%pathatm')
        call assert_equal_char('./crop', trim(state%crop%pathcrop), 'state%crop%pathcrop')
        call assert_equal_char('./drain', trim(state%drain%pathdrain), 'state%drain%pathdrain')
        call assert_equal_char('unit_out', trim(state%time%outfil), 'state%time%outfil')
        call assert_equal_int(1, state%time%nprintday, 'state%time%nprintday')
        call assert_equal_char('meteo.bbc', trim(state%atm%metfil), 'state%atm%metfil')
        call assert_equal_char('rain.bbc', trim(state%atm%rainfil), 'state%atm%rainfil')
        call assert_equal_int(2, state%ncrop, 'state%ncrop')
        call assert_equal_char('maizes.toml', trim(state%crop%cropfil(1)), 'state%crop%cropfil(1)')
        call assert_equal_char('potatod.toml', trim(state%crop%cropfil(2)), 'state%crop%cropfil(2)')
        call assert_true(state%time%tstart > 0.0d0, 'state%time%tstart > 0')
        call assert_equal_real(3.0d0, state%time%tend - state%time%tstart, 1.0d-9, 'tend - tstart')
        call assert_equal_char('swap.dra.toml', trim(state%drain%drfil), 'state%drain%drfil')
        call assert_equal_int(6, state%boundary%swbotb, 'state%boundary%swbotb')
        call assert_equal_int(2, state%boundary%swqhbot, 'state%boundary%swqhbot')
        call assert_equal_real(-150.0d0, state%boundary%deepgw, 1.0d-12, 'state%boundary%deepgw')
        call assert_equal_real(-120.0d0, state%boundary%hbot, 1.0d-12, 'state%boundary%hbot')
        call assert_equal_real(-0.01d0, state%boundary%qbot, 1.0d-12, 'state%boundary%qbot')
        call assert_equal_int(1, state%boundary%swpondmx, 'state%boundary%swpondmx')
        call assert_equal_real(2.5d0, state%boundary%pondmx, 1.0d-12, 'state%boundary%pondmx')
        call assert_equal_real(0.8d0, state%boundary%rsro, 1.0d-12, 'state%boundary%rsro')
        call assert_equal_real(0.12d0, state%boundary%runon, 1.0d-12, 'state%boundary%runon')
        call assert_equal_real(0.03d0, state%boundary%qtop, 1.0d-12, 'state%boundary%qtop')
        call assert_equal_real(-75.0d0, state%soil%gwli, 1.0d-12, 'state%soil%gwli')
        call assert_equal_real(-75.0d0, state%soil%gwl, 1.0d-12, 'state%soil%gwl')
        call assert_equal_real(1.5d0, state%soil%pondmx, 1.0d-12, 'state%soil%pondmx')
        call assert_equal_real(0.35d0, state%soil%cofred, 1.0d-12, 'state%soil%cofred')
        call assert_equal_int(2, state%soil%swhyst, 'state%soil%swhyst')
        call assert_equal_real(0.2d0, state%soil%tau, 1.0d-12, 'state%soil%tau')
        call assert_equal_real(88.0d0, state%soil%gwlconv, 1.0d-12, 'state%soil%gwlconv')
        call assert_equal_real(0.02d0, state%soil%CritDevh1Cp, 1.0d-12, 'state%soil%CritDevh1Cp')
        call assert_equal_real(0.12d0, state%soil%CritDevh2Cp, 1.0d-12, 'state%soil%CritDevh2Cp')
        call assert_equal_int(45, state%soil%msteps, 'state%soil%msteps')
        call assert_true(allocated(state%surfwater%owltab), 'state%surfwater%owltab allocated')
        call assert_equal_int(1, size(state%surfwater%owltab, 1), 'size(state%surfwater%owltab,1)')
        call assert_equal_int(10, state%surfwater%nmper, 'state%surfwater%nmper default')

        call swap_state_finalize(state)

        call assert_equal_char('', trim(state%atm%metfil), 'state%atm%metfil after finalize')
        call assert_equal_char('', trim(state%atm%rainfil), 'state%atm%rainfil after finalize')
        call assert_equal_char('', trim(state%atm%pathatm), 'state%atm%pathatm after finalize')
        call assert_equal_real(0.0d0, state%atm%lat, 1.0d-12, 'state%atm%lat after finalize')
        call assert_equal_int(0, state%boundary%swbotb, 'state%boundary%swbotb after finalize')
        call assert_equal_real(0.0d0, state%boundary%qbot, 1.0d-12, 'state%boundary%qbot after finalize')
        call assert_equal_real(0.0d0, state%boundary%pondmx, 1.0d-12, 'state%boundary%pondmx after finalize')
        call assert_equal_real(0.0d0, state%soil%gwl, 1.0d-12, 'state%soil%gwl after finalize')
        call assert_equal_real(0.0d0, state%soil%pondmx, 1.0d-12, 'state%soil%pondmx after finalize')
        call assert_equal_int(0, state%soil%swhyst, 'state%soil%swhyst after finalize')
        call assert_true(.not. allocated(state%surfwater%owltab), 'state%surfwater%owltab deallocated after finalize')
        call assert_equal_int(0, state%surfwater%nmper, 'state%surfwater%nmper reset after finalize')

        call remove_file_if_exists(cfg)
    end subroutine test_readswaptoml_maps_core_fields

!> Validate that drainage TOML file is read into `state%drain`.
    subroutine test_readswaptoml_reads_drainage_file()
        type(swap_state_t) :: state
        integer :: iu
        character(len=*), parameter :: cfg = 'unit_readswaptoml_drain_config.toml'
        character(len=*), parameter :: drcfg = 'unit_swap.dra.toml'

        call swap_state_init(state, numnod=3, numlay=1, nrlevs=5, ncrop=1)

        open(newunit=iu, file=cfg, status='replace', action='write')
        write(iu,'(A)') '[general.paths]'
        write(iu,'(A)') 'drain = "./"'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage]'
        write(iu,'(A)') 'drainage_file = "unit_swap"'
        close(iu)

        open(newunit=iu, file=drcfg, status='replace', action='write')
        write(iu,'(A)') '[drainage.basic.general]'
        write(iu,'(A)') 'dramet = 3'
        write(iu,'(A)') 'swdivd = 1'
        write(iu,'(A)') 'swdislay = 2'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage.basic.general.dislay_table]'
        write(iu,'(A)') 'swtopdislay = [1, 0]'
        write(iu,'(A)') 'ztopdislay = [-120.0, -80.0]'
        write(iu,'(A)') 'ftopdislay = [0.4, 0.6]'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage.basic.hooghoudt_ernst]'
        write(iu,'(A)') 'shape = 0.8'
        write(iu,'(A)') 'entres = 20.0'
        write(iu,'(A)') 'basegw = -200.0'
        write(iu,'(A)') 'zbotdr = -90.0'
        write(iu,'(A)') 'wetper = 30.0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage.basic.resistance]'
        write(iu,'(A)') 'nrlevs = 2'
        write(iu,'(A)') 'swintfl = 1'
        write(iu,'(A)') 'cofintflb = 0.5'
        write(iu,'(A)') 'expintflb = 0.7'
        write(iu,'(A)') 'swtopnrsrf = 1'
        write(iu,'(A)') ''
        write(iu,'(A)') '[[drainage.basic.resistance.levels]]'
        write(iu,'(A)') 'level = 1'
        write(iu,'(A)') 'drares = 100.0'
        write(iu,'(A)') 'infres = 200.0'
        write(iu,'(A)') 'swallo = 1'
        write(iu,'(A)') 'l = 20.0'
        write(iu,'(A)') 'zbotdr = -100.0'
        write(iu,'(A)') 'swdtyp = 2'
        write(iu,'(A)') ''
        write(iu,'(A)') '[[drainage.basic.resistance.levels]]'
        write(iu,'(A)') 'level = 2'
        write(iu,'(A)') 'drares = 110.0'
        write(iu,'(A)') 'infres = 210.0'
        write(iu,'(A)') 'swallo = 2'
        write(iu,'(A)') 'l = 25.0'
        write(iu,'(A)') 'zbotdr = -110.0'
        write(iu,'(A)') 'swdtyp = 1'
        write(iu,'(A)') ''
        write(iu,'(A)') '[drainage.extended.characteristics]'
        write(iu,'(A)') 'nrsrf = 2'
        write(iu,'(A)') 'swnrsrf = 1'
        write(iu,'(A)') 'rsurfdeep = 30.0'
        write(iu,'(A)') 'rsurfshallow = 10.0'
        write(iu,'(A)') ''
        write(iu,'(A)') '[[drainage.extended.systems]]'
        write(iu,'(A)') 'lev = 1'
        write(iu,'(A)') 'swdtyp = 0'
        write(iu,'(A)') 'l = 250.0'
        write(iu,'(A)') 'zbotdre = -1093.0'
        write(iu,'(A)') 'gwlinf = -350.0'
        write(iu,'(A)') 'rdrain = 150.0'
        write(iu,'(A)') 'rinfi = 4000.0'
        write(iu,'(A)') 'rentry = 0.8'
        write(iu,'(A)') 'rexit = 0.9'
        write(iu,'(A)') 'widthr = 100.0'
        write(iu,'(A)') 'taludr = 0.66'
        close(iu)

        call ReadSwapToml_state(state, cfg)

        call assert_equal_char('unit_swap', trim(state%drain%drfil), 'state%drain%drfil')
        call assert_equal_int(2, state%drain%nrlevs, 'state%drain%nrlevs')
        call assert_equal_int(3, state%drain%dramet, 'state%drain%dramet')
        call assert_equal_int(1, state%drain%swdivd, 'state%drain%swdivd')
        call assert_equal_int(2, state%drain%swdislay, 'state%drain%swdislay')
        call assert_equal_real(0.8d0, state%drain%shape, 1.0d-12, 'state%drain%shape')
        call assert_equal_real(20.0d0, state%drain%entres, 1.0d-12, 'state%drain%entres')
        call assert_equal_real(-200.0d0, state%drain%basegw, 1.0d-12, 'state%drain%basegw')
        call assert_equal_real(-100.0d0, state%drain%zbotdr(1), 1.0d-12, 'state%drain%zbotdr(1)')
        call assert_equal_real(110.0d0, state%drain%drares(2), 1.0d-12, 'state%drain%drares(2)')
        call assert_equal_real(210.0d0, state%drain%infres(2), 1.0d-12, 'state%drain%infres(2)')
        call assert_equal_int(1, state%drain%swallo(1), 'state%drain%swallo(1)')
        call assert_equal_real(250.0d0, state%drain%L(1), 1.0d-12, 'state%drain%L(1)')
        call assert_equal_real(-1093.0d0, state%drain%zbotdr(1), 1.0d-12, 'state%drain%zbotdr(1) extended')
        call assert_equal_real(-350.0d0, state%drain%gwlinf(1), 1.0d-12, 'state%drain%gwlinf(1)')
        call assert_equal_real(150.0d0, state%drain%rdrain(1), 1.0d-12, 'state%drain%rdrain(1)')
        call assert_equal_real(4000.0d0, state%drain%rinfi(1), 1.0d-12, 'state%drain%rinfi(1)')
        call assert_equal_real(0.8d0, state%drain%rentry(1), 1.0d-12, 'state%drain%rentry(1)')
        call assert_equal_real(0.9d0, state%drain%rexit(1), 1.0d-12, 'state%drain%rexit(1)')
        call assert_equal_real(100.0d0, state%drain%widthr(1), 1.0d-12, 'state%drain%widthr(1)')
        call assert_equal_real(0.66d0, state%drain%taludr(1), 1.0d-12, 'state%drain%taludr(1)')
        call assert_equal_int(1, state%drain%swnrsrf, 'state%drain%swnrsrf')
        call assert_equal_real(30.0d0, state%drain%rsurfdeep, 1.0d-12, 'state%drain%rsurfdeep')
        call assert_equal_real(10.0d0, state%drain%rsurfshallow, 1.0d-12, 'state%drain%rsurfshallow')
        call assert_equal_int(1, state%drain%swtopdislay(1), 'state%drain%swtopdislay(1)')
        call assert_equal_real(-120.0d0, state%drain%zTopDisLay(1), 1.0d-12, 'state%drain%zTopDisLay(1)')
        call assert_equal_real(0.6d0, state%drain%fTopDisLay(2), 1.0d-12, 'state%drain%fTopDisLay(2)')
        call assert_true(allocated(state%surfwater%qqhtab), 'state%surfwater%qqhtab allocated')
        call assert_equal_int(5, size(state%surfwater%owltab, 1), 'size(state%surfwater%owltab,1) for nrlevs=5')

        call swap_state_finalize(state)

        call remove_file_if_exists(cfg)
        call remove_file_if_exists(drcfg)
    end subroutine test_readswaptoml_reads_drainage_file

!> Assert integer equality.
    subroutine assert_equal_int(expected, actual, label)
        integer, intent(in) :: expected, actual
        character(len=*), intent(in) :: label

        if (actual /= expected) then
            write(*,'(A,1X,A,1X,I0,1X,A,1X,I0)') 'Assertion failed:', trim(label), actual, '/=', expected
            error stop 1
        end if
    end subroutine assert_equal_int

!> Assert real equality with absolute tolerance.
    subroutine assert_equal_real(expected, actual, tol, label)
        real(8), intent(in) :: expected, actual, tol
        character(len=*), intent(in) :: label

        if (abs(actual - expected) > tol) then
            write(*,'(A,1X,A,1X,ES16.8,1X,A,1X,ES16.8)') 'Assertion failed:', trim(label), actual, '/=', expected
            error stop 1
        end if
    end subroutine assert_equal_real

!> Assert character equality.
    subroutine assert_equal_char(expected, actual, label)
        character(len=*), intent(in) :: expected, actual
        character(len=*), intent(in) :: label

        if (trim(actual) /= trim(expected)) then
            write(*,'(A,1X,A,1X,A,1X,A,1X,A)') 'Assertion failed:', trim(label), '"' // trim(actual) // '"', '/=', '"' // trim(expected) // '"'
            error stop 1
        end if
    end subroutine assert_equal_char

!> Assert logical truth.
    subroutine assert_true(condition, label)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: label

        if (.not. condition) then
            write(*,'(A,1X,A)') 'Assertion failed:', trim(label)
            error stop 1
        end if
    end subroutine assert_true

!> Remove a file when it exists.
    subroutine remove_file_if_exists(path)
        character(len=*), intent(in) :: path
        integer :: iu
        logical :: exists

        inquire(file=path, exist=exists)
        if (exists) then
            open(newunit=iu, file=path, status='old', action='readwrite')
            close(iu, status='delete')
        end if
    end subroutine remove_file_if_exists

end program test_readswaptoml_state
