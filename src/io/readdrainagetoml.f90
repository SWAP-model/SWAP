!> TOML-based drainage configuration reader.
!!
!! Loads drainage configuration from `swap.dra.toml`-style files directly
!! into `drainage_state_t`.
!!
!! @note
!! Legacy drainage parsing is implemented in `ReadSwap/rddrb` and related
!! routines. This module is a state-first TOML reader used by
!! `ReadSwapToml_state`.
module readdrainagetoml_mod
    use tomlf, only: toml_table, toml_array, toml_error, toml_load, get_value, len
    use drainage_state_mod, only: drainage_state_t
    use swap_log, only: log_info
    implicit none
    private

    public :: ReadDrainageToml_state

contains

    !> Read drainage TOML file into explicit drainage state.
    !!
    !! @param[inout] drain Drainage state container
    !! @param[in] drainage_toml_path Path to drainage TOML configuration file
    subroutine ReadDrainageToml_state(drain, drainage_toml_path)
        type(drainage_state_t), intent(inout) :: drain
        character(len=*),       intent(in)    :: drainage_toml_path

        type(toml_table), allocatable :: doc
        type(toml_table), pointer     :: tab, subtab, item
        type(toml_array), pointer     :: arr
        type(toml_error), allocatable :: err

        integer :: istat, i, n, idx

        call toml_load(doc, trim(drainage_toml_path), error=err)
        if (allocated(err)) then
            call fatalerr('ReadDrainageToml_state', trim(err%message))
        end if

        call get_value(doc, 'drainage', tab, requested=.false., stat=istat)
        if (istat /= 0 .or. .not. associated(tab)) then
            call log_info('ReadDrainageToml', 'No [drainage] section found in: ' // trim(drainage_toml_path))
            return
        end if

        call get_value(tab, 'basic', subtab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(subtab)) call read_drainage_basic(subtab, drain)

        call get_value(tab, 'extended', subtab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(subtab)) call read_drainage_extended(subtab, drain)

        call get_value(tab, 'basic', subtab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(subtab)) then
            call get_value(subtab, 'resistance', item, requested=.false., stat=istat)
            if (istat == 0 .and. associated(item)) then
                call get_value(item, 'levels', arr, requested=.false., stat=istat)
                if (istat == 0 .and. associated(arr)) then
                    n = min(len(arr), max(0, drain%nrlevs))
                    do i = 1, n
                        call get_value(arr, i, item, stat=istat)
                        if (istat /= 0) cycle

                        idx = i
                        call get_value(item, 'level', idx, default=idx)
                        if (.not. allocated(drain%drares)) cycle
                        if (idx < 1 .or. idx > size(drain%drares)) cycle

                        call get_value(item, 'drares', drain%drares(idx), default=drain%drares(idx))
                        call get_value(item, 'infres', drain%infres(idx), default=drain%infres(idx))
                        call get_value(item, 'swallo', drain%swallo(idx), default=drain%swallo(idx))
                        call get_value(item, 'l',      drain%L(idx),      default=drain%L(idx))
                        call get_value(item, 'zbotdr', drain%zbotdr(idx), default=drain%zbotdr(idx))
                        call get_value(item, 'swdtyp', drain%swdtyp(idx), default=drain%swdtyp(idx))
                    end do
                end if
            end if
        end if

        call get_value(tab, 'extended', subtab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(subtab)) then
            call get_value(subtab, 'systems', arr, requested=.false., stat=istat)
            if (istat == 0 .and. associated(arr)) then
                n = min(len(arr), max(0, drain%nrlevs))
                do i = 1, n
                    call get_value(arr, i, item, stat=istat)
                    if (istat /= 0) cycle

                    idx = i
                    call get_value(item, 'lev', idx, default=idx)
                    if (.not. allocated(drain%rdrain)) cycle
                    if (idx < 1 .or. idx > size(drain%rdrain)) cycle

                    call get_value(item, 'swdtyp',  drain%swdtyp(idx), default=drain%swdtyp(idx))
                    call get_value(item, 'l',       drain%L(idx),      default=drain%L(idx))
                    call get_value(item, 'zbotdre', drain%zbotdr(idx), default=drain%zbotdr(idx))
                    call get_value(item, 'gwlinf',  drain%gwlinf(idx), default=drain%gwlinf(idx))
                    call get_value(item, 'rdrain',  drain%rdrain(idx), default=drain%rdrain(idx))
                    call get_value(item, 'rinfi',   drain%rinfi(idx),  default=drain%rinfi(idx))
                    call get_value(item, 'rentry',  drain%rentry(idx), default=drain%rentry(idx))
                    call get_value(item, 'rexit',   drain%rexit(idx),  default=drain%rexit(idx))
                    call get_value(item, 'widthr',  drain%widthr(idx), default=drain%widthr(idx))
                    call get_value(item, 'taludr',  drain%taludr(idx), default=drain%taludr(idx))
                end do
            end if
        end if

        call log_info('ReadDrainageToml', 'Loaded drainage TOML: ' // trim(drainage_toml_path))

    end subroutine ReadDrainageToml_state


    !> Read `[drainage.basic]` section.
    !!
    !! @param[in] basic Basic drainage TOML table
    !! @param[inout] drain Drainage state container
    subroutine read_drainage_basic(basic, drain)
        type(toml_table), intent(inout)        :: basic
        type(drainage_state_t), intent(inout)  :: drain

        type(toml_table), pointer :: tab, subtab
        type(toml_array), pointer :: arr
        integer :: istat, n

        call get_value(basic, 'general', tab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(tab)) then
            call get_value(tab, 'dramet',   drain%dramet,   default=drain%dramet)
            call get_value(tab, 'swdivd',   drain%swdivd,   default=drain%swdivd)
            call get_value(tab, 'swdislay', drain%swdislay, default=drain%swdislay)

            call get_value(tab, 'dislay_table', subtab, requested=.false., stat=istat)
            if (istat == 0 .and. associated(subtab)) then
                call get_value(subtab, 'swtopdislay', arr, requested=.false., stat=istat)
                if (istat == 0 .and. associated(arr) .and. allocated(drain%swtopdislay)) then
                    n = min(len(arr), size(drain%swtopdislay))
                    call read_int_array(arr, drain%swtopdislay, n)
                end if

                call get_value(subtab, 'ztopdislay', arr, requested=.false., stat=istat)
                if (istat == 0 .and. associated(arr) .and. allocated(drain%zTopDisLay)) then
                    n = min(len(arr), size(drain%zTopDisLay))
                    call read_real_array(arr, drain%zTopDisLay, n)
                end if

                call get_value(subtab, 'ftopdislay', arr, requested=.false., stat=istat)
                if (istat == 0 .and. associated(arr) .and. allocated(drain%fTopDisLay)) then
                    n = min(len(arr), size(drain%fTopDisLay))
                    call read_real_array(arr, drain%fTopDisLay, n)
                end if
            end if
        end if

        call get_value(basic, 'hooghoudt_ernst', tab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(tab)) then
            call get_value(tab, 'shape',  drain%shape,  default=drain%shape)
            call get_value(tab, 'entres', drain%entres, default=drain%entres)
            call get_value(tab, 'basegw', drain%basegw, default=drain%basegw)

            if (allocated(drain%zbotdr)) then
                if (size(drain%zbotdr) >= 1) call get_value(tab, 'zbotdr', drain%zbotdr(1), default=drain%zbotdr(1))
            end if
            if (allocated(drain%wetper)) then
                if (size(drain%wetper) >= 1) call get_value(tab, 'wetper', drain%wetper(1), default=drain%wetper(1))
            end if
        end if

        call get_value(basic, 'resistance', tab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(tab)) then
            call get_value(tab, 'nrlevs', drain%nrlevs, default=drain%nrlevs)
            call get_value(tab, 'swintfl', drain%swnrsrf, default=drain%swnrsrf)
            call get_value(tab, 'cofintflb', drain%cofintfl, default=drain%cofintfl)
            call get_value(tab, 'expintflb', drain%expintfl, default=drain%expintfl)
            call get_value(tab, 'swtopnrsrf', drain%SwTopnrsrf, default=drain%SwTopnrsrf)
        end if

    end subroutine read_drainage_basic


    !> Read `[drainage.extended]` section.
    !!
    !! @param[in] ext Extended drainage TOML table
    !! @param[inout] drain Drainage state container
    subroutine read_drainage_extended(ext, drain)
        type(toml_table), intent(inout)        :: ext
        type(drainage_state_t), intent(inout)  :: drain

        type(toml_table), pointer :: tab
        integer :: istat

        call get_value(ext, 'characteristics', tab, requested=.false., stat=istat)
        if (istat == 0 .and. associated(tab)) then
            call get_value(tab, 'nrsrf', drain%nrlevs, default=drain%nrlevs)
            call get_value(tab, 'swnrsrf', drain%swnrsrf, default=drain%swnrsrf)
            call get_value(tab, 'rsurfdeep', drain%rsurfdeep, default=drain%rsurfdeep)
            call get_value(tab, 'rsurfshallow', drain%rsurfshallow, default=drain%rsurfshallow)
        end if

    end subroutine read_drainage_extended


    !> Read a TOML integer array into a state array.
    !!
    !! @param[in] arr TOML array
    !! @param[inout] target Target Fortran array
    !! @param[in] n Number of values to copy
    subroutine read_int_array(arr, target, n)
        type(toml_array), intent(inout) :: arr
        integer, intent(inout)          :: target(:)
        integer, intent(in)             :: n

        integer :: i, istat

        do i = 1, n
            call get_value(arr, i, target(i), stat=istat)
        end do
    end subroutine read_int_array


    !> Read a TOML real array into a state array.
    !!
    !! @param[in] arr TOML array
    !! @param[inout] target Target Fortran array
    !! @param[in] n Number of values to copy
    subroutine read_real_array(arr, target, n)
        type(toml_array), intent(inout) :: arr
        real(8), intent(inout)          :: target(:)
        integer, intent(in)             :: n

        integer :: i, istat

        do i = 1, n
            call get_value(arr, i, target(i), stat=istat)
        end do
    end subroutine read_real_array

end module readdrainagetoml_mod
