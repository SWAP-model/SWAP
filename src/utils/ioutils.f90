module io_utils_mod
    implicit none
    private
    public :: parse_output_extensions

contains

    subroutine parse_output_extensions(extensions, n_ext, &
                                    swwba, swend, swvap, swbal, swblc, &
                                    swsba, swate, swbma, swdrf, swswb, &
                                    swini, swinc, swcrp, swstr, swirg, &
                                    csv, csv_tz)
        implicit none
        
        ! Input
        integer, intent(in) :: n_ext
        character(len=3), dimension(n_ext), intent(in) :: extensions
        
        ! Output - individual switches
        integer, intent(out) :: swwba, swend, swvap, swbal, swblc
        integer, intent(out) :: swsba, swate, swbma, swdrf, swswb
        integer, intent(out) :: swini, swinc, swcrp, swstr, swirg
        integer, intent(out) :: csv, csv_tz
        
        ! Local
        integer :: i
        character(len=3) :: ext_lower
        
        ! Initialize all switches to 0
        swwba = 0; swend = 0; swvap = 0; swbal = 0; swblc = 0
        swsba = 0; swate = 0; swbma = 0; swdrf = 0; swswb = 0
        swini = 0; swinc = 0; swcrp = 0; swstr = 0; swirg = 0
        csv = 0; csv_tz = 0
        
        ! Parse extensions and set corresponding switches
        do i = 1, n_ext
            ext_lower = to_lower(extensions(i))
            
            select case (trim(ext_lower))
                case ('wba')
                    swwba = 1
                case ('end')
                    swend = 1
                case ('vap')
                    swvap = 1
                case ('bal')
                    swbal = 1
                case ('blc')
                    swblc = 1
                case ('sba')
                    swsba = 1
                case ('ate')
                    swate = 1
                case ('bma')
                    swbma = 1
                case ('drf')
                    swdrf = 1
                case ('swb')
                    swswb = 1
                case ('ini')
                    swini = 1
                case ('inc')
                    swinc = 1
                case ('crp')
                    swcrp = 1
                case ('str')
                    swstr = 1
                case ('irg')
                    swirg = 1
                case ('csv')
                    swcsv = 1
                case ('csv_tz')
                    swcsv_tz = 1
                case default
                    ! Optionally handle unknown extensions
                    continue
            end select
        end do
        
    end subroutine parse_output_extensions

    !> Helper function to convert to lowercase
    !!
    !! Even though fortran is case insensitive, in this kind of comparison
    !! we are checking if character data (not identifier) are the same, and caps and lower case
    !! letters have a different value.
    !!
    function to_lower(str) result(lower_str)
        character(len=*), intent(in) :: str
        character(len=len(str)) :: lower_str
        integer :: i, ic
        
        lower_str = str
        do i = 1, len_trim(str)
            ic = iachar(str(i:i))
            if (ic >= iachar('A') .and. ic <= iachar('Z')) then
                lower_str(i:i) = achar(ic + 32)
            end if
        end do
    end function to_lower

end module io_utils_mod