!=====================================================================
! text_parser
!---------------------------------------------------------------------
! Low-level parser for structured tabular text files.
!
! Supported formats:
!   - CSV files (.csv)
!   - Whitespace-delimited text files (.txt, .dat)
!
! Responsibilities:
!   - Detect file format from extension
!   - Read and validate column headers
!   - Skip blank and comment lines
!   - Validate consistent column counts across rows
!   - Store raw table contents as strings
!   - Preserve original file line numbers for diagnostics
!
! Design notes:
!   - CSV parsing preserves empty fields
!   - Whitespace parsing treats repeated spaces/tabs as a
!     single separator
!   - All values are stored as fixed-length character strings
!   - No support for quotation marks in csv files. 
!
! Example:
!
!   type(TextTable) :: table
!
!   call parse_text_table('forcing.txt', table)
!
!   ! Access:
!   ! table%header(:)
!   ! table%values(irow,icol)    string valued (because of timestamps)
!   ! table%line_numbers(:)
!
!   call table%free()
!
!=====================================================================
module text_parser
    use str_utils, only: to_lower, to_upper, ends_with, inttostr
    implicit none
    private

    public :: TextTable
    public :: parse_text_table

    integer, parameter :: FORMAT_UNKNOWN    = 0
    integer, parameter :: FORMAT_CSV        = 1
    integer, parameter :: FORMAT_WHITESPACE = 2
    
    integer, parameter :: STR_LEN = 512

    type :: TextTable
        character(len=STR_LEN), allocatable :: header(:)
        character(len=STR_LEN), allocatable :: values(:,:)   ! (nrow, ncol)
        integer, allocatable :: line_numbers(:)
        integer :: nrow = 0
        integer :: ncol = 0
        integer :: format = FORMAT_UNKNOWN

    contains
        procedure :: find_column => text_table_find_column
        procedure :: free        => text_table_free
    end type TextTable

contains

    subroutine parse_text_table(filename, table)
        character(*), intent(in) :: filename
        type(TextTable), intent(out) :: table

        integer :: unit, ios, header_line
        character(len=STR_LEN), allocatable :: header(:)
        integer :: ncol
        logical :: ok
        character(len=1024) :: errmsg

        errmsg = ''
        ok = .false.

        call table%free()

        table%format = infer_format(filename)

        if (table%format == FORMAT_UNKNOWN) then
            errmsg = 'Unsupported text file extension: '//trim(filename)// &
                    '. Expected .csv, .txt, or .dat.'
            write(*,'(A)') trim(errmsg)
            stop 1
        end if

        open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
        if (ios /= 0) then
            errmsg = 'Could not open text file: '//trim(filename)
            write(*,'(A)') trim(errmsg)
            stop 1
        end if

        call read_header(unit, filename, table%format, header, ncol, header_line, ok, errmsg)

        if (.not. ok) then
            close(unit)
            write(*,'(A)') trim(errmsg)
            stop 1
        end if

        allocate(table%header(ncol))
        table%header = header
        table%ncol = ncol

        call count_data_rows(unit, filename, table%format, ncol, header_line, &
                     table%nrow, ok, errmsg)

        if (.not. ok) then
            close(unit)
            write(*,'(A)') trim(errmsg)
            stop 1
        end if

        allocate(table%values(table%nrow, table%ncol))
        allocate(table%line_numbers(table%nrow))

        rewind(unit)

        call read_table_data(unit, filename, table%format, header_line, &
                            table%values, table%line_numbers, ok, errmsg)
        close(unit)

        if (.not. ok) then
            write(*,'(A)') trim(errmsg)
            stop 1
        end if

    end subroutine parse_text_table

    integer function text_table_find_column(self, name) result(idx)
        class(TextTable), intent(in) :: self
        character(*),     intent(in) :: name

        integer :: j

        idx = 0

        do j = 1, self%ncol
            if (to_upper(trim(self%header(j))) == to_upper(trim(name))) then
                idx = j
                return
            end if
        end do
    end function text_table_find_column


    subroutine text_table_free(self)
        class(TextTable), intent(inout) :: self

        if (allocated(self%header)) deallocate(self%header)
        if (allocated(self%values)) deallocate(self%values)
        if (allocated(self%line_numbers)) deallocate(self%line_numbers)

        self%nrow = 0
        self%ncol = 0
        self%format = FORMAT_UNKNOWN
    end subroutine text_table_free

    !=========================================================
    !   Local Subroutines and functions
    !=========================================================

    ! Find out the format of the file
    pure integer function infer_format(filename) result(mode)
        character(*), intent(in) :: filename
        character(len=STR_LEN) :: lower

        lower = to_lower(trim(filename))

        if (ends_with(lower, '.csv')) then
            mode = FORMAT_CSV
        else if (ends_with(lower, '.txt') .or. ends_with(lower, '.dat')) then
            mode = FORMAT_WHITESPACE
        else
            mode = FORMAT_UNKNOWN
        end if
    end function infer_format

    subroutine read_header(unit, filename, format, header, ncol, header_line, ok, errmsg)
        integer, intent(in) :: unit
        character(*), intent(in) :: filename
        integer, intent(in) :: format
        character(len=STR_LEN), allocatable, intent(out) :: header(:)
        integer, intent(out) :: ncol, header_line
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        character(len=STR_LEN) :: line
        integer :: ios, iline

        ok = .false.
        errmsg = ''
        ncol = 0
        header_line = 0
        iline = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            iline = iline + 1

            if (is_blank_line(line) .or. is_comment_line(line)) cycle

            call split_line(line, format, header, ncol)
            header_line = iline

            if (ncol <= 0) then
                errmsg = 'Invalid header in '//trim(filename)//' at line '// &
                        trim(inttostr(iline))//': no columns found.'
                return
            end if

            call validate_header(header, filename, iline, ok, errmsg)
            return
        end do

        errmsg = 'No header found in '//trim(filename)//'.'
    end subroutine read_header

    subroutine count_data_rows(unit, filename, format, ncol, start_line, nrow, ok, errmsg)
        integer, intent(in) :: unit
        character(*), intent(in) :: filename
        integer, intent(in) :: format
        integer, intent(in) :: ncol
        integer, intent(in) :: start_line
        integer, intent(out) :: nrow
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        character(len=STR_LEN) :: line
        integer :: ios, iline, nf

        ok = .false.
        errmsg = ''
        nrow = 0
        iline = start_line

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            iline = iline + 1

            if (is_blank_line(line) .or. is_comment_line(line)) cycle

            nf = count_fields(line, format)

            if (nf /= ncol) then
                errmsg = 'Invalid number of columns in '//trim(filename)// &
                        ' at line '//trim(inttostr(iline))//': expected '// &
                        trim(inttostr(ncol))//' but found '//trim(inttostr(nf))//'.'
                return
            end if

            nrow = nrow + 1
        end do

        if (nrow <= 0) then
            errmsg = 'No data rows found in '//trim(filename)//'.'
            return
        end if

        ok = .true.
    end subroutine count_data_rows

    subroutine read_table_data(unit, filename, format, header_line, values, line_numbers, ok, errmsg)
        integer, intent(in) :: unit
        character(*), intent(in) :: filename
        integer, intent(in) :: format
        integer, intent(in) :: header_line
        character(len=STR_LEN), intent(out) :: values(:,:)
        integer, intent(out) :: line_numbers(:)
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        character(len=STR_LEN) :: line
        character(len=STR_LEN), allocatable :: fields(:)
        integer :: ios, iline, irow, nf

        ok = .false.
        errmsg = ''
        values = ''
        line_numbers = 0

        iline = 0
        irow = 0

        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            iline = iline + 1

            if (iline <= header_line) cycle
            if (is_blank_line(line) .or. is_comment_line(line)) cycle

            call split_line(line, format, fields, nf)

            if (nf /= size(values, 2)) then
                errmsg = 'Invalid number of columns in '//trim(filename)// &
                         ' at line '//trim(inttostr(iline))//': expected '// &
                         trim(inttostr(size(values, 2)))//' but found '// &
                         trim(inttostr(nf))//'.'
                return
            end if

            irow = irow + 1

            if (irow > size(values, 1)) then
                errmsg = 'Internal parser error while reading '//trim(filename)// &
                         ': more data rows found than expected.'
                return
            end if

            values(irow, :) = fields(:)
            line_numbers(irow) = iline
        end do

        if (irow /= size(values, 1)) then
            errmsg = 'Internal parser error while reading '//trim(filename)// &
                     ': expected '//trim(inttostr(size(values, 1)))// &
                     ' data rows but read '//trim(inttostr(irow))//'.'
            return
        end if

        ok = .true.
    end subroutine read_table_data


    pure integer function count_fields(line, format) result(n)
        character(*), intent(in) :: line
        integer,      intent(in) :: format

        select case (format)
        case (FORMAT_CSV)
            n = count_csv_fields(line)

        case (FORMAT_WHITESPACE)
            n = count_whitespace_fields(line)

        case default
            n = 0
        end select
    end function count_fields

    subroutine split_line(line, format, fields, nfields)
        character(*), intent(in) :: line
        integer,      intent(in) :: format
        character(len=STR_LEN), allocatable, intent(out) :: fields(:)
        integer, intent(out) :: nfields

        nfields = count_fields(line, format)

        if (nfields <= 0) then
            allocate(fields(0))
            return
        end if

        allocate(fields(nfields))

        select case (format)
        case (FORMAT_CSV)
            call split_csv_line(line, fields)
        case (FORMAT_WHITESPACE)
            call split_whitespace_line(line, fields)
        end select
    end subroutine split_line


    pure logical function is_blank_line(line) result(blank)
        character(*), intent(in) :: line
        blank = (len_trim(line) == 0)
    end function is_blank_line


    pure logical function is_comment_line(line) result(is_comment)
        character(*), intent(in) :: line
        character(len=len(line)) :: s

        s = adjustl(line)
        is_comment = len_trim(s) > 0 .and. s(1:1) == '#'
    end function is_comment_line

    pure logical function is_space(ch) result(space)
        character(len=1), intent(in) :: ch
        integer :: ia

        ia = iachar(ch)
        space = ia == 32 .or. ia == 9
    end function is_space

    pure integer function count_csv_fields(line) result(n)
        character(*), intent(in) :: line
        integer :: i, last

        last = len_trim(line)

        if (last == 0) then
            n = 0
            return
        end if

        n = 1
        do i = 1, last
            if (line(i:i) == ',') n = n + 1
        end do
    end function count_csv_fields

    pure integer function count_whitespace_fields(line) result(n)
        character(*), intent(in) :: line
        integer :: i, last
        logical :: in_field

        last = len_trim(line)
        n = 0
        in_field = .false.

        do i = 1, last
            if (is_space(line(i:i))) then
                in_field = .false.
            else
                if (.not. in_field) n = n + 1
                in_field = .true.
            end if
        end do
    end function count_whitespace_fields

    subroutine split_whitespace_line(line, fields)
        character(*), intent(in) :: line
        character(len=STR_LEN), intent(out) :: fields(:)

        integer :: i, last, start_pos, nf
        logical :: in_field

        last = len_trim(line)
        nf = 0
        in_field = .false.
        start_pos = 0
        fields = ''

        do i = 1, last
            if (is_space(line(i:i))) then
                if (in_field) then
                    nf = nf + 1
                    fields(nf) = adjustl(line(start_pos:i-1))
                    in_field = .false.
                end if
            else
                if (.not. in_field) then
                    start_pos = i
                    in_field = .true.
                end if
            end if
        end do

        if (in_field) then
            nf = nf + 1
            fields(nf) = adjustl(line(start_pos:last))
        end if
    end subroutine split_whitespace_line

    subroutine split_csv_line(line, fields)
        character(*), intent(in) :: line
        character(len=STR_LEN), intent(out) :: fields(:)

        integer :: i, last, start_pos, nf

        last = len_trim(line)
        start_pos = 1
        nf = 0
        fields = ''

        do i = 1, last
            if (line(i:i) == ',') then
                nf = nf + 1
                fields(nf) = adjustl(line(start_pos:i-1))
                start_pos = i + 1
            end if
        end do

        nf = nf + 1
        fields(nf) = adjustl(line(start_pos:last))
    end subroutine split_csv_line

    subroutine validate_header(header, filename, line_number, ok, errmsg)
        character(len=STR_LEN), intent(in) :: header(:)
        character(*), intent(in) :: filename
        integer, intent(in) :: line_number
        logical, intent(out) :: ok
        character(*), intent(out) :: errmsg

        integer :: i, j

        ok = .false.
        errmsg = ''

        do i = 1, size(header)
            if (len_trim(header(i)) == 0) then
                errmsg = 'Empty column name in '//trim(filename)//' at line '// &
                        trim(inttostr(line_number))//'.'
                return
            end if
        end do

        do i = 1, size(header) - 1
            do j = i + 1, size(header)
                if (trim(header(i)) == trim(header(j))) then
                    errmsg = 'Duplicate column name "'//trim(header(i))//'" in '// &
                            trim(filename)//' at line '//trim(inttostr(line_number))//'.'
                    return
                end if
            end do
        end do

        ok = .true.
    end subroutine validate_header

end module text_parser