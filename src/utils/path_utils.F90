!> Small utilities for manipulating file paths.
!! Pure string-based routines; they do not check whether files/directories exist.
module path_utils

    use str_utils, only: ends_with, to_lower

    implicit none
    private

    public :: path_dirname
    public :: path_basename
    public :: path_stem
    public :: path_extension
    public :: path_join
    public :: path_has_extension
    public :: path_replace_extension

contains

    pure integer function last_separator(path) result(pos)
        character(*), intent(in) :: path
        integer :: p1, p2

        p1 = index(path, '/',  back=.true.)
        p2 = index(path, '\',  back=.true.)
        pos = max(p1, p2)
    end function last_separator


    !> Return directory part of a path.
    !! Examples:
    !!   "output/test.nc" -> "output"
    !!   "test.nc"        -> ""
    pure function path_dirname(path, include_separator) result(dirname)
        character(*), intent(in) :: path
        logical, intent(in), optional :: include_separator
        character(:), allocatable :: dirname

        logical :: keep_sep
        integer :: p, n

        keep_sep = .false.
        if (present(include_separator)) keep_sep = include_separator

        n = len_trim(path)
        p = last_separator(path(:n))

        if (p == 0) then
            dirname = ''
        else
            if (keep_sep) then
                dirname = path(:p)
            else
                dirname = path(:p-1)
            end if
        end if
    end function path_dirname


    !> Return file name part of a path.
    !! Examples:
    !!   "output/test.nc" -> "test.nc"
    !!   "test.nc"        -> "test.nc"
    pure function path_basename(path) result(basename)
        character(*), intent(in) :: path
        character(:), allocatable :: basename
        integer :: p, n

        n = len_trim(path)
        p = last_separator(path(:n))

        if (p > 0 .and. p < n) then
            basename = path(p+1:n)
        else if (p == n) then
            basename = ''
        else
            basename = path(:n)
        end if
    end function path_basename


    !> Return basename without final extension.
    !! Examples:
    !!   "output/test.nc" -> "test"
    !!   "archive.tar.gz" -> "archive.tar"
    pure function path_stem(path) result(stem)
        character(*), intent(in) :: path
        character(:), allocatable :: stem
        character(:), allocatable :: base
        integer :: p

        base = path_basename(path)
        p = index(base, '.', back=.true.)

        if (p > 1) then
            stem = base(:p-1)
        else
            stem = base
        end if
    end function path_stem


    !> Return final extension, including the dot.
    !! Examples:
    !!   "test.nc"        -> ".nc"
    !!   "archive.tar.gz" -> ".gz"
    !!   "README"         -> ""
    pure function path_extension(path) result(ext)
        character(*), intent(in) :: path
        character(:), allocatable :: ext
        character(:), allocatable :: base
        integer :: p, n

        base = path_basename(path)
        n = len_trim(base)
        p = index(base, '.', back=.true.)

        if (p > 1 .and. p < n) then
            ext = base(p:n)
        else
            ext = ''
        end if
    end function path_extension

    pure function path_replace_extension(path, new_ext) result(new_path)
        character(*), intent(in) :: path
        character(*), intent(in) :: new_ext
        character(:), allocatable :: new_path

        character(:), allocatable :: ext

        ! Ensure extension starts with '.'
        if (len_trim(new_ext) == 0) then
            ext = ''
        else if (new_ext(1:1) == '.') then
            ext = trim(new_ext)
        else
            ext = '.' // trim(new_ext)
        end if

        if (len_trim(path_extension(path)) > 0) then
            new_path = path_dirname(path, include_separator=.true.) // &
                    path_stem(path) // ext
        else
            new_path = trim(path) // ext
        end if
    end function path_replace_extension


    !> Join directory and file name using '/'.
    !! If dirname is empty, returns basename.
    pure function path_join(dirname, basename) result(path)
        character(*), intent(in) :: dirname
        character(*), intent(in) :: basename
        character(:), allocatable :: path

        integer :: nd
        character :: sep

        nd = len_trim(dirname)

        ! Empty dirname -> just return basename.
        if (nd == 0) then
            path = trim(basename)
            return
        end if

        ! Already ends in a separator.
        if (dirname(nd:nd) == '/' .or. dirname(nd:nd) == '\') then
            path = trim(dirname) // trim(basename)
            return
        end if

        ! Reuse existing separator style if present.
        if (index(dirname, '\') > 0) then
            sep = '\'
        else
            sep = '/'
        end if

        path = trim(dirname) // sep // trim(basename)

    end function path_join


    !> True if path has the requested extension.
    !! The extension can be given with or without leading dot.
    !! Examples:
    !!   path_has_extension("a.nc", "nc")  -> true
    !!   path_has_extension("a.nc", ".nc") -> true
    pure logical function path_has_extension(path, ext, match_case) result(has_ext)
        character(*), intent(in) :: path
        character(*), intent(in) :: ext
        logical, intent(in), optional :: match_case

        logical :: lmatch_case
        character(:), allocatable :: expected
        character(:), allocatable :: lhs, rhs

        lmatch_case = .false.
        if (present(match_case)) lmatch_case = match_case

        if (len_trim(ext) == 0) then
            has_ext = .false.
            return
        end if

        if (ext(1:1) == '.') then
            expected = trim(ext)
        else
            expected = '.' // trim(ext)
        end if

        if (lmatch_case) then
            has_ext = ends_with(trim(path), expected)
        else
            lhs = to_lower(trim(path))
            rhs = to_lower(expected)
            has_ext = ends_with(lhs, rhs)
        end if

    end function path_has_extension

end module path_utils