!=====================================================================================
!>                      NetCDF Input/Output utilities 
!! This module provides low-level wrappers for the NetCDF Fortran-90 API
!! for safely opening, creating, reading, and writing NetCDF files.  
!! Includes helpers for:
!!   - File management (open, create, close, sync, redefine)
!!   - Dimension and variable definitions
!!   - Reading and writing real-valued data (1D–3D slabs or full variables)
!!   - Attribute handling (get/put string, integer, and real attributes)
!!   - Simple metadata queries (variable type, rank, dimensions)
!! All routines call nc_check for consistent error handling.
!! Some routines and documentation were refined with assistance from ChatGPT (OpenAI).
!======================================================================================

module netcdf_io
    use precision_types, only: rk, ik
    use netcdf
    use find_utils,      only: find_name

    implicit none
    private 

    type :: NcFile
       ! Stores the state of an open Netcdf file
       integer :: ncid   = -1              ! Handle ID assigned by NetCDF when file is open
       character(:), allocatable :: path   ! path to the file on disk
    end type NcFile
    public :: NcFile
    ! Readers
    public :: nc_open, nc_close, nc_check
    public :: nc_has_var, nc_dimlen, nc_var_ndims, nc_var_dims, nc_var_type, nc_len_of_var_1d
    public :: nc_get_att_str, nc_get_att_int, nc_get_att_real
    public :: nc_read_real, nc_read_real_1d, nc_read_real_2d, nc_read_real_3d, nc_read_real_1d_slice
    public :: nc_read_calendar_info
    ! Writers
    public :: nc_open_rw, nc_create, nc_sync, nc_redef, nc_enddef
    public :: nc_def_dim, nc_def_var_real, nc_def_var_double, nc_put_att_str, nc_put_att_real, nc_put_att_int
    public :: nc_write_real, nc_write_int

    interface nc_read_real
    !> Interface for reading real-valued variables from a NetCDF file.
    !! Automatically dispatches to the appropriate dimension routine (1D, 2D, or 3D)
    !! based on the dimensionality of the target array.
    !! Example: call nc_read_real(db, 'sst', start=[t,1,1], count=[1,nlat,nlon], data_array=sst2d)
        module procedure nc_read_real_aux_1d
        module procedure nc_read_real_aux_2d
        module procedure nc_read_real_aux_3d
    end interface

    interface nc_write_real
    ! Interface for writing real variables into a NetCDF file.
    !   (Using interfaces to keep compiler-friendly)
        module procedure nc_write_real_aux_1d
        module procedure nc_write_real_aux_2d
        module procedure nc_write_real_aux_3d
    end interface

    interface nc_write_int
    !> Interface for writing integer variables into a NetCDF file.
        module procedure nc_write_int_aux_1d
        module procedure nc_write_int_aux_2d
        module procedure nc_write_int_aux_3d
    end interface

   
  contains

    subroutine nc_check(status, where)
        ! Helper to test the success of any call to a Netcdf function
        integer, intent(in) :: status
        character(*), intent(in) :: where
        if (status /= nf90_noerr) then
            write(*,'(A,1X,A,1X,A)') 'FATAL(netcdf_io):', trim(where), trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine nc_check

    subroutine nc_open(db, path)
        ! Opens a Netcdf file
        ! db is the NcFile data structure and path is the filepath string
        type(NcFile), intent(inout) :: db
        character(*), intent(in)    :: path
        call nc_close(db)
        call nc_check(nf90_open(trim(path), nf90_nowrite, db%ncid), 'open '//trim(path))
        db%path = trim(path)
    end subroutine nc_open

    subroutine nc_close(db)
        !> Close an open NetCDF file associated with the given NcFile handle.
        !! Calls nf90_close, resets the ncid to -1, and deallocates the stored path.
        type(NcFile), intent(inout) :: db
        if (db%ncid /= -1) call nc_check(nf90_close(db%ncid), 'close '//merge(db%path,'(unknown)',allocated(db%path)))
        db%ncid = -1
        if (allocated(db%path)) deallocate(db%path)
    end subroutine nc_close

    logical function nc_has_var(db, varname)
        !> Checks if a variable with the given name exists in the NetCDF file.
        !! Returns .true. if the variable is found, .false. otherwise.
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer :: vid, stat
        stat = nf90_inq_varid(db%ncid, trim(varname), vid)
        nc_has_var = (stat == nf90_noerr)
    end function nc_has_var

    integer function nc_dimlen(db, dimname)
        !> Returns the length of a named dimension in the NetCDF file.
        !! Retrieves the dimension ID and returns its size.
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: dimname
        integer :: did, len
        call nc_check(nf90_inq_dimid(db%ncid, trim(dimname), did), 'inq_dimid('//trim(dimname)//')')
        call nc_check( nf90_inquire_dimension(db%ncid, did, len=len), 'inquire_dimension('//trim(dimname)//')' )
        nc_dimlen = len
    end function nc_dimlen


    integer function nc_var_ndims(db, varname)
        !> Returns the number of dimensions of a variable in the NetCDF file.
        ! db is the NcFile data structure and varname is the variable name. 
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer :: vid, nd, xtype, dimids(nf90_max_var_dims)

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        call nc_check( nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=nd, dimids=dimids), &
                    'inquire_variable('//trim(varname)//')' )
        nc_var_ndims = nd
    end function nc_var_ndims

    ! Length of a 1D coordinate variable  
    integer function nc_len_of_var_1d(db, varname) result(n)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        character(len=:), allocatable :: dimnames(:)
        integer, allocatable :: dimlens(:)

        call nc_var_dims(db, varname, dimnames, dimlens)
        if (size(dimlens) /= 1) then
            n = -1     ! signal error to caller
            return
        end if
        n = dimlens(1)
    end function nc_len_of_var_1d

    subroutine nc_var_dims(db, varname, dimnames, dimlens)
        !> Retrieve the names and lengths of all dimensions used by a variable in the NetCDF file.
        !! @param[in]  db        NcFile handle for the open NetCDF file.
        !! @param[in]  varname   Name of the variable to inspect.
        !! @param[out] dimnames  Array of dimension names (in the variable’s order).
        !! @param[out] dimlens   Array of corresponding dimension lengths.
        
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        character(:), allocatable, intent(out) :: dimnames(:)
        integer,               allocatable, intent(out) :: dimlens(:)
        integer :: vid, nd, xtype, dimids(nf90_max_var_dims), i, did, len
        character(len=nf90_max_name) :: dname

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        
        call nc_check( nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=nd, dimids=dimids), &
                        'inquire_variable('//trim(varname)//')' )

        allocate(character(len=nf90_max_name) :: dimnames(nd))
        allocate(dimlens(nd))

        do i=1, nd
            did = dimids(i)
            call nc_check( nf90_inquire_dimension(db%ncid, did, name=dname, len=len), 'inquire_dimension' )
            dimnames(i) = trim(dname)
            dimlens(i)  = len
        end do
    end subroutine nc_var_dims


    integer function nc_var_type(db, varname)
        !> Returns the NetCDF data type (xtype) of a variable in the open file.
        !! @param[in]  db       NcFile handle for the open NetCDF file.
        !! @param[in]  varname  Name of the variable to inspect.
        !! @return     Integer NetCDF type code (e.g., nf90_float, nf90_double, nf90_int).
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer :: vid, xtype, nd, dimids(nf90_max_var_dims)
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=nd, dimids=dimids), 'inq_var('//trim(varname)//')')
        nc_var_type = xtype
    end function nc_var_type

    subroutine nc_get_att_str(db, varname, attname, value, present)
        !> Retrieve a string attribute from a NetCDF variable.
        !! @param[in]  db        NcFile handle for the open file.
        !! @param[in]  varname   Name of the variable that owns the attribute.
        !! @param[in]  attname   Name of the attribute to read.
        !! @param[out] value     Allocatable string containing the attribute value (trimmed).
        !! @param[out] present   .true. if the attribute exists and was read successfully.
        type(NcFile), intent(in)  :: db
        character(*), intent(in)  :: varname, attname
        character(:), allocatable, intent(out) :: value
        logical,      intent(out) :: present
        integer :: vid, stat, lenp           ! <<< lenp declared here
        character(len=:), allocatable :: data_array

        stat = nf90_inq_varid(db%ncid, trim(varname), vid)
        if (stat /= nf90_noerr) then
            present = .false.; return
        end if

        stat = nf90_inquire_attribute(db%ncid, vid, trim(attname), len=lenp)
        if (stat /= nf90_noerr) then
            present = .false.; return
        end if

        allocate(character(len=lenp) :: data_array)
        stat = nf90_get_att(db%ncid, vid, trim(attname), data_array)
        if (stat == nf90_noerr) then
            value = trim(data_array); present = .true.
        else
            present = .false.
        end if
    end subroutine nc_get_att_str


    subroutine nc_get_att_int(db, varname, attname, value, present)
        !> Retrieves an integer attribute from a NetCDF variable.
        !! @param[in]  db        NcFile handle for the open file.
        !! @param[in]  varname   Name of the variable that owns the attribute.
        !! @param[in]  attname   Name of the attribute to read.
        !! @param[out] value     Integer variable receiving the attribute value.
        !! @param[out] present   .true. if the attribute exists and was read successfully.

        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, attname
        integer, intent(out)     :: value
        logical, intent(out)     :: present
        integer :: vid, stat
        stat = nf90_inq_varid(db%ncid, trim(varname), vid)
        if (stat /= nf90_noerr) then; present=.false.; return; end if
        stat = nf90_get_att(db%ncid, vid, trim(attname), value)
        present = (stat == nf90_noerr)
    end subroutine nc_get_att_int

    subroutine nc_get_att_real(db, varname, attname, value, present)
        !> Retrieves a real attribute from a NetCDF variable.
        !! @param[in]  db        NcFile handle for the open file.
        !! @param[in]  varname   Name of the variable that owns the attribute.
        !! @param[in]  attname   Name of the attribute to read.
        !! @param[out] value     Real variable receiving the attribute value (kind=rk).
        !! @param[out] present   .true. if the attribute exists and was read successfully.

        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, attname
        real(rk), intent(out)    :: value
        logical, intent(out)     :: present
        integer :: vid, stat
        stat = nf90_inq_varid(db%ncid, trim(varname), vid)
        if (stat /= nf90_noerr) then; present=.false.; return; end if
        stat = nf90_get_att(db%ncid, vid, trim(attname), value)
        present = (stat == nf90_noerr)
    end subroutine nc_get_att_real

    ! ---------- Calendar/units helper ----------
    subroutine nc_read_calendar_info(db, time_var, calendar, units, present_calendar, present_units)
        !> Retrieve calendar and units attributes from a time coordinate variable in a NetCDF file.
        !! @param[in]  db                NcFile handle for the open NetCDF file.
        !! @param[in]  time_var          Name of the time coordinate variable (e.g., "time").
        !! @param[out] calendar          Allocatable string receiving the calendar attribute value.
        !! @param[out] units             Allocatable string receiving the units attribute value.
        !! @param[out] present_calendar  .true. if the "calendar" attribute exists and was read successfully.
        !! @param[out] present_units     .true. if the "units" attribute exists and was read successfully.
        !! Example: call nc_read_calendar_info(db, 'time', cal, units, has_cal, has_units)
        type(NcFile), intent(in)           :: db
        character(*), intent(in)           :: time_var   !! name of your time coordinate variable
        character(:), allocatable, intent(out) :: calendar, units
        logical, intent(out)               :: present_calendar, present_units
        call nc_get_att_str(db, trim(time_var), 'calendar', calendar, present_calendar)
        call nc_get_att_str(db, trim(time_var), 'units',    units,    present_units)
    end subroutine nc_read_calendar_info


  ! ---------- Core reads (real) -----------------------------------------------------
  !-------------------------------------------------------------------------------------------
    subroutine nc_read_real_1d(db, varname, data_vec)
        !> Read an entire one-dimensional real variable from an open NetCDF file.
        !! @param[in]  db        NcFile handle for the open file.
        !! @param[in]  varname   Name of the variable to read.
        !! @param[out] data_vec  Output 1-D real vector containing all values of the variable.
        !! Example: call nc_read_real_1d(db, 'latitude', lat)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        real(rk),     intent(out):: data_vec(:)

        integer :: vid, xtype, ndims, dimids(NF90_MAX_VAR_DIMS), natts

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')

        call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, &
                                            dimids=dimids, nAtts=natts), 'inquire_variable('//trim(varname)//')')

        if (ndims /= 1) then
            call nc_check(NF90_EINVALCOORDS, 'read_real_1d: variable '//trim(varname)//' must be 1-D')
        end if

        ! Treat anything except CHAR as numeric; this keeps us portable across NetCDF-Fortran versions.
        if (xtype == NF90_CHAR) then
            call nc_check(NF90_EBADTYPE, 'read_real_1d: variable '//trim(varname)//' is character, not numeric')
        end if
        ! NetCDF will auto-convert integer/float types to real(rk).
        call nc_check(nf90_get_var(db%ncid, vid, data_vec), 'get_var('//trim(varname)//')')
    end subroutine nc_read_real_1d

    subroutine nc_read_real_2d(db, varname, data_array)
        !> Read an entire two-dimensional real variable from an open NetCDF file.
        !! @param[in]  db          NcFile handle for the open file.
        !! @param[in]  varname     Name of the variable to read.
        !! @param[out] data_array  Output 2-D real array containing all values of the variable.
        !! Example: call nc_read_real_2d(db, 'mask', mask)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        real(rk), intent(out)    :: data_array(:,:)
        integer :: vid, xtype, ndims, dimids(NF90_MAX_VAR_DIMS), natts
        integer :: j, lenj
        integer :: shp(2)
        character(len=64) :: want_s, got_s

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                        'inquire_variable('//trim(varname)//')')

        if (ndims /= 2) call nc_check(NF90_EINVALCOORDS, 'read_real_2d: variable '//trim(varname)//' must be 2-D')
        if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'read_real_2d: variable '//trim(varname)//' is character')

        shp = shape(data_array)
        do j = 1, 2
            call nc_check(nf90_inquire_dimension(db%ncid, dimids(j), len=lenj), 'inquire_dimension('//trim(varname)//')')
            if (lenj /= shp(j)) then
                write(want_s,'(I0,",",I0)') shp(1), shp(2)
                write(got_s ,'(I0,",",I0)') merge(lenj, -1, j==1), merge(lenj, -1, j==2)  ! quick hint; or query both lens
                call nc_check(NF90_EEDGE, 'read_real_2d: size mismatch for '//trim(varname)//' want('//trim(want_s)// &
                                        ') vs file dims(~'//trim(got_s)//')')
            end if
        end do

        call nc_check(nf90_get_var(db%ncid, vid, data_array), 'get_var('//trim(varname)//')')
    end subroutine nc_read_real_2d

    subroutine nc_read_real_3d(db, varname, data_array)
        !> Read an entire three-dimensional real variable from an open NetCDF file.
        !! @param[in]  db          NcFile handle for the open file.
        !! @param[in]  varname     Name of the variable to read.
        !! @param[out] data_array  Output 3-D real array containing all values of the variable.
        !! Example: call nc_read_real_3d(db, 'temperature', temp3d)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        real(rk), intent(out)    :: data_array(:,:,:)
        integer :: vid, xtype, ndims, dimids(NF90_MAX_VAR_DIMS), natts
        integer :: j, lenj
        integer :: shp(3)
        character(len=96) :: want_s, got_s

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                        'inquire_variable('//trim(varname)//')')

        if (ndims /= 3) call nc_check(NF90_EINVALCOORDS, 'read_real_3d: variable '//trim(varname)//' must be 3-D')
        if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'read_real_3d: variable '//trim(varname)//' is character')

        shp = shape(data_array)
        do j = 1, 3
            call nc_check(nf90_inquire_dimension(db%ncid, dimids(j), len=lenj), 'inquire_dimension('//trim(varname)//')')
            if (lenj /= shp(j)) then
                write(want_s,'(I0,",",I0,",",I0)') shp(1), shp(2), shp(3)
                write(got_s ,'(I0,",",I0,",",I0)') merge(lenj, -1, j==1), merge(lenj, -1, j==2), merge(lenj, -1, j==3)
                call nc_check(NF90_EEDGE, 'read_real_3d: size mismatch for '//trim(varname)//' want('//trim(want_s)// &
                                        ') vs file dims(~'//trim(got_s)//')')
            end if
        end do

        call nc_check(nf90_get_var(db%ncid, vid, data_array), 'get_var('//trim(varname)//')')
    end subroutine nc_read_real_3d

    !--------Auxiliar subroutines used when nc_read_real is called -------------

    subroutine nc_read_real_aux_1d(db, varname, start, count, step, data_array)
        !> Reads a one-dimensional real variable from an open NetCDF file.
        !! @param[in]  db       NcFile handle for the open file.
        !! @param[in]  varname  Name of the variable to read.
        !! @param[in]  start    Starting index in each dimension (size = 1).
        !! @param[in]  count    Number of elements to read (size = 1).
        !! @param[in]  step   Optional step between elements (size = 1).
        !! @param[out] data_array      Output real vector receiving the data.
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(out)    :: data_array(:)
        integer :: vid
        if (size(start)/=1 .or. size(count)/=1) stop 'nc_read_real(1D): start/count size must be 1'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=1) stop 'nc_read_real(1D): step size must be 1'
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'get_var('//trim(varname)//')')
        else
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count), 'get_var('//trim(varname)//')')
        end if
    end subroutine

    subroutine nc_read_real_aux_2d(db, varname, start, count, step, data_array)
        !> Reads a two-dimensional real variable from an open NetCDF file.
        !! @param[in]  db       NcFile handle for the open file.
        !! @param[in]  varname  Name of the variable to read.
        !! @param[in]  start    Starting index in each dimension (size = 1).
        !! @param[in]  count    Number of elements to read (size = 1).
        !! @param[in]  step   Optional step between elements (size = 1).
        !! @param[out] data_array      Output real vector receiving the data.
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(out)    :: data_array(:,:)
        integer :: vid
        if (size(start)/=2 .or. size(count)/=2) stop 'nc_read_real(2D): start/count size must be 2'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=2) stop 'nc_read_real(2D): step size must be 2'
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'get_var('//trim(varname)//')')
        else
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count), 'get_var('//trim(varname)//')')
        end if
    end subroutine

    subroutine nc_read_real_aux_3d(db, varname, start, count, step, data_array)
        !> Reads a 3-dimensional real variable from an open NetCDF file.
        !! @param[in]  db       NcFile handle for the open file.
        !! @param[in]  varname  Name of the variable to read.
        !! @param[in]  start    Starting index in each dimension (size = 1).
        !! @param[in]  count    Number of elements to read (size = 1).
        !! @param[in]  step   Optional step between elements (size = 1).
        !! @param[out] data_array      Output real vector receiving the data.
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(out)    :: data_array(:,:,:)
        integer :: vid
        if (size(start)/=3 .or. size(count)/=3) stop 'nc_read_real(3D): start/count size must be 3'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=3) stop 'nc_read_real(3D): step size must be 3'
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'get_var('//trim(varname)//')')
        else
            call nc_check(nf90_get_var(db%ncid, vid, data_array, start=start, count=count), 'get_var('//trim(varname)//')')
        end if
    end subroutine  


    ! ------------- Read a slice of a 1D array
    subroutine nc_read_real_1d_slice(db, varname, dim_name, i0, i1, out)
        ! Read a contiguous slice from a 1-D numeric variable along its sole dimension.
        ! Requires that the variable is 1-D and its dimension name matches dim_name.        
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, dim_name
        integer,      intent(in) :: i0, i1
        real(rk),     intent(out):: out(:)

        integer :: vid, ndims, dimids(NF90_MAX_VAR_DIMS), xtype, natts
        character(len=:), allocatable :: dnames(:)
        integer, allocatable :: dlens(:)
        integer :: len, start, count

        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid),        'inq_varid('//trim(varname)//')')
        call nc_check(nf90_inquire_variable(db%ncid, vid, xtype=xtype, ndims=ndims, dimids=dimids, nAtts=natts), &
                        'inquire_variable('//trim(varname)//')')
        if (ndims /= 1) call nc_check(NF90_EINVALCOORDS, 'read_1d_slice: '//trim(varname)//' must be 1-D')
        if (xtype == NF90_CHAR) call nc_check(NF90_EBADTYPE, 'read_1d_slice: '//trim(varname)//' is character')

        call nc_var_dims(db, varname, dnames, dlens)
        if (find_name(dnames, trim(dim_name)) /= 1) then
            call nc_check(NF90_EBADDIM, 'read_1d_slice: dim name mismatch for '//trim(varname))
        end if

        len = dlens(1)
        if (i0 < 1 .or. i1 < i0 .or. i1 > len) call nc_check(NF90_EEDGE, 'read_1d_slice: indices out of range')
        if (size(out) /= (i1 - i0 + 1))         call nc_check(NF90_EEDGE, 'read_1d_slice: output size mismatch')

        start = i0
        count = i1 - i0 + 1

        call nc_check(nf90_get_var(db%ncid, vid, out, start=(/start/), count=(/count/)), 'get_var slice '//trim(varname))
    end subroutine nc_read_real_1d_slice



  !-------------------------------------------------------------------------------------------------
  ! ---------------------------------- Netcdf writers ----------------------------------------------

    subroutine nc_create(db, path, overwrite)
        !> Create a new NetCDF file on disk.
        !! Closes any previously open file in the NcFile handle and creates a new NetCDF-4 file,
        !! optionally overwriting existing files if overwrite=.true. (default) or preserving them if .false.
        !! @param[inout] db       NcFile handle to initialize or reuse.
        !! @param[in]    path     Path for the new NetCDF file.
        !! @param[in]    overwrite  Optional flag: .true. to overwrite existing file, .false. to protect it.
        !! Example: call nc_create(db, 'output.nc', overwrite=.true.)

        type(NcFile), intent(inout) :: db
        character(*), intent(in)    :: path
        logical,      intent(in), optional :: overwrite
        integer :: mode

        call nc_close(db)
        mode = nf90_netcdf4
        if (present(overwrite)) then
            if (overwrite) then
                mode = ior(mode, nf90_clobber)      ! overwrite existing file
            else
                mode = ior(mode, nf90_noclobber)    ! fail if file exists
            end if
        else
            mode = ior(mode, nf90_clobber)        ! default: overwrite
        end if
        call nc_check(nf90_create(trim(path), mode, db%ncid), 'create '//trim(path))
        db%path = trim(path)
    end subroutine nc_create

    subroutine nc_open_rw(db, path)
        !> Open an existing NetCDF file in read–write mode.
        !! Closes any previously open file in the NcFile handle, then reopens the new one in write mode.
        !! @param[inout] db    NcFile handle to initialize or reuse.
        !! @param[in]    path  Path to the existing NetCDF file to open for modification.
        !! Example: call nc_open_rw(db, 'results.nc')
        type(NcFile), intent(inout) :: db
        character(*), intent(in)    :: path
        call nc_close(db)
        call nc_check(nf90_open(trim(path), nf90_write, db%ncid), 'open_rw '//trim(path))
        db%path = trim(path)
    end subroutine nc_open_rw

    subroutine nc_sync(db)
        !> Flushes all data_arrayfered changes for an open NetCDF file to disk.
        !! Ensures that any data or metadata written so far are synchronized with the physical file.
        !! @param[in] db  NcFile handle for an open NetCDF file in write mode.
        !! Example: call nc_sync(db)
        type(NcFile), intent(in) :: db
        call nc_check(nf90_sync(db%ncid), 'sync '//merge(db%path,'(unknown)',allocated(db%path)))
    end subroutine nc_sync

    subroutine nc_redef(db)
        type(NcFile), intent(in) :: db
        integer :: status

        status = nf90_redef(db%ncid)

        if (status == NF90_EINDEFINE) then
            ! Already in define mode – nothing to do, treat as success
            return
        end if

        call nc_check(status, 'redef')
    end subroutine nc_redef


    !> Exit NetCDF redefinition mode and resume data writing.
    !! Completes structural changes started after nc_redef and finalizes the file’s header definition.
    !! @param[in] db  NcFile handle for the open NetCDF file.
    !! Example: call nc_enddef(db)
    subroutine nc_enddef(db)     
        type(NcFile), intent(in) :: db
        call nc_check(nf90_enddef(db%ncid), 'enddef')
    end subroutine nc_enddef

    !------------ Define dimension and variables-----------------------------------
    !------------------------------------------------------------------------------

    !> Defines a new dimension in the open NetCDF file.
    !! Using NF90_UNLIMITED for an editable dimension.
    !! @param[in]  db      NcFile handle for the open file (in define mode).
    !! @param[in]  name    Name of the dimension to create.
    !! @param[in]  length  Length of the dimension (or NF90_UNLIMITED).
    !! @return     dimid   The NetCDF dimension ID of the newly defined dimension.
    !! Example: dim_t = nc_def_dim(db, 'time', NF90_UNLIMITED)
    integer function nc_def_dim(db, name, length) result(dimid)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: name
        integer,      intent(in) :: length  ! use NF90_UNLIMITED for unlimited
        call nc_check(nf90_def_dim(db%ncid, trim(name), length, dimid), 'def_dim('//trim(name)//')')
    end function nc_def_dim


    !> Defines a new real-valued variable with the given dimensions.
    !! @param[in]  db           NcFile handle for the open file (in define mode).
    !! @param[in]  name         Variable name.
    !! @param[in]  dimids       Dimension IDs in the variable’s storage order (fastest-varying last).
    !! @param[in]  contiguous   Optional flag: if .true., set contiguous layout (no chunking).
    !! @return     varid        The NetCDF variable ID of the newly defined variable.
    !! Example: var_sst = nc_def_var_real(db, 'sst', [dim_t, dim_lat, dim_lon])
    integer function nc_def_var_real(db, name, dimids, contiguous) result(varid)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: name
        integer,      intent(in) :: dimids(:)
        logical,      intent(in), optional :: contiguous
        integer, allocatable :: cs(:)

        call nc_check(nf90_def_var(db%ncid, trim(name), nf90_real, dimids, varid), &
                        'def_var('//trim(name)//')')

        if (present(contiguous)) then
            if (contiguous) then
            ! Fortran requires a non-empty array even though value is ignored for CONTIGUOUS.
            allocate(cs(1)); cs = 1
            call nc_check(nf90_def_var_chunking(db%ncid, varid, nf90_contiguous, cs), 'chunking(contig)')
            deallocate(cs)
            end if
        end if
    end function nc_def_var_real

    ! For real64 data
    integer function nc_def_var_double(db, name, dimids) result(varid)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: name
        integer, intent(in)      :: dimids(:)
        call nc_check(nf90_def_var(db%ncid, trim(name), nf90_double, dimids, varid), &
                        'def_var_double('//trim(name)//')')
    end function nc_def_var_double


    !> Defines a new integer variable with the given dimensions.
    !! @param[in]  db         NcFile handle for the open file (in define mode).
    !! @param[in]  name       Variable name.
    !! @param[in]  dimids     Dimension IDs for the variable (fastest-varying last).
    !! @return     varid      NetCDF variable ID of the newly defined integer variable.
    !! Example: var_days = nc_def_var_int(db, 'days_in_year', [dim_t])
    integer function nc_def_var_int(db, name, dimids) result(varid)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: name
        integer,      intent(in) :: dimids(:)
        call nc_check(nf90_def_var(db%ncid, trim(name), nf90_int, dimids, varid), &
                        'def_var('//trim(name)//')')
    end function nc_def_var_int



    ! --- Attribute setters (file-level: varname='' or global) --------------


    !> Writes a string attribute to a variable or as a global attribute in the open NetCDF file.
    !! Automatically distinguishes between variable-specific and global attributes via the optional flag.
    !! @param[in]  db       NcFile handle for the open file (in define mode or data mode).
    !! @param[in]  varname  Name of the variable to attach the attribute to (ignored if global=.true.).
    !! @param[in]  attname  Name of the attribute to write.
    !! @param[in]  value    Character value of the attribute.
    !! @param[in]  global   Optional flag: .true. to write a global attribute instead of a variable attribute.
    !! Example: 
    !!   call nc_put_att_str(db, 'time', 'units', 'days since 1900-01-01')
    !!   call nc_put_att_str(db, '', 'title', 'ShelfModel output', global=.true.)
    subroutine nc_put_att_str(db, varname, attname, value, global)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, attname, value
        logical,      intent(in), optional :: global
        integer :: vid
        if (present(global) .and. global) then
            vid = nf90_global
        else
            call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        end if
        call nc_check(nf90_put_att(db%ncid, vid, trim(attname), trim(value)), 'put_att_str('//trim(attname)//')')
    end subroutine nc_put_att_str

    !> Writes a real attribute to a variable or as a global attribute in the open NetCDF file.
    !! Automatically distinguishes between variable-specific and global attributes via the optional flag.
    subroutine nc_put_att_real(db, varname, attname, value, global)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, attname
        real(rk),     intent(in) :: value
        logical,      intent(in), optional :: global
        integer :: vid
        if (present(global) .and. global) then
            vid = nf90_global
        else
            call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        end if
        call nc_check(nf90_put_att(db%ncid, vid, trim(attname), value), 'put_att_real('//trim(attname)//')')
    end subroutine nc_put_att_real

    !> Writes an integer attribute to a variable or as a global attribute in the open NetCDF file.
    !! Automatically distinguishes between variable-specific and global attributes via the optional flag.
    subroutine nc_put_att_int(db, varname, attname, value, global)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname, attname
        integer,      intent(in) :: value
        logical,      intent(in), optional :: global
        integer :: vid
        if (present(global) .and. global) then
            vid = nf90_global
        else
            call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        end if
        call nc_check(nf90_put_att(db%ncid, vid, trim(attname), value), 'put_att_int('//trim(attname)//')')
    end subroutine nc_put_att_int

    !------------------------------ Data writers -------------------------------------------------------

    !> Write a 1-D real array slice to a NetCDF variable.
    !! @param[in] db      NcFile handle (open for write)
    !! @param[in] varname Variable name
    !! @param[in] start   Start indices (size=1)
    !! @param[in] count   Counts (size=1)
    !! @param[in] step    Optional step/stride (size=1, defaults to 1)
    !! @param[in] data_array     Real vector with data to write
    subroutine nc_write_real_aux_1d(db, varname, start, count, step, data_array)    
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(in)     :: data_array(:)
        integer :: vid
        if (size(start)/=1 .or. size(count)/=1) stop 'nc_write_real(1D): start/count size must be 1'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=1 .or. any(step<1)) stop 'nc_write_real(1D): step size=1, values>=1'
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

    !> Write a 2-D real array slice to a NetCDF variable.
    subroutine nc_write_real_aux_2d(db, varname, start, count, step, data_array)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(in)     :: data_array(:,:)
        integer :: vid
        if (size(start)/=2 .or. size(count)/=2) stop 'nc_write_real(2D): start/count size must be 2'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=2 .or. any(step<1)) stop 'nc_write_real(2D): step size=2, values>=1'
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

    !> Write a 3-D real array slice into a NetCDF variable (e.g., time × lat × lon).
    subroutine nc_write_real_aux_3d(db, varname, start, count, step, data_array)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        real(rk), intent(in)     :: data_array(:,:,:)
        integer :: vid
        if (size(start)/=3 .or. size(count)/=3) stop 'nc_write_real(3D): start/count size must be 3'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=3 .or. any(step<1)) stop 'nc_write_real(3D): step size=3, values>=1'
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

    !-------- Integer Data -----------

    !> Write a 1-D integer array slice to a NetCDF variable.
    !! called by the interface nc_write_int
    !! @param[in]  db      NcFile handle (open for write)
    !! @param[in]  varname Variable name
    !! @param[in]  start   Start indices (size=1)
    !! @param[in]  count   Counts (size=1)
    !! @param[in]  step    Optional step between elements (defaults to 1)
    !! @param[in]  data_array     Integer vector containing data to write
    subroutine nc_write_int_aux_1d(db, varname, start, count, step, data_array)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        integer, intent(in)      :: data_array(:)
        integer :: vid
        if (size(start)/=1 .or. size(count)/=1) stop 'nc_write_int(1D): start/count size must be 1'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

    !> Write a 2-D integer array slice to a NetCDF variable.
    ! called by the interface nc_write_int
    subroutine nc_write_int_aux_2d(db, varname, start, count, step, data_array)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        integer, intent(in)      :: data_array(:,:)
        integer :: vid
        if (size(start)/=2 .or. size(count)/=2) stop 'nc_write_int(2D): start/count size must be 2'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=2 .or. any(step<1)) stop 'nc_write_real(2D): step size=2, values>=1'
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

    !> Write a 3-D integer array slice to a NetCDF variable.
    ! called by the interface nc_write_int
    subroutine nc_write_int_aux_3d(db, varname, start, count, step, data_array)
        type(NcFile), intent(in) :: db
        character(*), intent(in) :: varname
        integer, intent(in)      :: start(:), count(:)
        integer, intent(in), optional :: step(:)
        integer, intent(in)      :: data_array(:,:,:)
        integer :: vid
        if (size(start)/=3 .or. size(count)/=3) stop 'nc_write_int(3D): start/count size must be 3'
        call nc_check(nf90_inq_varid(db%ncid, trim(varname), vid), 'inq_varid('//trim(varname)//')')
        if (present(step)) then
            if (size(step)/=3 .or. any(step<1)) stop 'nc_write_int(3D): step size=3, values>=1'
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count, stride=step), 'put_var('//trim(varname)//')')
        else
            call nc_check(nf90_put_var(db%ncid, vid, data_array, start=start, count=count), 'put_var('//trim(varname)//')')
        end if
    end subroutine

end module netcdf_io
