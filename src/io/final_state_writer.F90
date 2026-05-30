!> Write final model state to a CSV restart-style file.
!! The file is written in long format so variables on layer centres,
!! interfaces, surfaces, bottoms, and scalars can coexist safely.
module final_state_writer

   use precision_types,   only: rk
   use grids,             only: VerticalGrid
   use variable_registry, only: VarMetadata

   implicit none
   private

   public :: write_final_state_csv

contains

    subroutine write_final_state_csv(path, grid, vars)
        character(*),       intent(in) :: path
        type(VerticalGrid), intent(in) :: grid
        type(VarMetadata),  intent(in) :: vars(:)

        integer :: unit, ios, j
        character(len=512) :: iomsg

        open(newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios, iomsg=iomsg)

        if (ios /= 0) then
            write(*,*) 'ERROR: could not open final state file: ', trim(path)
            write(*,*) '       iostat=', ios, ' iomsg=', trim(iomsg)
            stop 1
        end if

        call write_header(unit)

        do j = 1, size(vars)
            if (.not. vars(j)%state_var) cycle

            select case (vars(j)%n_space_dims)
            case (1)
                call write_profile(unit, grid, vars(j))
            case (0)
                call write_scalar(unit, vars(j))
            case default
                write(*,*) 'ERROR: unsupported n_space_dims for final state variable: ', trim(vars(j)%name)
                stop 1
            end select
        end do

        close(unit, iostat=ios, iomsg=iomsg)
        if (ios /= 0) then
            write(*,*) 'ERROR: could not close final state file: ', trim(path)
            write(*,*) '       iostat=', ios, ' iomsg=', trim(iomsg)
            stop 1
        end if
    end subroutine write_final_state_csv


    subroutine write_profile(unit, grid, var)
        integer,           intent(in) :: unit
        type(VerticalGrid),intent(in) :: grid
        type(VarMetadata), intent(in) :: var

        integer :: k, n
        real(rk), allocatable :: zcoord(:)

        if (.not. associated(var%data_1d)) then
            write(*,*) 'ERROR: final state profile has no associated data_1d: ', trim(var%name)
            stop 1
        end if

        n = size(var%data_1d)

        select case (trim(var%vert_coord))
        case ('centre')
            if (n /= size(grid%z)) then
                write(*,*) 'ERROR: centre final state variable has wrong length: ', trim(var%name)
                write(*,*) '       size(data)=', n, ' size(grid%z)=', size(grid%z)
                stop 1
            end if
            allocate(zcoord(n))
            zcoord = grid%z

        case ('interface')
            if (n /= size(grid%z_w)) then
                write(*,*) 'ERROR: interface final state variable has wrong length: ', trim(var%name)
                write(*,*) '       size(data)=', n, ' size(grid%z_w)=', size(grid%z_w)
                stop 1
            end if
            allocate(zcoord(n))
            zcoord = grid%z_w

        case default
            write(*,*) 'ERROR: unsupported profile vert_coord in final state file: ', &
                        trim(var%name), ' vert_coord=', trim(var%vert_coord)
            stop 1
        end select

        do k = 1, n
            call write_row(unit, var, k, zcoord(k), var%data_1d(k))
        end do

        deallocate(zcoord)
    end subroutine write_profile


    subroutine write_scalar(unit, var)
        integer,           intent(in) :: unit
        type(VarMetadata), intent(in) :: var

        if (.not. associated(var%data_0d)) then
            write(*,*) 'ERROR: final state scalar has no associated data_0d: ', trim(var%name)
            stop 1
        end if

        call write_row(unit, var, 0, 0.0_rk, var%data_0d)
    end subroutine write_scalar


    subroutine write_row(unit, var, k, z, value)
        integer,           intent(in) :: unit
        type(VarMetadata), intent(in) :: var
        integer,           intent(in) :: k
        real(rk),          intent(in) :: z, value

        character(:), allocatable :: name, units, long_name, vert_coord
        integer :: ios
        character(len=512) :: iomsg

        name       = csv_escape(var%name)
        units      = csv_escape(var%units)
        long_name  = csv_escape(var%long_name)
        vert_coord = csv_escape(trim(var%vert_coord))

        write(unit,'(A,",",A,",",I0,",",ES24.16,",",ES24.16,",",A,",",A)', &
            iostat=ios, iomsg=iomsg) &
            name, vert_coord, k, z, value, units, long_name

        if (ios /= 0) then
            write(*,*) 'ERROR: failed writing final state row for variable: ', trim(var%name)
            write(*,*) '       k=', k
            write(*,*) '       iostat=', ios, ' iomsg=', trim(iomsg)
            stop 1
        end if
    end subroutine write_row


    pure function csv_escape(text) result(out)
        character(*), intent(in) :: text
        character(:), allocatable :: out

        integer :: i
        logical :: must_quote

        must_quote = index(text, ',')      > 0 .or. &
                        index(text, '"')      > 0 .or. &
                        index(text, char(10)) > 0 .or. &
                        index(text, char(13)) > 0

        if (.not. must_quote) then
            out = trim(text)
            return
        end if

        out = '"'
        do i = 1, len_trim(text)
            if (text(i:i) == '"') then
                out = out // '""'
            else
                out = out // text(i:i)
            end if
        end do
        out = out // '"'
    end function csv_escape

    subroutine write_header(unit)
        integer, intent(in) :: unit

        integer :: ios
        character(len=512) :: iomsg

        write(unit,'(A)', iostat=ios, iomsg=iomsg) &
            'name,vert_coord,depth_index,depth,value,units,long_name'

        if (ios /= 0) then
            write(*,*) 'ERROR: failed writing final state CSV header.'
            write(*,*) '       iostat=', ios, ' iomsg=', trim(iomsg)
            stop 1
        end if
    end subroutine write_header

end module final_state_writer