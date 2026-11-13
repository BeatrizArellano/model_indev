! Regression test for explicit exponent formatting and legacy restart parsing.
program test_restart_csv_format
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    real(real64) :: values(6), actual
    character(len=64) :: formatted, legacy
    integer :: i, ios, output_unit

    values = [ 1.2345678901234567_real64, tiny(1.0_real64), &
               nearest(0.0_real64, 1.0_real64), 0.0_real64, &
              -123.45678901234567_real64, -nearest(0.0_real64, 1.0_real64) ]

    open(newunit=output_unit, file='restart_format_test.csv', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot create CSV fixture'
    write(output_unit,'(A)') 'name,vert_coord,index,depth,value,units,long_name'

    do i = 1, size(values)
        write(formatted,'(ES24.16E3)',iostat=ios) values(i)
        if (ios /= 0) error stop 'formatting failed'
        if (index(formatted, 'E') == 0) error stop 'missing exponent marker'
        if (index(formatted, '*') /= 0) error stop 'field overflow'
        read(formatted,*,iostat=ios) actual
        if (ios /= 0) error stop 'Fortran round-trip parse failed'
        if (.not. ieee_is_finite(actual)) error stop 'nonfinite round-trip result'
        if (actual /= values(i)) error stop 'non-identical round-trip value'
        write(output_unit,'(A,",",A,",",I0,",",ES24.16E3,",",ES24.16E3,",",A,",",A)') &
             'test', 'centre', i, real(i,real64), values(i), '1', 'format regression'
    end do
    close(output_unit)

    legacy = '1.0671817950170925-321'
    read(legacy,*,iostat=ios) actual
    if (ios /= 0) error stop 'legacy exponentless format was rejected'
    if (.not. ieee_is_finite(actual)) error stop 'legacy parse was nonfinite'
    if (actual <= 0.0_real64) error stop 'legacy subnormal value was lost'
    print '(A)', 'Fortran restart format checks passed'
end program test_restart_csv_format
