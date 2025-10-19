! Useful simple functions for statistics or sorting vectors
module stats_utils
  use precision_types, only: rk    ! Importing real64 and int64

  implicit none
  private

  public :: median_rk, sort_real_inplace

contains

    ! Median of a real vector (simple version)
    ! - Empty input -> 0.0_rk
    ! - Even n -> mean of the two middle values after sorting
    ! - Complexity: O(n log n) due to sort_real_inplace
    ! - NOTE: NaNs will propagate to the result if present.
    !------------------------------------------------------------
    pure real(rk) function median_rk(v) result(med)
        real(rk), intent(in) :: v(:)
        real(rk), allocatable :: w(:)
        integer :: n, k
        n = size(v)
        if (n == 0) then
            med = 0.0_rk; return
        end if
        allocate(w(n)); w = v
        call sort_real_inplace(w)
        k = n/2
        if (mod(n,2)==0) then
            med = 0.5_rk*(w(k) + w(k+1))
        else
            med = w(k+1)
        end if
    end function median_rk

    ! Simple function to sort an array of reals in ascending order
    ! Good only for small arrays (e.g., n < ~5e3) - Use carefully
    pure subroutine sort_real_inplace(a)
        real(rk), intent(inout) :: a(:)
        integer :: i, j
        real(rk) :: key
        do i=2, size(a)
            key = a(i); j = i-1
            do while (j>=1 .and. a(j)>key)
                a(j+1) = a(j); j = j-1
            end do
            a(j+1) = key
        end do
    end subroutine sort_real_inplace

end module stats_utils
