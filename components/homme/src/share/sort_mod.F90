#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module sort_mod

  implicit none
  private

  public :: &
       ! sortints(a): a(r,c) is sorted along a(1,:), and a(2:r,:) are carried
       ! along. An example application is to set a(2,:) = 1:c. After sortints
       ! returns, a(2,:) contains the permutation induced by the sort of a(1,:).
       sortints, &
       ! sortints3(a): a(3,c) is sorted lexicographically on a(1,:), then
       ! a(2,:), then a(3,:). If the columns are distinct (e.g., a(3,:) = 1:c),
       ! the result does not depend on the qsort implementation.
       sortints3

contains
  
  subroutine sortints(a)
    integer, intent(inout) :: a(:,:)

    call sortints_impl(a, size(a,1), size(a,2))
  end subroutine sortints

  subroutine sortints_impl(a, nrow, ncol)
    use iso_c_binding, only: c_int, c_size_t, c_ptr, c_sizeof, c_loc

    interface
       subroutine qsort(base, nmemb, size, cmp) bind(c)
         use iso_c_binding, only: c_ptr, c_size_t, c_int
         type (c_ptr), value, intent(in) :: base
         integer (kind=c_size_t), value, intent(in) :: nmemb, size
         interface
            function cmp(a, b) result(o) bind(c)
              use iso_c_binding, only: c_int
              integer (kind=c_int), intent(in) :: a, b
              integer (kind=c_int) :: o
            end function cmp
         end interface
       end subroutine qsort
    end interface

    integer (kind=c_int), target, intent(inout) :: a(nrow,ncol)
    integer, intent(in) :: nrow, ncol
    type (c_ptr) :: a_ptr
    integer (kind=c_size_t) :: r_s, c_s

    a_ptr = c_loc(a)
    r_s = nrow
    c_s = ncol
    call qsort(a_ptr, c_s, r_s*c_sizeof(c_int), cmp_ints)
  end subroutine sortints_impl

  subroutine sortints3(a)
    use parallel_mod, only: abortmp
    integer, intent(inout) :: a(:,:)

    if (size(a,1) /= 3) call abortmp('sortints3: a must have 3 rows')
    call sortints3_impl(a, size(a,2))
  end subroutine sortints3

  subroutine sortints3_impl(a, ncol)
    use iso_c_binding, only: c_int, c_size_t, c_ptr, c_sizeof, c_loc

    interface
       subroutine qsort3(base, nmemb, size, cmp) bind(c, name='qsort')
         use iso_c_binding, only: c_ptr, c_size_t, c_int
         type (c_ptr), value, intent(in) :: base
         integer (kind=c_size_t), value, intent(in) :: nmemb, size
         interface
            function cmp(a, b) result(o) bind(c)
              use iso_c_binding, only: c_int
              integer (kind=c_int), intent(in) :: a(3), b(3)
              integer (kind=c_int) :: o
            end function cmp
         end interface
       end subroutine qsort3
    end interface

    integer (kind=c_int), target, intent(inout) :: a(3,ncol)
    integer, intent(in) :: ncol
    type (c_ptr) :: a_ptr
    integer (kind=c_size_t) :: c_s

    a_ptr = c_loc(a)
    c_s = ncol
    call qsort3(a_ptr, c_s, 3*c_sizeof(a(1,1)), cmp_ints3)
  end subroutine sortints3_impl

  function cmp_ints3(a, b) result(o) bind(c)
    use iso_c_binding, only: c_int

    integer (kind=c_int), intent(in) :: a(3), b(3)
    integer (kind=c_int) :: o
    integer :: i

    o = 0
    do i = 1, 3
       if (a(i) < b(i)) then
          o = -1
          return
       elseif (a(i) > b(i)) then
          o = 1
          return
       end if
    end do
  end function cmp_ints3

  function cmp_ints(a, b) result(o) bind(c)
    use iso_c_binding, only: c_int

    integer (kind=c_int), intent(in) :: a, b
    integer (kind=c_int) :: o

    if (a < b) then
       o = -1
    elseif (a == b) then
       o = 0
    else
       o = 1
    end if
  end function cmp_ints

end module sort_mod
