
! --------------------------------------------------------------------
module func_sort
contains
! --------------------------------------------------------------------
! INTEGER FUNCTION  FindMinimum():
!    This function returns the location of the minimum in the section
! between Start and End.
! --------------------------------------------------------------------

   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      REAL(4), DIMENSION(1:), INTENT(IN) :: x
      INTEGER, INTENT(IN)                :: Start, End
      REAL(4)                            :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)		! assume the first is the min
      Location = Start			! record its position
      DO i = Start+1, End		! start with next elements
         IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location        	! return the position
   END FUNCTION  FindMinimum
end module func_sort


! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! -------------------------
! -------------------------

SUBROUTINE  rSort2_ex(x, y, z, Size)
use func_sort
  IMPLICIT  NONE
      INTEGER, INTENT(IN)                   :: Size
      REAL(4), DIMENSION(Size), INTENT(INOUT) :: x
      REAL(4), DIMENSION(Size), INTENT(INOUT) :: y
      REAL(4), DIMENSION(Size,size), INTENT(INOUT) :: z
      INTEGER                               :: i,j
      INTEGER                               :: Location

      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
         CALL  Swap(y(i), y(Location))	! swap this and the minimum
         do j=1,size
!!!            CALL  Swap(z(i,j), z(Location,j))	! swap this and the minimum
               CALL  Swap(z(j,i), z(j,Location))	! swap this and the minimum
         end do
      END DO

    END SUBROUTINE  RSort2_Ex

!----------------------------!
   SUBROUTINE  rSort2(x, y, Size)
     use func_sort
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                   :: Size
      REAL(4), DIMENSION(Size), INTENT(INOUT) :: x
      REAL(4), DIMENSION(Size), INTENT(INOUT) :: y
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
         CALL  Swap(y(i), y(Location))	! swap this and the minimum
      END DO


    END SUBROUTINE  RSort2

! -------------------------

   SUBROUTINE  rSort(x, Size)
     use func_sort
      IMPLICIT  NONE
      REAL(4), DIMENSION(1:), INTENT(INOUT) :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location

      DO i = 1, Size-1			! except for the last
         Location = FindMinimum(x, i, Size)	! find min from this to last
         CALL  Swap(x(i), x(Location))	! swap this and the minimum
      END DO

    END SUBROUTINE  RSort
   
! -------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL(4), INTENT(INOUT) :: a, b
      REAL(4)                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap
