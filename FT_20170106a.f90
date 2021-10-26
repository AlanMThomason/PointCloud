
subroutine print_matrix(matrix, matrix_d)
implicit none
integer, intent(inout) :: matrix_d
real(8), intent(inout), dimension(0:matrix_d - 1, 0:matrix_d ) :: matrix
!f2py intent(in, out) :: matrix

integer :: i, j

DO i = 0, matrix_d - 1
  DO j = 0, matrix_d
    write (*, 100, ADVANCE = 'NO')  matrix(i, j)
    100 FORMAT(F9.4)
  ENDDO
  print *, "  "
ENDDO

print *, "  "
end subroutine

subroutine gauss(matrix, matrix_d)
implicit none
integer, intent(inout) :: matrix_d
real(8), intent(inout), dimension(0:matrix_d - 1, 0:matrix_d) :: matrix
!f2py intent(in, out) :: matrix
real(8), dimension(0:matrix_d - 1) :: coeffs

!Written by Alan Michel Thomason, Plymouth Machine Integration, lLC
! 2014-08-28
! Based on information in Advanced Engineering Mathematics, 8th ed.
! Erwin Kreyszig
! special attention to p. 888

!This Gausian elimation program was specifically written to solve a 7th degree polynomial equation bridging
!  two pieces of camshaft profile with know s, v, a, and j.
!The coefficients should be available in an 8 by 8 matrix (i=7 x j=7), augmented to the right with the knowns in
! column 9 (j=8) and the row number in column 10 (j=9).

real(8) :: coeff, holding_cell, multiplier, normalizer, exponent_var, temp
integer :: loop_watchdog, loop_watchdog_max = 1000
real(8) :: max_thus_far
integer :: row_with_max
integer :: row, h, i, i_base, j, k, m, z
integer :: skip

!integer :: loop_watchdog_max = 1000
!integer :: k = 0

k = 0
skip = 0

!print *, "Gauss, Step 1"
!call print_matrix(matrix, matrix_d)
!goto 900

!Go through and reorder the matrix.
!print *, "Reorder matrix"
!WRITE(*,*) "Reorder matrix"

z = 0
i_base = 0
loop_watchdog = 0

100 loop_watchdog = loop_watchdog + 1
IF (loop_watchdog .GT. loop_watchdog_max) THEN
  print *, "Watchdog Value exceeded, exiting prematurely"
  GOTO 900 
ENDIF


max_thus_far = 0  ! Absolute value
row_with_max = i_base
!FOR i = i_base To matrix_d  ! We only need to work with the rows that have not yet been solved.
DO i = i_base, matrix_d - 1
  IF (ABS(matrix(i, z)) .GT. max_thus_far) THEN
    max_thus_far = ABS(matrix(i, z))
    row_with_max = i
  ENDIF
ENDDO

!print *, "max_thus_far ", max_thus_far

IF (row_with_max .GT. i_base) THEN
  DO j = 0, matrix_d 
    holding_cell = matrix(i_base, j)
    matrix(i_base, j) = matrix(row_with_max, j)
    matrix(row_with_max, j) = holding_cell
  ENDDO
ENDIF





i_base = i_base + 1
IF (i_base .LE. matrix_d-1) THEN 
  GOTO 100
ENDIF    

!print *, "Gauss Step 2, Reordered"
!print *, "i_base = ", i_base
!call print_matrix(matrix, matrix_d)

!CALL print_matrix(matrix, matrix_d)


IF (matrix(z, z) .NE. 0) THEN
  Multiplier = 1 / matrix(z, z)
ELSE
  print *, "Ending Matrix solver early because a leading coefficient is zero.  z = ", z
  GOTO 900
ENDIF

DO j = z, matrix_d
!FOR j = z TO matrix_d OR j = 8 Then matrix(z, j) = multiplier * matrix(z, j)
  IF (j .LE. matrix_d) THEN 
    matrix(z, j) = multiplier * matrix(z, j)
  ENDIF    
ENDDO

!print *, "Multiplier = ", multiplier
!call print_matrix(matrix, matrix_d)

!goto 900

Normalizer = 1

DO i = z + 1, matrix_d - 1
  IF (matrix(i, z) .NE. 0) THEN
    Multiplier = -1 / matrix(i, z)
    DO j = z,8
      IF (j .LE. matrix_d) THEN
        temp = (multiplier * matrix(i, j))
        matrix(i, j) = temp + matrix(z, j)
      ENDIF
    ENDDO
  ENDIF
  !print *, "z = ", z, " Multiplier = ", multiplier
  !call print_matrix(matrix, matrix_d)
ENDDO



!goto 900

z = z + 1
i_base = z
IF (z < matrix_d) THEN 
  GOTO 100
ENDIF    

!Part of Gaussian elimination is to reorder the rows of the matrix from greatest to least with regards to the
!  particular coefficient being evaluated.
!Now go back up through substituting values found into the rows one by one to get a diagonal of 1

!print *, "Cancelling"
DO i = matrix_d -1, 1, -1
  DO k = i, matrix_d -1, 1
    DO m = 1, 2, 1    
      multiplier = matrix(i - 1, k)
      DO j = 0, matrix_d
        IF (j .LE. matrix_d) THEN 
          matrix(i - 1, j) = matrix(i - 1, j) - multiplier * matrix(k, j)
        ENDIF        
        !call print_matrix(matrix, matrix_d)
      ENDDO
    ENDDO   ! m is here simply to get around a floating point problem.
  ENDDO
ENDDO

DO i = 0, matrix_d-1
  IF (i < matrix_d) THEN
    coeffs(i) = matrix(i, matrix_d+1)
  ELSE
    coeffs(i) = 0
  ENDIF
ENDDO

!print *, "Final Result"
!call print_matrix(matrix, matrix_d)


900 i = 0
!print *, "End"
END SUBROUTINE

subroutine fitplane(p,q,r,g, npoints)
implicit none
integer*4, intent(inout) :: npoints
real(8), intent(inout), dimension(0:npoints-1) :: p, q, r
!f2py intent(in, out) :: p, q, r
real(8), intent(inout), dimension(0:3) :: g
!f2py intent(in, out) :: g

real(8) :: r_sqd, r_bar, ss_res, ss_tot
real(8), dimension(0:2, 0:3) :: D
real(8), dimension(0:2) :: F
real(8), dimension(0:npoints - 1) :: res, tot
integer :: i, j, d_dim
!G contains the results is the three dimensions that form a plane and the r^2 value at the end.
! Create the necessary matrices to find the plane representing all of the points.

!write (*, "Test")
!print *, "Fitplane 01"
d_dim = 3
IF (npoints .EQ. 0) THEN
    print *, "No points to evaluate"        
        G = [0.0,0.0,0.0,0.0]    
        GOTO 900
ENDIF

!print *, "Fitplane 02"
DO i = 0, d_dim
  D(i,:) = 0.0
ENDDO    

DO i = 0, d_dim-1
  F(i) = 0.0
ENDDO 


!print *, "Fitplane 03"
D(0,0) = sum(p * p)
D(1,0) = sum(p * q)
D(2,0) = sum (p)

D(0,1) = sum(p * q)
D(1,1) = sum(q * q)
D(2,1) = sum(q)

D(0,2) = sum(p)
D(1,2) = sum(q)
D(2,2) = npoints

D(0,3) = sum(p*r)
D(1,3) = sum(q*r)
D(2,3) = sum(r)

!print *, "Fitplane 04"
call gauss(D, d_dim)
!print *, "Fitplane 05"

DO i = 0,2
  G(i) = D(i,3)
ENDDO    
!print *, "Fitplane 06"
    
!Determine R values    
    r_bar = sum(r) / npoints
    tot = r - r_bar
    tot = tot * tot
    res = (G(0)*p + G(1)*q + G(2) - r)
    res = res * res
    ss_res = sum(res)/npoints
    ss_tot = sum(tot)/npoints

    r_sqd = (ss_res/ss_tot)
    r_sqd = 1 - r_sqd
    
    G(3) = r_sqd
    !print *, "Fitplane 07"

900 i = 0
!print *, "End"
END SUBROUTINE

subroutine crunch(big, little, mask, bigd, littled)
implicit none
integer, intent(inout) :: bigd, littled
real(8), intent(inout), dimension(0:bigd-1) :: big
!f2py intent(in, out) :: big
real(8), intent(inout), dimension(0:littled-1) :: little
!f2py intent(in, out) :: Little
logical, intent(inout), dimension(0:bigd-1) :: mask
!f2py intent(in, out) :: Mask

! This subroutine takes all of the elements out of a 1-d array in which the 
!  boolean array (of the same size as array big) is true and puts them into
!  a smaller array little.

integer*4 :: i, j

    j = 0
    do i = 0, bigd-1 
        if (mask(i) .eqv. .TRUE.) then        
            little(j) = big(i)  
            j = j + 1
        end if
    end do    

end subroutine

subroutine crunch_all(big1, big2, big3, little1, little2, little3, mask, bigd, littled)
implicit none
integer, intent(inout) :: bigd, littled
real(8), intent(inout), dimension(0:bigd-1) :: big1, big2, big3
!f2py intent(in, out) :: big1, big2, big3
real(8), intent(inout), dimension(0:littled-1) :: little1, little2, little3
!f2py intent(in, out) :: little1, little2, little3
logical, intent(inout), dimension(0:bigd-1) :: mask
!f2py intent(in, out) :: Mask

! This subroutine takes all of the elements out of a 1-d array in which the 
!  boolean array (of the same size as array big) is true and puts them into
!  a smaller array little.

integer*4 :: i, j

    j = 0
    do i = 0, bigd-1 
        if (mask(i) .eqv. .TRUE.) then        
            little1(j) = big1(i)  
            little2(j) = big2(i) 
            little3(j) = big3(i) 
            j = j + 1
        end if
    end do    

end subroutine

subroutine crunch2(big1, big2, big3, little1, little2, little3, mask, bigd, littled)

use omp_lib
implicit none
integer, intent(inout) :: bigd, littled
real(8), intent(inout), dimension(0:bigd-1) :: big1, big2, big3
!f2py intent(in, out) :: big1, big2, big3
real(8), intent(inout), dimension(0:littled-1) :: little1, little2, little3
!f2py intent(in, out) :: little1, little2, little3
logical, intent(inout), dimension(0:bigd-1) :: mask
!f2py intent(in, out) :: Mask

! This subroutine takes all of the elements out of a 1-d array in which the 
!  boolean array (of the same size as array big) is true and puts them into
!  a smaller array little.

integer*4 :: i, j

! From the following website:
! https://people.sc.fsu.edu/~jburkardt/f_src/multitask_openmp/multitask_openmp.f90

!$omp parallel 

    !$omp sections

    !$omp section
    call crunch ( big1, little1, mask, bigd, littled )
    !$omp section
    call crunch ( big2, little2, mask, bigd, littled )
    !$omp section
    call crunch ( big3, little3, mask, bigd, littled )
    !$omp end sections

!$omp end parallel

end subroutine

subroutine reducearraysize(xbig, ybig, zbig, xlittle, ylittle, zlittle, bigd, littled)
implicit none
integer, intent(inout) :: bigd, littled
real(8), intent(inout), dimension(0:bigd-1) :: xbig, ybig, zbig 
!f2py intent(in, out) :: xbig, ybig, zbig
real(8), intent(inout), dimension(0:littled-1) :: xlittle, ylittle, zlittle
!f2py intent(in, out) :: xlittle, ylittle, zlittle

integer*4 :: i, j, factor
factor = int(bigd / littled)
write(*,*) "Factor = ", factor

    i = 0
    j = 0
    do while (i < bigd-1)
        xlittle(j) = xbig(i)
        ylittle(j) = ybig(i)
        zlittle(j) = zbig(i) 
        i = i + factor
        j = j + 1
    end do    

end subroutine

subroutine setonorigin(Xpos, Ypos, Zpos, Dim)
implicit none
integer*4, intent(in) :: Dim
real(8), intent(inout), dimension(0:Dim - 1) :: Xpos, Ypos, Zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos

    Xpos = Xpos - MINVAL(Xpos)
    Ypos = Ypos - MINVAL(Ypos)
    Zpos = Zpos - MINVAL(Zpos)
   
end subroutine

subroutine setonpoint2(p, q, r, xpos, ypos, zpos, index_pt, dimofarray)
implicit none
integer*4, intent(in) :: dimofarray, index_pt
real(8), intent(inout), dimension(0:dimofarray-1) :: p, q, r
!f2py intent(in, out) :: p, q, r 
real(8), intent(inout), dimension(0:dimofarray-1) :: xpos, ypos, zpos
!f2py intent(in, out) :: xpos, ypos, zpos
    
real(8) :: p_off, q_off, r_off
p_off = p(index_pt)
q_off = q(index_pt)
r_off = r(index_pt)

!            xpos = p - p(index_pt)
!            ypos = q - q(index_pt)
!            zpos = r - r(index_pt)

    xpos = p - p_off
    ypos = q - q_off
    zpos = r - r_off
      
   
end subroutine

subroutine setonpoint3(p, q, r, xpos, ypos, zpos, stilltosort, dimofarrays, index_pt)
implicit none
integer*4, intent(in) :: dimofarrays, index_pt
real(8), intent(inout), dimension(0:dimofarrays-1) :: p, q, r
!f2py intent(in, out) :: p, q, r
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos
!f2py intent(in, out) :: xpos, ypos, zpos
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
    
real(8) :: p_off, q_off, r_off
p_off = p(index_pt)
q_off = q(index_pt)
r_off = r(index_pt)

!            xpos = p - p(index_pt)
!            ypos = q - q(index_pt)
!            zpos = r - r(index_pt)

where(stilltosort)
    xpos = p - p_off
    ypos = q - q_off
    zpos = r - r_off
end where
      
   
end subroutine

subroutine setonpoint(Xpos, Ypos, Zpos, index_pt, dimofarrays)
implicit none
integer*4, intent(inout) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: Xpos, Ypos, Zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
integer*4, intent(inout) :: index_pt




!rint *, "Index = ", Index
!print *, Xpos(Index)
!print *, Ypos


    Xpos = Xpos - Xpos(index_pt)
    Ypos = Ypos - Ypos(index_pt)
    Zpos = Zpos - Zpos(index_pt)

    !    Where (StillToSort) Xpos = Xpos - Xpos(Index)
    !    Where (StillToSort) Ypos = Ypos - Ypos(Index)
    !    Where (StillToSort) Zpos = Zpos - Zpos(Index)
   
end subroutine

subroutine sort(Xpos, Ypos, Zpos, index_pts, Dim)
implicit none
integer*4, intent(in) :: Dim
real(8), intent(inout), dimension(0:Dim-1) :: Xpos, Ypos, Zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
integer*4, intent(inout), dimension(0:Dim-1) :: index_pts
!f2py intent(in,out) :: index_pts
integer*4 :: i, j, i_minloc, Total, ResetLevel, MaskFactor, MaskCount
real(8) :: Dist_Min
real(8), dimension(0:Dim-1) :: Dist
logical, dimension(0:Dim-1) :: StillToSort, SpeedUpMask, CombinedMask, AllTrue
integer*2, dimension(0:Dim-1) :: SortCount, AllZeros
external :: Distance
external :: SetOnOrigin

real(8) :: Xref, Yref, Zref, Threshold
real(8), dimension(0:Dim-1) :: Xnew, Ynew, Znew

REAL :: t1,t2,rate 
INTEGER*4 :: c1,c2,c3, c4, c5, c6, c7, cr,cm
real(8) :: MaskRatio


    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
      
    CALL CPU_TIME(t1)
    CALL SYSTEM_CLOCK(c1)
    
    rate = REAL(cr)
    

    do i = 0, Dim    
        StillToSort(i) = .TRUE.
        SpeedUpMask(i) = .TRUE.
        CombinedMask(i) = .TRUE.
        AllTrue(i) = .TRUE.
        SortCount(i) = 0 
        AllZeros(i) = 0
    end do
    
    ResetLevel = Max(Dim * 1 / 10000,1)
    Print *, "Reset Level = ", ResetLevel
    MaskFactor = Dim * 9/10
    !Where (SpeedUpMask) SortCount = 1
    !Print *, "Total SpeedUpMask", SUM(SortCount)
    !Total = Count(SpeedUpMask)
    !Print *, "Total = ", Count(SpeedUpMask)
    !Print *, "Dim = ", Dim
    
!do i = 0, 1e3    
!    StillToSort(i) = .TRUE.
!end do

!    do i = 1e3, Dim    
!        StillToSort(i) = .FALSE.
!    end do
    
    


    Print *, "In Sort subroutine, Step 1"
    CALL SetOnOrigin(Xpos, Ypos, Zpos, Dim)
    !CALL Distance(Xpos, Ypos, Zpos, Dist, Dim)
    !Print *, "Just called Distance in sort, X, Y, Z, Dist = "
    !Print *, Xpos
    !Print *, Ypos
    !Print *, Ypos
    !Print *, Dist
    !Dist = Xpos*Xpos + Ypos*Ypos + Zpos*Zpos
    !Dist_Min = MAXVAL(Dist)
    !Print *, "Hello There from sort"
    !Print *, Dist
    Print *, "Sort, 2"
    do i = 0,Dim
    !Print *, "Sort, Step 2"    
    !do i = 0, 1    
        Print *, "Sort, 3"
        !Print *, Xpos, Ypos, Zpos
        CALL SYSTEM_CLOCK(c2)
        !CALL Distance(Xpos, Ypos, Zpos, Dist, Dim, i, StillToSort)   
        !Xref = Xpos(i)
        !Yref = Ypos(i)
        !Zref = Zpos(i)
        !Where (StillToSort) Dist = (Xpos-Xref)*(Xpos-Xref) + (Ypos-Yref)*(Ypos-Yref) + (Zpos-Zref)*(Zpos-Zref)
        CombinedMask = StillToSort .AND. SpeedUpMask
        !Print *, "CombinedMask", Count(CombinedMask)
        Print *, "Sort, 4"
        Where (CombinedMask) Dist = (Xpos-Xpos(i))*(Xpos-Xpos(i)) + (Ypos-Ypos(i))*(Ypos-Ypos(i)) + (Zpos-Zpos(i))*(Zpos-Zpos(i))
        Print *, "Sort, 5"
        CALL SYSTEM_CLOCK(c3)
        !Print *, "Sort, Step 3"
        !Dist_Min = MINVAL(Dist, 1, StillToSort) 
        !Print *, "Dist = "
        !Print *, Dist
        i_minloc = MINLOC(Dist,1, CombinedMask) - 1
        Print *, "Sort, 6"
        CALL SYSTEM_CLOCK(c4)
        index_pts(i) = i_minloc
        CALL SYSTEM_CLOCK(c5)
        StillToSort(i_minloc) = .FALSE.
        Print *, "Sort, 7"
        CALL SYSTEM_CLOCK(c6)
        !Print *, StillToSort    
        Print *, "i, ResetLevel", i, ResetLevel
        if (MOD(i,ResetLevel) .LT. 1) then 
            Print *, "Sort, 8"
            Print *, "CombinedMask", Count(CombinedMask)
            Print *, "i , i_minloc = ", i, i_minloc
            WRITE(*,*) "Distance : ",(c3 - c2)/rate
            WRITE(*,*) "MINLOC : ",(c4 - c3)/rate
            WRITE(*,*) "Assignment : ",(c5 - c4)/rate
            !WRITE(*,*) "Still to sort : ",(c6 - c5)/rate
            !WRITE(*,*) "SetOnPoint : ",(c7 - c6)/rate
            Print *, "Loop Time = ", (c6-c2)/rate, " sec"
            
            MaskRatio = 0
            do while (MaskRatio < 3)
                SpeedUpMask = AllTrue                              
                Threshold = (MAXVAL(Dist,StillToSort)-MINVAL(Dist, StillToSort))/MaskFactor
                Where (Dist > Threshold) SpeedUpMask = .FALSE.
                CombinedMask = StillToSort .AND. SpeedUpMask
                MaskCount = Count(CombinedMask)
                MaskRatio = Real(MaskCount / ResetLevel)
                Print *, "MaskRatio = ", MaskRatio
                if (MaskRatio < 3) then 
                    MaskFactor = MaskFactor / 1e3
                end if                
            end do
            
            !SortCount = AllZeros    
            !Print *, "Total SpeedUpMask", Count(SpeedUpMask)
            !Print *, Dist
            Print *, "Maxval(Xpos), Maxval(Dist), Minval, Threshold ", MAXVAL(Xpos), MAXVAL(Dist), MINVAL(Dist), Threshold
            !Print *, SpeedUpMask
            !Where (SpeedUpMask) SortCount = 1
            
            Print *, "Total SpeedUpMask", Count(SpeedUpMask)
        end if        
        !do j = 0, i
            ! Have we already indentified this cell?         
        !end do
        !CALL SetOnPoint(Xpos, Ypos, Zpos, Dim, i_minloc, StillToSort)
        CALL SYSTEM_CLOCK(c7)
    end do   
    
    WRITE(*,*) "Distance : ",(c3 - c2)/rate
    WRITE(*,*) "MINLOC : ",(c4 - c3)/rate
    WRITE(*,*) "Assignment : ",(c5 - c4)/rate
    WRITE(*,*) "Still to sort : ",(c6 - c5)/rate
    WRITE(*,*) "SetOnPoint : ",(c7 - c6)/rate
    
end subroutine

subroutine distance(xpos, ypos, zpos, dist, stilltosort, dimofarrays)
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
real(8), intent(inout), dimension(0:dimofarrays-1) :: Dist
!f2py intent(in,out) :: Dist
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
    
    WHERE(stilltosort) dist = sqrt(Xpos*Xpos + Ypos*Ypos + Zpos*Zpos) 
end subroutine

subroutine distance2(xpos, ypos, zpos, dist, dimofarrays)
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
real(8), intent(inout), dimension(0:dimofarrays-1) :: Dist
!f2py intent(in,out) :: Dist
    
    dist = sqrt(Xpos*Xpos + Ypos*Ypos + Zpos*Zpos)
end subroutine

subroutine distance3(p, q, r, dist, stilltosort, index_pt, dimofarrays)
implicit none
integer*4, intent(in) :: dimofarrays, index_pt
real(8), intent(inout), dimension(0:dimofarrays-1) :: p, q, r
!f2py intent(in,out) :: Xpos, Ypos, Zpos
real(8), intent(inout), dimension(0:dimofarrays-1) :: Dist
!f2py intent(in,out) :: Dist
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort

real(8) :: p_off, q_off, r_off
p_off = p(index_pt)
q_off = q(index_pt)
r_off = r(index_pt)
    
    WHERE(stilltosort) 
        dist = sqrt((p-p_off)*(p-p_off) + (q-q_off)*(q-q_off) + (r-r_off)*(r-r_off))
    END WHERE    
end subroutine

subroutine dist_pt_to_line(xpos, ypos, zpos, dist, stilltosort, xa, ya, za, xn, yn, zn, dimofarrays)
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
real(8), intent(inout), dimension(0:dimofarrays-1) :: Dist
!f2py intent(in,out) :: Dist
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout) :: xa, ya, za, xn, yn, zn
real(8) , dimension(0: dimofarrays - 1 ) :: xproj, yproj, zproj, dist2, dist3, dist4, dotab

real dotn
dotn = xn*xn + yn*yn +zn*zn

WRITE(*,*) "xa, ya, za, xn, yn, zn ", xa, ya, za, xn, yn, zn

    WHERE(stilltosort)         
        dist = (xa - xpos)*(xa - xpos) + (ya - ypos)*(ya - ypos) + (za - zpos)*(za - zpos) - &
                (((xa - xpos) * xn + (ya - ypos) * yn + (za - zpos) * zn))**2 / dotn 
        dist = dist ** 0.5  
    END WHERE
  
end subroutine

!subroutine project_pts_on_plane(xpos, ypos, zpos, xproj, yproj, zproj, stilltosort, xa, ya, za, &
!    xn, yn, zn, xi, yi, zi, xj, yj, zj, dimofarrays)
subroutine project_pts_on_plane(xpos, ypos, zpos, xproj, yproj, zproj, stilltosort, xa, ya, za, &
    xn, yn, zn, coeffs, dimofarrays)
    
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos, xproj, yproj, zproj
!f2py intent(in,out) :: Xpos, Ypos, Zpos, xproj, yproj, zproj
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout) :: xa, ya, za, xn, yn, zn
!f2py intent(in, out) :: xa, ya, za, xn, yn, zn
real(8), intent(inout), dimension(0:5) :: coeffs
!f2py intent(in,out) :: coeffs
real(8) , dimension(0: dimofarrays - 1 ) :: x, y, z, dist2, dist3, dist4, dotab
integer :: i, ver, hor

real(8) :: dotn, distance, magnitude, xi0, yi0, zi0
real(8) :: xi, yi, zi, xj, yj, zj
!, xi, yi, zi, xj, yj, zj
real(8), dimension(0: dimofarrays - 1) :: xpl, ypl, zpl
real(8), dimension(0:2, 0:3) :: D
integer d_dim

d_dim = 3

! A nice visual explanation of the process can be found below.
!http://stackoverflow.com/questions/9605556/how-to-project-a-3d-point-to-a-3d-plane#9605695

! First create arrays of x, y, and z points zeroed to xa, ya, za
x = xpos - xa
y = ypos - ya
z = zpos - za

dotn = (xn*xn + yn*yn +zn*zn)**0.5

xn = xn / dotn
yn = yn / dotn
zn = zn / dotn

!! Find the most suitable 'up' position, then find the most suitable horizontal axis
!! 0 = x axis, 1 = y axis, 2 = z axis
!
ver = 0
hor = 0
!IF ((zn <= yn) .AND. (zn <= xn)) THEN
!    ver = 2
!    IF (yn <= xn) THEN
!        hor = 1
!    ENDIF
!ELSE
!    IF (yn <= xn) THEN
!        ver = 1
!        IF (zn <= xn) THEN
!            hor = 2
!        ENDIF
!    ELSE
!        IF (zn <= yn) THEN
!            hor = 2
!        ELSE
!            hor = 1
!        ENDIF
!    ENDIF
!ENDIF


IF ((abs(zn) >= abs(yn)) .AND. (abs(zn) >= abs(xn))) THEN
    ver = 1
ELSE
    IF ((abs(yn) >= abs(xn)) .AND. (abs(yn) >= abs(zn))) THEN
        ver = 2
    ELSE
        hor = 1
        ver = 2
    ENDIF
ENDIF

!WRITE(*,*) "Ver, hor", hor, ver

! First project all the points onto the plane.
WHERE(stilltosort)         
    xpl = x - (x*xn + y*yn + z*zn) * xn
    ypl = y - (x*xn + y*yn + z*zn) * yn
    zpl = z - (x*xn + y*yn + z*zn) * zn
END WHERE

! At this point we have projected the point on to the arbitrary plane, but we need to orient this to two dimensions
! We need vectors starting from xa, ya, za that are oriented with 

xi0 = 0.0
yi0 = 0.0
zi0 = 0.0
IF (hor == 0) THEN
    xi0 = 1.0
ELSE
    IF (hor == 1) THEN
        yi0 = 1.0
    ELSE
        zi0 = 1.0
    ENDIF
ENDIF

dotn = xi0*xn + yi0 * yn + zi0 * zn 
xi = xi0 - dotn * xn
yi = yi0 - dotn * yn
zi = zi0 - dotn * zn

magnitude = (xi*xi + yi*yi + zi*zi)**0.5
xi = xi/ magnitude
yi = yi/magnitude
zi = zi/magnitude

!WRITE(*,*) "xi, yi, zi = ", xi, yi, zi

!That gives us the horizontal unit vector.  We want a vertical unit vector j which is perpendicular to 
! both the i vector and such that n = i x j

!print *, "Fitplane 03"

! cross product of:
!x  y   z
!xi yi zi
!xj yj zj
! = 
!xn yn zn



!D(0,0) = 0.0
!D(1,0) = zi
!D(2,0) = -1*yi
!
!D(0,1) = -1*zi
!D(1,1) = 0.0
!D(2,1) = xi
!
!D(0,2) = yi
!D(1,2) = -1 * xi
!D(2,2) = 0.0
!
!D(0,3) = xn
!D(1,3) = yn
!D(2,3) = zn
!
!CALL print_matrix(D, 3)
!
!call gauss(D, d_dim)
!
!xj = D(0,3)
!yj = D(1, 3)
!zj = D(2, 3)
!

xj = yn * zi - zn * yi
yj = zn * xi - xn * zi
zj = xn * yi - yn * xi
magnitude = (xj*xj + yj*yj + zj*zj)**0.5
xj = xj / magnitude
yj = yj / magnitude
zj = zj / magnitude

coeffs(0) = xi
coeffs(1) = yi
coeffs(2) = zi
coeffs(3) = xj
coeffs(4) = yj
coeffs(5) = zj

!WRITE(*,*) "xj, yj, zj = ", xj, yj, zj

! At this point we have two unit vectors upon which we can project all of our points


WHERE(stilltosort)         
    xproj = xpl*xi + ypl*yi + zpl*zi
    yproj = xpl*xj + ypl*yj + zpl*zj
    zproj = x*xn + y*yn + z*zn
    !Note that zproj is the distance from the old point to the projection plane
END WHERE
  
end subroutine

subroutine project_pts_on_plane2(xpos, ypos, zpos, xproj, yproj, zproj, stilltosort, x0, y0, z0, &
    xi, yi, zi, xj, yj, zj, xk, yk, zk, dimofarrays)
    
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos, xproj, yproj, zproj
!f2py intent(in,out) :: Xpos, Ypos, Zpos, xproj, yproj, zproj
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout) :: x0, y0, z0, xi, yi, zi, xj, yj, zj, xk, yk, zk
!f2py intent(in, out) :: x0, y0, z0, xi, yi, zi, xj, yj, zj, xk, yk, zk
real(8) , dimension(0: dimofarrays - 1 ) :: x, y, z, dist2, dist3, dist4, dotab

real(8) :: dotn, distance, magnitude, xi0, yi0, zi0
real(8), dimension(0: dimofarrays - 1) :: xpl, ypl, zpl

! A nice visual explanation of the process can be found below.
!http://stackoverflow.com/questions/9605556/how-to-project-a-3d-point-to-a-3d-plane#9605695

! First create arrays of x, y, and z points zeroed to x0, y0, z0
    x = xpos - x0
    y = ypos - y0
    z = zpos - z0

    ! First project all the points onto the plane.
    WHERE(stilltosort)         
        xpl = x - (x*xk + y*yk + z*zk) * xk
        ypl = y - (x*xk + y*yk + z*zk) * yk
        zpl = z - (x*xk + y*yk + z*zk) * zk
    END WHERE

    WHERE(stilltosort)         
        xproj = xpl*xi + ypl*yi + zpl*zi
        yproj = xpl*xj + ypl*yj + zpl*zj
        zproj = x*xk + y*yk + z*zk
    END WHERE
  
end subroutine

subroutine dist_plane_to_line(xpos, ypos, zpos, dist, stilltosort, xc, yc, zc, xn, yn, zn, dimofarrays)
implicit none
integer, intent(in) :: dimofarrays
real(8), intent(inout), dimension(0:dimofarrays-1) :: xpos, ypos, zpos
!f2py intent(in,out) :: Xpos, Ypos, Zpos
real(8), intent(inout), dimension(0:dimofarrays-1) :: Dist
!f2py intent(in,out) :: Dist
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout) :: xc, yc, zc, xn, yn, zn

real linemag_sqd
linemag_sqd = xn*xn + yn*yn +zn*zn

    WHERE(stilltosort)         
        dist = (((xpos - xc)*xn+(ypos - yc) * yn + (zpos - zc)*zn / linemag_sqd ) **2 ) * &
            (xn * xn + yn * yn + zn * zn)
        !dist = abs((
    END WHERE
end subroutine

subroutine estimate_density(density_big, xin, yin, zin, stilltosort, littled, bigd)
implicit none
integer*4, intent(in) :: littled, bigd
real(8), intent(inout), dimension(0:bigd - 1) :: density_big
real(8), intent(inout), dimension(0:littled-1) :: xin, yin, zin
!f2py intent(in,out) :: xin, yin, zin
logical, intent(inout), dimension(0:bigd-1) :: stilltosort
!f2py intent(in, out) :: stilltosort

integer*4 :: i, j, select_len, samples
real(8), dimension(0: littled - 1) :: x, y, z, x_ctrd, y_ctrd, z_ctrd, i_circle_r, dist
integer*4, dimension(0: bigd - 1) :: i_orig
integer*4, dimension(0: littled - 1) :: i_circle
real(8) :: density_scalar, dist_limited
logical, dimension(0: littled - 1) :: stilltosort2
integer*4 :: i_maxloc

WRITE(*,*) "Fortran Estimate Density"

stilltosort2 = .TRUE.
select_len = 299  

x = xin
y = yin
z = zin

density_big = 0.0

!WRITE(*,*) "1"
j = 0
DO i = 0, bigd - 1
    i_orig(i) = i
    IF (stilltosort(i) .EQV. .TRUE.) THEN
        i_circle(j) = i
        j = j + 1
    END IF
END DO

!WRITE(*,*) "2"
! Now self.i_circle should represent all the indices from which data was taken.
    
samples = min(3, littled)
select_len = min(select_len, int(littled/2))

!samples is how many data points around each point we will evaluate
!select_len is the number of datapoints we are going to take unless we reach the density limits first.

!WRITE(*,*) "3"

DO i = 0, littled - 1
    stilltosort2 = .TRUE.
!DO i = 0, 1
    !WRITE(*,*) i
    x_ctrd = x - x(i)
    y_ctrd = y - y(i)
    z_ctrd = z - z(i)
    !WRITE(*,*) "x(i), y(i) ", x(i), y(i)
    
    !WRITE(*,*) "4"
    
    dist = 1 / (MAX((x_ctrd * x_ctrd + y_ctrd * y_ctrd),0.00001))**0.5
    
    !WRITE(*,*) "dist " , dist
    density_scalar = 0.0
    stilltosort2(i) = .FALSE.
    
    !WRITE(*,*) "5"
    
    DO j = 1, samples - 1
        !WRITE(*,*) "6: ", i, j
        !WRITE(*,*) "dist ", dist
        !WRITE(*,*) "littled , bigd", littled, bigd
        !WRITE(*,*) "dist ", dist(0), dist(1), dist(2), "stilltosort ", stilltosort2(0), stilltosort2(1), stilltosort2(2)
        !dist_limited = ABS(z_ctrd(i_minloc)-z_ctrd(i) / dist(i_minloc)
        i_maxloc = MAXLOC(dist, 1, stilltosort2) - 1
        !WRITE(*,*) "dist ", dist(0), dist(1), dist(2), "stilltosort ", stilltosort2(0), stilltosort2(1), stilltosort2(2)
        !WRITE(*,*) "Min of dist ", MINVAL(dist), MINVAL(dist, stilltosort2), "Location ", i_minloc, "dist(i_minloc ", dist(i_minloc)
        !WRITE(*,*) "dist ", dist(0), dist(1), dist(2), "stilltosort ", stilltosort2(0), stilltosort2(1), stilltosort2(2)
        !WRITE(*,*) dist
        !WRITE(*,*) "j = ", j
        i_maxloc = MAX(0, i_maxloc)
            
        !IF (dist(i_minloc) .LT. 0.0001) THEN
        !    dist_limited = 0.0001
        !    WRITE(*,*) "Limited distance, x_ctrd, y_ctrd, z_ctrd = ", x_ctrd(i_minloc), y_ctrd(i_minloc), z_ctrd(i_minloc)
        !ELSE
        !    dist_limited = dist(i_minloc)
        !ENDIF
        
!        density_scalar = density_scalar + 1 / MAX(0.0001, dist(i_minloc)**0.5)
        !density_scalar = density_scalar + ABS(z_ctrd(i_minloc)-z_ctrd(i)) / dist_limited
        density_scalar = density_scalar + dist(i_maxloc)
!        density_scalar = density_scalar + 1 / MAX(0.0001, dist(i_minloc))
        !WRITE(*,*) "i ", i, " MINLOC ", i_minloc, " dist ", dist(i_minloc), " density_scalar ", density_scalar
        !WRITE(*,*) density_scalar
        !WRITE(*,*) i_minloc, i_circle(i_minloc)
        stilltosort2(i_maxloc) = .FALSE.
    END DO
    !stilltosort2(i) = .TRUE.
    density_big(i_circle(i)) = density_scalar
    !stilltosort2 = .TRUE.
END DO    
    
END subroutine

subroutine estimate_density2(density_big, xin, yin, zin, stilltosort, littled, bigd)
implicit none
integer*4, intent(in) :: littled, bigd
real(8), intent(inout), dimension(0:bigd - 1) :: density_big
real(8), intent(inout), dimension(0:littled-1) :: xin, yin, zin
!f2py intent(in,out) :: xin, yin, zin
logical, intent(inout), dimension(0:bigd-1) :: stilltosort
!f2py intent(in, out) :: stilltosort

integer*4 :: i, j, select_len, samples
real(8), dimension(0: littled - 1) :: x, y, z, x_ctrd, y_ctrd, z_ctrd, i_circle_r, dist
integer*4, dimension(0: bigd - 1) :: i_orig
integer*4, dimension(0: littled - 1) :: i_circle
real(8) :: density_scalar, dist_limited, dist_ind, dist_sqd, d_max, d_min, d_avg
logical, dimension(0: littled - 1) :: stilltosort2
integer*4 :: i_minloc

!WRITE(*,*) "Fortran Estimate Density"

stilltosort2 = .TRUE.
select_len = 299  

x = xin
y = yin
z = zin

density_big = 0.0

j = 0
DO i = 0, bigd - 1
    i_orig(i) = i
    IF (stilltosort(i) .EQV. .TRUE.) THEN
        i_circle(j) = i
        j = j + 1
    END IF
END DO

! Now self.i_circle should represent all the indices from which data was taken.
    
samples = min(30, littled)
select_len = min(select_len, int(littled/2))

!samples is how many data points around each point we will evaluate
!select_len is the number of datapoints we are going to take unless we reach the density limits first.


stilltosort2 = .TRUE.
!i = 0
DO i = 0, littled - 1

!DO i = 0, 1
    !WRITE(*,*) i
    x_ctrd = x - x(i)
    y_ctrd = y - y(i)
    z_ctrd = z - z(i)
    !WRITE(*,*) "x(i), y(i) ", x(i), y(i)
    
!    dist = 1 / (MAX((x_ctrd * x_ctrd + y_ctrd * y_ctrd),0.00001))**0.5
!    dist = 1 / (MAX((x_ctrd * x_ctrd + y_ctrd * y_ctrd),0.00001))
!    dist = sqrt(x_ctrd * x_ctrd + y_ctrd * y_ctrd)
    dist = (x_ctrd * x_ctrd + y_ctrd * y_ctrd)
!    dist = x_ctrd * x_ctrd + y_ctrd * y_ctrd
    
    !WRITE(*,*) "dist " , dist
    density_scalar = 0.0
    stilltosort2(i) = .FALSE.
    
!    DO j = 1, samples - 1
!        i_minloc = MINLOC(dist, 1, stilltosort2) - 1
!        i_minloc = MAX(0, i_minloc)
!        density_scalar = density_scalar + dist(i_minloc)
!        stilltosort2(i_minloc) = .FALSE.
!    END DO


    DO j = 1, samples - 1
        i_minloc = MINLOC(dist, 1, stilltosort2) - 1
        i_minloc = MAX(0, i_minloc)
    !    density_scalar = density_scalar + dist(i_minloc)
        dist_ind = abs(sqrt(x(i_minloc) * x(i_minloc) + y(i_minloc) * y(i_minloc)) - &
            sqrt(x(i) * x(i) + y(i) * y(i)))
        dist_sqd = dist_ind * dist_ind
        density_scalar = density_scalar + dist_sqd
        stilltosort2(i_minloc) = .FALSE.
    END DO
    

    !stilltosort2(i) = .TRUE.
    density_big(i_circle(i)) = 1 / density_scalar / 1000000
!    i = i_minloc
    !stilltosort2 = .TRUE.
END DO  

!d_max = MAX(density_big, stilltosort)
!   d_min = MIN(density_big, stilltosort)
d_avg = SUM(density_big,stilltosort) / littled  

Do i = 0, littled -1
    density_big(i_circle(i)) = MAX(density_big(i_circle(i)), d_avg)
END DO
    
END subroutine

subroutine findmin(array, stilltosort, minimum_pos, minimum_val, gsize, dimofarrays)
implicit none
integer*4, intent(in) :: dimofarrays, gsize
real(8), intent(inout), dimension(0:dimofarrays-1) :: array
!f2py intent(in,out) :: Xpos, Ypos, Zpos
logical, intent(inout), dimension(0:dimofarrays-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
integer*4, intent(inout), dimension(0:gsize - 1) :: minimum_pos
real(8), intent(inout), dimension(0:gsize - 1) :: minimum_val

integer :: minimum_pos0, minimum_pos1, minimum_pos2, minimum_pos3
real(8) :: minimum_val0, minimum_val1, minimum_val2, minimum_val3

logical, dimension(0:dimofarrays-1) :: disttoogreat
real(8),  dimension(0:dimofarrays-1) :: roughdist
    
integer*4 :: i


!minimum_pos = -1
!minimum_val = 0

minimum_pos1 = -1
minimum_pos2 = -1
minimum_pos3 = -1
minimum_val1 = 0
minimum_val2 = 0
minimum_val3 = 0


DO i = 0, dimofarrays - 1
  IF (stilltosort(i)) THEN
    IF (minimum_pos1 < 0) THEN
      minimum_pos1 = i
      minimum_val1 = array(i)
      minimum_val2 = array(i)+1
      minimum_val3 = array(i)+2
    ELSE IF ((array(i) < minimum_val3) .OR. (minimum_pos3 < 0)) THEN
      IF ((array(i) < minimum_val2) .OR. (minimum_pos2 < 0)) THEN
            minimum_pos3 = minimum_pos2
            minimum_val3 = minimum_val2  
        IF (array(i) < minimum_val1) THEN
            minimum_pos2 = minimum_pos1
            minimum_val2 = minimum_val1           
            minimum_pos1 = i
            minimum_val1 = array(i)
        ELSE 
          minimum_pos2 = i
          minimum_val2 = array(i)
        ENDIF
      ELSE
        minimum_pos3 = i
        minimum_val3 = array(i)
      ENDIF
    ENDIF
  ENDIF
ENDDO

!WRITE(*,*) "Minimum_pos1", minimum_pos1, "array(minimum_pos)", array(minimum_pos1)
minimum_pos(1) = minimum_pos1
minimum_pos(2) = minimum_pos2
minimum_pos(3) = minimum_pos3

!WRITE(*,*) minimum_val1
!WRITE(*,*) minimum_val2
!WRITE(*,*) minimum_val3
!
!WRITE(*,*) minimum_pos1
!WRITE(*,*) minimum_pos2
!WRITE(*,*) minimum_pos3
!
minimum_val(1) = minimum_val1
minimum_val(2) = minimum_val2
minimum_val(3) = minimum_val3

end subroutine

subroutine group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                        stilltosort_r, g0_max, g1_max, r_min, max_dist, i_min, i_max, gsize, npoints)
implicit none
integer*4, intent(in) :: npoints
real(8), intent(in), dimension(0:npoints-1) :: p, q, r
real(8), intent(inout), dimension(0:npoints-1) :: g0, g1, g2, g3, i1, i2, i3, i4, stilltosort_r
logical, intent(inout), dimension(0:npoints-1) :: stilltosort
integer*4, intent(in) :: gsize
real(8), intent(in) :: g0_max, g1_max, r_min, max_dist
integer*4, intent(in) :: i_min, i_max

logical, dimension(0:npoints-1) :: CombinedMask, alltrue, r_too_low, j_loop
integer*4 :: i, j, i_minloc, gsize_wk, npoints_wk, total_points
real(8), dimension(0:gsize - 1) :: gx, gy, gz, gsize_zeros, group
real(8), dimension(0:npoints - 1) :: x, y, z, dist, group_r
real(8), dimension(0:3) :: gplane
integer*4, dimension(0:gsize - 1) :: group_i

real(8), dimension(0:npoints - 1) :: p_wk, q_wk, r_wk
logical :: error_flag

total_points = 0

p_wk = p
q_wk = q
r_wk = r

gsize_wk = gsize
npoints_wk = npoints
! the above is necessary due to a problem with using fortran code
! from python.  If I declare the array size in fitplane as 'IN'
! in fitplane and try to use it in python, it decodes improperly.
! Therefore, I leave all code to be used in python as 'INOUT' but
! if I also use it in FORTRAN code which properly indicates
! array sizes as read only, a warning is thrown.

!stilltosort_2_wk = stilltosort

alltrue = .TRUE.
error_flag = .FALSE.

group_r = 0
DO i = 0, gsize - 1
    group(i) = 0.0
    group_i(i) = -1
    gx(i) = 0.0
    gy(i) = 0.0
    gz(i) = 0.0
    gsize_zeros = 0.0
ENDDO

!WRITE(*,*) "g0_max, g1_max, r_min ", g0_max, g1_max, r_min

DO i = i_min, i_max
    IF (MOD(i,10000) == 0) THEN   
        PRINT *, "i = ", i, " / ", i_max, " / ", npoints
    ENDIF
    IF (stilltosort(i)) THEN
!    IF (alltrue(i) .EQV. .FALSE.) THEN
        stilltosort(i) = .FALSE.
        dist = max_dist
        CALL distance3(p_wk,q_wk,r_wk,dist,stilltosort, i, npoints)
        gx = gsize_zeros 
        gy = gsize_zeros
        gz = gsize_zeros
        !DO j = 0, gsize - 1
        !    group_i(j) = -1
        !ENDDO  
        group_i = -1
        group_i(0) = i

!        j_loop = stilltosort

        CALL findmin(dist, stilltosort, group_i, group, gsize, npoints)


        gx(0) = p(i)
        gy(0) = q(i)
        gz(0) = r(i)
        error_flag = .FALSE.
        DO j = 1, gsize - 1

            i_minloc = group_i(j)
            IF (i_minloc .EQ. group_i(j-1)) THEN
                error_flag = .TRUE.
            ENDIF
            gx(j) = p(i_minloc)
            gy(j) = q(i_minloc)
            gz(j) = r(i_minloc)

          SELECT CASE (j)
            case (1)
               i1(i) = REAL(i_minloc)
            case (2)
              i2(i) = REAL(i_minloc)        
            case (3)
              i3(i) = REAL(i_minloc)        
            end select        
        ENDDO     
        
! Good up til this point
        group(1) = dist(group_i(1))
        group(2) = dist(group_i(2))
        group(3) = dist(group_i(3))

        IF (error_flag .EQV. .FALSE.) THEN
!            WRITE (*,*) gx, gy, gz, gplane, gsize_wk
            CALL fitplane(gx,gy,gz,gplane, gsize_wk) 
            group_r(i) = gplane(3)
            g0(i) = gplane(0)
            g1(i) = gplane(1)
            g2(i) = gplane(2)
            g3(i) = gplane(3)
            
            IF (gplane(3) .GT.r_min) THEN
                IF (abs(gplane(0)) .LT. g0_max) THEN
                    IF (abs(gplane(1)) .LT. g1_max) THEN
                        !stilltosort_2_wk(i) = .TRUE.
                        stilltosort_r(i) = 1.0
                        stilltosort_r(INT(i1(i))) = 1.0
                        stilltosort_r(INT(i2(i))) = 1.0
                        stilltosort_r(INT(i3(i))) = 1.0
                        stilltosort(INT(i1(i))) = .FALSE.
                        stilltosort(INT(i2(i))) = .FALSE.
                        stilltosort(INT(i3(i))) = .FALSE.
                        
!                        total_points = total_points + 1
                    ENDIF                
                ENDIF            
            ENDIF 
        
        ELSE 
            WRITE(*,*) "Error found"
        ENDIF
        
! Problem is somewhere above this point

        
    ENDIF
ENDDO
!WRITE(*,*) "Total points ", total_points
!stilltosort_2 = stilltosort_2_wk

end subroutine



subroutine group_points_threaded2(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                        stilltosort_r, g0_max, g1_max, r_min, max_dist, i_min, i_max, &
                        index_little, index_little_bw, thread, gsize, littled, npoints)
implicit none
integer*4, intent(in) :: npoints, littled
real(8), intent(in), dimension(0:littled-1) :: p, q, r
real(8), intent(inout), dimension(0:npoints-1) :: g0, g1, g2, g3, i1, i2, i3, i4, stilltosort_r
logical, intent(inout), dimension(0:npoints-1) :: stilltosort
integer, intent(in) :: gsize
real(8), intent(in) :: g0_max, g1_max, r_min, max_dist
integer*4, intent(in) :: i_min, i_max, thread
integer*4, intent(inout), dimension(0:littled - 1) :: index_little
integer*4, intent(inout), dimension(0:npoints - 1) :: index_little_bw

logical, dimension(0:npoints-1) :: CombinedMask, alltrue, r_too_low, j_loop
logical, dimension(0:littled - 1) :: stilltosort_l
integer*4 :: i, j, i_minloc, i_minloc_last, gsize_wk, npoints_wk, total_points, index_l, i_l
real(8), dimension(0:gsize - 1) :: gx, gy, gz, gsize_zeros, group, group_l
real(8), dimension(0:littled - 1) :: x, y, z, dist
real(8), dimension(0:npoints - 1) :: group_r
real(8), dimension(0:3) :: gplane
integer*4, dimension(0:gsize - 1) :: group_i, group_i_l


real(8), dimension(0:npoints - 1) :: p_wk, q_wk, r_wk
integer :: per_complete

logical :: error_flag

total_points = 0

p_wk = p
q_wk = q
r_wk = r

gsize_wk = gsize

stilltosort_l = .TRUE.

!stilltosort_2_wk = stilltosort

alltrue = .TRUE.
error_flag = .FALSE.

group_r = 0
DO i = 0, gsize - 1
    group(i) = 0.0
    group_i(i) = -1
    gx(i) = 0.0
    gy(i) = 0.0
    gz(i) = 0.0
    gsize_zeros = 0.0
ENDDO

!WRITE(*,*) "g0_max, g1_max, r_min ", g0_max, g1_max, r_min
WRITE(*,*) "thread ", thread, " started. npoints ", npoints, " littled = ", littled, " i_min = ", i_min, " i_max = ", i_max 

DO i = i_min, i_max
    IF (MOD(i,10000) == 0) THEN 
        per_complete = INT(100*(i - i_min)/(i_max-i_min))
        PRINT *, "thread ", thread, "   ", per_complete, " % complete"
    ENDIF
    IF (stilltosort(i)) THEN
!    IF (alltrue(i) .EQV. .FALSE.) THEN
        stilltosort(i) = .FALSE.
        index_l = index_little_bw(i)
        stilltosort_l(index_l) = .FALSE.
        dist = max_dist
        CALL distance3(p_wk,q_wk,r_wk,dist,stilltosort_l, index_l, littled)
        gx = gsize_zeros 
        gy = gsize_zeros
        gz = gsize_zeros
        !DO j = 0, gsize - 1
        !    group_i(j) = -1
        !ENDDO  
        group_i_l = -1
        group_i_l(0) = index_l

!        j_loop = stilltosort

        CALL findmin(dist, stilltosort_l, group_i_l, group, gsize, littled)
        error_flag = .FALSE.
    
    
        group_i(0) = i
!        do i_l = 1, gsize - 1
!
!                group_i(i_l) = index_little(group_i_l(i_l))
!                IF (group_i(i_l) .LT. 0) THEN
!                    error_flag = .TRUE.
!                ENDIF
!            ENDIF
!        ENDDO
        

 !       gx(0) = p_wk(i)
 !       gy(0) = q_wk(i)
 !       gz(0) = r_wk(i)


        

        i_minloc_last = -1
        error_flag = .FALSE.
        DO j = 0, gsize - 1
            IF (error_flag .EQV. .FALSE.) THEN
                i_l = group_i_l(j)
                IF (i_l .LT. 0) THEN
                    WRITE(*,*) "Thread ", thread, "group_i_l(",j,") < 0 "
                    error_flag = .TRUE.
                ELSE
                    IF(i_l .GT. (littled - 1)) THEN
                        WRITE(*,*) "Thread ", thread, "group_i_l(",j,") > littled - 1 ", "i_l = ",i_l
                        error_flag = .TRUE.
                    ELSE
                        i_minloc = index_little(i_l)
                        if (i_minloc .LT. 0) THEN
                            WRITE(*,*) "Thread ", thread, "i_minloc",j,") < 0 " 
                            error_flag = .TRUE.
                        ELSE
                            IF (i_minloc .GT. (npoints - 1)) THEN
                                WRITE(*,*) "Thread ", thread, "i_minloc",j,") > npoints - 1 "
                                error_flag = .TRUE.
                            ELSE
                                IF (i_minloc .EQ. i_minloc_last) THEN
                                    WRITE(*,*) "Thread ", thread, "non-unique minimum points"
                                    error_flag = .TRUE.
                                ELSE
                                    i_minloc_last = i_minloc
                                ENDIF
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF
            ENDIF
            gx(j) = p(i_l)
            gy(j) = q(i_l)
            gz(j) = r(i_l)

            SELECT CASE (j)
                case (1)
                    i1(i) = REAL(i_minloc)
                case (2)
                    i2(i) = REAL(i_minloc)        
                case (3)
                    i3(i) = REAL(i_minloc)        
            end select        
        ENDDO     
        
        IF (error_flag .EQV. .FALSE.) THEN
        
! Good up til this point
 !           group(1) = dist(group_i(1))
 !           group(2) = dist(group_i(2))
 !           group(3) = dist(group_i(3))

!            WRITE (*,*) gx, gy, gz, gplane, gsize_wk
            gplane = 0.0
            CALL fitplane(gx,gy,gz,gplane, gsize_wk) 
            group_r(i) = gplane(3)
            g0(i) = gplane(0)
            g1(i) = gplane(1)
            g2(i) = gplane(2)
            g3(i) = gplane(3)
            

!            stilltosort(i1(1)) = .FALSE.
!            stilltosort(i1(2)) = .FALSE.
!            stilltosort(i1(3)) = .FALSE.
!            stilltosort_l(index_little_bw(INT(i1(i)))) = .FALSE.
!            stilltosort_l(index_little_bw(INT(i2(i)))) = .FALSE.
!            stilltosort_l(index_little_bw(INT(i3(i)))) = .FALSE.
            
            IF (gplane(3) .GT.r_min) THEN
                IF (abs(gplane(0)) .LT. g0_max) THEN
                    IF (abs(gplane(1)) .LT. g1_max) THEN
                        !stilltosort_2_wk(i) = .TRUE.
                        stilltosort_r(i) = 1.0
                        stilltosort_r(INT(i1(i))) = 1.0
                        stilltosort_r(INT(i2(i))) = 1.0
                        stilltosort_r(INT(i3(i))) = 1.0
                        stilltosort(INT(i1(1))) = .FALSE.
                        stilltosort(INT(i1(2))) = .FALSE.
                        stilltosort(INT(i1(3))) = .FALSE.
                        stilltosort_l(index_little_bw(INT(i1(i)))) = .FALSE.
                        stilltosort_l(index_little_bw(INT(i2(i)))) = .FALSE.
                        stilltosort_l(index_little_bw(INT(i3(i)))) = .FALSE.

                        
                        
!                        total_points = total_points + 1
                    ENDIF                
                ENDIF            
            ENDIF 
        
        ELSE 
            WRITE(*,*) "Error found"
        ENDIF
        
! Problem is somewhere above this point

        
    ENDIF
ENDDO
!WRITE(*,*) "Total points ", total_points
!stilltosort_2 = stilltosort_2_wk

end subroutine

subroutine group_points(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, stilltosort_r, &
                        g0_max, g1_max, r_min, gsize, npoints)
use omp_lib
implicit none
integer*4, intent(inout) :: npoints
real(8), intent(inout), dimension(0:npoints-1) :: p, q, r, g0, g1, g2, g3, i1, i2, i3, i4
!f2py intent(in,out) :: p, q, r, g0, g1, g2, g3, i1, i2, i3, i4
logical, intent(inout), dimension(0:npoints-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout), dimension(0:npoints-1) :: stilltosort_r
!f2py intent(in, out) :: stilltosort_r
integer, intent(inout) :: gsize
real(8), intent(inout) :: g0_max, g1_max, r_min


logical, dimension(0:npoints-1) :: stilltosort_2, stilltosort_2a, stilltosort_2b, stilltosort_2c, &
                                    CombinedMask, alltrue, r_too_low, j_loop

! This routine takes in points defined by p, q, and r and arranges them into 
! groups of gsize.  The points taken are the closest point.

real(8), dimension(0:npoints - 1) :: x, y, z, dist
real(8) :: max_x, max_y, max_z, max_dist
integer*4 :: i, j, k, i_minloc, threads, i_min, i_max
real(8), dimension(0:gsize - 1) :: gx, gy, gz, gsize_zeros, group
! gplane will contain the results of a planar regression, with element
! 3 containing the R^2 value.
real(8), dimension(0:3) :: gplane
integer*4, dimension(0:gsize - 1) :: group_i



REAL :: t1,t2,rate, ta0, ta1, ta, tb0, tb1, tb, tc0, tc1, tc, td0, td1, td, te0, te1, te, tf0
REAL :: tf1, tf, tg0, tg1, tg, th0, th1, th, t_tot, t_ml0, t_ml1, t_ml
INTEGER :: c1,c2,cr,cm
REAL(8) :: minimum_val
INTEGER :: minimum_pos
integer*4, dimension(0 : 2) :: breaks
integer*4 :: instances

threads = 3
j = 0
k = 0

stilltosort_r = 0.0
WHERE(stilltosort) stilltosort_r = 1.0
instances = int(sum(stilltosort_r))
write(*,*) "instances ", instances

IF (instances > 100 * threads) THEN
    breaks(0) = 0
    do i = 0, npoints - 1
        if (stilltosort(i)) then
            j = j + 1
            if (j > (instances / threads + 10)) then
                j = 0
                k = k + 1
                breaks(k) = i
            end if
        endif    
    enddo
ELSE
    breaks(0) = 0
    DO j = 1, threads - 1
        breaks(j) = npoints - 1
    END DO
    IF(instances < 3) THEN
        WRITE(*,*) "Insufficient points to evaluate"
        GOTO 900
    ENDIF
ENDIF

write(*,*) "Breaks ", breaks

minimum_val = 0
minimum_pos = 0


CALL system_clock(count_rate=cr)
CALL system_clock(count_max=cm)
      
CALL CPU_TIME(t1)
CALL SYSTEM_CLOCK(c1)

rate = REAL(cr)
!WRITE(*,*) "system_clock rate ",rate    
      

x = p
y = q
z = r
g0 = 0.0
g1 = 0.0
g2 = 0.0
g3 = 0.0
i1 = 0.0
i2 = 0.0
i3 = 0.0
i4 = 0.0

max_x = MAXVAL(x, 1)
max_y = MAXVAL(y, 1)
max_z = MAXVAL(z, 1)

max_dist = max_x * max_x + max_y * max_y + max_z * max_z
max_dist = max_dist*10 + 1.0

! Set entire dist array to the maximum possible distance
dist = max_dist
!print *, "Hello, entered calculation group_points"
!print *, "Maximum distance", max_dist

alltrue = .TRUE.


DO i = 0, gsize - 1
    group(i) = 0.0
    group_i(i) = -1
    gx(i) = 0.0
    gy(i) = 0.0
    gz(i) = 0.0
    gsize_zeros = 0.0
ENDDO
! Should make the above just assignments to 0

!print *, "stilltosort", stilltosort
stilltosort_2 = stilltosort
stilltosort_r = 0.0
!stilltosort_2a = stilltosort
!stilltosort_2b = stilltosort
!stilltosort_2c = stilltosort
!print *, "stilltosort_2 ", stilltosort_2
ta = 0
tb = 0
tc = 0
td = 0
te = 0
tf = 0
tg = 0
t_ml = 0

i_min = 0
i_max = npoints - 1
!First loop in group
write(*,*) "about to enter parallel section"
IF (breaks(1) .NE. (i_max)) THEN
    !$omp parallel shared (stilltosort)
    !private (p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, g0_max, g1_max, r_min, max_dist, gsize, npoints)

    !$omp sections

        !$omp section
        call group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist, 0,breaks(1), gsize, npoints)

        !$omp section
        call group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist, breaks(1) + 1,breaks(2), gsize, npoints)
        !$omp section

        call group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist,breaks(2) + 1, npoints - 1, gsize, npoints)
    !$omp end sections

    !$omp end parallel




    print *, "End of first group in sort"

    !stilltosort = stilltosort_2
    !print *, stilltosort

    !stilltosort_r = 0.0
    
    !WHERE(stilltosort_2a) stilltosort_r = 1.0
    !WHERE(stilltosort_2b) stilltosort_r = 1.0
    !WHERE(stilltosort_2c) stilltosort_r = 1.0
    
    !WHERE(stilltosort) stilltosort_r = 1.0
    write(*,*) sum(stilltosort_r)
ELSE
   call group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, stilltosort_r, &
                                g0_max, g1_max, r_min, max_dist,0, npoints - 1, gsize, npoints)
ENDIF

!CALL CPU_TIME(t2)
CALL SYSTEM_CLOCK(c2)
      
WRITE(*,*) "system_clock : ",(c2 - c1)/rate
!WRITE(*,*) "cpu_time     : ",(t2-t1)
t_tot = t2-t1

!WRITE(*,*) "time a: ", ta
!WRITE(*,*) "time b: ", tb
!WRITE(*,*) "time c: ", tc
!WRITE(*,*) "time d: ", td
!WRITE(*,*) "time e: ", te
!WRITE(*,*) "time f: ", tf
!WRITE(*,*) "time g: ", tg
!WRITE(*,*) "time ml : ", t_ml

!WRITE(*,*) "time a+...+g", ta + tb + tc + td + te + tf + tg
!WRITE(*,*) "time_tot ", t_tot


900 print *, "group_points_4"            
end subroutine  

subroutine group_points2(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, stilltosort_r, &
                        g0_max, g1_max, r_min, gsize, littled, npoints)
use omp_lib
implicit none
integer*4, intent(inout) :: npoints, littled
real(8), intent(inout), dimension(0:npoints-1) :: p, q, r, g0, g1, g2, g3, i1, i2, i3, i4
!f2py intent(in,out) :: p, q, r, g0, g1, g2, g3, i1, i2, i3, i4
logical, intent(inout), dimension(0:npoints-1) :: stilltosort
!f2py intent(in, out) :: stilltosort
real(8), intent(inout), dimension(0:npoints-1) :: stilltosort_r
!f2py intent(in, out) :: stilltosort_r
integer, intent(inout) :: gsize
real(8), intent(inout) :: g0_max, g1_max, r_min


logical, dimension(0:npoints-1) :: stilltosort_2, CombinedMask, alltrue, r_too_low, j_loop

! This routine takes in points defined by p, q, and r and arranges them into 
! groups of gsize.  The points taken are the closest point.

real(8), dimension(0:littled - 1) :: x, y, z, dist
integer*4, dimension(0:littled - 1) :: index_little
logical*4, dimension(0:littled - 1) :: stilltosort_l
integer*4, dimension(0:npoints - 1) :: index_little_bw
real(8) :: max_x, max_y, max_z, max_dist
integer*4 :: i, j, k, i_minloc, threads, i_min, i_max
real(8), dimension(0:gsize - 1) :: gx, gy, gz, gsize_zeros, group
! gplane will contain the results of a planar regression, with element
! 3 containing the R^2 value.
real(8), dimension(0:3) :: gplane
integer*4, dimension(0:gsize - 1) :: group_i



REAL :: t1,t2,rate, ta0, ta1, ta, tb0, tb1, tb, tc0, tc1, tc, td0, td1, td, te0, te1, te, tf0
REAL :: tf1, tf, tg0, tg1, tg, th0, th1, th, t_tot, t_ml0, t_ml1, t_ml
INTEGER :: c1,c2,cr,cm
REAL(8) :: minimum_val
INTEGER :: minimum_pos
integer*4, dimension(0 : 2) :: breaks
integer*4 :: instances

WRITE(*,*) "npoints = ",npoints, " littled = ", littled, " gsize = ", gsize
stilltosort_l = .TRUE.
j = 0
index_little_bw = -1
DO i = 0, npoints - 1
    IF (stilltosort(i)) THEN
        index_little(j) = i
        index_little_bw(i) = j
        x(j) = p(i)
        y(j) = q(i)
        z(j) = r(i)
        j = j + 1
    ENDIF
END DO



threads = 3
j = 0
k = 0

stilltosort_r = 0.0
WHERE(stilltosort) stilltosort_r = 1.0
instances = int(sum(stilltosort_r))
write(*,*) "instances ", instances

IF (instances > 100 * threads) THEN
    breaks(0) = 0
    do i = 0, npoints - 1
        if (stilltosort(i)) then
            j = j + 1
            if (j > (instances / threads + 10)) then
                j = 0
                k = k + 1
                breaks(k) = i
            end if
        endif    
    enddo
ELSE
    breaks(0) = 0
    DO j = 1, threads - 1
        breaks(j) = npoints - 1
    END DO
    IF(instances < 3) THEN
        WRITE(*,*) "Insufficient points to evaluate"
        GOTO 900
    ENDIF
ENDIF

write(*,*) "Breaks ", breaks

minimum_val = 0
minimum_pos = 0


CALL system_clock(count_rate=cr)
CALL system_clock(count_max=cm)
      
CALL CPU_TIME(t1)
CALL SYSTEM_CLOCK(c1)

rate = REAL(cr)
!WRITE(*,*) "system_clock rate ",rate    
      

!x = p
!y = q
!z = r
g0 = 0.0
g1 = 0.0
g2 = 0.0
g3 = 0.0
i1 = 0.0
i2 = 0.0
i3 = 0.0
i4 = 0.0

max_x = MAXVAL(x, 1)
max_y = MAXVAL(y, 1)
max_z = MAXVAL(z, 1)

max_dist = max_x * max_x + max_y * max_y + max_z * max_z
max_dist = max_dist*10 + 1.0

! Set entire dist array to the maximum possible distance
dist = max_dist
!print *, "Hello, entered calculation group_points"
!print *, "Maximum distance", max_dist

alltrue = .TRUE.

DO i = 0, gsize - 1
    group(i) = 0.0
    group_i(i) = -1
    gx(i) = 0.0
    gy(i) = 0.0
    gz(i) = 0.0
    gsize_zeros = 0.0
ENDDO
! Should make the above just assignments to 0

!print *, "stilltosort", stilltosort
stilltosort_2 = stilltosort
stilltosort_r = 0.0
!stilltosort_2a = stilltosort
!stilltosort_2b = stilltosort
!stilltosort_2c = stilltosort
!print *, "stilltosort_2 ", stilltosort_2
ta = 0
tb = 0
tc = 0
td = 0
te = 0
tf = 0
tg = 0
t_ml = 0

i_min = 0
i_max = npoints - 1
!First loop in group
write(*,*) "about to enter parallel section"
IF (breaks(1) .NE. (i_max)) THEN
    !$omp parallel shared (stilltosort)
    !private (p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, g0_max, g1_max, r_min, max_dist, gsize, npoints)

    !$omp sections

        !$omp section
        call group_points_threaded2(x, y, z, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist, &
                                0,breaks(1), index_little, index_little_bw, 1, gsize, littled, npoints)

        !$omp section
        call group_points_threaded2(x, y, z, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist, &
                                breaks(1) + 1, breaks(2), index_little, index_little_bw, 2, gsize, littled, npoints)
        !$omp section

        call group_points_threaded2(x, y, z, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, &
                                stilltosort_r, g0_max, g1_max, r_min, max_dist, &
                                breaks(2) + 1, npoints - 1, index_little, index_little_bw, 3,  gsize, littled, npoints)
    !$omp end sections

    !$omp end parallel




    print *, "End of first group in sort"

    !stilltosort = stilltosort_2
    !print *, stilltosort

    !stilltosort_r = 0.0
    
    !WHERE(stilltosort_2a) stilltosort_r = 1.0
    !WHERE(stilltosort_2b) stilltosort_r = 1.0
    !WHERE(stilltosort_2c) stilltosort_r = 1.0
    
    !WHERE(stilltosort) stilltosort_r = 1.0
    write(*,*) sum(stilltosort_r)
ELSE
   !call group_points_threaded(p, q, r, g0, g1, g2, g3, i1, i2, i3, i4, stilltosort, stilltosort_r, &
   !                             g0_max, g1_max, r_min, max_dist,0, npoints - 1, gsize, npoints)
ENDIF

!CALL CPU_TIME(t2)
CALL SYSTEM_CLOCK(c2)
      
WRITE(*,*) "system_clock : ",(c2 - c1)/rate
!WRITE(*,*) "cpu_time     : ",(t2-t1)
t_tot = t2-t1

!WRITE(*,*) "time a: ", ta
!WRITE(*,*) "time b: ", tb
!WRITE(*,*) "time c: ", tc
!WRITE(*,*) "time d: ", td
!WRITE(*,*) "time e: ", te
!WRITE(*,*) "time f: ", tf
!WRITE(*,*) "time g: ", tg
!WRITE(*,*) "time ml : ", t_ml

!WRITE(*,*) "time a+...+g", ta + tb + tc + td + te + tf + tg
!WRITE(*,*) "time_tot ", t_tot


900 print *, "group_points_4"            
end subroutine  



subroutine findlongeststring(p, q, r, g0, g1, g2, g3, g0avg, g1avg, g2avg, next_index, pts_total, &
                                stilltosort_in, stilltosort_r, z_err_lim1, z_err_lim2, littled, bigd)
implicit none
integer*4, intent(in) :: littled, bigd
real(8), intent(inout), dimension(0:bigd-1) :: p, q, r, g0, g1, g2, g3, g0avg, g1avg, g2avg, &
                                                pts_total, stilltosort_r, next_index
!f2py intent(in,out) :: p, q, r, g0, g1, g2, g3, g0avg, g1avg, g2avg, pts_total, stilltosort_r, next_index
logical, intent(inout), dimension(0:bigd-1) :: stilltosort_in
!f2py intent(in, out) :: stilltosort_in
real(8), intent(inout) :: z_err_lim1, z_err_lim2

logical, dimension(0:bigd-1) :: stilltosort
real(8) :: z_pred, g0_tot, g1_tot, g2_tot, temp, max_envelope, z_err1, z_err2
real(8), dimension(0:bigd - 1) :: index_big, empty_array, dist
integer*4 :: i, j, j_start, k, pts_max_index, pts_counter
logical :: found

stilltosort = stilltosort_in
!max_envelope = max(max(p), max(q), max(r), -1*min(p), -1*min(q),-1*min(r))
empty_array = 0.0
CALL distance(p, q, empty_array, dist, stilltosort, bigd)
i = 0

!index_big = real(stilltosort, 8)
WHERE (stilltosort)
    index_big = 1.0
ELSEWHERE 
    index_big = 0.0
END WHERE
pts_total = -1.0
next_index = -1.0

WRITE(*,*) "Limits ", z_err_lim1, z_err_lim2

pts_max_index = -1
DO i = 0, BigD - 1
     
!    Go through the entire set of points
    IF (stilltosort(i)) THEN
        ! Found a point which is in the index.  Reset found variable to FALSE"
        found = .FALSE.
        DO j = 0, i
                    !starting at 0, go through all prior points to see if there are examples
                    ! in which the current point is within some margin of the plane described by a previous point            
            IF (stilltosort(j)) THEN
                IF (pts_total(j) > 0.0) THEN
!                    IF (ABS(g2(i) - g2avg(j)) < 0.1) THEN
                    z_err1 = ABS(g2(i) - g2avg(j))
                    z_err2 = ABS(r(i) - (g0avg(j) * p(i) + g1avg(j) * q(i) + g2avg(j)))
!                    z_err2 = 0
!                    IF (ABS(g2(i) - g2avg(j)) < (dist_scale * dist(i))) THEN
                    IF ((z_err1 < (z_err_lim1 * dist(i))) .AND. (z_err2 < (z_err_lim2))) THEN
                        ! Found a candidate point, add this point to the chain:
                        j_start = j
                        k = j
                        pts_counter = 1
                        g0_tot = g0(k)
                        g1_tot = g1(k)
                        g2_tot = g2(k)
                        DO WHILE (next_index(k) .GT. -1 )
                            pts_counter = pts_counter + 1
                            k = INT(next_index(k))
                            g0_tot = g0_tot + g0(k)
                            g1_tot = g1_tot + g1(k)
                            g2_tot = g2_tot + g2(k)
                        ENDDO 
                        ! Done looking through chain
                        IF (k .NE. i) THEN
                            temp =  REAL(i)
                            next_index(k) = temp
                            ! In other words, our search through previous points for one that was already part of a 
                            !  was not successful, so don't create a circular link.
                        ENDIF      
                        !WRITE(*,*) "Calculating average plane values, pts_counter = ", pts_counter
                        g0avg(j_start) = g0_tot / pts_counter
                        g1avg(j_start) = g1_tot / pts_counter
                        g2avg(j_start) = g2_tot / pts_counter
                        pts_total(j_start) = pts_counter
                        IF (pts_max_index > -1) THEN
                            IF (pts_counter > pts_total(pts_max_index)) THEN
                                pts_max_index = j_start
                            ENDIF
                            stilltosort(i) = .FALSE.
                        ELSE  ! If pts_max_index = -1, then it has not yet been assigned a value.
                            pts_max_index = j_start
                        ENDIF
                        found = .TRUE.
                        EXIT
                    ENDIF   
                ENDIF                
            ENDIF  
            
       
        ENDDO   
        IF (.NOT. found) THEN
            ! No suitable string exists, so treat the current candidate as a new string.
            g0avg(i) = g0(i)
            g1avg(i) = g1(i)
            g2avg(i) = g2(i)
            pts_total(i) = 1
        ENDIF  
            
    ENDIF

ENDDO


    stilltosort_r(i) = 0.0
   
stilltosort_r = 0.0
i = pts_max_index
DO WHILE (i > -1)
    stilltosort_r(i) = 1.0
    i = int(next_index(i))            
ENDDO        

temp = SUM(stilltosort_r)            
WRITE(*,*) "From FindLongest string FORTRAN routine" , temp
   
END SUBROUTINE    
        
subroutine rotatepoints(p,q,angle,npoints)
implicit none
integer*4, intent(inout) :: npoints
real(8), intent(inout), dimension(0:npoints - 1) :: p, q
!f2py intent(in,out) :: p, q
real(8) :: angle

REAL(8), dimension(0:npoints - 1) :: p2, q2


    p2 = p * COS(angle) + q * SIN(angle)
    q2 =  -1.0*p * SIN(angle) + q * COS(angle) 
    
    p = p2
    q = q2

end subroutine

subroutine rotatepoints2(p, q, p2, q2, p_off, q_off, angle, npoints)
implicit none
integer*4, intent(inout) :: npoints
real(8), intent(inout), dimension(0:npoints - 1) :: p, q, p2, q2
!f2py intent(in,out) :: p, q, p2, q2
real(8), intent(inout) :: p_off, q_off, angle



    p2 = (p + p_off) * COS(angle) + (q + q_off) * SIN(angle)
    q2 =  -1.0*(p + p_off) * SIN(angle) + (q + q_off) * COS(angle) 
    
    !p = p2
    !q = q2

end subroutine

subroutine unitize(p,q,r, p2, q2, r2)
implicit none
real(8), intent(inout) :: p, q, r
real(8), intent(out) :: p2, q2, r2

real(8) :: length

length = sqrt(p*p + q*q + r*r)

p2 = p / length
q2 = q / length
r2 = r / length

end subroutine

subroutine crossproduct(x0, y0, z0, x1, y1, z1, x2, y2, z2)
implicit none
real(8), intent(inout) :: x0, y0, z0, x1, y1, z1
real(8), intent(out) :: x2, y2, z2

x2 = y0 * z1 - y1 * z0
y2 = x1 * z0 - x0 * z1
z2 = x0 * y1 - x1 * y0

end subroutine