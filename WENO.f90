
PROGRAM LINEAR1
IMPLICIT NONE
! ------------------------------------------------------------------ !
! DECLARATION OF VARIABLES
! ------------------------------------------------------------------ !
REAL::LB,RB,LENGTH     		          ! DOMAIN VARIABLES
REAL::T,TFINAL			 	  ! TIME VARIABLES
REAL::DX,DT      			  ! STEP SIZES
REAL::U	     				  ! PROPAGATION VELOCITY
REAL::CFL    				  ! CFL CONDITION
REAL::WENO    				  ! WENO DISCRETIZATION
INTEGER::POINTS				  ! GRID POINTS
INTEGER::I,J,K      			  ! LOOP COUNTERS
REAL, ALLOCATABLE, DIMENSION(:)::X 	  ! GRID COORDINATES
REAL, ALLOCATABLE, DIMENSION (:,:)::AN    ! ANALYTICAL SOLUTION
REAL, ALLOCATABLE, DIMENSION (:,:,:)::FRK ! RUNGE KUTTA OPERATOR
REAL, ALLOCATABLE, DIMENSION (:,:,:,:)::F ! NUMERICAL SOLUTION
CHARACTER(LEN=40)::F1,F2,FNAME1,FNAME2    ! FILE NAMES
! ------------------------------------------------------------------ !
! READ THE SIMULATION PARAMETERS FROM INPUT FILE
! ------------------------------------------------------------------ !
OPEN(0,FILE="Input.dat", FORM="formatted",ACTION="read")
READ(0,*)
READ(0,*)LB
READ(0,*)RB
READ(0,*)U
READ(0,*)TFINAL
READ(0,*)POINTS
READ(0,*)CFL
 CLOSE(0)
LENGTH=RB-LB
DX=LENGTH/(POINTS-1)
DT=CFL*DX/U
! ------------------------------------------------------------------ !
! MEMORY ALLOCATION
! ------------------------------------------------------------------ !
ALLOCATE(X(0:POINTS-1))		     ! X(COORDS)
ALLOCATE(AN(0:POINTS-1,0:1))         ! AN(COORDS,F0-F1)
ALLOCATE(FRK(-2:POINTS+1,0:2,0:1))   ! FRK(COORDS,RK-ITERATION,F0-F1)
ALLOCATE(F(-2:POINTS+1,0:1,0:4,0:1)) ! F(COORDS,T-STEP,SCHEME,F0-F1)
! ------------------------------------------------------------------ !
! GRID GENERATION
! ------------------------------------------------------------------ !
DO I=0,POINTS-1
X(I)=LB+DX*I
END DO
! ------------------------------------------------------------------ !
! INITIAL CONDITIONS & ANALYTICAL SOLUTIONS
! ------------------------------------------------------------------ !
DO I=0,POINTS-1
F(I,0,:,0)=SIN(2*4.0*(atan(1.0d0))*X(I))
F(I,0,:,1)=SIGN(1.0,X(I))
AN(I,0)=SIN(2*4.0*(atan(1.0d0))*(X(I)-U*TFINAL))
AN(I,1)=SIGN(1.0,X(I)-U*TFINAL)
END DO

! GHOST CELLS FOR THE WENO SCHEMES (F0) (EVOLVE IN TIME)
F(-1,0,:,0)=F(POINTS-2,0,:,0)
F(-2,0,:,0)=F(POINTS-3,0,:,0)
F(POINTS,0,:,0)=F(1,0,:,0)
F(POINTS+1,0,:,0)=F(2,0,:,0)

! GHOST CELLS FOR THE WENO SCHEMES (F1) (FIXED IN TIME)
F(-2:0,:,:,1)=F(0,0,0,1)
F(POINTS-1:POINTS+1,:,:,1)=F(POINTS-1,0,0,1)

! GHOST CELLS FOR THE RUNGE-KUTTA OPERATOR (FRK-F1) (FIXED IN TIME)
FRK(-2:0,:,1)=F(0,0,0,1)
FRK(POINTS-1:POINTS+1,:,1)=F(POINTS-1,0,0,1)
! ------------------------------------------------------------------ !
! NUMERICAL ENGINE
! ------------------------------------------------------------------ !
DO K=1,FLOOR((TFINAL/DT)+0.000000001)! TIME LOOP (CONSIDERING A E-9 TOLERANCE)
DO I=1,POINTS-1 ! SPACE LOOP
! FORWARD IN TIME - CENTERED IN SPACE (FTCS)
F(I,1,0,:)=F(I,0,0,:)-(0.5*U*(F(I+1,0,0,:)-F(I-1,0,0,:))*DT/DX)
! FIRST-ORDER UPWIND
F(I,1,1,:)=F(I,0,1,:)-(U*(F(I,0,1,:)-F(I-1,0,1,:))*DT/DX)
! LAX-FRIEDRICHS
F(I,1,2,:)=0.5*(F(I+1,0,2,:)+F(I-1,0,2,:))-(0.5*U*(F(I+1,0,2,:)-F(I-1,0,2,:))*DT/DX)
! LAX-WENDROFF
F(I,1,3,:)=F(I,0,3,:)-(0.5*U*(F(I+1,0,3,:)-F(I-1,0,3,:))*DT/DX)+0.5*((U*DT/DX)**2)*(F(I+1,0,3,:)-2*F(I,0,3,:)+F(I-1,0,3,:))
! WENO
F(I,1,4,0)=F(I,0,4,0)+DT*WENO(F(I-3:I+2,0,4,0),U,DX)
F(I,1,4,1)=F(I,0,4,1)+DT*WENO(F(I-3:I+2,0,4,1),U,DX)
END DO

! BOUNDARY CONDITIONS (F1) (FIXED IN TIME)
F(POINTS-1,1,:,1)=F(POINTS,1,:,1)

! GHOST CELLS FOR THE WENO SCHEMES (F0) (EVOLVE IN TIME)
F(-2:0,1,:,0)=F(POINTS-3:POINTS-1,1,:,0)
F(POINTS:POINTS+1,1,:,0)=F(1:2,1,:,0)

! IMPLEMENTING RUNGE-KUTTA (3TH-ORDER)
FRK(:,0,:)=F(:,1,4,:)
DO J=0,1 ! RUNGE-KUTTA LOOP
DO I=1,POINTS-1 ! SPACE LOOP
FRK(I,J+1,0)=(3/4.0)*((J*4/9.0)+1-J)*F(I,0,4,0)+(1/4.0)*((J*8/3.0)+1-J)*FRK(I,J,0)+(1/4.0)*((J*8/3.0)+1-J)*DT*WENO(FRK(I-3:I+2,J,0),U,DX)
FRK(I,J+1,1)=(3/4.0)*((J*4/9.0)+1-J)*F(I,0,4,1)+(1/4.0)*((J*8/3.0)+1-J)*FRK(I,J,1)+(1/4.0)*((J*8/3.0)+1-J)*DT*WENO(FRK(I-3:I+2,J,1),U,DX)
END DO

! BOUNDARY CONDITIONS (F1) (FIXED IN TIME)
FRK(POINTS-1,J+1,1)=FRK(POINTS,J+1,1)

! GHOST CELLS FOR THE WENO SCHEMES (F0) (EVOLVE IN TIME)
FRK(-2:0,J+1,0)=FRK(POINTS-3:POINTS-1,J+1,0)
FRK(POINTS:POINTS+1,J+1,0)=FRK(1:2,J+1,0)
END DO
F(:,1,4,:)=FRK(:,2,:)

! ITERATION UPDATE
F(:,0,:,:)=F(:,1,:,:)

END DO
! ------------------------------------------------------------------ !
! WRITE THE SIMULATION DATA TO OUTPUT FILE
! ------------------------------------------------------------------ !
WRITE(F1,FMT='(I10)') FLOOR(TFINAL)

! WRITES THE NUMERICAL SOLUTION FILES
DO J=0,4 ! NUMERICAL SCHEME LOOP
WRITE(F2,FMT='(I10)') J
FNAME1="Results_F0/T"//TRIM(ADJUSTL(F1))//"_NS"//TRIM(ADJUSTL(F2))//"_S.dat"
FNAME2="Results_F1/T"//TRIM(ADJUSTL(F1))//"_NS"//TRIM(ADJUSTL(F2))//"_S.dat"
OPEN(0,FILE=FNAME1, FORM="formatted",STATUS="replace")
OPEN(1,FILE=FNAME2, FORM="formatted",STATUS="replace")
DO I=0,POINTS-1 ! SPACE LOOP
WRITE(0,*)X(I), F(I,1,J,0)
WRITE(1,*)X(I), F(I,1,J,1)
END DO
 CLOSE(0)
 CLOSE(1)
END DO

! WRITES THE ANALYTICAL SOLUTION FILES
OPEN(0,FILE="Results_F0/T"//TRIM(ADJUSTL(F1))//"_AN.dat", FORM="formatted",STATUS="replace")
OPEN(1,FILE="Results_F1/T"//TRIM(ADJUSTL(F1))//"_AN.dat", FORM="formatted",STATUS="replace")
DO I=0,POINTS-1 ! SPACE LOOP
WRITE(0,*)X(I), AN(I,0)
WRITE(1,*)X(I), AN(I,1)
END DO
 CLOSE(0)
 CLOSE(1)

END PROGRAM
! ------------------------------------------------------------------ !
! WENO FUNCTION
! ------------------------------------------------------------------ !
REAL FUNCTION WENO(WF,U,DX)
IMPLICIT NONE
! ------------------------------------------------------------------ !
REAL, DIMENSION(6), INTENT(IN) :: WF	   ! GLOBAL STENCIL
REAL, INTENT(IN) :: U		 	   ! PROPAGATION VELOCITY
REAL, INTENT(IN) :: DX			   ! SPATIAL STEP
INTEGER :: I,J				   ! LOOP COUNTERS
REAL,PARAMETER :: e=0.000001		   ! NON-ZERO PARAMETER
REAL, ALLOCATABLE, DIMENSION (:)::R	   ! RECONSTRUCTION POLYNOMIAL
REAL, ALLOCATABLE, DIMENSION (:,:)::A,P,IS ! WEIGHTS, STENCIL POLYNOMIALS, AND SMOOTHNES INDICATORS
! ------------------------------------------------------------------ !
ALLOCATE(R(0:1))                           ! R(Xj+0.5,Xj-0.5)=R(LEFT-RIGHT)
ALLOCATE(A(0:2,0:1))                       ! A(SUBSTENCIL,LEFT-RIGHT)
ALLOCATE(P(0:2,0:1))                       ! P(SUBSTENCIL,LEFT-RIGHT)
ALLOCATE(IS(0:2,0:1))                      ! IS(SUBSTENCIL,LEFT-RIGHT)
! ------------------------------------------------------------------ !
DO J=0,1 ! LEFT-RIGHT LOOP
DO I=0,2 ! SUBSTENCIL LOOP
P(I,J)=0.5*(((WF(I+4-J)-2*WF(I+3-J)+WF(I+2-J))*(1.5-I)**2)+((WF(I+4-J)-WF(I+2-J))*(1.5-I))+2*WF(I+3-J)-((WF(I+4-J)-2*WF(I+3-J)+WF(I+2-J))/12))
IS(I,J)=0.5*(((WF(I+3-J)-WF(I+2-J))**2)+((WF(I+4-J)-WF(I+3-J))**2))+((WF(I+4-J)-2*WF(I+3-J)+WF(I+2-J))**2)
END DO
A(0,J)=1/(12.0*((e+IS(0,J))**3))
A(1,J)=1/(2.0*((e+IS(1,J))**3))
A(2,J)=1/(4.0*((e+IS(2,J))**3))
R(J)=(A(0,J)*P(0,J)+A(1,J)*P(1,J)+A(2,J)*P(2,J))/(A(0,J)+A(1,J)+A(2,J))
END DO
WENO=-U*(R(0)-R(1))/DX

END FUNCTION WENO
