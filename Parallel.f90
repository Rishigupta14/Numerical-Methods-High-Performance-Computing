PROGRAM LINEAR1
USE MPI
IMPLICIT NONE
! ------------------------------------------------------------------ !
! DECLARATION OF VARIABLES !
! ------------------------------------------------------------------ !
REAL::LB,RB,LENGTH     		  	 ! DOMAIN VARIABLES
REAL::T,TFINAL  			 ! TIME VARIABLES
REAL::DX,DT      			 ! STEP SIZES
REAL::U	     				 ! PROPAGATION VELOCITY
REAL::CFL    				 ! CFL CONDITION
REAL::WENO    				 ! WENO DISCRETIZATION
INTEGER::POINTS				 ! GRID POINTS
INTEGER::I,J,K     			 ! LOOP COUNTERS
REAL, ALLOCATABLE, DIMENSION(:)::WX,X 	 ! WRITEABLE/GRID COORDINATES
REAL, ALLOCATABLE, DIMENSION(:,:)::WA    ! WRITEABLE ANALYTICAL SOLUTION
REAL, ALLOCATABLE, DIMENSION(:,:,:)::WF  ! WRITEABLE SOLUTION
REAL, ALLOCATABLE, DIMENSION(:,:,:,:)::F ! NUMERICAL SOLUTION
REAL, ALLOCATABLE, DIMENSION(:,:,:)::FRK ! RUNGE KUTTA OPERATOR
CHARACTER(LEN=40)::F1,F2,FNAME1,FNAME2   ! FILE NAME
INTEGER::B,E      			 ! EDGE BOOLEANS
INTEGER::PPART				 ! POINTS PER PARTITON
INTEGER::iSTART,iEND			 ! PARTITION INDEXS
INTEGER::iERR				 ! INTEGER ERROR
INTEGER::STATUS(MPI_STATUS_SIZE)	 ! INTEGER STATUS
INTEGER::RANK	        		 ! PROCESS ID
INTEGER::SRANK	        		 ! SENDING PROCESS ID
INTEGER::NP				 ! NUMBER OF PROCESSES
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
! MPI INITIALIZATION
! ------------------------------------------------------------------ !
CALL MPI_INIT(iERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, iERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NP, iERR)
CALL MPI_BARRIER(MPI_COMM_WORLD, iERR)
PPART=POINTS/NP
B=0 ! FIRST PROCESSOR INDICATOR (B)
E=0 ! LAST PROCESSOR INDICATOR (E)
IF (RANK.EQ.NP-1) THEN
E=1 
iEND=POINTS-1
ELSE
IF (RANK.EQ.0) THEN
! ONLY ALLOCATED FOR THE MASTER PROCESS		   
ALLOCATE(WX(0:POINTS-1))	   ! WX(COORDS)
ALLOCATE(WA(0:POINTS-1,0:1))       ! WA(COORDS,F0-F1) 
ALLOCATE(WF(0:POINTS-1,0:4,0:1))   ! WF(COORDS,SCHEME,F0-F1)
B=1
END IF
iEND=(RANK+1)*PPART-1
END IF
iSTART=RANK*PPART
! ------------------------------------------------------------------ !
! MEMORY ALLOCATION
! ------------------------------------------------------------------ !
ALLOCATE(X(iSTART-3+B:iEND+2))	  	   ! X(COORDS)
ALLOCATE(FRK(iSTART-3+B:iEND+2,0:2,0:1))   ! FRK(COORDS,RK-ITERATION,F0-F1)
ALLOCATE(F(iSTART-3+B:iEND+2,0:1,0:4,0:1)) ! F(COORDS,T-STEP,SCHEME,F0-F1)
! ------------------------------------------------------------------ !
! GRID GENERATION
! ------------------------------------------------------------------ !
DO I=iSTART-3*(1-B),iEND+2*(1-E)
X(I)=LB+DX*(I)
END DO

! GHOST CELLS FOR THE WENO SCHEMES
X(iSTART-2)=(1-B)*X(iSTART-2)+B*(RB-2*DX)
X(iSTART-1)=(1-B)*X(iSTART-1)+B*(RB-DX)
X(iEND+1)=(1-E)*X(iEND+1)+E*(LB+DX)
X(iEND+2)=(1-E)*X(iEND+2)+E*(LB+2*DX)
! ------------------------------------------------------------------ !
! INITIAL CONDITIONS
! ------------------------------------------------------------------ !
DO I=iSTART-3+B,iEND+2
F(I,0,:,0)=SIN(2*4.0*(atan(1.0d0))*X(I))
F(I,0,:,1)=SIGN(1.0,X(I))
END DO

! GHOST CELLS FOR THE WENO SCHEMES (F1) (FIXED IN TIME)
F(iSTART-2*B:iSTART,:,:,1)=(1-B)*F(iSTART,0,0,1)+B*SIGN(1.0,LB)
F(iEND:iEND+2*E,:,:,1)=(1-E)*F(iEND,0,0,1)+E*SIGN(1.0,RB)

! GHOST CELLS FOR THE RUNGE-KUTTA OPERATOR (FRK-F1) (FIXED IN TIME)
FRK(iSTART-2*B:iSTART,:,1)=F(iSTART,0,0,1)
FRK(iEND:iEND+2*E,:,1)=F(iEND,0,0,1)
! ------------------------------------------------------------------ !
! NUMERICAL ENGINE
! ------------------------------------------------------------------ !
DO K=1,FLOOR((TFINAL/DT)+0.000000001)! TIME LOOP (CONSIDERING A E-9 TOLERANCE)
DO I=iSTART+B,iEND ! SPACE LOOP
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
F(iEND,1,:,1)=F(iEND+E,1,:,1)

! BOUNDARY COMMUNICATIONS
CALL MPI_SEND(F(iEND-2:iEND,1,:,0:(1-E)),3*5*(2-E),MPI_DOUBLE_PRECISION,(1-E)*(RANK+1),1,MPI_COMM_WORLD,iERR)
CALL MPI_SEND(F(iSTART+B:iSTART+1+B,1,:,0:(1-B)),2*5*(2-B),MPI_DOUBLE_PRECISION,(1-B)*(RANK-1)+B*(NP-1),2,MPI_COMM_WORLD,iERR)
CALL MPI_RECV(F(iSTART-3+B:iSTART-1+B,1,:,0:(1-B)),3*5*(2-B),MPI_DOUBLE_PRECISION,(1-B)*(RANK-1)+B*(NP-1),1,MPI_COMM_WORLD,STATUS,iERR)
CALL MPI_RECV(F(iEND+1:iEND+2,1,:,0:(1-E)),2*5*(2-E),MPI_DOUBLE_PRECISION,(1-E)*(RANK+1),2,MPI_COMM_WORLD,STATUS,iERR)

! IMPLEMENTING RUNGE-KUTTA (3TH-ORDER)
FRK(:,0,:)=F(:,1,4,:)
DO J=0,1 ! RUNGE-KUTTA LOOP
DO I=iSTART+B,iEND ! SPACE LOOP
FRK(I,J+1,0)=(3/4.0)*((J*4/9.0)+1-J)*F(I,0,4,0)+(1/4.0)*((J*8/3.0)+1-J)*FRK(I,J,0)+(1/4.0)*((J*8/3.0)+1-J)*DT*WENO(FRK(I-3:I+2,J,0),U,DX) 
FRK(I,J+1,1)=(3/4.0)*((J*4/9.0)+1-J)*F(I,0,4,1)+(1/4.0)*((J*8/3.0)+1-J)*FRK(I,J,1)+(1/4.0)*((J*8/3.0)+1-J)*DT*WENO(FRK(I-3:I+2,J,1),U,DX)
END DO

! BOUNDARY CONDITIONS (F1) (FIXED IN TIME)
FRK(iEND,J+1,1)=FRK(iEND+E,J+1,1)

! BOUNDARY COMMUNICATIONS
CALL MPI_SEND(FRK(iEND-2:iEND,J+1,0:(1-E)),3*(2-E),MPI_DOUBLE_PRECISION,(1-E)*(RANK+1),3,MPI_COMM_WORLD,iERR)
CALL MPI_SEND(FRK(iSTART+B:iSTART+1+B,J+1,0:(1-B)),2*(2-B),MPI_DOUBLE_PRECISION,(1-B)*(RANK-1)+B*(NP-1),4,MPI_COMM_WORLD,iERR)
CALL MPI_RECV(FRK(iSTART-3+B:iSTART-1+B,J+1,0:(1-B)),3*(2-B),MPI_DOUBLE_PRECISION,(1-B)*(RANK-1)+B*(NP-1),3,MPI_COMM_WORLD,STATUS,iERR)
CALL MPI_RECV(FRK(iEND+1:iEND+2,J+1,0:(1-E)),2*(2-E),MPI_DOUBLE_PRECISION,(1-E)*(RANK+1),4,MPI_COMM_WORLD,STATUS,iERR)
END DO
F(:,1,4,:)=FRK(:,2,:)

! ITERATION UPDATE
F(:,0,:,:)=F(:,1,:,:)

END DO
! ------------------------------------------------------------------ !
! GATHER THE SIMULATION DATA FROM THE PROCESSES (MASTER PROCESS: 0 )
! ------------------------------------------------------------------ !
DO I=0,1 ! FUNCTIONS LOOP
DO J=0,4 ! SCHEMES LOOP
CALL MPI_GATHER(F(iSTART:iEND,1,J,I),(iEND-iSTART+1),MPI_DOUBLE_PRECISION,WF(:,J,I),(iEND-iSTART+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iERR)
END DO
END DO
CALL MPI_GATHER(X(iSTART:iEND),(iEND-iSTART+1),MPI_DOUBLE_PRECISION,WX,(iEND-iSTART+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iERR)
CALL MPI_BARRIER(MPI_COMM_WORLD, iERR)
! ------------------------------------------------------------------ !
! WRITE THE SIMULATION DATA TO OUTPUT FILE (ONLY MASTER PROCESS)
! ------------------------------------------------------------------ !
IF (B.EQ.1) THEN
WRITE(F1,FMT='(I10)') FLOOR(TFINAL)

! WRITES THE NUMERICAL SOLUTION FILES
DO J=0,4 ! NUMERICAL SCHEME LOOP
WRITE(F2,FMT='(I10)') J
FNAME1="Results_F0/T"//TRIM(ADJUSTL(F1))//"_NS"//TRIM(ADJUSTL(F2))//"_P.dat"
FNAME2="Results_F1/T"//TRIM(ADJUSTL(F1))//"_NS"//TRIM(ADJUSTL(F2))//"_P.dat"
OPEN(0,FILE=FNAME1, FORM="formatted",STATUS="replace")
OPEN(1,FILE=FNAME2, FORM="formatted",STATUS="replace")
DO I=0,POINTS-1 ! SPACE LOOP
WRITE(0,*)WX(I), WF(I,J,0)
WRITE(1,*)WX(I), WF(I,J,1)
END DO
 CLOSE(0)
 CLOSE(1)
END DO

! WRITES THE ANALYTICAL SOLUTION FILES
OPEN(0,FILE="Results_F0/T"//TRIM(ADJUSTL(F1))//"_AN.dat", FORM="formatted",STATUS="replace")
OPEN(1,FILE="Results_F1/T"//TRIM(ADJUSTL(F1))//"_AN.dat", FORM="formatted",STATUS="replace")
DO I=0,POINTS-1 ! SPACE LOOP
WRITE(0,*)WX(I), SIN(2*4.0*(atan(1.0d0))*(WX(I)-U*TFINAL))
WRITE(1,*)WX(I), SIGN(1.0,WX(I)-U*TFINAL)
END DO
 CLOSE(0)
 CLOSE(1)
END IF
! ------------------------------------------------------------------ !
! CLOSE THE MPI ENVIRONMENT
! ------------------------------------------------------------------ !
CALL MPI_FINALIZE(iERR)
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
ALLOCATE(R(0:1))      			   ! R(Xj+0.5,Xj-0.5)=R(LEFT-RIGHT)
ALLOCATE(A(0:2,0:1)) 			   ! A(SUBSTENCIL,LEFT-RIGHT)
ALLOCATE(P(0:2,0:1)) 		           ! P(SUBSTENCIL,LEFT-RIGHT)
ALLOCATE(IS(0:2,0:1)) 			   ! IS(SUBSTENCIL,LEFT-RIGHT)
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
