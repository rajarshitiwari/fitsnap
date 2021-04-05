MODULE RANDOM_NUMBERS_CLASS
  IMPLICIT NONE

CONTAINS

  SUBROUTINE INIT_RANDOM_SEED()
    USE ISO_FORTRAN_ENV, ONLY: INT64
    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: SEED(:)
    INTEGER :: I, N, UN, ISTAT, DT(8), PID
    INTEGER(INT64) :: T

    CALL RANDOM_SEED(SIZE = N)
    ALLOCATE(SEED(N))
    ! FIRST TRY IF THE OS PROVIDES A RANDOM NUMBER GENERATOR
    OPEN(NEWUNIT=UN, FILE="/DEV/URANDOM", ACCESS="STREAM", &
         FORM="UNFORMATTED", ACTION="READ", STATUS="OLD", IOSTAT=ISTAT)
    IF (ISTAT == 0) THEN
       READ(UN) SEED
       CLOSE(UN)
    ELSE
       ! FALLBACK TO XOR:ING THE CURRENT TIME AND PID. THE PID IS
       ! USEFUL IN CASE ONE LAUNCHES MULTIPLE INSTANCES OF THE SAME
       ! PROGRAM IN PARALLEL.
       CALL SYSTEM_CLOCK(T)
       IF (T == 0) THEN
          CALL DATE_AND_TIME(VALUES=DT)
          T = (DT(1) - 1970) * 365_INT64 * 24 * 60 * 60 * 1000 &
               + DT(2) * 31_INT64 * 24 * 60 * 60 * 1000 &
               + DT(3) * 24_INT64 * 60 * 60 * 1000 &
               + DT(5) * 60 * 60 * 1000 &
               + DT(6) * 60 * 1000 + DT(7) * 1000 &
               + DT(8)
       END IF
       PID = GETPID()
       T = IEOR(T, INT(PID, KIND(T)))
       DO I = 1, N
          SEED(I) = LCG(T)
       END DO
    END IF
    CALL RANDOM_SEED(PUT=SEED)
  CONTAINS
    ! THIS SIMPLE PRNG MIGHT NOT BE GOOD ENOUGH FOR REAL WORK, BUT IS
    ! SUFFICIENT FOR SEEDING A BETTER PRNG.
    FUNCTION LCG(S)
      INTEGER :: LCG
      INTEGER(INT64) :: S
      IF (S == 0) THEN
         S = 104729
      ELSE
         S = MOD(S, 4294967296_INT64)
      END IF
      S = MOD(S * 279470273_INT64, 4294967291_INT64)
      LCG = INT(MOD(S, INT(HUGE(0), INT64)), KIND(0))
    END FUNCTION LCG
  END SUBROUTINE INIT_RANDOM_SEED

END MODULE RANDOM_NUMBERS_CLASS

MODULE FIT_LAMMPS_CLASS
  IMPLICIT NONE

  TYPE DATA_FILE
     INTEGER                  :: FRAMES
     INTEGER                  :: NATS
     CHARACTER (LEN=100)      :: INP_DATA
     CHARACTER (LEN=100)      :: INP_ENER
     CHARACTER (LEN=100)      :: INP_FORCES
     CHARACTER (LEN=100)      :: INP_TRAJ
     REAL(8), ALLOCATABLE     :: ENER(:)
     CHARACTER(LEN=5), ALLOCATABLE :: LABEL(:)
     REAL(8), ALLOCATABLE     :: X(:,:)
     REAL(8), ALLOCATABLE     :: FX(:,:)
     REAL(8), ALLOCATABLE     :: FY(:,:)
     REAL(8), ALLOCATABLE     :: FZ(:,:)
     REAL(8)                  :: WEIGHT
  END TYPE DATA_FILE

  TYPE SYSTEM
     TYPE(DATA_FILE), ALLOCATABLE   :: DATA(:)
     INTEGER                        :: NDATA
     CHARACTER (LEN=100)            :: INP
     CHARACTER (LEN=100)            :: INP_FIX
     CHARACTER (LEN=100)            :: INP_FIT
     INTEGER                        :: TOT_FRAMES
     INTEGER                        :: NPAR2FIT
   CONTAINS
     PROCEDURE                      :: READ_SYS
     !PROCEDURE                      :: DELETE
  END TYPE SYSTEM

  TYPE KERNEL_KIND
     INTEGER :: NENVS=0
     REAL(8) :: SIGMA=0
     REAL(8), ALLOCATABLE  :: B(:,:)
   CONTAINS
     PROCEDURE :: DEALLOC => DEALLOC_KERNEL_KIND
  END TYPE KERNEL_KIND

  TYPE KERNEL_GLOBAL
     INTEGER :: NKINDS=0
     TYPE(KERNEL_KIND), ALLOCATABLE :: K(:)
   CONTAINS
     PROCEDURE :: DEALLOC => DEALLOC_KERNEL_GLOBAL
  END TYPE KERNEL_GLOBAL

CONTAINS

  SUBROUTINE DEALLOC_KERNEL_GLOBAL(THIS)
    IMPLICIT NONE
    CLASS(KERNEL_GLOBAL) :: THIS
    INTEGER              :: I
    IF ( ALLOCATED(THIS % K) ) THEN
       DO I = 1, THIS % NKINDS
          CALL THIS % K(I) % DEALLOC()
       ENDDO
       DEALLOCATE(THIS % K)
    END IF
    THIS % NKINDS = 0
    RETURN
  END SUBROUTINE DEALLOC_KERNEL_GLOBAL

  SUBROUTINE DEALLOC_KERNEL_KIND(THIS)
    IMPLICIT NONE
    CLASS(KERNEL_KIND)    :: THIS
    IF( ALLOCATED(THIS%B) ) DEALLOCATE(THIS%B)
    THIS%NENVS = 0
    THIS%SIGMA = 0
    RETURN
  END SUBROUTINE DEALLOC_KERNEL_KIND

  SUBROUTINE READ_SYS(SYS, DATA_FILE, FIT_ENER, FIT_FORCES)
    IMPLICIT NONE
    CLASS(SYSTEM)           :: SYS
    INTEGER                 :: I,L,V,K,M,N
    CHARACTER(LEN=100)      :: LABEL,DATA_FILE
    LOGICAL                 :: FIT_FORCES,FIT_ENER

    IF (ALLOCATED(SYS%DATA)) THEN
       DO I = 1, SYS%NDATA
          IF (ALLOCATED(SYS%DATA(I)%X)   ) DEALLOCATE(SYS%DATA(I)%X   )
          IF (ALLOCATED(SYS%DATA(I)%ENER)) DEALLOCATE(SYS%DATA(I)%ENER)
          IF (ALLOCATED(SYS%DATA(I)%FX)  ) DEALLOCATE(SYS%DATA(I)%FX  )
          IF (ALLOCATED(SYS%DATA(I)%FY)  ) DEALLOCATE(SYS%DATA(I)%FY  )
          IF (ALLOCATED(SYS%DATA(I)%FZ)  ) DEALLOCATE(SYS%DATA(I)%FZ  )
       END DO
       DEALLOCATE(SYS%DATA)
    END IF

    OPEN(8,FILE=TRIM(DATA_FILE))

    SYS%TOT_FRAMES = 0
    SYS%NPAR2FIT = 0

    READ(8, *) SYS%NDATA

    ALLOCATE(SYS%DATA(SYS%NDATA))

    DO I = 1, SYS%NDATA
       
       IF (FIT_FORCES .AND. FIT_ENER) THEN
          READ(8, *) SYS%DATA(I)%INP_DATA, SYS%DATA(I)%FRAMES, SYS%DATA(I)%INP_TRAJ, SYS%DATA(I)%INP_ENER,&
               SYS%DATA(I)%INP_FORCES,SYS%DATA(I)%WEIGHT
       END IF
       IF (FIT_ENER .AND. (.NOT. FIT_FORCES) ) THEN
          READ(8, *) SYS%DATA(I)%INP_DATA, SYS%DATA(I)%FRAMES, SYS%DATA(I)%INP_TRAJ, &
               SYS%DATA(I)%INP_ENER, SYS%DATA(I)%WEIGHT
       END IF
       IF ((.NOT. FIT_ENER) .AND. FIT_FORCES ) THEN
          READ(8, *) SYS%DATA(I)%INP_DATA, SYS%DATA(I)%FRAMES, SYS%DATA(I)%INP_TRAJ, &
               SYS%DATA(I)%INP_FORCES, SYS%DATA(I)%WEIGHT
       END IF

       OPEN(12, FILE=SYS%DATA(I)%INP_TRAJ)

       IF(FIT_ENER)THEN
          ALLOCATE(SYS%DATA(I)%ENER(SYS%DATA(I)%FRAMES))
          OPEN(13, FILE=SYS%DATA(I)%INP_ENER)
       END IF

       IF (FIT_FORCES) THEN
          OPEN(14, FILE=SYS%DATA(I)%INP_FORCES)
       END IF

       DO L = 1, SYS%DATA(I)%FRAMES

          SYS%TOT_FRAMES = SYS%TOT_FRAMES + 1

          READ(12, *) SYS%DATA(I)%NATS
          READ(12, *)
          IF (FIT_FORCES) THEN
             READ(14, *)
             READ(14, *)
          END IF

          IF(.NOT. ALLOCATED(SYS%DATA(I)%LABEL)) THEN
             ALLOCATE(SYS%DATA(I)%LABEL(SYS%DATA(I)%NATS))
          END IF
          IF (.NOT. ALLOCATED(SYS%DATA(I)%X)) THEN
             ALLOCATE(SYS%DATA(I)%X(SYS%DATA(I)%FRAMES, 3*SYS%DATA(I)%NATS))
          END IF
          IF(.NOT. ALLOCATED(SYS%DATA(I)%FX) .AND. FIT_FORCES)THEN
             ALLOCATE(SYS%DATA(I)%FX(SYS%DATA(I)%FRAMES, SYS%DATA(I)%NATS))
          END IF
          IF(.NOT. ALLOCATED(SYS%DATA(I)%FY) .AND. FIT_FORCES)THEN
             ALLOCATE(SYS%DATA(I)%FY(SYS%DATA(I)%FRAMES, SYS%DATA(I)%NATS))
          END IF
          IF(.NOT. ALLOCATED(SYS%DATA(I)%FZ) .AND. FIT_FORCES)THEN
             ALLOCATE(SYS%DATA(I)%FZ(SYS%DATA(I)%FRAMES, SYS%DATA(I)%NATS))
          END IF

          V = 1
          DO K = 1, SYS%DATA(I)%NATS
             !
             READ(12, *) SYS%DATA(I)%LABEL(K), SYS%DATA(I)%X(L,V), SYS%DATA(I)%X(L,V+1), SYS%DATA(I)%X(L,V+2)
             !
             IF(FIT_FORCES)THEN
                READ(14, *) SYS%DATA(I)%FX(L,K), SYS%DATA(I)%FY(L,K), SYS%DATA(I)%FZ(L,K)
                SYS%NPAR2FIT = SYS%NPAR2FIT + 3
             END IF
             !
             V = V + 3
             !
          END DO

          IF (FIT_ENER) THEN
             READ(13, *) SYS%DATA(I)%ENER(L)
             SYS%NPAR2FIT = SYS%NPAR2FIT + 1
          END IF

       END DO

       CLOSE(12)
       IF(FIT_ENER) CLOSE(13)
       IF(FIT_FORCES) CLOSE(14)

    END DO

    CLOSE(8)

    RETURN
  END SUBROUTINE READ_SYS


END MODULE FIT_LAMMPS_CLASS

MODULE COMMON_VAR
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_PTR, C_INT
  USE FIT_LAMMPS_CLASS
  !
  TYPE(SYSTEM)              :: SYS
  TYPE (C_PTR), ALLOCATABLE :: LMP(:)
  TYPE(KERNEL_GLOBAL)       :: KERNEL
  LOGICAL                   :: FIT_FORCES = .FALSE., PCA = .FALSE., CS=.FALSE., FIT_ENER = .FALSE.
  LOGICAL                   :: FIT_TENSORS = .TRUE.
  LOGICAL                   :: SKIP_FIT = .FALSE., PRINT_BI = .FALSE., REFINE= .FALSE., METRIX = .TRUE.
  REAL(8)                   :: CM_VAL = 1.0D0, REFINE_TEMP, THR_KERNEL = 0.5D0
  INTEGER                   :: REFINE_MAXITER, ITER, TENS_ORDER, E0CS = 0
  !
END MODULE COMMON_VAR

MODULE LAPACK_INVERSE
  IMPLICIT NONE

CONTAINS

  SUBROUTINE MAT_INV(X, N)
    IMPLICIT NONE
    INTEGER                              :: INFO,LWORK,N,IALLOC
    INTEGER, DIMENSION(N)                :: IPIV
    REAL(8), DIMENSION(N,N)              :: X_INV
    COMPLEX(8), DIMENSION(N,N)           :: CX_INV
    CLASS(*), INTENT(IN), DIMENSION(N,N) :: X
    REAL(8), ALLOCATABLE, DIMENSION(:)   :: WORK

    SELECT TYPE(X)

    TYPE IS (COMPLEX(8))

       CX_INV = X
       
       LWORK = (N)**2
       ALLOCATE(WORK(LWORK), STAT=IALLOC)
       CALL ZGETRF(N, N, CX_INV, N, IPIV, INFO)
       CALL ZGETRI(N, CX_INV, N, IPIV, WORK, LWORK, INFO)

       X = CX_INV

    TYPE IS (REAL(8))

       X_INV = X

       LWORK = (N)**2
       ALLOCATE(WORK(LWORK), STAT=IALLOC)
       CALL DGETRF(N, N, X_INV, N, IPIV, INFO)
       CALL DGETRI(N, X_INV, N, IPIV, WORK, LWORK, INFO)

       X = X_INV

    END SELECT

    RETURN
  END SUBROUTINE MAT_INV

END MODULE LAPACK_INVERSE

MODULE LAPACK_DIAG_SIMM
  IMPLICIT NONE

CONTAINS

  SUBROUTINE NEW_DIAG(N, A, W)
    IMPLICIT NONE
    INTEGER                               :: I, J, INFO, S, K, N, IALLOC
    INTEGER                               :: L, INF, INFR
    CLASS(*), DIMENSION(N,N)              :: A
    COMPLEX(8), DIMENSION(N*(N+1)/2)      :: AP
    REAL(8),  DIMENSION(N)                :: W
    REAL(8), ALLOCATABLE, DIMENSION(:)    :: WORK
    COMPLEX(8), ALLOCATABLE, DIMENSION(:) :: CWORK


    SELECT TYPE(A)
       !
    TYPE IS (REAL(8))
       L = (3 * N - 1)
       ALLOCATE(WORK(L), STAT=IALLOC)
       CALL DSYEV('V', 'U', N, A, N, W, WORK, L, INFR)
       DEALLOCATE(WORK)
       IF (INFR .NE. 0) THEN
          WRITE(*, *) 'dsyev diagonalization failed'
          FLUSH(6)
          STOP
       END IF
       !
    TYPE IS (COMPLEX(8))
       !
       S = 1
       DO J = 1, N
          DO I = 1, J
             AP(S) = A(I,J)
             S = S + 1
          END DO
       END DO
       !
       L = (2 * N - 1) + 1000
       ALLOCATE(CWORK(L), STAT=IALLOC)
       K = (3 * N - 2) + 1000
       ALLOCATE(WORK(K), STAT=IALLOC)
       CALL ZHPEV('V', 'U', N, AP, W, A, N, CWORK, WORK, INF)
       DEALLOCATE(CWORK)
       DEALLOCATE(WORK)
       IF (INF .NE. 0) THEN
          WRITE(*,*) 'zhpev diagonalization failed'
          STOP
       END IF
       !
    END SELECT

    RETURN
  END SUBROUTINE NEW_DIAG

END MODULE LAPACK_DIAG_SIMM

MODULE LAPACK_DIAG_ASIMM
  IMPLICIT NONE

CONTAINS

  SUBROUTINE NEW_DIAG2(N, A, W) ! NO INTENT?
    IMPLICIT NONE
    INTEGER                             :: I, J, INFO, S, K, N, IALLOC
    INTEGER                             :: L, INF, INFR
    CLASS(*), DIMENSION(N,N)            :: A
    COMPLEX(8), DIMENSION(N*(N+1)/2)    :: AP
    COMPLEX(8), DIMENSION(N,N)          :: VRI, VLI
    COMPLEX(8), DIMENSION(N)            :: W
    REAL(8), DIMENSION(N)               :: WR, WI
    REAL(8), DIMENSION(N,N)             :: VL, VR
    REAL(8), ALLOCATABLE,  DIMENSION(:) :: WORK
    COMPLEX(8), ALLOCATABLE,  DIMENSION(:) :: CWORK

    SELECT TYPE(A)
       !
    TYPE IS (REAL(8))
       !
       L = 4 * N
       ALLOCATE(WORK(L), STAT=IALLOC)
       CALL DGEEV('V','V', N, A, N, WR, WI, VL, N, VR, N, WORK, L, INFR)
       DEALLOCATE(WORK)
       IF (INFR .NE. 0) THEN
          WRITE(*,*) 'dgeev diagonalization failed'
          FLUSH(6)
          STOP
       END IF
       !
       DO J = 1, N
          W(J) = CMPLX(WR(J), WI(J), 8)
       END DO
       !
       A = VR
       !
    TYPE IS (COMPLEX(8))

       L = 6 * N
       ALLOCATE(CWORK(L), STAT=IALLOC)
       ALLOCATE(WORK(L), STAT=IALLOC)
       CALL ZGEEV('V', 'V', N, A, N, W, VLI, N, VRI, N, CWORK, L, WORK, INFR)
       DEALLOCATE(CWORK)
       DEALLOCATE(WORK)
       IF (INF .NE. 0) THEN
          WRITE(*,*) 'zgeev diagonalization failed'
       END IF
       A = VLI
    END SELECT

    RETURN
  END SUBROUTINE NEW_DIAG2

END MODULE LAPACK_DIAG_ASIMM

MODULE LISTS_CLASS
  IMPLICIT NONE

  TYPE :: LIST_NODE
     CLASS(*), POINTER      :: KEY => NULL()
     CLASS(LIST_NODE), POINTER  :: NEXT => NULL()
     CLASS(LIST_NODE), POINTER  :: PREV => NULL()
  END TYPE LIST_NODE
  
  TYPE LIST
     CLASS(LIST_NODE), POINTER :: HEAD => NULL()
     CLASS(LIST_NODE), POINTER :: TAIL => NULL()
     CLASS(LIST_NODE), POINTER :: NODE => NULL()
     INTEGER                   :: NELEM
   CONTAINS
     PROCEDURE                 :: INIT => INIT_LIST
     PROCEDURE                 :: SKIP => SKIP_NODE
     PROCEDURE                 :: REW  => REWIND_NODE
     PROCEDURE                 :: RM => REMOVE_NODE
     PROCEDURE                 :: REBOOT => REBOOT_LIST
     PROCEDURE                 :: DELETE => DESTROY_LIST
     PROCEDURE                 :: ADD_NODE
     PROCEDURE                 :: RD_DBL_NODE
     PROCEDURE                 :: RD_CMPLX_NODE
     PROCEDURE                 :: RD_INT_NODE
     GENERIC                   :: RD_VAL => RD_INT_NODE, RD_DBL_NODE, RD_CMPLX_NODE
  END TYPE LIST

CONTAINS



!!!!!   LISTS GENERAL FUNCTIONS

  SUBROUTINE SKIP_NODE(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)  :: THIS_LIST
    IF( ASSOCIATED(THIS_LIST%NODE%NEXT) ) THEN
       THIS_LIST%NODE => THIS_LIST%NODE%NEXT
    ELSE
       THIS_LIST%NODE=>THIS_LIST%TAIL
    END IF
    RETURN
  END SUBROUTINE SKIP_NODE

  SUBROUTINE DESTROY_LIST(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST) :: THIS_LIST
    DO WHILE (THIS_LIST%NELEM .GT. 0)
       THIS_LIST%NODE => THIS_LIST%HEAD
       CALL THIS_LIST%RM()
    END DO
    RETURN
  END SUBROUTINE DESTROY_LIST

  SUBROUTINE REMOVE_NODE(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)  :: THIS_LIST
    CLASS(LIST_NODE), POINTER :: TMP_NODE

    IF (THIS_LIST%NELEM .EQ. 0) RETURN

    IF (THIS_LIST%NELEM .EQ. 1)THEN

       IF (ASSOCIATED(THIS_LIST%NODE%KEY)) THEN
          DEALLOCATE(THIS_LIST%NODE%KEY)
       END IF
       IF (ASSOCIATED(THIS_LIST%NODE)) THEN
          DEALLOCATE(THIS_LIST%NODE)
       END IF
       THIS_LIST%NODE => NULL()
       THIS_LIST%HEAD => NULL()
       THIS_LIST%TAIL => NULL()

       THIS_LIST%NELEM = THIS_LIST%NELEM - 1
       
    ELSE

       IF( .NOT. ASSOCIATED(THIS_LIST%NODE%NEXT) .AND. &
            ASSOCIATED(THIS_LIST%NODE%PREV) ) THEN
          THIS_LIST%TAIL => THIS_LIST%NODE%PREV
          TMP_NODE       => THIS_LIST%NODE
          THIS_LIST%NODE => THIS_LIST%TAIL
          THIS_LIST%TAIL%NEXT => NULL()
          IF (ASSOCIATED(TMP_NODE%KEY)) THEN
             DEALLOCATE(TMP_NODE%KEY)
             TMP_NODE%KEY => NULL()
          END IF
          IF (ASSOCIATED(TMP_NODE)) THEN
             DEALLOCATE(TMP_NODE)
             TMP_NODE => NULL()
          END IF
          THIS_LIST%NELEM = THIS_LIST%NELEM - 1
       END IF

       IF( .NOT. ASSOCIATED(THIS_LIST%NODE%PREV) .AND. &
            ASSOCIATED(THIS_LIST%NODE%NEXT) ) THEN
          THIS_LIST%HEAD => THIS_LIST%NODE%NEXT
          TMP_NODE       => THIS_LIST%NODE
          THIS_LIST%NODE => THIS_LIST%HEAD
          THIS_LIST%HEAD%PREV => NULL()
          IF (ASSOCIATED(TMP_NODE%KEY)) THEN
             DEALLOCATE(TMP_NODE%KEY)
             TMP_NODE%KEY => NULL()
          END IF
          IF (ASSOCIATED(TMP_NODE)) THEN
             DEALLOCATE(TMP_NODE)
             TMP_NODE => NULL()
          END IF
          THIS_LIST%NELEM = THIS_LIST%NELEM - 1
       END IF

       IF( ASSOCIATED(THIS_LIST%NODE%NEXT) .AND. &
            ASSOCIATED(THIS_LIST%NODE%PREV) ) THEN
          THIS_LIST%NODE%PREV%NEXT => THIS_LIST%NODE%NEXT
          THIS_LIST%NODE%NEXT%PREV => THIS_LIST%NODE%PREV
          TMP_NODE => THIS_LIST%NODE
          TMP_NODE%KEY => THIS_LIST%NODE%KEY
          THIS_LIST%NODE => THIS_LIST%NODE%NEXT
          IF (ASSOCIATED(TMP_NODE%KEY)) THEN
             DEALLOCATE(TMP_NODE%KEY)
             TMP_NODE%KEY => NULL()
          END IF
          IF(ASSOCIATED(TMP_NODE)) THEN
             DEALLOCATE(TMP_NODE)
             TMP_NODE => NULL()
          END IF
          THIS_LIST%NELEM = THIS_LIST%NELEM - 1
       END IF

    END IF

    RETURN
  END SUBROUTINE REMOVE_NODE

  SUBROUTINE REWIND_NODE(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)  :: THIS_LIST
    IF (ASSOCIATED(THIS_LIST%NODE%PREV)) THEN
       THIS_LIST%NODE => THIS_LIST%NODE%PREV
    ELSE
       THIS_LIST%NODE => THIS_LIST%HEAD
    END IF
    RETURN
  END SUBROUTINE REWIND_NODE

  SUBROUTINE REBOOT_LIST(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)   :: THIS_LIST
    IF ( ASSOCIATED(THIS_LIST%NODE) &
         .AND. ASSOCIATED(THIS_LIST%HEAD) ) THEN
       THIS_LIST%NODE => THIS_LIST%HEAD
    END IF
    RETURN
  END SUBROUTINE REBOOT_LIST

  SUBROUTINE LAST_NODE_LIST(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)   :: THIS_LIST
    THIS_LIST%NODE => THIS_LIST%TAIL
    RETURN
  END SUBROUTINE LAST_NODE_LIST

  SUBROUTINE INIT_LIST(THIS_LIST)
    IMPLICIT NONE
    CLASS(LIST)    :: THIS_LIST
    THIS_LIST%NELEM = 0
    RETURN
  END SUBROUTINE INIT_LIST

  SUBROUTINE ADD_NODE(THIS_LIST, VAL)
    IMPLICIT NONE
    CLASS(LIST)               :: THIS_LIST
    CLASS(*), POINTER         :: ARROW
    CLASS(*), OPTIONAL        :: VAL
    CLASS(LIST_NODE), POINTER :: TMP_NODE

    IF (THIS_LIST%NELEM .EQ. 0) THEN
       ALLOCATE(THIS_LIST%HEAD)
       ALLOCATE(THIS_LIST%NODE)
       IF (PRESENT(VAL)) THEN
          SELECT TYPE (VAL)
          TYPE IS (INTEGER)
             ALLOCATE(INTEGER::THIS_LIST%HEAD%KEY)
             ALLOCATE(INTEGER::THIS_LIST%NODE%KEY)
          TYPE IS (REAL(8))
             ALLOCATE(REAL(8)::THIS_LIST%HEAD%KEY)
             ALLOCATE(REAL(8)::THIS_LIST%NODE%KEY)
          TYPE IS (COMPLEX(8))
             ALLOCATE(COMPLEX(8)::THIS_LIST%HEAD%KEY)
             ALLOCATE(COMPLEX(8)::THIS_LIST%NODE%KEY)
          END SELECT
       END IF
       THIS_LIST%NODE => THIS_LIST%HEAD
       THIS_LIST%TAIL => THIS_LIST%HEAD
       THIS_LIST%NELEM = THIS_LIST%NELEM + 1
       TMP_NODE => THIS_LIST%HEAD
    ELSE
       ALLOCATE(TMP_NODE)
       IF (PRESENT(VAL)) THEN
          SELECT TYPE (VAL)
          TYPE IS (INTEGER)
             ALLOCATE(INTEGER::TMP_NODE%KEY)
          TYPE IS (REAL(8))
             ALLOCATE(REAL(8)::TMP_NODE%KEY)
          TYPE IS (COMPLEX(8))
             ALLOCATE(COMPLEX(8)::TMP_NODE%KEY)
          END SELECT
       END IF
       THIS_LIST%TAIL%NEXT => TMP_NODE
       TMP_NODE%PREV => THIS_LIST%TAIL
       THIS_LIST%TAIL => TMP_NODE
       THIS_LIST%NELEM = THIS_LIST%NELEM + 1
    END IF
    
    IF ( PRESENT(VAL) ) THEN
       SELECT TYPE (VAL)
          
       TYPE IS (INTEGER)
          SELECT TYPE (ARROW => TMP_NODE%KEY)
          TYPE IS (INTEGER)
             ARROW = VAL
          END SELECT
       TYPE IS (REAL(8))
          SELECT TYPE (ARROW => TMP_NODE%KEY)
          TYPE IS (REAL(8))
             ARROW = VAL
          END SELECT
       TYPE IS (COMPLEX(8))
          SELECT TYPE (ARROW => TMP_NODE%KEY)
          TYPE IS (COMPLEX(8))
             ARROW = VAL
          END SELECT
          
       END SELECT
    END IF
    
    TMP_NODE => NULL()
    
    RETURN
  END SUBROUTINE ADD_NODE

  SUBROUTINE RD_INT_NODE(THIS, VAL)
    IMPLICIT NONE
    CLASS(LIST)       :: THIS
    INTEGER           :: VAL
    CLASS(*), POINTER :: BHO

    SELECT TYPE (BHO => THIS%NODE%KEY)
    TYPE IS (INTEGER)
       VAL = BHO
    END SELECT

    RETURN
  END SUBROUTINE RD_INT_NODE

  SUBROUTINE RD_DBL_NODE(THIS, VAL)
    IMPLICIT NONE
    CLASS(LIST)       :: THIS
    REAL(8)           :: VAL
    CLASS(*), POINTER :: BHO

    SELECT TYPE (BHO => THIS%NODE%KEY)
    TYPE IS (REAL(8))
       VAL = BHO
    END SELECT

    RETURN
  END SUBROUTINE RD_DBL_NODE

  SUBROUTINE RD_CMPLX_NODE(THIS, VAL)
    IMPLICIT NONE
    CLASS(LIST)       :: THIS
    COMPLEX(8)        :: VAL
    CLASS(*), POINTER :: BHO

    SELECT TYPE (BHO => THIS%NODE%key)
    TYPE IS (COMPLEX(8))
       VAL = BHO
    END SELECT

    RETURN
  END SUBROUTINE RD_CMPLX_NODE


END MODULE LISTS_CLASS

MODULE META_CLASS
  USE LISTS_CLASS
  IMPLICIT NONE

  TYPE MEMORY_POT
     TYPE(LIST) :: GAUSS
     REAL(8)    :: WEIGHT = 1.0D0
     REAL(8)    :: SIGMA  = 0.2D0
  END TYPE MEMORY_POT

  TYPE(MEMORY_POT) :: META1, META2
  LOGICAL          :: DO_META = .TRUE.
  REAL(8)          :: DIST1, DIST2, HEIGHT1, HEIGHT2
  REAL(8)          :: F1(3), F2(3), F3(3)

END MODULE META_CLASS

SUBROUTINE GET_LSMF_SNAP
  USE FIT_LAMMPS_CLASS
  USE LAPACK_DIAG_SIMM
  USE LAPACK_INVERSE
  USE COMMON_VAR
  USE LAMMPS
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_PTR, C_INT
  IMPLICIT NONE
  INTEGER                :: AAA, BBB, L, I, K, N, V, S, T, LL, TT
  INTEGER                :: SS, VV, INF, LWORK, DIMA, DIMB
  REAL(8), ALLOCATABLE   :: X(:), B(:), A(:,:), WORK(:), AA(:,:)
  REAL(8), ALLOCATABLE   :: W(:), TVAL(:,:), YY(:,:), Y(:,:), A2(:,:)
  REAL(8), ALLOCATABLE   :: ACM(:,:), BCM(:), BB(:)
  REAL(8)                :: DUMP, AVE

  REAL (C_DOUBLE), POINTER   :: ENER(:,:) => NULL()
  REAL (C_DOUBLE), POINTER   :: F(:,:) => NULL()
  REAL (C_DOUBLE), POINTER   :: KIND_NAT(:) => NULL()
  REAL (C_DOUBLE), POINTER   :: EREF => NULL()
  REAL (C_DOUBLE), POINTER   :: FX_REF(:) => NULL()
  REAL (C_DOUBLE), POINTER   :: FY_REF(:) => NULL()
  REAL (C_DOUBLE), POINTER   :: FZ_REF(:) => NULL()
  REAL (C_DOUBLE), POINTER   :: ID_DBL(:)=> NULL()
  INTEGER, ALLOCATABLE       :: MAP(:), ID(:)

  INTEGER                         :: NKINDS, BI_ORDER, NPAR, NPAR2FIT, TOT_FRAMES
  INTEGER                         :: NATS2FIT, QUADFLAG
  REAL(8)                         :: GEN_CUTOFF
  CHARACTER(LEN=10000)            :: SNAP_STRING, SNAP_STRING2
  REAL(8), ALLOCATABLE            :: CUTOFF(:), RADII(:), SIGMA(:)
  INTEGER, ALLOCATABLE            :: KIND_COUNT(:), TYPE2FIT(:)
  CHARACTER(LEN=2), ALLOCATABLE   :: LABEL(:)

  TOT_FRAMES = SYS%TOT_FRAMES
  NPAR2FIT   = SYS%NPAR2FIT

  IF(PRINT_BI) OPEN(121, FILE='Bi_compoenents.dat')

  OPEN(16, FILE = SYS%INP_FIT)
  READ(16, *) GEN_CUTOFF, BI_ORDER, NPAR, QUADFLAG
  READ(16, *) NKINDS
  ALLOCATE(LABEL(NKINDS))
  ALLOCATE(TYPE2FIT(NKINDS))
  ALLOCATE(CUTOFF(NKINDS))
  ALLOCATE(SIGMA(NKINDS))
  ALLOCATE(RADII(NKINDS))
  ALLOCATE(KIND_COUNT(NKINDS))
  KIND_COUNT = 0
  DO I = 1, NKINDS
     READ(16,*) LABEL(I), TYPE2FIT(I), RADII(I), CUTOFF(I), SIGMA(I)
     WRITE(*,*) LABEL(I), TYPE2FIT(I), RADII(I), CUTOFF(I), SIGMA(I)
  END DO
  CLOSE(16)
  WRITE(SNAP_STRING, *) GEN_CUTOFF, '1.0000', BI_ORDER, (RADII(I), I = 1, NKINDS),&
       & (CUTOFF(I), I = 1, NKINDS), 'quadraticflag ', QUADFLAG
  WRITE(SNAP_STRING2, *) (TYPE2FIT(I), I = 1, NKINDS)
  WRITE(*,*) TRIM(SNAP_STRING)
  WRITE(*,*) 'Fitting types: ',TRIM(SNAP_STRING2)

  IF ( .NOT. SKIP_FIT ) THEN
     OPEN(222, FILE = 'snapcoeff')
     WRITE(222, *) NKINDS, NPAR
     DO I = 1, NKINDS
        WRITE(222,*) LABEL(I), RADII(I), CUTOFF(I)
        DO N = 1, NPAR
           WRITE(222, *) 0.0D0
        END DO
     END DO
     FLUSH(222)
     CLOSE(222)
  END IF

  ALLOCATE(LMP(SYS%NDATA))

  PRINT*,'INITIALIZING LAMMPS POINTERS: TOTAL = ', SYS%NDATA
  DO I = 1, SYS%NDATA

     CALL LAMMPS_OPEN_NO_MPI ('lmp -log none -screen none', LMP(I))
     CALL LAMMPS_FILE (LMP(I), SYS%INP)
     CALL LAMMPS_COMMAND (LMP(I), 'read_data ' // tRIM(SYS%DATA(I)%INP_DATA))
     CALL LAMMPS_COMMAND (LMP(I), 'group fitsnap type ' // TRIM(SNAP_STRING2))
     CALL LAMMPS_FILE (LMP(I), SYS%INP_FIX)
     CALL LAMMPS_COMMAND (LMP(I), &
          'compute sna_e all sna/atom ' // TRIM(SNAP_STRING) //&
          ' diagonal 3 rmin0 0 switchflag 1')
     CALL LAMMPS_COMMAND (LMP(I), 'compute type all property/atom type')
     CALL LAMMPS_COMMAND (LMP(I), 'compute id all property/atom id')
     CALL LAMMPS_COMMAND (LMP(I), 'compute pe_ener all pe')

     IF (FIT_FORCES) THEN
        CALL LAMMPS_COMMAND (LMP(I), &
             'compute sna_f all snad/atom ' // TRIM(SNAP_STRING)//&
             ' diagonal 3 rmin0 0 switchflag 1')
        CALL LAMMPS_COMMAND (LMP(I), 'compute f_x all property/atom fx')
        CALL LAMMPS_COMMAND (LMP(I), 'compute f_y all property/atom fy')
        CALL LAMMPS_COMMAND (LMP(I), 'compute f_z all property/atom fz')
     END IF

  END DO


  !! DO KERNEL !!

  IF (REFINE) THEN

     CALL KERNEL%DEALLOC()

     KERNEL%NKINDS = NKINDS
     ALLOCATE(KERNEL%K(NKINDS))
     DO I = 1, NKINDS
        KERNEL%K(I)%SIGMA = SIGMA(I)
        KERNEL%K(I)%NENVS = 0
     END DO

     DO I = 1, SYS%NDATA
        CALL LAMMPS_COMMAND (LMP(I), 'run 0')
        CALL LAMMPS_EXTRACT_COMPUTE (KIND_NAT, LMP(I), 'type', 1, 1)
        CALL LAMMPS_EXTRACT_COMPUTE (ENER, LMP(I), 'sna_e', 1, 2)
        KIND_COUNT = 0
        DO T = 1, SYS%DATA(I)%NATS
           KIND_COUNT(NINT(KIND_NAT(T))) = KIND_COUNT(NINT(KIND_NAT(T))) + 1
        END DO
        DO T = 1, NKINDS
           KERNEL%K(T)%NENVS = KERNEL%K(T)%NENVS + KIND_COUNT(T) * SYS%DATA(I)%FRAMES
        END DO
     END DO
     DO T = 1, NKINDS
        ALLOCATE(KERNEL%K(T)%B(KERNEL%K(T)%NENVS, SIZE(ENER,1)))
     END DO

     DO I = 1, NKINDS
        KERNEL%K(I)%NENVS = 0
     END DO

     DO I = 1, SYS%NDATA
        DO L = 1, SYS%DATA(I)%FRAMES

           CALL LAMMPS_SCATTER_ATOMS (LMP(I), 'x', SYS%DATA(I)%X(L, 1:3*SYS%DATA(I)%NATS))
           CALL LAMMPS_COMMAND (LMP(I), 'run 0')
           CALL LAMMPS_EXTRACT_COMPUTE (KIND_NAT, LMP(I), 'type', 1, 1)
           CALL LAMMPS_EXTRACT_COMPUTE (ENER, LMP(I), 'sna_e', 1, 2)
           CALL LAMMPS_EXTRACT_COMPUTE (EREF, LMP(I), 'pe_ener', 0, 0)
           IF (.NOT. ALLOCATED(ID))  ALLOCATE(ID(SYS%DATA(I)%NATS))
           IF (.NOT. ALLOCATED(MAP)) ALLOCATE(MAP(SYS%DATA(I)%NATS))
           CALL LAMMPS_EXTRACT_COMPUTE (ID_DBL, LMP(I), 'id', 1, 1)

           ID = INT(ID_DBL)
           ID_DBL => NULL()

           DO K = 1, SYS%DATA(I)%NATS
              MAP(ID(K)) = K
           END DO

           DO T = 1, SYS%DATA(I)%NATS
              V = KERNEL%K(NINT(KIND_NAT(MAP(T))))%NENVS + 1
              DO K = 1, SIZE(ENER, 1)
                 KERNEL%K(NINT(KIND_NAT(MAP(T))))%B(V, K) = ENER(K, MAP(T))
              END DO
              KERNEL%K(NINT(KIND_NAT(MAP(T))))%NENVS = KERNEL%K(NINT(KIND_NAT(MAP(T))))%NENVS + 1
           END DO

           ID_DBL   => NULL()
           ENER     => NULL()
           EREF     => NULL()
           KIND_NAT => NULL()

        END DO

        DEALLOCATE(ID)
        DEALLOCATE(MAP)

     END DO

  END IF ! REFINE

  !! DO FITTING !!

  IF (.NOT. SKIP_FIT) THEN

     V  = 1
     VV = 1

     DO I = 1, SYS%NDATA

        DO L = 1, SYS%DATA(I)%FRAMES

           CALL LAMMPS_SCATTER_ATOMS (LMP(I), 'x', SYS%DATA(I)%X(L, 1:3*SYS%DATA(I)%NATS))
           CALL LAMMPS_COMMAND (LMP(I), 'run 0')
           CALL LAMMPS_EXTRACT_COMPUTE (KIND_NAT, LMP(I), 'type', 1, 1)
           CALL LAMMPS_EXTRACT_COMPUTE (ENER, LMP(I), 'sna_e', 1, 2)
           CALL LAMMPS_EXTRACT_COMPUTE (EREF, LMP(I), 'pe_ener', 0, 0)
           IF (.NOT. ALLOCATED(ID) ) ALLOCATE(ID(SYS%DATA(I)%NATS))
           IF (.NOT. ALLOCATED(MAP)) ALLOCATE(MAP(SYS%DATA(I)%NATS))
           CALL LAMMPS_EXTRACT_COMPUTE (ID_DBL, LMP(I), 'id', 1, 1)

           ID = INT(ID_DBL)
           ID_DBL => NULL()

           DO K = 1, SYS%DATA(I)%NATS
              MAP(ID(K)) = K
           END DO

           IF (L .EQ. 1) THEN

              NPAR = SIZE(ENER, 1) + 1
              KIND_COUNT = 0

              DO T=1,SYS%DATA(I)%NATS
                 KIND_COUNT(NINT(KIND_NAT(T))) = KIND_COUNT(NINT(KIND_NAT(T))) + 1
              ENDDO
              
              DIMA = NPAR * NKINDS
              DIMB = NPAR2FIT

              IF (CS) DIMB = DIMB + (NPAR-1) * NKINDS
              IF (E0CS .EQ. 1) DIMB = DIMB + NKINDS - 1
              IF (E0CS .EQ. 2) DIMB = DIMB + NKINDS

              IF (.NOT. ALLOCATED(B)) THEN
                 ALLOCATE(B(DIMB))
                 B = 0.0D0
              END IF
              IF (.NOT. ALLOCATED(A)) THEN
                 ALLOCATE(A(DIMB, DIMA))
                 A = 0.0D0
              END IF

              !            if(.not.cs)then
              !             if(.not.allocated(B)) allocate(B(npar2fit+nkinds-1))
              !             if(.not.allocated(B)) allocate(B(npar2fit+nkinds))
              !             if(.not.allocated(B)) allocate(B(npar2fit))
              !             if(.not.allocated(A)) then
              !              allocate(A(npar2fit+nkinds-1,npar*nkinds))
              !              allocate(A(npar2fit+nkinds,npar*nkinds))
              !              allocate(A(npar2fit,npar*nkinds))
              !              A=0.0d0
              !             endif
              !            else    ! compressive sensing
              !             if(.not.allocated(B)) allocate(B(npar2fit+nkinds-1+(npar-1)*nkinds))
              !             if(.not.allocated(B)) allocate(B(npar2fit+nkinds+(npar-1)*nkinds))
              !             if(.not.allocated(B)) allocate(B(npar2fit+(npar-1)*nkinds))
              !             if(.not.allocated(A)) then
              !              allocate(A(npar2fit+(npar-1)*nkinds+nkinds-1,npar*nkinds))
              !              allocate(A(npar2fit+(npar-1)*nkinds+nkinds,npar*nkinds))
              !              allocate(A(npar2fit+(npar-1)*nkinds,npar*nkinds))
              !              A=0.0d0
              !             endif
              !            endif

              IF (PCA) THEN
                 OPEN(1313, FILE = 'PCA.dat')
                 !             if(.not.allocated(A2)) allocate(A2(tot_frames,(npar-1)*nkinds))
                 !             if(.not.allocated(BB)) allocate(BB(tot_frames))
                 !             if(.not.allocated(W)) allocate(W((npar-1)*nkinds))
                 !             if(.not.allocated(Y)) allocate(Y(tot_frames,(npar-1)*nkinds))
                 !             if(.not.allocated(YY)) allocate(YY((npar-1)*nkinds,(npar-1)*nkinds))
                 !             if(.not.allocated(Tval)) allocate(Tval(tot_frames,(npar-1)*nkinds))
                 IF (.NOT. ALLOCATED(A2)  ) ALLOCATE(A2(TOT_FRAMES, NPAR-1))
                 IF (.NOT. ALLOCATED(BB)  ) ALLOCATE(BB(TOT_FRAMES))
                 IF (.NOT. ALLOCATED(W)   ) ALLOCATE(W(NPAR-1))
                 IF (.NOT. ALLOCATED(Y)   ) ALLOCATE(Y(TOT_FRAMES, NPAR-1))
                 IF (.NOT. ALLOCATED(YY)  ) ALLOCATE(YY(NPAR-1, NPAR-1))
                 IF (.NOT. ALLOCATED(TVAL)) ALLOCATE(TVAL(TOT_FRAMES, NPAR-1))
              END IF

           END IF
           
           IF (FIT_ENER) THEN
              
              S = 1
              DO K = 1, NKINDS
                 A(V, S) = KIND_COUNT(K)
                 S = S + NPAR
              END DO
              
              IF(PRINT_BI) WRITE(121, *) ENER(:, 1), ENER(:, 2), ENER(:, 3)

              DO K = 2, NPAR
                 DO T = 1, SYS%DATA(I)%NATS
                    
                    S = K + (NINT(KIND_NAT(T)) - 1) * NPAR
                    A(V, S) = A(V, S) + ENER(K-1, MAP(T)) * SYS%DATA(I)%WEIGHT
                    IF (PCA) THEN
                       IF (KIND_NAT(T) .EQ. 1) THEN
                          SS = K - 1 + (NINT(KIND_NAT(T)) - 1) * (NPAR - 1)
                          A2(VV, SS) = A(VV, SS) + ENER(K-1, MAP(T))
                       END IF
                    END IF

                 END DO
              END DO
              
              B(V) = (SYS%DATA(I)%ENER(L) - EREF) * SYS%DATA(I)%WEIGHT
              V  = V  + 1
              VV = VV + 1

              ENER => NULL()
              EREF => NULL()

           END IF ! FITENER

           IF (FIT_FORCES) THEN

              CALL LAMMPS_EXTRACT_COMPUTE (F, LMP(I), 'sna_f', 1, 2)
              CALL LAMMPS_EXTRACT_COMPUTE (FX_REF, LMP(I), 'f_x', 1, 1)
              CALL LAMMPS_EXTRACT_COMPUTE (FY_REF, LMP(I), 'f_y', 1, 1)
              CALL LAMMPS_EXTRACT_COMPUTE (FZ_REF, LMP(I), 'f_z', 1, 1)

              DO T = 1, SYS%DATA(I)%NATS

                 S = 1
                 DO K = 1, NKINDS
                    A(V, S) = 0.0D0
                    S = S + NPAR
                 END DO

                 DO N = 1, NKINDS
                    S = (N-1) * (3 * (NPAR-1)) + 1  ! NPAR-1?
                    DO K = 2, NPAR
                       A(V, (N-1) * NPAR + K) = F(S, MAP(T)) * SYS%DATA(I)%WEIGHT
                       S = S + 1
                    END DO
                 END DO
                 B(V) = (SYS%DATA(I)%FX(L, T) - FX_REF(MAP(T))) * SYS%DATA(I)%WEIGHT
                 V = V + 1

                 S = 1
                 DO K = 1, NKINDS
                    A(V, S) = 0.0D0
                    S = S + NPAR
                 END DO

                 DO N = 1, NKINDS
                    S = (N-1) * (3 * (NPAR-1)) + 1 + NPAR - 1
                    DO K = 2, NPAR
                       A(V, (N-1) * NPAR + K) = F(S, MAP(T)) * SYS%DATA(I)%WEIGHT
                       S = S + 1
                    END DO
                 END DO
                 B(V) = (SYS%DATA(I)%FY(L,T)-FY_REF(MAP(T))) * SYS%DATA(I)%WEIGHT
                 V = V + 1

                 S = 1
                 DO K = 1, NKINDS
                    A(V, S) = 0.0D0
                    S = S + NPAR
                 END DO

                 DO N = 1, NKINDS
                    S = (N-1) * (3 * (NPAR-1)) + 1 + 2 * (NPAR-1)
                    DO K = 2, NPAR
                       A(V, (N-1) * NPAR + K) = F(S,MAP(T)) * SYS%DATA(I)%WEIGHT
                       S = S + 1
                    END DO
                 END DO
                 B(V) = (SYS%DATA(I)%FZ(L, T) - FZ_REF(MAP(T))) * SYS%DATA(I)%WEIGHT
                 V = V + 1

              END DO

              F      => NULL()
              FX_REF => NULL()
              FY_REF => NULL()
              FZ_REF => NULL()

           END IF  ! END IF ON FORSE
        END DO  ! CICLO ON FRAMES

        IF (ALLOCATED(MAP)) DEALLOCATE(MAP)
        IF (ALLOCATED(ID) ) DEALLOCATE(ID)
        PRINT*,'COMPUTING BISPECTRUM: LINE => ', I, "FRAMES IN LINE: ", SYS%DATA(I)%FRAMES
     END DO   ! CICLO SU DATA

     IF (E0CS .EQ. 1) THEN
        S = 1
        DO K = 1, NKINDS - 1
           B(V)    = 0.0D0
           A(V, S) = 1.0D0
           S = S + NPAR
           V = V + 1
        END DO
     END IF

     IF (E0CS .EQ. 2) THEN
        S = 1
        DO K = 1, NKINDS
           B(V)    = 0.0D0
           A(V, S) = 1.0D2 !. WHY?
           S = S + NPAR
           V = V + 1
        END DO
     END IF

     ! COMPRESSIVE SENSING

     IF (CS) THEN
        DO K = 1, NKINDS
           DO L = 2, NPAR
              B(V) = 0.0D0
              S = (K-1) * NPAR + L
              A(V, S) = CM_VAL
              V = V + 1
           END DO
        END DO
     END IF

     ! PRINCIPAL COMPONENT ANALYSIS

     IF (PCA) THEN

        Y = A2
        DO L = 1, SIZE(Y, 2)
           AVE = 0.0D0
           DO S = 1, SIZE(Y, 1)
              AVE = AVE + Y(S, L)
           END DO
           Y(1:SIZE(Y, 1), L) = Y(1:SIZE(Y, 1), L) - AVE / SIZE(Y, 1)
        END DO
        YY = MATMUL(TRANSPOSE(Y), Y)
        CALL NEW_DIAG(SIZE(YY, 1), YY, W)

        IF (.NOT.ALLOCATED(AA)) ALLOCATE(AA(TOT_FRAMES,4))

        AA = MATMUL(A2, YY(:, (SIZE(W)-4):SIZE(W)))

        W = SQRT(ABS(W))
        W = W / SUM(W)
        WRITE(1313, *) 'Principal Components Values'
        DO K = SIZE(W), 1, -1
           WRITE(1313, *) W(K)
        END DO
        DO L = 1, SIZE(AA, 1)
           WRITE(1313, *) (AA(L, K), K = 1, 4)
        END DO

        !!
        !         aaa=70
        !         bbb=tot_frames
        !         lwork=bbb+64*bbb+1000
        !         allocate(work(lwork))
        !         call dgels('N',bbb,aaa,1,AA,bbb,BB,bbb,WORK,LWORK,inf)
        !         deallocate(work)
        !         if(inf.ne.0)then
        !          write(*,*) 'zgels failed',inf
        !           stop
        !         else

        !         YY=transpose(YY)
        !         BB=matmul(YY(:,1:70),BB(1:70))

        !         open(222,file='snapcoeff')

        !         l=1
        !         write(222,*) nkinds,npar
        !         do i =1,nkinds
        !          write(222,*) label(i),radii(i),cutoff(i)
        !          write(222,*) 0.000000000
        !          do n=1,npar-1
        !           write(222,*) BB(l)
        !           l=l+1
        !          enddo
        !         enddo
        !         close(222)
        
        !        endif
        !!

     END IF

     !        aaa=npar*nkinds
     !        if(cs)then
     !         bbb=npar2fit+(nkinds*(npar-1))!-1
     !         bbb=npar2fit+(nkinds*npar)-1
     !        else
     !         bbb=npar2fit!+nkinds!-1
     !         bbb=npar2fit
     !        endif

     LWORK = DIMB + 64 * DIMB + 1000 ! VERY DISTURBING !
     ALLOCATE(WORK(LWORK))
     CALL DGELS('N', DIMB, DIMA, 1, A, DIMB, B, DIMB, WORK, LWORK, INF)
     DEALLOCATE(WORK)
     IF (INF .NE. 0) THEN
        WRITE(*,*) 'zgels failed', INF
        STOP
     ELSE
        OPEN(222, FILE = 'snapcoeff')
        L = 1
        WRITE(222, *) NKINDS, NPAR
        DO I = 1, NKINDS
           WRITE(222, *) LABEL(I), RADII(I), CUTOFF(I)
           DO N = 1, NPAR
              WRITE(222, *) B(L)
              L = L + 1
           END DO
        END DO

        !        dump=0.0d0
        !        do i=npar*nkinds+1,size(B)
        !         dump=dump+B(i)**2
        !        enddo

        !        write(222,*) 'dgels residulal: ',dump
        CLOSE(222)

        OPEN(333, FILE = 'snapparam')
        WRITE(333, *) 'rcutfac ', GEN_CUTOFF
        WRITE(333, *) 'twojmax ', BI_ORDER
        WRITE(333, *) 'quadraticflag ', QUADFLAG
        WRITE(333, *) 'rfac0 1.00000'
        WRITE(333, *) 'rmin0 0'
        WRITE(333, *) 'diagonalstyle 3'
        WRITE(333, *) 'switchflag 1'
        CLOSE(333)

        DO I=1,SYS%NDATA
           CALL LAMMPS_FILE (LMP(I),SYS%INP_FIX)
        END DO

     END IF

  END IF ! SKIP_FIT

  RETURN
END SUBROUTINE GET_LSMF_SNAP

SUBROUTINE REFINE_SNAP
  USE FIT_LAMMPS_CLASS
  USE LAPACK_DIAG_SIMM
  USE LAPACK_INVERSE
  USE RANDOM_NUMBERS_CLASS
  USE COMMON_VAR
  USE LAMMPS
  USE LISTS_CLASS
  USE META_CLASS
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_PTR, C_INT
  IMPLICIT NONE

  INTEGER                         :: L, I, J, V, K
  REAL(8), ALLOCATABLE            :: R(:)
  !
  CHARACTER(LEN=10000)            :: SNAP_STRING
  REAL(8), ALLOCATABLE            :: CUTOFF(:), RADII(:)
  INTEGER, ALLOCATABLE            :: KIND_COUNT(:), TYPE2FIT(:)
  CHARACTER (LEN=2),ALLOCATABLE   :: LABEL(:)
  INTEGER                         :: NKINDS,BI_ORDER, NPAR, NPAR2FIT, TOT_FRAMES
  INTEGER                         :: NATS2FIT, QUADFLAG
  REAL(8)                         :: GEN_CUTOFF, PI

  REAL (C_DOUBLE), POINTER        :: ID_DBL(:)=> NULL()
  REAL (C_DOUBLE), POINTER        :: KIND_NAT(:) => NULL()
  INTEGER, ALLOCATABLE            :: MAP(:), ID(:)

  REAL(8), ALLOCATABLE   :: MAX_KERNEL(:)
  REAL(8)                :: DFT_ENER,VAL,RAND_SEED

  REAL (C_DOUBLE), POINTER   :: B(:,:) => NULL()
  REAL (C_DOUBLE), POINTER   :: FF_ENER => NULL()

  ALLOCATE(LMP(1))

  PI=ACOS(-1.0D0)

  CALL LAMMPS_OPEN_NO_MPI ('lmp -log none', LMP(1))
  CALL LAMMPS_FILE (LMP(1), SYS%INP)
  CALL LAMMPS_COMMAND (LMP(1), 'read_data ' // TRIM(SYS%DATA(1)%INP_DATA))
  CALL LAMMPS_FILE (LMP(1), SYS%INP_FIX)
  CALL LAMMPS_COMMAND (LMP(1), 'compute pe_ener all pe')

  OPEN(16, FILE = SYS%INP_FIT)
  READ(16, *) GEN_CUTOFF, BI_ORDER, NPAR, QUADFLAG
  READ(16, *) NKINDS
  ALLOCATE(LABEL(NKINDS))
  ALLOCATE(TYPE2FIT(NKINDS))
  ALLOCATE(CUTOFF(NKINDS))
  ALLOCATE(RADII(NKINDS))
  ALLOCATE(KIND_COUNT(NKINDS))
  KIND_COUNT = 0
  DO I = 1, NKINDS
     READ(16, *) LABEL(I), TYPE2FIT(I), RADII(I), CUTOFF(I)
     WRITE(*, *) LABEL(I), TYPE2FIT(I), RADII(I), CUTOFF(I)
  END DO
  CLOSE(16)
  WRITE(SNAP_STRING, *) GEN_CUTOFF, '1.0000', BI_ORDER, (RADII(I), I = 1, NKINDS),&
       & (CUTOFF(I), I = 1, NKINDS), 'quadraticflag ',QUADFLAG

  CALL LAMMPS_COMMAND (LMP(1), &
       'compute sna_e all sna/atom ' // TRIM(SNAP_STRING)//&
       ' diagonal 3 rmin0 0 switchflag 1')
  CALL LAMMPS_COMMAND (LMP(1), 'compute type all property/atom type')
  CALL LAMMPS_COMMAND (LMP(1), 'compute id all property/atom id')

  ! MINIMIZE ENERGY

  CALL LAMMPS_COMMAND (LMP(1), 'thermo 1')
  CALL LAMMPS_COMMAND (LMP(1), 'thermo_style custom step time temp pe etotal ')
  !         CALL LAMMPS_COMMAND (LMP(1), 'minimize 1.0e-8 1.0e-8 1000 100000')

  WRITE(SNAP_STRING, *) (LABEL(I) // ' ', I = 1, SIZE(LABEL))

  CALL LAMMPS_COMMAND (LMP(1), 'dump xyz_dump all xyz 5 geo_opt.xyz')

  CALL LAMMPS_COMMAND (LMP(1), 'dump_modify            xyz_dump element '//trim(snap_string))

  ! SET VELOCITIES

  WRITE(SNAP_STRING, *) REFINE_TEMP

  CALL LAMMPS_COMMAND (LMP(1), 'timestep 0.25')
  CALL LAMMPS_COMMAND (LMP(1), 'variable t equal' // TRIM(SNAP_STRING))

  CALL RANDOM_NUMBER(RAND_SEED)
  WRITE(SNAP_STRING, *) NINT(RAND_SEED * 10000.0D0)

  CALL LAMMPS_COMMAND (LMP(1), 'velocity all create $t ' // TRIM(SNAP_STRING) // ' dist gaussian')
  CALL LAMMPS_COMMAND (LMP(1), 'velocity all zero linear')
  CALL LAMMPS_COMMAND (LMP(1), 'velocity all zero angular')

  CALL LAMMPS_COMMAND (LMP(1), 'group atom1 id 1')
  CALL LAMMPS_COMMAND (LMP(1), 'group atom2 id 4')
  !         CALL LAMMPS_COMMAND (LMP(1), 'group atom3 id 7')
  CALL META1%GAUSS%INIT()
  !         CALL META2%GAUSS%INIT()

  CALL LAMMPS_COMMAND (LMP(1), 'fix 1 all nvt temp $t $t 100.0 tchain 3')

  OPEN(13, FILE='new_geo.xyz',    ACCESS='APPEND')
  OPEN(14, FILE='kernel_max.dat', ACCESS='APPEND')
  OPEN(15, FILE='new_geo.ener',   ACCESS='APPEND')
  OPEN(16, FILE='meta.dat',       ACCESS='APPEND')

  IF (.NOT. ALLOCATED(MAX_KERNEL)) ALLOCATE(MAX_KERNEL(SYS%DATA(1)%NATS))
  ! IF (.NOT. ALLOCATED(MAX_KERNEL)) ALLOCATE(MAX_KERNEL(KERNEL%NKINDS))

  WRITE(14, *) '## nkinds: ', KERNEL%NKINDS, 'nenvs: ', (KERNEL%K(I)%NENVS, I = 1, KERNEL%NKINDS)
  FLUSH(14)

  DO I = 1, 400000              ! RAJARSHI: WHY THIS PARTICULAR NUMBER?

     CALL LAMMPS_COMMAND (LMP(1), 'run 4')
     CALL LAMMPS_EXTRACT_COMPUTE (B, LMP(1), 'sna_e', 1, 2)
     CALL LAMMPS_EXTRACT_COMPUTE (FF_ENER, LMP(1), 'pe_ener', 0, 0)
     CALL LAMMPS_EXTRACT_COMPUTE (KIND_NAT, LMP(1), 'type', 1, 1)
     IF (.NOT. ALLOCATED(ID)) ALLOCATE(ID(SYS%DATA(1)%NATS))
     IF (.NOT. ALLOCATED(MAP)) ALLOCATE(MAP(SYS%DATA(1)%NATS))
     CALL LAMMPS_EXTRACT_COMPUTE (ID_DBL, LMP(1), 'id', 1, 1)

     IF (DO_META) THEN

        IF (ALLOCATED(R)) DEALLOCATE(R)
        CALL LAMMPS_GATHER_ATOMS (LMP(1), 'x', 3, R)

        DIST1 = (R(1) - R(10))**2
        DIST1 = DIST1 + (R(2) - R(11))**2
        DIST1 = DIST1 + (R(3) - R(12))**2
        DIST1 = SQRT(DIST1)

        !           DIST2=(R(4)-R(7))**2
        !           DIST2=DIST2+(R(5)-R(8))**2
        !           DIST2=DIST2+(R(6)-R(9))**2
        !           DIST2=SQRT(DIST2)

        F1 = 0.0D0
        F2 = 0.0D0
        !           F3=0.0D0

        IF (META1%GAUSS%NELEM .GE. 1) THEN

           CALL META1%GAUSS%REBOOT()
           !            CALL META2%GAUSS%REBOOT()

           DO J = 1, META1%GAUSS%NELEM

              CALL META1%GAUSS%RD_VAL(HEIGHT1)
              !   call meta2%gauss%rd_val(height2)
              !
              !f1(1)=f1(1)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !f1(2)=f1(2)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !f1(3)=f1(3)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1
              !
              !f2(1)=f2(1)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !f2(2)=f2(2)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !f2(3)=f2(3)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2)*(dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1
              !
              !f1(1)=f1(1)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !f1(2)=f1(2)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !f1(3)=f1(3)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1
              !
              !f2(1)=f2(1)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !f2(2)=f2(2)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !f2(3)=f2(3)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1
              !
              !f1(1)=f1(1)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(4)-r(19))/meta1%sigma**2/dist2
              !f1(2)=f1(2)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(5)-r(20))/meta1%sigma**2/dist2
              !f1(3)=f1(3)+2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(6)-r(21))/meta1%sigma**2/dist2
              !
              !f3(1)=f3(1)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(4)-r(19))/meta1%sigma**2/dist2
              !f3(2)=f3(2)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(5)-r(20))/meta1%sigma**2/dist2
              !f3(3)=f3(3)-2*meta1%weight*exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*(dist2-height2)*(r(6)-r(21))/meta1%sigma**2/dist2
              !
              CALL META1%GAUSS%SKIP()
              !   CALL META2%GAUSS%SKIP()

           END DO

        END IF

        !! CONTRAINT DISTANCES
        IF (DIST1 .GT. 8.0D0 .AND. DIST1 .LT. 16.0D0) THEN
           ! WHY MULTIPLY BY 2000 ? !
           F1(1) = F1(1) + 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(1) - R(10)) / DIST1
           F1(2) = F1(2) + 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(2) - R(11)) / DIST1
           F1(3) = F1(3) + 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(3) - R(12)) / DIST1
           F2(1) = F2(1) - 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(1) - R(10)) / DIST1
           F2(2) = F2(2) - 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(2) - R(11)) / DIST1
           F2(3) = F2(3) - 2000.D0 * SIN(PI * DIST1 / 8.0D0) * PI/8.0D0 * (R(3) - R(12)) / DIST1
        END IF

        IF (I .GE. 2) THEN
           CALL LAMMPS_COMMAND (LMP(1), 'unfix add1')
           CALL LAMMPS_COMMAND (LMP(1), 'unfix add2')
           !   CALL LAMMPS_COMMAND (LMP(1), 'unfix add3')
        END IF

        WRITE(SNAP_STRING, "(3F12.6)") F1(1), F1(2), F1(3)
        CALL LAMMPS_COMMAND (LMP(1), 'fix add1 atom1 addforce ' // TRIM(SNAP_STRING))
        WRITE(SNAP_STRING, "(3F12.6)") F2(1), F2(2), F2(3)
        CALL LAMMPS_COMMAND (LMP(1), 'fix add2 atom2 addforce ' // TRIM(SNAP_STRING))
        !   WRITE(SNAP_STRING,"(3F12.6)") F3(1), F3(2), F3(3)
        !   CALL LAMMPS_COMMAND (LMP(1), 'fix add3 atom3 addforce ' // TRIM(SNAP_STRING))

        IF (MOD(DBLE(I), 10.0D0) .LT. 1.0E-6) THEN
           CALL META1%GAUSS%ADD_NODE(DIST1)
           !            CALL META2%GAUSS%ADD_NODE(DIST2)
           WRITE(16, *) I, DIST1!, DIST2
        END IF

     END IF

     ID = INT(ID_DBL)
     ID_DBL => NULL()

     DO K = 1, SYS%DATA(1)%NATS
        MAP(ID(K)) = K
     END DO

     MAX_KERNEL = 0.0D0

     DO K = 1, SIZE(B, 2)
        L = NINT(KIND_NAT(MAP(K)))
        DO V = 1, KERNEL%K(L)%NENVS
           VAL = 0.0D0
           DO J = 1, SIZE(B, 1)
              VAL = VAL - (B(J, MAP(K)) - KERNEL%K(L)%B(V, J))**2
           END DO
           VAL = EXP(VAL / 2 * KERNEL%K(L)%SIGMA**2)
           IF (VAL .GT. MAX_KERNEL(MAP(K))) MAX_KERNEL(MAP(K)) = VAL
        END DO
     END DO

     WRITE(14, *) I, MAX_KERNEL
     FLUSH(14)

     IF (ANY(MAX_KERNEL .LT. THR_KERNEL)) THEN
        IF (ALLOCATED(R)) DEALLOCATE(R)
        CALL LAMMPS_GATHER_ATOMS (LMP(1),'x', 3, R)
        WRITE(13, *) SYS%DATA(1)%NATS
        WRITE(13, *)
        L = 1
        DO J = 1, SYS%DATA(1)%NATS
           WRITE(13, *) SYS%DATA(1)%LABEL(J), R(L), R(L+1), R(L+2)
           L = L + 3
        END DO
        FLUSH(13)
        CALL EXECUTE_COMMAND_LINE('./run_DFT_scf.x')
        OPEN(16, FILE = 'last_dft_ener.dat'); READ(16, *) DFT_ENER; CLOSE(16)
        WRITE(15, *) DFT_ENER, FF_ENER
        FLUSH(15)
        CLOSE(13); CLOSE(14); CLOSE(15)
        CALL LAMMPS_CLOSE (LMP(1))
        DEALLOCATE(LMP)
        RETURN                  ! RAJARSHI: NEXT TWO STATEMENTS WON'T RUN!
        WRITE(*, *) 'PUPPA'
        FLUSH(6)
     END IF

  END DO

  CALL LAMMPS_CLOSE (LMP(1))
  DEALLOCATE(LMP)
  REFINE = .FALSE.
  CLOSE(13)
  CLOSE(14)
  CLOSE(15)

  RETURN
END SUBROUTINE REFINE_SNAP

PROGRAM FITSNAP
  USE FIT_LAMMPS_CLASS
  USE COMMON_VAR
  IMPLICIT NONE
  INTEGER              :: L
  CHARACTER(LEN=100)   :: COMMAND, INPUT, DATAS, OUTPUT, NEW_DATAS

  IF (IARGC() .EQ. 0) THEN
     WRITE(*, *) 'FitSnap Usage:'
     WRITE(*, *) '-datas   : list of files to be fitted'
     WRITE(*, *) '-inp     : Lammps like input to set put calculation'
     WRITE(*, *) '-pot_fit : Lammps like input for the potentials to be fitted '
     WRITE(*, *) '-pot_fix : Lammps like input to be appended, not touched by the optimization'
     WRITE(*, *) '-ener    : switch on the fitting of energies'
     WRITE(*, *) '-forces  : switch on the fitting of forces'
     WRITE(*, *) '-tensor  : switch on the fitting of a tensor'
     WRITE(*, *) '-compress <val> : activates ridge-regression'
     WRITE(*, *) '-pca            : activate principal components analysys.'
     WRITE(*, *) '-print_bi       : print bispectrum components'
     WRITE(*, *) '-out'
     WRITE(*, *) '-refine <iter> <temp> <thr>'
     STOP
  END IF

  DO L = 1, IARGC()

     CALL GETARG(L, COMMAND)

     IF (TRIM(COMMAND) .EQ. '-datas') THEN
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) DATAS
     END IF
     IF (TRIM(COMMAND) .EQ. '-inp') THEN
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) SYS%INP
     END IF
     IF (TRIM(COMMAND) .EQ. '-pot_fit') THEN
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) SYS%INP_FIT
     END IF
     IF (TRIM(COMMAND) .EQ. '-pot_fix') THEN
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) SYS%INP_FIX
     END IF
     IF (TRIM(COMMAND) .EQ. '-forces')   FIT_FORCES = .TRUE.
     IF (TRIM(COMMAND) .EQ. '-ener')     FIT_ENER   = .TRUE.
     IF (TRIM(COMMAND) .EQ. '-skip_fit') SKIP_FIT   = .TRUE.
     IF (TRIM(COMMAND) .EQ. '-print_bi') PRINT_BI   = .TRUE.
     IF (TRIM(COMMAND) .EQ. '-pca')      PCA        = .TRUE.
     IF (TRIM(COMMAND) .EQ. '-refine') THEN
        REFINE = .TRUE.
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) REFINE_MAXITER
        CALL GETARG(L + 2, COMMAND);   READ(COMMAND, *) REFINE_TEMP
        CALL GETARG(L + 3, COMMAND);   READ(COMMAND, *) THR_KERNEL
     END IF
     IF (TRIM(COMMAND) .EQ. '-compress') THEN
        CS = .TRUE.
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) CM_VAL
     END IF
     IF (TRIM(COMMAND) .EQ. '-E0') THEN
        !! -E0 0 -E0 1:
        CALL GETARG(L + 1, COMMAND);   READ(COMMAND, *) E0CS
     END IF

  END DO

  ITER = 0
  NEW_DATAS = DATAS

  DO WHILE ( (ITER .LE. REFINE_MAXITER .AND. REFINE) .OR. ITER .EQ. 0 )

     IF (ITER .GT. 0) THEN

        CALL EXECUTE_COMMAND_LINE('tail -n $( head -n 1 ' // TRIM(DATAS) // ' ) ' // TRIM(DATAS) // ' > new_datas')
        CALL EXECUTE_COMMAND_LINE('sed -i "1i $(( 1+$( head -n 1 ' // TRIM(DATAS) // ') ))" new_datas')

        NEW_DATAS = 'new_datas'
        OPEN(9, FILE = NEW_DATAS, ACCESS = 'APPEND')

        IF (FIT_ENER .AND. (.NOT. FIT_FORCES)) THEN
           WRITE(9, *) TRIM(SYS%DATA(1)%INP_DATA) // ' ', ITER, 'new_geo.xyz', ' new_geo.ener 1.0'
           CLOSE(9)
        END IF
        IF (FIT_ENER .AND. FIT_FORCES) THEN
           WRITE(9, *) TRIM(SYS%DATA(1)%INP_DATA) // ' ', ITER, 'new_geo.xyz', ' new_geo.ener new_geo.force 1.0'
           CLOSE(9)
        END IF
        IF ((.NOT. FIT_ENER) .AND. FIT_FORCES) THEN
           WRITE(9, *) TRIM(SYS%DATA(1)%INP_DATA) // ' ', ITER, 'new_geo.xyz', ' new_geo.force 1.0'
           CLOSE(9)
        END IF

     END IF

     CALL SYS%READ_SYS(NEW_DATAS, FIT_ENER, FIT_FORCES)
     CALL GET_LSMF_SNAP
     PRINT*,'BEFORE GET_CHI2'
     CALL GET_CHI2
     IF (REFINE) CALL REFINE_SNAP
     ITER = ITER + 1

  END DO
  WRITE(*, '("fitsnap.x done")')
  
END PROGRAM FITSNAP

SUBROUTINE GET_CHI2
  USE FIT_LAMMPS_CLASS
  USE COMMON_VAR
  USE LAMMPS
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_DOUBLE, C_PTR, C_INT
  IMPLICIT NONE
  REAL(8)     :: CHI_VAL_ENER, CHI_VAL_FX, CHI_VAL_FY, CHI_VAL_FZ
  INTEGER     :: L, I, K, N, V
  REAL(C_DOUBLE), POINTER :: ENER      => NULL()
  REAL(C_DOUBLE), POINTER :: FX(:)     => NULL()
  REAL(C_DOUBLE), POINTER :: FY(:)     => NULL()
  REAL(C_DOUBLE), POINTER :: FZ(:)     => NULL()
  REAL(C_DOUBLE), POINTER :: ID_DBL(:) => NULL()
  INTEGER, ALLOCATABLE    :: MAP(:), ID(:)

  CHI_VAL_ENER = 0.0D0
  CHI_VAL_FX   = 0.0D0
  CHI_VAL_FY   = 0.0D0
  CHI_VAL_FZ   = 0.0D0

  OPEN(222, FILE = 'energy_rms.dat')
  OPEN(333, FILE = 'force_rms.dat')
  WRITE(222, *) 'RMS Energies'
  WRITE(333, *) 'RMS Forces'

  !!      CALCOLA CHI2
  DO I = 1, SYS%NDATA

     DO L = 1, SYS%DATA(I)%FRAMES

        CALL LAMMPS_SCATTER_ATOMS (LMP(I), 'x', SYS%DATA(I)%X(L, 1:3*SYS%DATA(I)%NATS))
        CALL LAMMPS_COMMAND (LMP(I), 'run 0')

        IF (FIT_ENER) THEN

           CALL LAMMPS_EXTRACT_COMPUTE (ENER, LMP(I), 'pe_ener', 0, 0)

           CHI_VAL_ENER = CHI_VAL_ENER + ( (ENER - SYS%DATA(I)%ENER(L)) / SYS%DATA(I)%NATS )**2

           !WRITE(222, *) I, L, ENER / SYS%DATA(I)%NATS, SYS%DATA(I)%ENER(L) / SYS%DATA(I)%NATS, (ENER-SYS%DATA(I)%ENER(L)) / SYS%DATA(I)%NATS
           WRITE(222, *) I, L, ENER, SYS%DATA(I)%ENER(L), (ENER - SYS%DATA(I)%ENER(L)), SYS%DATA(I)%NATS

           ENER => NULL()

        END IF

        IF (FIT_FORCES) THEN

           CALL LAMMPS_EXTRACT_COMPUTE (FX, LMP(I), 'f_x', 1, 1)
           CALL LAMMPS_EXTRACT_COMPUTE (FY, LMP(I), 'f_y', 1, 1)
           CALL LAMMPS_EXTRACT_COMPUTE (FZ, LMP(I), 'f_z', 1, 1)
           IF (.NOT. ALLOCATED(ID))  ALLOCATE(ID(SYS%DATA(I)%NATS))
           IF (.NOT. ALLOCATED(MAP)) ALLOCATE(MAP(SYS%DATA(I)%NATS))
           CALL LAMMPS_EXTRACT_COMPUTE (ID_DBL, LMP(I), 'id', 1, 1)

           ID = INT(ID_DBL)
           ID_DBL => NULL()

           DO K = 1, SYS%DATA(I)%NATS
              MAP(ID(K)) = K
           END DO

           DO K = 1, SYS%DATA(I)%NATS
              CHI_VAL_FX = CHI_VAL_FX + ( FX(MAP(K)) - SYS%DATA(I)%FX(L, K) )**2
              CHI_VAL_FY = CHI_VAL_FY + ( FY(MAP(K)) - SYS%DATA(I)%FY(L, K) )**2
              CHI_VAL_FZ = CHI_VAL_FZ + ( FZ(MAP(K)) - SYS%DATA(I)%FZ(L, K) )**2
              WRITE(333, *) K, FX(MAP(K)), FY(MAP(K)), FZ(MAP(K)), SYS%DATA(I)%FX(L,K), SYS%DATA(I)%FY(L,K), SYS%DATA(I)%FZ(L,K)
           END DO

           FX => NULL()
           FY => NULL()
           FZ => NULL()

        END IF

     END DO      ! CICLO SU FRAMES

     CALL LAMMPS_CLOSE (LMP(I))
     IF (ALLOCATED(MAP)) DEALLOCATE(MAP)
     IF (ALLOCATED(ID))  DEALLOCATE(ID)
     PRINT*,'CALCULATING CHI2: LINE => ', I
  END DO   ! CICLO SU DATA

  IF (ALLOCATED(LMP)) DEALLOCATE(LMP)

  CHI_VAL_ENER = SQRT(CHI_VAL_ENER / SYS%TOT_FRAMES)
  !WRITE(222, *) 'Total RMS Energies (Kcal/mol/atom):', SQRT(CHI_VAL_ENER / SYS%TOT_FRAMES)
  WRITE(222, *) 'Total RMS Energies (Kcal/mol/atom):', CHI_VAL_ENER
  CLOSE(222)
  CHI_VAL_FX = SQRT(CHI_VAL_FX/SYS%TOT_FRAMES)
  CHI_VAL_FY = SQRT(CHI_VAL_FY/SYS%TOT_FRAMES)
  CHI_VAL_FZ = SQRT(CHI_VAL_FZ/SYS%TOT_FRAMES)
  !WRITE(333, *) 'Total RMS Forces (Kcal/mol/Ang): ', SQRT(CHI_VAL_FX/SYS%TOT_FRAMES), SQRT(CHI_VAL_FY/SYS%TOT_FRAMES), SQRT(CHI_VAL_FZ/SYS%TOT_FRAMES)
  WRITE(333, *) 'Total RMS Forces (Kcal/mol/Ang): ', CHI_VAL_FX, CHI_VAL_FY, CHI_VAL_FZ
  CLOSE(333)

  RETURN
END SUBROUTINE GET_CHI2
