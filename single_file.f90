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

     V = 1
     VV = 1

     do i=1,sys%ndata

        do l=1,sys%data(i)%frames

           call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
           call lammps_command (lmp(i), 'run 0')
           call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
           call lammps_extract_compute (ener, lmp(i), 'sna_e', 1, 2)
           call lammps_extract_compute (Eref, lmp(i), 'pe_ener', 0, 0)
           if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
           if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
           call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)

           id=INT(id_dbl)
           id_dbl=>null()

           do k=1,sys%data(i)%nats
              map(id(k))=k
           enddo

           if(l.eq.1)then

              npar=size(ener,1)+1
              kind_count=0

              do t=1,sys%data(i)%nats
                 kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
              enddo

              dimA=npar*nkinds
              dimB=npar2fit

              if(cs) dimB=dimB+(npar-1)*nkinds
              if(e0cs.eq.1) dimB=dimB+nkinds-1
              if(e0cs.eq.2) dimB=dimB+nkinds

              if(.not.allocated(B))then
                 allocate(B(dimB))
                 B=0.0d0
              endif
              if(.not.allocated(A))then
                 allocate(A(dimB,dimA))
                 A=0.0d0
              endif

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

              if(pca)then
                 open(1313,file='PCA.dat')
                 !             if(.not.allocated(A2)) allocate(A2(tot_frames,(npar-1)*nkinds))
                 !             if(.not.allocated(BB)) allocate(BB(tot_frames))
                 !             if(.not.allocated(W)) allocate(W((npar-1)*nkinds))
                 !             if(.not.allocated(Y)) allocate(Y(tot_frames,(npar-1)*nkinds))
                 !             if(.not.allocated(YY)) allocate(YY((npar-1)*nkinds,(npar-1)*nkinds))
                 !             if(.not.allocated(Tval)) allocate(Tval(tot_frames,(npar-1)*nkinds))
                 if(.not.allocated(A2)) allocate(A2(tot_frames,(npar-1)))
                 if(.not.allocated(BB)) allocate(BB(tot_frames))
                 if(.not.allocated(W)) allocate(W((npar-1)))
                 if(.not.allocated(Y)) allocate(Y(tot_frames,(npar-1)))
                 if(.not.allocated(YY)) allocate(YY((npar-1),(npar-1)))
                 if(.not.allocated(Tval)) allocate(Tval(tot_frames,(npar-1)))
              endif

           endif

           if(fit_ener)then

              s=1
              do k=1,nkinds
                 A(v,s)=kind_count(k)
                 s=s+npar
              enddo

              if(print_bi) write(121,*) ener(:,1),ener(:,2),ener(:,3)
              do k=2,npar
                 do t=1,sys%data(i)%nats


                    s=k+(nint(kind_nat(t))-1)*npar
                    A(v,s)=A(v,s)+ener(k-1,map(t))*sys%data(i)%weight
                    if(pca)then
                       if(kind_nat(t).eq.1)then
                          ss=k-1+(nint(kind_nat(t))-1)*(npar-1)
                          A2(vv,ss)=A(vv,ss)+ener(k-1,map(t))
                       endif
                    endif

                 enddo
              enddo

              B(v)=(sys%data(i)%ener(l)-Eref)*sys%data(i)%weight
              v=v+1
              vv=vv+1

              ener=>null()
              Eref=>null()

           endif ! fitener

           if(fit_forces )then

              call lammps_extract_compute (f, lmp(i), 'sna_f', 1, 2)
              call lammps_extract_compute (fx_ref, lmp(i), 'f_x', 1, 1)
              call lammps_extract_compute (fy_ref, lmp(i), 'f_y', 1, 1)
              call lammps_extract_compute (fz_ref, lmp(i), 'f_z', 1, 1)

              do t=1,sys%data(i)%nats

                 s=1
                 do k=1,nkinds
                    A(v,s)=0.0
                    s=s+npar
                 enddo

                 do n=1,nkinds
                    s=(n-1)*(3*(npar-1))+1  ! npar-1?
                    do k=2,npar
                       A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
                       s=s+1
                    enddo
                 enddo
                 B(v)=(sys%data(i)%fx(l,t)-fx_ref(map(t)))*sys%data(i)%weight
                 v=v+1

                 s=1
                 do k=1,nkinds
                    A(v,s)=0.0
                    s=s+npar
                 enddo

                 do n=1,nkinds
                    s=(n-1)*(3*(npar-1))+1+npar-1
                    do k=2,npar
                       A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
                       s=s+1
                    enddo
                 enddo
                 B(v)=(sys%data(i)%fy(l,t)-fy_ref(map(t)))*sys%data(i)%weight
                 v=v+1

                 s=1
                 do k=1,nkinds
                    A(v,s)=0.0
                    s=s+npar
                 enddo

                 do n=1,nkinds
                    s=(n-1)*(3*(npar-1))+1+2*(npar-1)
                    do k=2,npar
                       A(v,((n-1)*npar+k))=f(s,map(t))*sys%data(i)%weight
                       s=s+1
                    enddo
                 enddo
                 B(v)=(sys%data(i)%fz(l,t)-fz_ref(map(t)))*sys%data(i)%weight
                 v=v+1

              enddo

              f=>null()
              fx_ref=>null()
              fy_ref=>null()
              fz_ref=>null()

           endif  ! end if on forse
        enddo  ! ciclo on frames

        if(allocated(map)) deallocate(map)
        if(allocated(id)) deallocate(id)
        PRINT*,'COMPUTING BISPECTRUM: LINE => ', I, "FRAMES IN LINE: ", SYS%DATA(I)%FRAMES
     enddo   ! ciclo su data

     if(e0cs.eq.1)then
        s=1
        do k=1,nkinds-1
           B(v)=0.0d0
           A(v,s)=1.0d0
           s=s+npar
           v=v+1
        enddo
     endif

     if(e0cs.eq.2)then
        s=1
        do k=1,nkinds
           B(v)=0.0d0
           A(v,s)=1.0d2
           s=s+npar
           v=v+1
        enddo
     endif

     ! compressive sensing

     if(cs)then
        do k=1,nkinds
           do l=2,npar
              B(v)=0.0d0
              s=(k-1)*npar+l
              A(v,s)=cm_val
              v=v+1
           enddo
        enddo
     endif


     ! principal component analysis

     if(pca)then

        Y=A2
        do l=1,size(Y,2)
           ave=0.0d0
           do s=1,size(Y,1)
              ave=ave+Y(s,l)
           enddo
           Y(1:size(Y,1),l)=Y(1:size(Y,1),l)-ave/size(Y,1)
        enddo
        YY=matmul(transpose(Y),Y)
        call new_diag(size(YY,1),YY,W)

        if(.not.allocated(AA)) allocate(AA(tot_frames,4))

        AA=matmul(A2,YY(:,(size(W)-4):size(W)))

        W=sqrt(abs(W))
        W=W/sum(W)
        write(1313,*) 'Principal Components Values'
        do k=size(W),1,-1
           write(1313,*) W(k)
        enddo
        do l=1,size(AA,1)
           write(1313,*) (AA(l,k),k=1,4)
        enddo

!!!
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
!!!

     endif

     !        aaa=npar*nkinds
     !        if(cs)then
     !         bbb=npar2fit+(nkinds*(npar-1))!-1
     !         bbb=npar2fit+(nkinds*npar)-1
     !        else
     !         bbb=npar2fit!+nkinds!-1
     !         bbb=npar2fit
     !        endif
     lwork=dimB+64*dimB+1000
     allocate(work(lwork))
     call dgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
     deallocate(work)
     if(inf.ne.0)then
        write(*,*) 'zgels failed',inf
        stop
     else


        open(222,file='snapcoeff')

        l=1
        write(222,*) nkinds,npar
        do i =1,nkinds
           write(222,*) label(i),radii(i),cutoff(i)
           do n=1,npar
              write(222,*) B(l)
              l=l+1
           enddo
        enddo

        !        dump=0.0d0
        !        do i=npar*nkinds+1,size(B)
        !         dump=dump+B(i)**2
        !        enddo

        !        write(222,*) 'dgels residulal: ',dump
        close(222)

        open(333,file='snapparam')
        write(333,*) 'rcutfac ',gen_cutoff
        write(333,*) 'twojmax ',bi_order
        write(333,*) 'quadraticflag ',quadflag
        write(333,*) 'rfac0 1.00000'
        write(333,*) 'rmin0 0'
        write(333,*) 'diagonalstyle 3'
        write(333,*) 'switchflag 1'

        close(333)

        do i=1,sys%ndata
           call lammps_file (lmp(i),sys%inp_fix)
        enddo

     endif

  endif ! skip_fit

  return
end subroutine get_lsmf_snap

subroutine refine_snap
  use fit_lammps_class
  use lapack_diag_simm
  use lapack_inverse
  use random_numbers_class
  use common_var
  use LAMMPS
  use lists_class
  use meta_class
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
  implicit none

  integer                         :: l,i,j,v,k
  real(8), allocatable   :: r(:)

  character (len=10000)           :: snap_string
  real(8), allocatable   :: cutoff(:),radii(:)
  integer, allocatable            :: kind_count(:),type2fit(:)
  character (len=2),allocatable   :: label(:)
  integer                         :: nkinds,bi_order,npar,npar2fit,tot_frames
  integer                         :: nats2fit,quadflag
  real(8)                :: gen_cutoff,pi

  real (C_double), pointer        :: id_dbl(:)=> null()
  real (c_double), pointer        :: kind_nat(:) => null()
  integer, allocatable            :: map(:),id(:)

  real(8), allocatable   :: max_kernel(:)
  real(8)                :: dft_ener,val,rand_seed

  real (C_double), pointer   :: B(:,:) => null()
  real (c_double), pointer   :: ff_ener => null()

  allocate(lmp(1))

  pi=acos(-1.0d0)

  call lammps_open_no_mpi ('lmp -log none', lmp(1))
  call lammps_file (lmp(1),sys%inp)
  call lammps_command (lmp(1),'read_data '//trim(sys%data(1)%inp_data))
  call lammps_file (lmp(1),sys%inp_fix)
  call lammps_command (lmp(1), 'compute pe_ener all pe')

  open(16,file=sys%inp_fit)
  read(16,*) gen_cutoff,bi_order,npar,quadflag
  read(16,*) nkinds
  allocate(label(nkinds))
  allocate(type2fit(nkinds))
  allocate(cutoff(nkinds))
  allocate(radii(nkinds))
  allocate(kind_count(nkinds))
  kind_count=0
  do i=1,nkinds
     read(16,*) label(i),type2fit(i),radii(i),cutoff(i)
     write(*,*) label(i),type2fit(i),radii(i),cutoff(i)
  enddo
  close(16)
  write(snap_string,*) gen_cutoff,'1.0000',bi_order,(radii(i),i=1,nkinds),(cutoff(i),i=1,nkinds),&
       'quadraticflag ',quadflag

  call lammps_command (lmp(1), &
       'compute sna_e all sna/atom '//trim(snap_string)//&
       ' diagonal 3 rmin0 0 switchflag 1')
  call lammps_command (lmp(1),&
       'compute type all property/atom type')
  call lammps_command (lmp(1),&
       'compute id all property/atom id')

  ! minimize energy

  call lammps_command (lmp(1),'thermo 1')
  call lammps_command (lmp(1),&
       'thermo_style custom step time temp pe etotal ')
  !         call lammps_command (lmp(1), &
  !                'minimize 1.0e-8 1.0e-8 1000 100000')

  write(snap_string,*) (label(i)//' ',i=1,size(label))

  call lammps_command (lmp(1),&
       'dump xyz_dump all xyz 5 geo_opt.xyz')

  call lammps_command (lmp(1),&
       'dump_modify            xyz_dump element '//trim(snap_string))

  ! set velocities

  write(snap_string,*) refine_temp

  call lammps_command (lmp(1),'timestep 0.25')
  call lammps_command (lmp(1),'variable t equal'//trim(snap_string))

  call random_number(rand_seed)
  write(snap_string,*) nint(rand_seed*10000.0d0)

  call lammps_command (lmp(1), &
       'velocity all create $t '//trim(snap_string)//&
       ' dist gaussian')
  call lammps_command (lmp(1),'velocity all zero linear')
  call lammps_command (lmp(1),'velocity all zero angular')


  call lammps_command (lmp(1),'group atom1 id 1')
  call lammps_command (lmp(1),'group atom2 id 4')
  !         call lammps_command (lmp(1),'group atom3 id 7')
  call meta1%gauss%init()
  !         call meta2%gauss%init()

  call lammps_command (lmp(1), &
       'fix 1 all nvt temp $t $t 100.0 tchain 3')

  open(13,file='new_geo.xyz',access='append')
  open(14,file='kernel_max.dat',access='append')
  open(15,file='new_geo.ener',access='append')
  open(16,file='meta.dat',access='append')

  if(.not.allocated(max_kernel)) &
       allocate(max_kernel(sys%data(1)%nats))
  !         if(.not.allocated(max_kernel)) &
  !                allocate(max_kernel(kernel%nkinds))

  write(14,*) '## nkinds: ',kernel%nkinds,'nenvs: ',(kernel%K(i)%nenvs,i=1,kernel%nkinds)
  flush(14)

  do i=1,400000

     call lammps_command (lmp(1), 'run 4')
     call lammps_extract_compute (B, lmp(1), 'sna_e', 1, 2)
     call lammps_extract_compute (ff_ener, lmp(1), 'pe_ener', 0, 0)
     call lammps_extract_compute (kind_nat, lmp(1), 'type', 1, 1)
     if(.not. allocated(id)) allocate(id(sys%data(1)%nats))
     if(.not. allocated(map)) allocate(map(sys%data(1)%nats))
     call lammps_extract_compute (id_dbl, lmp(1), 'id', 1, 1)

     if(do_meta)then

        if(allocated(r)) deallocate(r)
        call lammps_gather_atoms (lmp(1),'x', 3, r)

        dist1=(r(1)-r(10))**2
        dist1=dist1+(r(2)-r(11))**2
        dist1=dist1+(r(3)-r(12))**2
        dist1=sqrt(dist1)

        !           dist2=(r(4)-r(7))**2
        !           dist2=dist2+(r(5)-r(8))**2
        !           dist2=dist2+(r(6)-r(9))**2
        !           dist2=sqrt(dist2)

        f1=0.0d0
        f2=0.0d0
        !           f3=0.0d0

        if(meta1%gauss%nelem.ge.1)then

           call meta1%gauss%reboot()
           !            call meta2%gauss%reboot()

           do j=1,meta1%gauss%nelem

              call meta1%gauss%rd_val(height1)
              !             call meta2%gauss%rd_val(height2)

              !             f1(1)=f1(1)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !             f1(2)=f1(2)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !             f1(3)=f1(3)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1

              !             f2(1)=f2(1)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !             f2(2)=f2(2)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !             f2(3)=f2(3)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2)*&
              !                 (dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1

              !             f1(1)=f1(1)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !             f1(2)=f1(2)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !             f1(3)=f1(3)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1

              !             f2(1)=f2(1)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(4)-r(7))/meta1%sigma**2/dist1
              !             f2(2)=f2(2)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(5)-r(8))/meta1%sigma**2/dist1
              !             f2(3)=f2(3)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist1-height1)*(r(6)-r(9))/meta1%sigma**2/dist1

              !             f1(1)=f1(1)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(4)-r(19))/meta1%sigma**2/dist2
              !             f1(2)=f1(2)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(5)-r(20))/meta1%sigma**2/dist2
              !             f1(3)=f1(3)+2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(6)-r(21))/meta1%sigma**2/dist2

              !             f3(1)=f3(1)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(4)-r(19))/meta1%sigma**2/dist2
              !             f3(2)=f3(2)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(5)-r(20))/meta1%sigma**2/dist2
              !             f3(3)=f3(3)-2*meta1%weight*&
              !                 exp(-(height1-dist1)**2/meta1%sigma**2-(height2-dist2)**2/meta2%sigma**2)*&
              !                 (dist2-height2)*(r(6)-r(21))/meta1%sigma**2/dist2

              call meta1%gauss%skip()
              !             call meta2%gauss%skip()

           enddo

        endif

        !! contraint distances
        if(dist1.gt.8.0d0 .and. dist1.lt.16.0d0)then
           f1(1)=f1(1)+2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(1)-r(10))/dist1
           f1(2)=f1(2)+2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(2)-r(11))/dist1
           f1(3)=f1(3)+2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(3)-r(12))/dist1
           f2(1)=f2(1)-2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(1)-r(10))/dist1
           f2(2)=f2(2)-2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(2)-r(11))/dist1
           f2(3)=f2(3)-2000.d0*sin(pi*dist1/8.0d0)*pi/8.0d0*(r(3)-r(12))/dist1
        endif

        if(i.ge.2)then
           call lammps_command (lmp(1), 'unfix add1')
           call lammps_command (lmp(1), 'unfix add2')
           !             call lammps_command (lmp(1), 'unfix add3')
        endif

        write(snap_string,"(3F12.6)") f1(1),f1(2),f1(3)
        call lammps_command (lmp(1), 'fix add1 atom1 addforce '&
             //trim(snap_string))
        write(snap_string,"(3F12.6)") f2(1),f2(2),f2(3)
        call lammps_command (lmp(1), 'fix add2 atom2 addforce '&
             //trim(snap_string))
        !            write(snap_string,"(3F12.6)") f3(1),f3(2),f3(3)
        !            call lammps_command (lmp(1), 'fix add3 atom3 addforce '&
        !                                 //trim(snap_string))

        if(mod(dble(i),10.0d0).lt.1.0e-6)then
           call meta1%gauss%add_node(dist1)
           !            call meta2%gauss%add_node(dist2)
           write(16,*) i,dist1!,dist2
        endif

     endif

     id=INT(id_dbl)
     id_dbl=>null()

     do k=1,sys%data(1)%nats
        map(id(k))=k
     enddo

     max_kernel=0.0d0

     do k=1,size(B,2)
        l=nint(kind_nat(map(k)))
        do v=1,kernel%K(l)%nenvs
           val=0.0d0
           do j=1,size(B,1)
              val=val-(B(j,map(k))-kernel%K(l)%B(v,j))**2
           enddo
           val=exp(val/2*kernel%K(l)%sigma**2)
           if(val.gt.max_kernel(map(k))) max_kernel(map(k))=val
        enddo
     enddo

     write(14,*) i,max_kernel
     flush(14)

     if(any(max_kernel.lt.thr_kernel))then
        if(allocated(r)) deallocate(r)
        call lammps_gather_atoms (lmp(1),'x', 3, r)
        write(13,*) sys%data(1)%nats
        write(13,*)
        l=1
        do j=1,sys%data(1)%nats
           write(13,*) sys%data(1)%label(j),r(l),r(l+1),r(l+2)
           l=l+3
        enddo
        flush(13)
        call execute_command_line('./run_DFT_scf.x')
        open(16,file='last_dft_ener.dat')
        read(16,*) dft_ener
        close(16)
        write(15,*) dft_ener,ff_ener
        flush(15)
        close(13)
        close(14)
        close(15)
        call lammps_close (lmp(1))
        deallocate(lmp)
        return
        write(*,*) 'PUPPA'
        flush(6)
     endif

  enddo

  call lammps_close (lmp(1))
  deallocate(lmp)
  refine=.false.
  close(13)
  close(14)
  close(15)

  return
end subroutine refine_snap

program fitsnap
  use fit_lammps_class
  use common_var
  implicit none
  integer                        :: l
  character (len=100)            :: command,input,datas,output,new_datas

  if(iargc().eq.0)then
     write(*,*) 'FitSnap Usage:'
     write(*,*) '-datas : list of files to be fitted'
     write(*,*) '-inp : Lammps like input to set put calculation'
     write(*,*) '-pot_fit : Lammps like input for the potentials to be &
          fitted '
     write(*,*) '-pot_fix : Lammps like input to be appended, not &
          touched by the optimization'
     write(*,*) '-ener   : switch on the fitting of energies'
     write(*,*) '-forces : switch on the fitting of forces'
     write(*,*) '-tensor : switch on the fitting of a tensor'
     write(*,*) '-compress <val> : activates ridge-regression'
     write(*,*) '-pca            : activate principal components &
          analysys.'
     write(*,*) '-print_bi       : print bispectrum components'
     write(*,*) '-out'
     write(*,*) '-refine <iter> <temp> <thr>'
     stop
  endif

  do l=1,iargc()

     call getarg(l,command)

     if(trim(command).eq.'-datas')then
        call getarg(l+1,command)
        read(command,*) datas
     endif
     if(trim(command).eq.'-inp')then
        call getarg(l+1,command)
        read(command,*) sys%inp
     endif
     if(trim(command).eq.'-pot_fit')then
        call getarg(l+1,command)
        read(command,*) sys%inp_fit
     endif
     if(trim(command).eq.'-pot_fix')then
        call getarg(l+1,command)
        read(command,*) sys%inp_fix
     endif
     if(trim(command).eq.'-forces')then
        fit_forces=.true.
     endif
     if(trim(command).eq.'-ener')then
        fit_ener=.true.
     endif
     if(trim(command).eq.'-skip_fit')then
        skip_fit=.true.
     endif
     if(trim(command).eq.'-print_bi')then
        print_bi=.true.
     endif
     if(trim(command).eq.'-pca')then
        pca=.true.
     endif
     if(trim(command).eq.'-refine')then
        refine=.true.
        call getarg(l+1,command)
        read(command,*) refine_maxiter
        call getarg(l+2,command)
        read(command,*) refine_temp
        call getarg(l+3,command)
        read(command,*) thr_kernel
     endif
     if(trim(command).eq.'-compress')then
        cs=.true.
        call getarg(l+1,command)
        read(command,*) cm_val
     endif
     if(trim(command).eq.'-E0')then
!!!! -E0 0 -E0 1:
        call getarg(l+1,command)
        read(command,*) e0cs
     endif

  enddo

  iter=0
  new_datas=datas

  do while ( (iter.le.refine_maxiter .and. refine) .or. iter.eq.0 )

     if(iter.gt.0)then

        call execute_command_line('tail -n $( head -n 1 '&
             //trim(datas)//' ) '//trim(datas)//' > new_datas')
        call execute_command_line('sed -i "1i $(( 1+$( head -n 1 '&
             //trim(datas)//') ))" new_datas'  )

        new_datas='new_datas'
        open(9,file=new_datas,access='append')

        if(fit_ener .and. (.not. fit_forces))then
           write(9,*) trim(sys%data(1)%inp_data)//' ',&
                iter,'new_geo.xyz',' new_geo.ener 1.0'
           close(9)
        endif
        if(fit_ener .and. fit_forces)then
           write(9,*) trim(sys%data(1)%inp_data)//' ',&
                iter,'new_geo.xyz',&
                ' new_geo.ener new_geo.force 1.0'
           close(9)
        endif
        if((.not. fit_ener) .and. fit_forces)then
           write(9,*) trim(sys%data(1)%inp_data)//' ',&
                iter,'new_geo.xyz',&
                ' new_geo.force 1.0'
           close(9)
        endif

     endif

     call sys%read_sys(new_datas,fit_ener,fit_forces)
     call get_lsmf_snap
     PRINT*,'BEFORE GET_CHI2'
     call get_chi2
     if(refine) call refine_snap
     iter=iter+1

  enddo

  return
end program fitsnap

subroutine get_chi2
  use fit_lammps_class
  use common_var
  use LAMMPS
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
  implicit none
  real(8)     :: chi_val_ener,chi_val_fx,chi_val_fy,chi_val_fz
  integer              :: l,i,k,n,v
  real (C_double), pointer :: ener => null()
  real (C_double), pointer :: fx(:) => null()
  real (C_double), pointer :: fy(:) => null()
  real (C_double), pointer :: fz(:) => null()
  real (C_double), pointer   :: id_dbl(:)=> null()
  integer, allocatable       :: map(:),id(:)

  chi_val_ener=0.0d0
  chi_val_fx=0.0d0
  chi_val_fy=0.0d0
  chi_val_fz=0.0d0

  open(222,file='energy_rms.dat')
  open(333,file='force_rms.dat')
  write(222,*) 'RMS Energies'
  write(333,*) 'RMS Forces'

  !!      calcola chi2
  do i=1,sys%ndata

     do l=1,sys%data(i)%frames

        call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
        call lammps_command (lmp(i), 'run 0')

        if(fit_ener)then

           call lammps_extract_compute (ener, lmp(i), 'pe_ener', 0, 0)

           chi_val_ener=chi_val_ener+( (ener-sys%data(i)%ener(l))/sys%data(i)%nats )**2

           !            write(222,*) i,l,ener/sys%data(i)%nats,sys%data(i)%ener(l)/sys%data(i)%nats,&
           !                 (ener-sys%data(i)%ener(l))/sys%data(i)%nats
           write(222,*) i, l, ener, sys%data(i)%ener(l), (ener-sys%data(i)%ener(l)), sys%data(i)%nats

           ener=>null()

        endif

        if(fit_forces)then

           call lammps_extract_compute (fx, lmp(i), 'f_x', 1, 1)
           call lammps_extract_compute (fy, lmp(i), 'f_y', 1, 1)
           call lammps_extract_compute (fz, lmp(i), 'f_z', 1, 1)
           if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
           if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
           call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)

           id=INT(id_dbl)
           id_dbl=>null()

           do k=1,sys%data(i)%nats
              map(id(k))=k
           enddo

           do k=1,sys%data(i)%nats
              chi_val_fx=chi_val_fx+( fx(map(k))-sys%data(i)%fx(l,k) )**2
              chi_val_fy=chi_val_fy+( fy(map(k))-sys%data(i)%fy(l,k) )**2
              chi_val_fz=chi_val_fz+( fz(map(k))-sys%data(i)%fz(l,k) )**2
              write(333,*) k,fx(map(k)),fy(map(k)),fz(map(k)),sys%data(i)%fx(l,k),sys%data(i)%fy(l,k),sys%data(i)%fz(l,k)
           enddo

           fx=>null()
           fy=>null()
           fz=>null()

        endif

     enddo      ! ciclo su frames

     call lammps_close (lmp(i))
     if(allocated(map)) deallocate(map)
     if(allocated(id)) deallocate(id)
     PRINT*,'CALCULATING CHI2: LINE => ', I
  enddo   ! ciclo su data

  if(allocated(lmp)) deallocate(lmp)

  write(222,*) 'Total RMS Energies (Kcal/mol/atom):',sqrt(chi_val_ener/sys%tot_frames)
  write(333,*) 'Total RMS Forces (Kcal/mol/Ang): ',sqrt(chi_val_fx/sys%tot_frames),&
       sqrt(chi_val_fy/sys%tot_frames),&
       sqrt(chi_val_fz/sys%tot_frames)
  close(222)
  close(333)

  return
end subroutine get_chi2
