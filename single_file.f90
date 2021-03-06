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

  subroutine dealloc_kernel_kind(this)
    implicit none
    class(kernel_kind)    :: this
    if( allocated(this%B) ) deallocate(this%B)
    this%nenvs=0
    this%sigma=0
    return
  end subroutine dealloc_kernel_kind

  subroutine read_sys(sys,data_file,fit_ener,fit_forces)
    implicit none
    class(system)           :: sys
    integer                 :: i,l,v,k,m,n
    character(len=100)      :: label,data_file
    logical                 :: fit_forces,fit_ener

    if(allocated(sys%data))then
       do i=1,sys%ndata
          if(allocated(sys%data(i)%x)) deallocate(sys%data(i)%x)
          if(allocated(sys%data(i)%ener)) deallocate(sys%data(i)%ener)
          if(allocated(sys%data(i)%fx)) deallocate(sys%data(i)%fx)
          if(allocated(sys%data(i)%fy)) deallocate(sys%data(i)%fy)
          if(allocated(sys%data(i)%fz)) deallocate(sys%data(i)%fz)
       enddo
       deallocate(sys%data)
    endif

    open(8,file=trim(data_file))

    sys%tot_frames=0
    sys%npar2fit=0

    read(8,*) sys%ndata

    allocate(sys%data(sys%ndata))

    do i=1,sys%ndata

       if(fit_forces .and. fit_ener)then
          read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj,sys%data(i)%inp_ener,&
               sys%data(i)%inp_forces,sys%data(i)%weight
       endif

       if(fit_ener .and.  (.not. fit_forces) )then
          read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj, &
               sys%data(i)%inp_ener,sys%data(i)%weight
       endif

       if((.not. fit_ener) .and. fit_forces )then
          read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj, &
               sys%data(i)%inp_forces,sys%data(i)%weight
       endif

       open(12,file=sys%data(i)%inp_traj)

       if(fit_ener)then
          allocate(sys%data(i)%ener(sys%data(i)%frames))
          open(13,file=sys%data(i)%inp_ener)
       endif

       if(fit_forces)then
          open(14,file=sys%data(i)%inp_forces)
       endif

       do l=1,sys%data(i)%frames

          sys%tot_frames=sys%tot_frames+1

          read(12,*) sys%data(i)%nats
          read(12,*)
          if(fit_forces)then
             read(14,*)
             read(14,*)
          endif

          if(.not.allocated(sys%data(i)%label))then
             allocate(sys%data(i)%label(sys%data(i)%nats))
          endif
          if(.not.allocated(sys%data(i)%x))then
             allocate(sys%data(i)%x(sys%data(i)%frames,3*sys%data(i)%nats))
          endif
          if(.not.allocated(sys%data(i)%fx) .and. fit_forces)then
             allocate(sys%data(i)%fx(sys%data(i)%frames,sys%data(i)%nats))
          endif
          if(.not.allocated(sys%data(i)%fy) .and. fit_forces)then
             allocate(sys%data(i)%fy(sys%data(i)%frames,sys%data(i)%nats))
          endif
          if(.not.allocated(sys%data(i)%fz) .and. fit_forces)then
             allocate(sys%data(i)%fz(sys%data(i)%frames,sys%data(i)%nats))
          endif

          v=1
          do k=1,sys%data(i)%nats

             read(12,*) sys%data(i)%label(k),sys%data(i)%x(l,v),sys%data(i)%x(l,v+1),sys%data(i)%x(l,v+2)

             if(fit_forces)then
                read(14,*) sys%data(i)%fx(l,k),sys%data(i)%fy(l,k),sys%data(i)%fz(l,k)
                sys%npar2fit=sys%npar2fit+3
             endif
             v=v+3

          enddo

          if(fit_ener)then
             read(13,*) sys%data(i)%ener(l)
             sys%npar2fit=sys%npar2fit+1
          endif

       enddo

       close(12)
       if(fit_ener) close(13)
       if(fit_forces) close(14)

    enddo

    close(8)

    return
  end subroutine read_sys


end module fit_lammps_class

module common_var
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
  use fit_lammps_class
  !
  TYPE(system)                   :: sys
  type (C_ptr), allocatable      :: lmp(:)
  type(kernel_global)            :: kernel
  logical                        :: fit_forces=.false.,pca=.false.,cs=.false.,fit_ener=.false.
  logical                        :: fit_tensors=.true.
  logical                        :: skip_fit=.false.,print_bi=.false.,refine=.false.,metrix=.true.
  real(8)               :: cm_val=1.0d0,refine_temp,thr_kernel=0.5d0
  integer                        :: refine_maxiter,iter,tens_order,e0cs=0
  !
end module common_var

module lapack_inverse
  implicit none

contains

  subroutine mat_inv(X,N)
    implicit none
    integer                                         :: info,lwork,N,ialloc
    integer, dimension(N)                           :: ipiv
    real(8),  dimension(N,N)               :: X_inv
    double complex,  dimension(N,N)                 :: cX_inv
    CLASS(*), intent(in), dimension(N,N)            :: X
    real(8), allocatable, dimension(:)     :: work


    select type(X)

    type is (complex(8))

       cX_inv=X

       lwork=(N)**2
       allocate(work(lwork),stat=ialloc)
       call zgetrf(N,N,cX_inv,N,IPIV,info)
       call zgetri(N,cX_inv,N,IPIV,work,lwork,info)

       X=cX_inv

    type is (real(8))

       X_inv=X

       lwork=(N)**2
       allocate(work(lwork),stat=ialloc)
       call dgetrf(N,N,X_inv,N,IPIV,info)
       call dgetri(N,X_inv,N,IPIV,work,lwork,info)

       X=X_inv

    end select

    return
  end subroutine mat_inv

end module lapack_inverse

module lapack_diag_simm
  implicit none

contains

  subroutine new_diag(N,A,W)
    implicit none
    integer                                      :: i,j,INFO,s,k,N,ialloc
    integer                                      :: l,inf,infr
    CLASS(*),  dimension(N,N)                    :: A
    DOUBLE COMPLEX,    DIMENSION(N*(N+1)/2)      :: AP
    REAL(8),  DIMENSION(N)              :: W
    REAL(8), ALLOCATABLE,  DIMENSION(:) :: work
    DOUBLE COMPLEX,   ALLOCATABLE,  DIMENSION(:) :: cwork


    select type(A)

    type is (real(8))


       l=(3*N-1)
       allocate(work(l),stat=ialloc)
       call dsyev('V','U',N,A,N,W,work,l,infr)
       deallocate(work)
       if(infr.ne.0)then
          write(*,*) 'dsyev diagonalization failed'
          FLUSH(6)
          stop
       endif

    type is (complex(8))

       s=1
       do j=1,N
          do i=1,j
             AP(s)=A(i,j)
             s=s+1
          enddo
       enddo

       l=(2*N-1)+1000
       allocate(cwork(l),stat=ialloc)
       k=(3*N-2)+1000
       allocate(work(k),stat=ialloc)
       call zhpev('V','U',N,AP,W,A,N,cwork,work,inf)
       deallocate(cwork)
       deallocate(work)
       if(inf.ne.0)then
          write(*,*) 'zhpev diagonalization failed'
          stop
       endif

    end select

    return
  end subroutine new_diag

end module lapack_diag_simm

module lapack_diag_asimm
  implicit none

contains

  subroutine new_diag2(N,A,W)
    implicit none
    integer                                      :: i,j,INFO,s,k,N,ialloc
    integer                                      :: l,inf,infr
    CLASS(*),  dimension(N,N)                    :: A
    DOUBLE COMPLEX,    DIMENSION(N*(N+1)/2)      :: AP
    DOUBLE COMPLEX,    DIMENSION(N,N)            :: VRI,VLI
    DOUBLE COMPLEX,    DIMENSION(N)              :: W
    REAL(8),  DIMENSION(N)              :: WR,WI
    REAL(8),  DIMENSION(N,N)            :: VL,VR
    REAL(8), ALLOCATABLE,  DIMENSION(:) :: work
    DOUBLE COMPLEX,   ALLOCATABLE,  DIMENSION(:) :: cwork


    select type(A)

    type is (real(8))


       l=4*N
       allocate(work(l),stat=ialloc)
       call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,l,infr)
       deallocate(work)
       if(infr.ne.0)then
          write(*,*) 'dgeev diagonalization failed'
          FLUSH(6)
          stop
       endif

       do j=1,N
          W(j)=CMPLX(WR(j),WI(j),8)
       enddo

       A=VR

    type is (complex(8))

       l=6*N
       allocate(cwork(l),stat=ialloc)
       allocate(work(l),stat=ialloc)
       call zgeev('V','V',N,A,N,W,VLI,N,VRI,N,cwork,l,work,infr)
       deallocate(cwork)
       deallocate(work)
       if(inf.ne.0)then
          write(*,*) 'zgeev diagonalization failed'
       endif
       A=VLI

    end select

    return
  end subroutine new_diag2

end module lapack_diag_asimm

module lists_class
  implicit none


  type :: list_node
     class(*), pointer      :: key => null()
     class(list_node), pointer  :: next => null()
     class(list_node), pointer  :: prev => null()
  end type list_node

  type list
     class(list_node), pointer :: head => null()
     class(list_node), pointer :: tail => null()
     class(list_node), pointer :: node => null()
     integer                   :: nelem
   contains
     procedure                 :: init => init_list
     procedure                 :: skip => skip_node
     procedure                 :: rew  => rewind_node
     procedure                 :: rm => remove_node
     procedure                 :: reboot => reboot_list
     procedure                 :: delete => destroy_list
     procedure                 :: add_node
     procedure                 :: rd_dbl_node
     procedure                 :: rd_cmplx_node
     procedure                 :: rd_int_node
     generic                   :: rd_val => rd_int_node,rd_dbl_node,rd_cmplx_node
  end type list


contains



!!!!!   lists general functions

  subroutine skip_node(this_list)
    implicit none
    class(list)  :: this_list
    if( associated(this_list%node%next) ) then
       this_list%node=>this_list%node%next
    else
       this_list%node=>this_list%tail
    endif
    return
  end subroutine skip_node

  subroutine destroy_list(this_list)
    implicit none
    class(list) :: this_list
    do while (this_list%nelem.gt.0)
       this_list%node=>this_list%head
       call this_list%rm()
    enddo
    return
  end subroutine destroy_list

  subroutine remove_node(this_list)
    implicit none
    class(list)  :: this_list
    class(list_node), pointer :: tmp_node

    if (this_list%nelem .eq. 0) return

    if (this_list%nelem .eq. 1)then

       if(associated(this_list%node%key)) then
          deallocate(this_list%node%key)
       endif
       if(associated(this_list%node)) then
          deallocate(this_list%node)
       endif
       this_list%node=>null()
       this_list%head=>null()
       this_list%tail=>null()

       this_list%nelem=this_list%nelem-1

    else

       if( .not. associated(this_list%node%next) .and. &
            associated(this_list%node%prev) ) then
          this_list%tail=>this_list%node%prev
          tmp_node=>this_list%node
          this_list%node=>this_list%tail
          this_list%tail%next=>null()
          if(associated(tmp_node%key)) then
             deallocate(tmp_node%key)
             tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
             deallocate(tmp_node)
             tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
       endif


       if( .not. associated(this_list%node%prev) .and. &
            associated(this_list%node%next) ) then
          this_list%head=>this_list%node%next
          tmp_node=>this_list%node
          this_list%node=>this_list%head
          this_list%head%prev=>null()
          if(associated(tmp_node%key)) then
             deallocate(tmp_node%key)
             tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
             deallocate(tmp_node)
             tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
       endif

       if( associated(this_list%node%next) .and. &
            associated(this_list%node%prev) ) then
          this_list%node%prev%next=>this_list%node%next
          this_list%node%next%prev=>this_list%node%prev
          tmp_node=>this_list%node
          tmp_node%key=>this_list%node%key
          this_list%node=>this_list%node%next
          if(associated(tmp_node%key)) then
             deallocate(tmp_node%key)
             tmp_node%key=>null()
          endif
          if(associated(tmp_node)) then
             deallocate(tmp_node)
             tmp_node=>null()
          endif
          this_list%nelem=this_list%nelem-1
       endif

    endif

    return
  end subroutine remove_node


  subroutine rewind_node(this_list)
    implicit none
    class(list)  :: this_list
    if( associated(this_list%node%prev) ) then
       this_list%node=>this_list%node%prev
    else
       this_list%node=>this_list%head
    endif
    return
  end subroutine rewind_node


  subroutine reboot_list(this_list)
    implicit none
    class(list)   :: this_list
    if ( associated(this_list%node) &
         .and. associated(this_list%head )  ) then
       this_list%node=>this_list%head
    endif
    return
  end subroutine reboot_list


  subroutine last_node_list(this_list)
    implicit none
    class(list)   :: this_list
    this_list%node=>this_list%tail
    return
  end subroutine last_node_list


  subroutine init_list(this_list)
    implicit none
    class(list)    :: this_list
    this_list%nelem=0
    return
  end subroutine init_list


  subroutine add_node(this_list,val)
    implicit none
    class(list)                    :: this_list
    class(*),pointer               :: arrow
    class(*),optional              :: val
    class(list_node),pointer       :: tmp_node

    if(this_list%nelem.eq.0)then
       allocate(this_list%head)
       allocate(this_list%node)
       if(present(val))then
          select type (val)
          type is (integer)
             allocate(integer::this_list%head%key)
             allocate(integer::this_list%node%key)
          type is (real(8))
             allocate(real(8)::this_list%head%key)
             allocate(real(8)::this_list%node%key)
          type is (complex(8))
             allocate(complex(8)::this_list%head%key)
             allocate(complex(8)::this_list%node%key)
          end select
       endif
       this_list%node=>this_list%head
       this_list%tail=>this_list%head
       this_list%nelem=this_list%nelem+1
       tmp_node=>this_list%head
    else
       allocate(tmp_node)
       if(present(val))then
          select type (val)
          type is (integer)
             allocate(integer::tmp_node%key)
          type is (real(8))
             allocate(real(8)::tmp_node%key)
          type is (complex(8))
             allocate(complex(8)::tmp_node%key)
          end select
       endif
       this_list%tail%next=>tmp_node
       tmp_node%prev=>this_list%tail
       this_list%tail=>tmp_node
       this_list%nelem=this_list%nelem+1
    endif



    if ( present(val) ) then
       select type (val)

       type is (integer)
          select type (arrow=>tmp_node%key)
          type is (integer)
             arrow=val
          end select
       type is (real(8))
          select type (arrow=>tmp_node%key)
          type is (real(8))
             arrow=val
          end select
       type is (complex(8))
          select type (arrow=>tmp_node%key)
          type is (complex(8))
             arrow=val
          end select

       end select
    endif

    tmp_node=>null()

    return
  end subroutine add_node


  subroutine rd_int_node(this,val)
    implicit none
    class(list)           :: this
    integer               :: val
    class(*),pointer      :: bho

    select type (bho=>this%node%key)
    type is (integer)
       val=bho
    end select

    return
  end subroutine rd_int_node

  subroutine rd_dbl_node(this,val)
    implicit none
    class(list)           :: this
    real(8)      :: val
    class(*),pointer      :: bho

    select type (bho=>this%node%key)
    type is (real(8))
       val=bho
    end select

    return
  end subroutine rd_dbl_node

  subroutine rd_cmplx_node(this,val)
    implicit none
    class(list)           :: this
    complex(8)            :: val
    class(*),pointer      :: bho

    select type (bho=>this%node%key)
    type is (complex(8))
       val=bho
    end select

    return
  end subroutine rd_cmplx_node


end module lists_class

module meta_class
  use lists_class
  implicit none

  type memory_pot
     type(list)           :: gauss
     real(8)     :: weight=1.0
     real(8)     :: sigma=0.2
  end type memory_pot

  type(memory_pot)                :: meta1,meta2
  logical                         :: do_meta=.true.
  real(8)                :: dist1,dist2,height1,height2
  real(8)                :: f1(3),f2(3),f3(3)

end module meta_class

subroutine get_lsmf_snap
  use fit_lammps_class
  use lapack_diag_simm
  use lapack_inverse
  use common_var
  use LAMMPS
  use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
  implicit none
  integer                         :: aaa,bbb,l,i,k,n,v,s,t,ll,tt
  integer                         :: ss,vv,inf,LWORK,dimA,dimB
  real(8), allocatable   :: x(:),B(:),A(:,:),work(:),AA(:,:)
  real(8), allocatable   :: W(:),Tval(:,:),YY(:,:),Y(:,:),A2(:,:)
  real(8), allocatable   :: ACM(:,:),BCM(:),BB(:)
  real(8)                :: dump,ave

  real (C_double), pointer   :: ener(:,:) => null()
  real (C_double), pointer   :: f(:,:) => null()
  real (c_double), pointer   :: kind_nat(:) => null()
  real (c_double), pointer   :: Eref => null()
  real (C_double), pointer   :: fx_ref(:) => null()
  real (C_double), pointer   :: fy_ref(:) => null()
  real (C_double), pointer   :: fz_ref(:) => null()
  real (C_double), pointer   :: id_dbl(:)=> null()
  integer, allocatable       :: map(:),id(:)

  integer                         :: nkinds,bi_order,npar,npar2fit,tot_frames
  integer                         :: nats2fit,quadflag
  real(8)                :: gen_cutoff
  character (len=10000)           :: snap_string,snap_string2
  real(8), allocatable   :: cutoff(:),radii(:),sigma(:)
  integer, allocatable            :: kind_count(:),type2fit(:)
  character (len=2),allocatable   :: label(:)


  tot_frames=sys%tot_frames
  npar2fit=sys%npar2fit

  if(print_bi) open(121,file='Bi_compoenents.dat')

  open(16,file=sys%inp_fit)
  read(16,*) gen_cutoff,bi_order,npar,quadflag
  read(16,*) nkinds
  allocate(label(nkinds))
  allocate(type2fit(nkinds))
  allocate(cutoff(nkinds))
  allocate(sigma(nkinds))
  allocate(radii(nkinds))
  allocate(kind_count(nkinds))
  kind_count=0
  do i=1,nkinds
     read(16,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
     write(*,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
  enddo
  close(16)
  write(snap_string,*) gen_cutoff,'1.0000',bi_order,(radii(i),i=1,nkinds),(cutoff(i),i=1,nkinds),&
       'quadraticflag ',quadflag
  write(snap_string2,*) (type2fit(i),i=1,nkinds)
  write(*,*) trim(snap_string)
  write(*,*) 'Fitting types: ',trim(snap_string2)

  if(.not.skip_fit)then

     open(222,file='snapcoeff')

     write(222,*) nkinds,npar
     do i =1,nkinds
        write(222,*) label(i),radii(i),cutoff(i)
        do n=1,npar
           write(222,*) 0.0000000
        enddo
     enddo
     flush(222)

     close(222)

  endif


  allocate(lmp(sys%ndata))

  do i=1,sys%ndata

     call lammps_open_no_mpi ('lmp -log none', lmp(i))
     call lammps_file (lmp(i),sys%inp)
     call lammps_command (lmp(i),'read_data '//trim(sys%data(i)%inp_data))
     call lammps_command (lmp(i),'group fitsnap type '//trim(snap_string2))
     call lammps_file (lmp(i),sys%inp_fix)
     call lammps_command (lmp(i), &
          'compute sna_e all sna/atom '//trim(snap_string)//&
          ' diagonal 3 rmin0 0 switchflag 1')
     call lammps_command (lmp(i),&
          'compute type all property/atom type')
     call lammps_command (lmp(i),&
          'compute id all property/atom id')
     call lammps_command (lmp(i), 'compute pe_ener all pe')

     if(fit_forces)then
        call lammps_command (lmp(i), &
             'compute sna_f all snad/atom '//trim(snap_string)//&
             ' diagonal 3 rmin0 0 switchflag 1')
        call lammps_command (lmp(i), &
             'compute f_x all property/atom fx')
        call lammps_command (lmp(i), &
             'compute f_y all property/atom fy')
        call lammps_command (lmp(i), &
             'compute f_z all property/atom fz')
     endif

  enddo


!!! do kernel

  if(refine)then

     call kernel%dealloc()

     kernel%nkinds=nkinds
     allocate(kernel%K(nkinds))
     do i=1,nkinds
        kernel%K(i)%sigma=sigma(i)
        kernel%K(i)%nenvs=0
     enddo

     do i=1,sys%ndata
        call lammps_command (lmp(i), 'run 0')
        call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
        call lammps_extract_compute (ener, lmp(i), 'sna_e', 1, 2)
        kind_count=0
        do t=1,sys%data(i)%nats
           kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
        enddo
        do t=1,nkinds
           kernel%K(t)%nenvs=kernel%K(t)%nenvs+kind_count(t)*sys%data(i)%frames
        enddo
     enddo
     do t=1,nkinds
        allocate(kernel%K(t)%B(kernel%K(t)%nenvs,size(ener,1)))
     enddo

     do i=1,nkinds
        kernel%K(i)%nenvs=0
     enddo

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

           do t=1,sys%data(i)%nats
              v=kernel%K(nint(kind_nat(map(t))))%nenvs+1
              do k=1,size(ener,1)
                 kernel%K(nint(kind_nat(map(t))))%B(v,k)=ener(k,map(t))
              enddo
              kernel%K(nint(kind_nat(map(t))))%nenvs=&
                   kernel%K(nint(kind_nat(map(t))))%nenvs+1
           enddo

           id_dbl=>null()
           ener=>null()
           Eref=>null()
           kind_nat=>null()

        enddo

        deallocate(id)
        deallocate(map)

     enddo

  endif ! refine

!!! do fitting

  if(.not.skip_fit)then

     v=1
     vv=1

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
           write(222,*) i,l,ener,sys%data(i)%ener(l),&
                (ener-sys%data(i)%ener(l))

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
