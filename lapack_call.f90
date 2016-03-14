!###################################################################################################################################
!                                                     LAPACK callers
!###################################################################################################################################
subroutine lapack_zgeev(JOBVR, N, A, Wc)  ! Calls LAPACK's eigensolver for general complex matrices
! ZGEEV computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will 
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray 
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
!character(1) :: JOBVL  !  JOBVL is CHARACTER*1; = 'N': left eigenvectors of A are not computed; = 'V': left eigenvectors of are computed.
character(1) :: JOBVR  !  JOBVR is CHARACTER*1; = 'N': right eigenvectors of A are not computed;  = 'V': right eigenvectors of A are computed.
character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA = N  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(1:N,1:N)  ! A is COMPLEX*16 array, dimension (LDA,N); On entry, the N-by-N matrix A. On exit, A has been overwritten.
complex(8) :: Wc(1:N)  ! Wc is COMPLEX*16 array, dimension (N). Wc contains the computed eigenvalues.
!integer :: LDVL = N  !LDVL is INTEGER. The leading dimension of the array VL.  LDVL >= 1; if JOBVL = 'V', LDVL >= N.
complex(8) :: VL(1:N,1:N)  ! VL is COMPLEX*16 array, dimension (LDVL,N); If JOBVL = 'V', the left eigenvectors u(j) are stored one
                        ! after another in the columns of VL, in the same order as their eigenvalues.
                        ! If JOBVL = 'N', VL is not referenced. u(j) = VL(:,j), the j-th column of VL.
!integer :: LDVR = N  ! LDVR is INTEGER. The leading dimension of the array VR.  LDVR >= 1; if JOBVR = 'V', LDVR >= N.
complex(8) :: VR(1:N,1:N)  ! VR is COMPLEX*16 array, dimension (LDVR,N). If JOBVR = 'V', the right eigenvectors v(j) are stored one
                                     ! after another in the columns of VR, in the same order their eigenvalues.
                                     ! If JOBVR = 'N', VR is not referenced. v(j) = VR(:,j), the j-th column of VR.
                                     
!integer :: LWORK = 2*N  ! LWORK is INTEGER; The dimension of the array WORK.  LWORK >= max(1,2*N). For good performance, LWORK must generally be larger.
                  ! If LWORK = -1, then a workspace query is assumed; the routine only calculates the optimal size of the WORK array, returns
                  ! this value as the first entry of the WORK array, and no error related to LWORK is issued by XERBLA.                                     
complex(8) :: WORK(1:2*N)  ! WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)). On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
real(8) :: RWORK(1:2*N)  ! RWORK is DOUBLE PRECISION array, dimension (2*N)                    
integer :: INFO   ! INFO is INTEGER
                  ! = 0:  successful exit
                  ! < 0:  if INFO = -i, the i-th argument had an illegal value.
                  ! > 0:  if INFO = i, the QR algorithm failed to compute all the
                  ! eigenvalues, and no eigenvectors have been computed; elements and i+1:N of W contain eigenvalues which have converged.
                
 call zgeev ('N',   JOBVR, N, A, N,   Wc, VL, N,    VR, N,    WORK, 2*N,   RWORK, INFO)
!call zgeev (JOBVL, JOBVR, N, A, LDA, Wc, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zheevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for Hermitian complex matrices
! ZHEEVD computes all eigenvalues and, optionally, eigenvectors of a complex Hermitian matrix A.  If eigenvectors are desired, it uses a
! divide and conquer algorithm. The divide and conquer algorithm makes very mild assumptions about floating point arithmetic. It will 
! work on machines with a guard digit in add/subtract, or on those binary machines without guard digits which subtract like the Cray 
! X-MP, Cray Y-MP, Cray C-90, or Cray-2. It could conceivably fail on hexadecimal or decimal machines without guard digits, but we know of none.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
complex(8) :: A(N,N)  ! A is COMPLEX array, dimension (LDA, N). On entry, the Hermitian matrix A.  If UPLO = 'U', the
                      ! leading N-by-N upper triangular part of A contains the upper triangular part of the matrix A.  
                      ! If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower triangular part of the 
                      ! matrix A. On exit, if JOBZ = 'V', then if INFO = 0, A contains the orthonormal eigenvectors of the 
                      ! matrix A. If JOBZ = 'N', then on exit the lower triangle (if UPLO='L') or the upper triangle 
                      ! (if UPLO='U') of A, including the diagonal, is destroyed.
real(8) :: W(N)  ! W is REAL array, dimension (N). If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
          !The length of the array WORK.
          !If N <= 1,                LWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
          !If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.

          !If LWORK = -1, then a workspace query is assumed; the routine
          !only calculates the optimal sizes of the WORK, RWORK and
          !IWORK arrays, returns these values as the first entries of
          !the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
complex(8), allocatable :: WORK(:)  !WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LRWORK   ! LRWORK is INTEGER
          !The dimension of the array RWORK.
          !If N <= 1,                LRWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
          !If JOBZ  = 'V' and N > 1, LRWORK must be at least
           !              1 + 5*N + 2*N**2.

          !If LRWORK = -1, then a workspace query is assumed; the
          !routine only calculates the optimal sizes of the WORK, RWORK
          !and IWORK arrays, returns these values as the first entries
          !of the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
real(8), allocatable :: RWORK(:)    ! RWORK is DOUBLE PRECISION array,
                                        ! dimension (LRWORK)
         ! On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
integer :: LIWORK   ! LIWORK is INTEGER
          !The dimension of the array IWORK.
          !If N <= 1,                LIWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
          !If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.

          !If LIWORK = -1, then a workspace query is assumed; the
          !routine only calculates the optimal sizes of the WORK, RWORK
          !and IWORK arrays, returns these values as the first entries
          !of the WORK, RWORK and IWORK arrays, and no error message
          !related to LWORK or LRWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
          !On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
                  ! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements 
                  ! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the 
                  ! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns 
                  ! INFO/(N+1) through mod(INFO,N+1).
                  
     if (JOBZ == 'N') then ;   LWORK = N + 1 ;   LRWORK = N ;   LIWORK = 1
else if (JOBZ == 'V') then ;   LWORK = 2*N + N**2 ;   LRWORK = 1 + 5*N + 2*N**2 ;   LIWORK = 3 + 5*N
     endif
allocate( WORK(1:LWORK), RWORK(1:LRWORK), IWORK(1:LIWORK) )

call zheevd(JOBZ,'U',N,A,N,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
!call zheevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
  
deallocate( WORK, RWORK, IWORK )
  
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_dsyevd(JOBZ, N, A, W)  ! Calls LAPACK's eigensolver for symmetric real matrices
! DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
! real symmetric matrix A. If eigenvectors are desired, it uses a
! divide and conquer algorithm.!

! The divide and conquer algorithm makes very mild assumptions about
! floating point arithmetic. It will work on machines with a guard
! digit in add/subtract, or on those binary machines without guard
! digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
! Cray-2. It could conceivably fail on hexadecimal or decimal machines
! without guard digits, but we know of none.

! Because of large use of BLAS of level 3, DSYEVD needs N**2 more
! workspace than DSYEVX.
character(1) :: JOBZ  ! JOBZ is CHARACTER*1; = 'N':  Compute eigenvalues only; = 'V':  Compute eigenvalues and eigenvectors.
!character(1) :: UPLO = 'U'  ! UPLO is CHARACTER*1; = 'U':  Upper triangle of A is stored; = 'L':  Lower triangle of A is stored.
integer :: N  ! N is INTEGER; The order of the matrix A.  N >= 0.
!integer :: LDA  ! LDA is INTEGER. The leading dimension of the array A.  LDA >= max(1,N).
real(8) :: A(1:N,1:N)  ! A is DOUBLE PRECISION array, dimension (LDA, N)
          !On entry, the symmetric matrix A.  If UPLO = 'U', the
          !leading N-by-N upper triangular part of A contains the
          !upper triangular part of the matrix A.  If UPLO = 'L',
          !the leading N-by-N lower triangular part of A contains
          !the lower triangular part of the matrix A.
          !On exit, if JOBZ = 'V', then if INFO = 0, A contains the
          !orthonormal eigenvectors of the matrix A.
          !If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
          !or the upper triangle (if UPLO='U') of A, including the
          !diagonal, is destroyed.
real(8) :: W(1:N)  ! W is DOUBLE PRECISION array, dimension (N)
          !If INFO = 0, the eigenvalues in ascending order.
integer :: LWORK  ! LWORK is INTEGER
          !The dimension of the array WORK.
          !If N <= 1,               LWORK must be at least 1.
          !If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
          !If JOBZ = 'V' and N > 1, LWORK must be at least
          !                                      1 + 6*N + 2*N**2.

          !If LWORK = -1, then a workspace query is assumed; the routine
          !only calculates the optimal sizes of the WORK and IWORK
          !arrays, returns these values as the first entries of the WORK
          !and IWORK arrays, and no error message related to LWORK or
          !LIWORK is issued by XERBLA.
real(8), allocatable :: WORK(:)  !WORK is DOUBLE PRECISION array,
           !                              dimension (LWORK)
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: LIWORK   ! LIWORK is INTEGER
          !The dimension of the array IWORK.
          !If N <= 1,                LIWORK must be at least 1.
          !If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
          !If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.!

         ! If LIWORK = -1, then a workspace query is assumed; the
         ! routine only calculates the optimal sizes of the WORK and
         ! IWORK arrays, returns these values as the first entries of
         ! the WORK and IWORK arrays, and no error message related to
         ! LWORK or LIWORK is issued by XERBLA.
integer, allocatable :: IWORK(:)   ! IWORK is INTEGER array, dimension (MAX(1,LIWORK))
          !On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
integer :: INFO   ! INFO is INTEGER; = 0:  successful exit; < 0:  if INFO = -i, the i-th argument had an illegal value
                  ! > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed to converge; i off-diagonal elements 
                  ! of an intermediate tridiagonal form did not converge to zero; if INFO = i and JOBZ = 'V', then the 
                  ! algorithm failed to compute an eigenvalue while working on the submatrix lying in rows and columns 
                  ! INFO/(N+1) through mod(INFO,N+1).
                  
     if (JOBZ == 'N') then ;   LWORK = 2*N+1 ;   LIWORK = 1 
else if (JOBZ == 'V') then ;   LWORK =  1 + 6*N + 2*N**2 ;   LIWORK = 3 + 5*N
     endif
allocate( WORK(1:LWORK), IWORK(1:LIWORK) )

call dsyevd(JOBZ,'U',N,A,N,W,WORK,LWORK,IWORK,LIWORK,INFO)
!call dsyevd(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,IWORK,LIWORK,INFO)

deallocate(WORK, IWORK)
    
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine lapack_zgeqrfp(N, Z, Q, R)  ! Calls LAPACK QR factorization 
! For a complex square matrix A, it returns the orthogonalized Q and upper triangular R matrices
implicit none
!ZGEQRFP computes a QR factorization of a complex M-by-N matrix A:
! A = Q * R. The diagonal entries of R are real and nonnegative.

!integer :: M !is INTEGER
          !The number of rows of the matrix A.  M >= 0.
integer :: N !is INTEGER
          !The number of columns of the matrix A.  N >= 0.
!integer :: LDA !is INTEGER
!          !The leading dimension of the array A.  LDA >= max(1,M).
complex(8) :: A(1:N,1:N) !is COMPLEX*16 array, dimension (LDA,N)
          !On entry, the M-by-N matrix A.
          !On exit, the elements on and above the diagonal of the array
          !contain the min(M,N)-by-N upper trapezoidal matrix R (R is
          !upper triangular if m >= n); the elements below the diagonal,
          !with the array TAU, represent the orthogonal matrix Q as a
          !product of min(m,n) elementary reflectors (see Further
          !Details).
complex(8) :: TAU(1:N) !is COMPLEX*16 array, dimension (min(M,N))
          !The scalar factors of the elementary reflectors
          !Further Details
          !The matrix Q is represented as a product of elementary reflectors
          !  Q = H(1) H(2) . . . H(k), where k = min(m,n).
          !Each H(i) has the form
          !  H(i) = I - tau * v * v**H 
          !where tau is a real scalar, and v is a real vector with
          !v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
          !and tau in TAU(i).
!integer :: LWORK !is INTEGER
          !The dimension of the array WORK. The dimension can be divided into three parts.
          ! 1) The part for the triangular factor T. If the very last T is not bigger 
          !    than any of the rest, then this part is NB x ceiling(K/NB), otherwise, 
          !   NB x (K-NT), where K = min(M,N) and NT is the dimension of the very last T              
          !2) The part for the very last T when T is bigger than any of the rest T. 
          !   The size of this part is NT x NT, where NT = K - ceiling ((K-NX)/NB) x NB,
          !   where K = min(M,N), NX is calculated by
           !        NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
          ! 3) The part for dlarfb is of size max((N-M)*K, (N-M)*NB, K*NB, NB*NB)
          ! So LWORK = part1 + part2 + part3
          ! If LWORK = -1, then a workspace query is assumed; the routine
          ! only calculates the optimal size of the WORK array, returns
          ! this value as the first entry of the WORK array, and no error
          ! message related to LWORK is issued by XERBLA.
complex(8) :: WORK(1:3*N*N) !is COMPLEX*16 array, dimension (MAX(1,LWORK))
          !On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
integer :: INFO !is INTEGER
          != 0:  successful exit
          !< 0:  if INFO = -i, the i-th argument had an illegal value
integer :: j, k  ! Auxiliary variables
complex(8) :: identity(1:N,1:N), proj(1:N,1:N), v(1:N)
complex(8) :: Z(1:N,1:N), Q(1:N,1:N), R(1:N,1:N)  ! Input and output matrices

A = Z
 call zgeqrf(N, N, A, N, TAU, WORK, 3*N*N, INFO) ! It did not find zgeqrfp in the last compilation
!call zgeqrf(M, N, A, LDA, TAU, WORK, LWORK, INFO)

R = 0.d0 ;   forall (j = 1:N, k = 1:N , k >= j ) R(j,k) = A(j,k)

! Computing the unitary
call identity_c(N, identity) ;   Q = identity
do j = 1, N
  if (j > 1) v(1:j-1) = 0.d0 ;   v(j) = 1.d0 ;   v(j+1:N) = A(j+1:N,j) ;   call projector(v, N, proj)
  Q = matmul(Q,(identity-TAU(j)*proj)) 
enddo

end
!###################################################################################################################################