!###################################################################################################################################
!                                                        Discord
!###################################################################################################################################
real(8) function discord_oz_bds(rho)  ! Returns the Ollivier-Zurek discord for 2-qubit Bell-diagonal states
! Ref: S. Luo, Quantum discord for two-qubit systems, PRA 77, 042303 (2008).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: MI_bds, ccorr_hv_bds  ! For the mutual information and classical correlation functions

discord_oz_bds = MI_bds(rho) - ccorr_hv_bds(rho)

end
!---------------------------------------
real(8) function ccorr_hv_bds(rho)  ! Returns the Henderson-Vedral classical correlation for 2-qubit Bell-diagonal states
! Ref: S. Luo, Quantum discord for two-qubit systems, PRA 77, 042303 (2008).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: c11, c22, c33, c  ! Elements of the correlation vector and the maximal one
real(8) :: log2   ! For the log base two function

! Elements of the correlation vector
 c33 = 4.d0*dble(rho(1,1)) - 1.d0 ;   c11 = 2.d0*(dble(rho(1,4)) + dble(rho(2,3))) ;   c22 = 2.d0*(dble(rho(2,3)) - dble(rho(1,4)))
! Classical correlation
 c = max( abs(c11), abs(c22), abs(c33) )
 ccorr_hv_bds = 0.5d0*((1.d0-c)*log2(1.d0-c) + (1.d0+c)*log2(1.d0+c))

end
!---------------------------------------
real(8) function MI_bds(rho)  ! Returns the mutual information for 2-qubit Bell-diagonal states
! Here S(A) = S(B) = 1 and I = 2 - S(A,B)
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: c11, c22, c33  ! Elements of the correlation vector
real(8) :: l00, l01, l10, l11 ! For mutual information and classical correlation, and related auxiliary variables
real(8) :: log2  ! For the log2 function

! Elements of the correlation vector
 c33 = 4.d0*dble(rho(1,1)) - 1.d0 ;   c11 = 2.d0*(dble(rho(1,4)) + dble(rho(2,3))) ;   c22 = 2.d0*(dble(rho(2,3)) - dble(rho(1,4)))
! Eigenvalues of the Bell-diagonal two-qubit density matrix
l00 = ( 1.d0 + c11 - c22 + c33 )/4.d0 ;   l01 = ( 1.d0 + c11 + c22 - c33 )/4.d0
l10 = ( 1.d0 - c11 + c22 + c33 )/4.d0 ;   l11 = ( 1.d0 - c11 - c22 - c33 )/4.d0
! Mutual information (total correlation)
MI_bds = l00*log2(4.d0*l00) + l01*log2(4.d0*l01) + l10*log2(4.d0*l10) + l11*log2(4.d0*l11)

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_tr_xs(ssys, rho)  ! Returns the trace distance discord for 2-qubits X states
! Ref: F. Ciccarello, T. Tufarelli, and V. Giovannetti, Toward computability of trace distance discord, NJP 16, 013038 (2014).
implicit none
complex(8) :: rho(1:4,1:4)  ! The X density matrix
real(8) :: c11, c22, c33, a3, b3  ! For the nonzero correlation matrix and Bloch vector elements
real(8) :: x1, x2, x3, x4  ! Auxiliary variables
character(1) :: ssys  ! Determines which subsystem is the classical one in the extremization

 c11 = 2.d0*(abs(rho(2,3)) + abs(rho(1,4))) ;   c22 = 2.d0*(abs(rho(2,3)) - abs(rho(1,4)))
 c33 = 1.d0 - 2.d0*(dble(rho(2,2)) + dble(rho(3,3)))
 x1 = max(c11**2.d0,c22**2.d0) ;   x2 = min(c11**2.d0,c22**2.d0) ;   x3 = min(c33**2.d0,x1)
 
if ( ssys == 'a' ) then
  a3 = 2.d0*(dble(rho(1,1)) + dble(rho(2,2))) - 1.d0
  x4 = max(c33**2.d0,a3**2.d0+x2)
  discord_tr_xs = sqrt((x4*x1-x3*x2)/(x4-x3+x1-x2))
else if ( ssys == 'b' ) then
  b3 = 2.d0*(dble(rho(1,1)) + dble(rho(3,3))) - 1.d0
  x4 = max(c33**2.d0,b3**2.d0+x2)
  discord_tr_xs = sqrt((x4*x1-x3*x2)/(x4-x3+x1-x2))
endif

end
!---------------------------------------
real(8) function ccorr_tr_xs(rho)  ! Returns the trace distance classical correlation for 2-qubit X states
! Ref: P. C. Obando, F. M. Paula, and M. S. Sarandy, Trace-distance correlations for X states and emergence of the pointer basis 
!      in Markovian and non-Markovian regimes, PRA 92, 032307 (2015).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: c11, c22, c33, a3, b3  ! For the nonzero correlation matrix and Bloch vector elements

 c11 = 2.d0*(abs(rho(2,3)) + abs(rho(1,4))) ;   c22 = 2.d0*(abs(rho(2,3)) - abs(rho(1,4)))
 c33 = 1.d0 - 2.d0*(dble(rho(2,2)) + dble(rho(3,3))) ;   a3 = 2.d0*(dble(rho(1,1)) + dble(rho(2,2))) - 1.d0
 b3 = 2.d0*(dble(rho(1,1)) + dble(rho(3,3))) - 1.d0
 
 ccorr_tr_xs = max(abs(c11),abs(c22),abs(c33-a3*b3))

end
!---------------------------------------
real(8) function tcorr_tr_xs(rho)  ! Returns the total correlation for 2-qubit X states
! Ref: P. C. Obando, F. M. Paula, and M. S. Sarandy, Trace-distance correlations for X states and emergence of the pointer basis 
!      in Markovian and non-Markovian regimes, PRA 92, 032307 (2015).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: c11, c22, c33, a3, b3  ! For the nonzero correlation matrix and Bloch vector elements
real(8) :: kmax, kint, kmin  ! Auxiliary varibles

 c11 = 2.d0*(abs(rho(2,3)) + abs(rho(1,4))) ;   c22 = 2.d0*(abs(rho(2,3)) - abs(rho(1,4)))
 c33 = 1.d0 - 2.d0*(dble(rho(2,2)) + dble(rho(3,3)))
 a3 = 2.d0*(dble(rho(1,1)) + dble(rho(2,2))) - 1.d0 ;   b3 = 2.d0*(dble(rho(1,1)) + dble(rho(3,3))) - 1.d0
 
 !kmax = max(abs(c11),abs(c22),abs(c33-a3*b3)) ;   kmin = min(abs(c11),abs(c22),abs(c33-a3*b3))
 !kint = abs(c11) + abs(c22) + abs(c33-a3*b3) - kmax - kmin
 !tcorr_tr_xs = 0.5d0*(kmax + max(kmax,kmin+kint))

 ! The result seems to be the same for both expressions
 tcorr_tr_xs = 0.25*(abs(c11+c22+c33-a3*b3) + abs(c11+c22-c33+a3*b3) + abs(c11-c22+c33-a3*b3) + abs(c11-c22-c33+a3*b3))
 
end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_hs_2qb(ssys, rho)  ! Returns the 2-norm (or Hilbert-Schmidt) discord for 2-qubits states
! Ref: B. Dakić, V. Vedral e Č. Brukner, Necessary and sufficient condition for nonzero quantum discord, PRL 105, 190502 (2010).
implicit none
character(1) :: ssys  ! Tells if sub-system a or b is classical (in the minimization)
complex(8) :: rho(1:4,1:4)  ! The density matrix
integer :: m, n  ! Auxiliary variables for counters
real(8) :: norm_r, norm_hs  ! For the vector norm and Hilbert-Schmidt norm functions
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: mK(1:3,1:3)  ! Auxiliary matrix
real(8) :: W(1:3)  ! For the eigenvalues of K
real(8) :: trace_he  ! For the trace function

call stokes_parameters_2qb(rho, ma, mb, corr)

if (ssys == 'a' ) then
  do m = 1, 3 ;   do n = 1, 3 ;   mK(m,n) = ma(m)*ma(n) ;   enddo ;   enddo ;   mK = mK + matmul(corr,transpose(corr))
  call lapack_dsyevd('N', 3, mK, W)
  discord_hs_2qb = 0.25d0*( (norm_r(3, ma))**2.d0 + (norm_hs(3, 3, corr))**2.d0 - maxval(W) )
else if ( ssys == 'b' ) then
  do m = 1, 3 ;   do n = 1, 3 ;   mK(m,n) = mb(m)*mb(n) ;   enddo ;   enddo ;   mK = mK + matmul(corr,transpose(corr))
  call lapack_dsyevd('N', 3, mK, W)
  discord_hs_2qb = 0.25d0*( (norm_r(3, mb))**2.d0 + (norm_hs(3, 3, corr))**2.d0 - maxval(W) )
endif

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_hs(ssys, da, db, rho)  ! Returns the expression, from ref. below, for the square 2-norm discord
! Ref: S. Luo and S. Fu, Geometric measure of quantum discord, PRA 82, 034302 (2010)
implicit none
character(1), intent(in) :: ssys  ! Tells if sub-system a or b is classical one, in the minimization
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! The reduced states of the subsystems
real(8), allocatable :: bv(:)  ! For the Bloch vector, of rho_a or of rho_b
real(8), allocatable :: proj_bv(:,:)  ! For the projector on the Bloch vector: x*x^t
real(8), allocatable :: corrmat(:,:)  ! For the correlation matrix
real(8), allocatable :: Xi(:,:)  ! For the S-correlation matrix, with S = A, B
real(8), allocatable :: W(:)  ! For the eigenvalues of e.g. Xi = (2/da*db)*( a*a^t + (2/db)*C*C^t)
integer :: dda, ddb  ! Auxiliary variable for the dimensions

dda = da**2 - 1 ;   ddb = db**2 - 1

if (ssys == 'a' ) then ! CQ states are 'classical'
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
  allocate( bv(1:dda) ) ;   call bloch_vector_gellmann(da, rho_a, bv) ;   deallocate( rho_a )  ! Computes the Bloch vector
  allocate( proj_bv(1:dda,1:dda) ) ;   call projector_re(bv, dda, proj_bv) ;   deallocate( bv )  
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann(da, db, rho, corrmat)  ! Computes the correlation matrix
  allocate( Xi(1:dda,1:dda) ) ;   Xi = (2.d0/dble(da*da*db))*( proj_bv + (2.d0/dble(db))*matmul(corrmat,transpose(corrmat)) ) ! A-correlation matrix
  deallocate( proj_bv, corrmat )
  allocate( W(1:dda) ) ;   call lapack_dsyevd('N', dda, Xi, W) ;   deallocate( Xi )
  discord_hs = sum(W(1:(da*(da-1)))) ;   deallocate( W )  ! for the ameliorated 2-norm discord
else if ( ssys == 'b' ) then  ! QC states are 'classical'
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
  allocate( bv(1:ddb) ) ;   call bloch_vector_gellmann(db, rho_b, bv) ;   deallocate( rho_b )  ! Computes the Bloch vector
  allocate( proj_bv(1:ddb,1:ddb) ) ;   call projector_re(bv, ddb, proj_bv) ;   deallocate( bv ) 
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann(da, db, rho, corrmat)  ! Computes the correlation matrix
  allocate( Xi(1:ddb,1:ddb) ) ;   Xi = (2.d0/dble(da*db*db))*( proj_bv + (2.d0/dble(da))*matmul(transpose(corrmat),corrmat) ) ! B-correlation matrix
  deallocate( proj_bv, corrmat )
  allocate( W(1:ddb) ) ;   call lapack_dsyevd('N', ddb, Xi, W) ;   deallocate( Xi )
  discord_hs = sum(W(1:(db*(db-1)))) ;   deallocate( W )  ! for the ameliorated 2-norm discord
endif

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_hsa(ssys, da, db, rho)  ! Returns the expression, from ref. below, for the amended square 2-norm discord
! Ref: S. J. Akhtarshenas, H. Mohammadi, S. Karimi, and Z. Azmi, Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
character(1), intent(in) :: ssys  ! Tells if sub-system a or b is classical one, in the minimization
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! The reduced states of the subsystems
real(8), allocatable :: bv(:)  ! For the Bloch vector, of rho_a or of rho_b
real(8), allocatable :: proj_bv(:,:)  ! For the projector on the Bloch vector: x*x^t
real(8), allocatable :: corrmat(:,:)  ! For the correlation matrix
real(8), allocatable :: Xi(:,:)  ! For the S-correlation matrix, with S = A, B
real(8), allocatable :: W(:)  ! For the eigenvalues of e.g. Xi = (2/da*db)*( a*a^t + (2/db)*C*C^t)
real(8) :: purity  ! For the purity function
integer :: dda, ddb  ! Auxiliary variable for the dimensions

dda = da**2 - 1 ;   ddb = db**2 - 1

if (ssys == 'a' ) then ! CQ states are 'classical'
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
  allocate( bv(1:dda) ) ;   call bloch_vector_gellmann(da, rho_a, bv) ;   deallocate( rho_a )  ! Computes the Bloch vector
  allocate( proj_bv(1:dda,1:dda) ) ;   call projector_re(bv, dda, proj_bv) ;   deallocate( bv )  
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann(da, db, rho, corrmat)  ! Computes the correlation matrix
  allocate( Xi(1:dda,1:dda) ) ;   Xi = (2.d0/dble(da*da*db))*( proj_bv + (2.d0/dble(db))*matmul(corrmat,transpose(corrmat)) ) ! A-correlation matrix
  deallocate( proj_bv, corrmat )
  allocate( W(1:dda) ) ;   call lapack_dsyevd('N', dda, Xi, W) ;   deallocate( Xi )
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
  discord_hsa = (1.d0/purity(db, rho_b))*sum(W(1:(da*(da-1)))) ;   deallocate( W, rho_b )  ! for the ameliorated 2-norm discord
else if ( ssys == 'b' ) then  ! QC states are 'classical'
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
  allocate( bv(1:ddb) ) ;   call bloch_vector_gellmann(db, rho_b, bv) ;   deallocate( rho_b )  ! Computes the Bloch vector
  allocate( proj_bv(1:ddb,1:ddb) ) ;   call projector_re(bv, ddb, proj_bv) ;   deallocate( bv ) 
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann(da, db, rho, corrmat)  ! Computes the correlation matrix
  allocate( Xi(1:ddb,1:ddb) ) ;   Xi = (2.d0/dble(da*db*db))*( proj_bv + (2.d0/dble(da))*matmul(transpose(corrmat),corrmat) ) ! B-correlation matrix
  deallocate( proj_bv, corrmat )
  allocate( W(1:ddb) ) ;   call lapack_dsyevd('N', ddb, Xi, W) ;   deallocate( Xi )
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
  discord_hsa = (1.d0/purity(da, rho_a))*sum(W(1:(db*(db-1)))) ;   deallocate( W, rho_a )  ! for the ameliorated 2-norm discord
endif

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_mid(da, db, rho)  ! Returns the expression, from ref. below, for the amended square 2-norm discord
! Ref: S. Luo, Using measurement-induced disturbance to characterize correlations as classical or quantum, PRA 77, 022301 (2008)
implicit none
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! The reduced states of the subsystems
complex(8), allocatable :: rhom(:,:)  ! For the measured density matrix
complex(8), allocatable :: proj_a(:,:), proj_b(:,:)  ! Auxiliary matrices for the projectors on the eigenvetors of rho_a and rho_b, respectively
complex(8), allocatable :: kp(:,:)  ! For the Kronecker product between the projectors above
real(8) :: mutual_information, mi_rho, mi_rhom  ! For the mutual information function
real(8) :: neumann  ! For the von Neumann entropy function
complex(8), allocatable :: MA(:,:), MB(:,:)  ! Auxiliary matrices for the eigenvetors of rho_a and rho_b, respectively
real(8), allocatable :: Wa(:), Wb(:)  ! Auxiliary vectors for the eigenvalues of rho_a and rho_b, respectively
integer :: d  ! For the total dimension
integer :: j, k  ! Auxiliary variables for counters

d = da*db

allocate( rho_a(1:da,1:da), rho_b(1:db,1:db) )
! The mutual information of rho is computed here to avoid calculating the reductions of rho two times
call partial_trace_a_he(rho, da, db, rho_b) ;   call partial_trace_b_he(rho, da, db, rho_a)
mi_rho = neumann(da, rho_a) + neumann(db, rho_b) - neumann(d, rho)

allocate( MA(1:da,1:da), Wa(1:da), MB(1:db,1:db), Wb(1:db) )
MA = rho_a ;   call lapack_zheevd('V', da, MA, Wa);   MB = rho_b ;   call lapack_zheevd('V', db, MB, Wb) 
deallocate( rho_a, rho_b, Wa, Wb ) ;   allocate( proj_a(1:da,1:da), proj_b(1:db,1:db), rhom(1:d,1:d), kp(1:d,1:d) )
rhom = 0.d0
do j = 1, da ;   do k = 1, db
  call projector(MA(:,j), da, proj_a) ;   call projector(MB(:,k), db, proj_b)
  call kronecker_product_c(proj_a, da, da, proj_b, db, db, kp)
  rhom = rhom + matmul(kp,matmul(rho,kp))
enddo ;   enddo
deallocate(MA, MB, proj_a, proj_b, kp)
mi_rhom = mutual_information(da, db, rhom) ;   deallocate(rhom)

discord_mid = mi_rho - mi_rhom

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_blocks(ssys, da, db, rho)  ! Returns the expression, from ref. below, for the 'easy' discord
! Ref: H. Cao, Z.-Q. Wu, L.-Y. Hu, X.-X. Xu, and J.-H. Huang, An easy measure of quantum correlation, QIP 14, 4103 (2015).
! REMARK. Not tested yet
implicit none
character(1), intent(in) :: ssys  ! Tells if sub-system a or b is classical one, in the minimization
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: MM(:,:), MM1(:,:), MM2(:,:)  ! For sum of commutators and for the sub-matrices
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm
integer :: j, l, m, n  ! Auxiliary variables for counters
integer :: ull, llm, lln  ! Auxiliary variables for the limits for the counters
real(8) :: disc  ! Auxiliary variable to compute discord

if (ssys == 'a' ) then ! CQ states are the 'classical' ones
  write(*,*) 'a'
else if ( ssys == 'b' ) then  ! QC states are the 'classical' ones
  disc = 0.d0
  allocate( MM(1:db,1:db), MM1(1:db,1:db),  MM2(1:db,1:db) ) ;   MM = 0.d0
  do j = 1, da 
    if (j < da) then ;  ull = da ;  else if (j == da) then ;   ull = da-1 ;   endif
    do l = 1, ull
      MM1 = rho((j-1)*db+1:j*db,(l-1)*db+1:l*db)
      if ( l < da ) then ;   llm = j ;   lln = l+1 ;   else if ( l == da ) then ;   llm = j+1 ;   lln = 1 ;   endif
      do m = llm, da
        do n = lln, da
          MM2 = rho((m-1)*db+1:m*db,(n-1)*db+1:n*db)
          MM = matmul(MM1,MM2) - matmul(MM2,MM1)
          disc = disc + 2.d0*((norm_hs(db, db, MM))**2.d0)
        enddo
        lln = 1
      enddo
    enddo
  enddo
  discord_blocks = sqrt(2.d0-2.d0*sqrt(1.d0-disc))
  deallocate( MM, MM1, MM2 )
endif

end
!###################################################################################################################################