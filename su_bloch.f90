!###################################################################################################################################
!                                                 Generators of SU(n)
!###################################################################################################################################
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)  ! Defines the three Pauli's matrices and the identity matrix
implicit none
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

sigma_0 = 0.d0 ;   sigma_0(1,1) = 1.d0 ;   sigma_0(2,2) = 1.d0
sigma_1 = 0.d0 ;   sigma_1(1,2) = 1.d0 ;   sigma_1(2,1) = 1.d0
sigma_2 = 0.d0 ;   sigma_2(1,2) = -(0.d0,1.d0) ;   sigma_2(2,1) = (0.d0,1.d0)
sigma_3 = 0.d0 ;   sigma_3(1,1) = 1.d0 ;   sigma_3(2,2) = -1.d0

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine gellmann(d, group, k, l, ggmm)  ! Returns one of the generalized Gell Mann matrices
!Ref. J. Maziero, "Computing partial traces and reduced density matrices", arXiv:1601.07458.
implicit none
integer :: d  ! Dimenson of the space, SU(d)
integer, intent(in) :: group  ! Sets from which of the three groups the Gell Mann matrix should be picked out (it's 1, 2, or 3)
integer, intent(in) :: k, l  ! Sets which matrix is to be picked up (if group == 1, then k == l and k shall be used for j)
complex(8), intent(out) :: ggmm(1:d,1:d)  ! The generalized Gell Mann matrices
real(8) :: mterm  ! Auxiliary variable for the elements of the diagonal generators
integer :: m  ! Auxiliary variable for counters
integer :: j  ! Auxiliary variable for the index of the diagonal generators

ggmm = 0.d0  

if ( group == 1 ) then ! Diagonal generators
  j = k ;   mterm = sqrt(2.d0/(dble(j*(j+1))))  ! G = sqrt(2/(j(j+1)))*(sum_{m=1...j)|m><m| -j|j+1><j+1|)
  do m = 1, j ;   ggmm(m,m) = mterm ;   enddo ;   ggmm(j+1,j+1) = -dble(j)*mterm
else if ( group == 2 ) then  ! Symmetric generators
  ggmm(k,l) = 1.d0 ;   ggmm(l,k) = 1.d0  ! G = |k><l| + |l><k|
else if ( group == 3 ) then  ! Asymmetric generators
  ggmm(k,l) = -(0.d0,1.d0) ;   ggmm(l,k) = (0.d0,1.d0)  ! G = -i(|k><l| - |l><k|)
endif

end
!###################################################################################################################################
!                                     BLOCH VECTOR & CORRELATION MATRIX (optimized)
!###################################################################################################################################
subroutine bloch_vector_gellmann(d, rho, bv)  ! Returns the BLOCH VECTOR using the generalized Gell Mann matrices as basis
! Care must be taken regarding the choice of indexes for the components, as shown below
! Ref. J. Maziero, "Computing partial traces and reduced density matrices", arXiv:1601.07458.
implicit none
integer :: d  ! Dimension of the Bloch vector and of the density matrix
complex(8), intent(in) :: rho(1:d,1:d)  ! Density matrix
real(8), intent(out) :: bv(1:(d**2-1))  ! Bloch vector
integer :: j, k, l  ! Auxiliary variables for counters

bv = 0.d0  ! Initializes the Bloch vector
! First d-1 components of the Bloch vector, corresponding to the digonal generators
do j = 1, d-1
  do k = 1, j ;   bv(j) = bv(j) + dble(rho(k,k)) ;   enddo
  bv(j) = ( bv(j) - dble(j)*dble(rho(j+1,j+1)) )*(dble(d)/sqrt(dble(2*j*(j+1))))
enddo
! Components of the Bloch vector corresponding to the non-diagonal symmetric generators
j = d-1
do k = 1, d-1 ;   do l = k+1, d  ! Here j starts in j = d and ends in j = d*(d+1)/2 - 1
  j = j + 1 ;  bv(j) = dble(d)*dble(rho(l,k))
enddo ;   enddo
! Components of the Bloch vector corresponding to the non-diagonal asymmetric generators
do k = 1, d-1 ;   do l = k+1, d  ! Here j starts in j = d*(d+1)/2 and ends in j = d**2 - 1
  j = j + 1 ;   bv(j) = dble(d)*aimag(rho(l,k))
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine corrmat_gellmann(da, db, rho, corrmat)  ! Returns the CORRELATION MATRIX associated to a density operator, using Gell Mann basis
implicit none
integer :: da, db  ! Dimensions of the sub-systems
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite density operator
real(8) :: corrmat(1:da**2-1,1:db**2-1)  ! The correlation matrix
! Auxiliary matrices (sub-matrices to form the correlation matrix)
real(8), allocatable :: cm11(:,:), cm12(:,:), cm13(:,:), cm21(:,:), cm22(:,:), cm23(:,:), cm31(:,:), cm32(:,:), cm33(:,:)
integer :: dda, ddb, dr  ! Auxiliary variables for the sub-matrices handling
integer :: j, k, m, n, p, q  ! Auxiliary variables for counters

dda = ((da*(da+1))/2) - 1 ;   ddb = ((db*(db+1))/2) - 1

!------------------------------------------------- Diagonal-Diagonal (1,1)
allocate( cm11(1:(da-1),1:(db-1)) ) ;   cm11 = 0.d0
do j = 1, da-1 ;   do k = 1, db-1
  do m = 1, (j+1) ;   do p = 1, (k+1) ;   dr = (m-1)*db + p
    if ( (m /= (j+1)) .and. (p /= (k+1)) ) then
      cm11(j,k) = cm11(j,k) + dble(rho(dr,dr))
    else if ( (m /= (j+1)) .and. (p == (k+1)) ) then
      cm11(j,k) = cm11(j,k) - dble(k)*dble(rho(dr,dr))
    else if ( (m == (j+1)) .and. (p /= (k+1)) ) then
      cm11(j,k) = cm11(j,k) - dble(j)*dble(rho(dr,dr))
    else if ( (m == (j+1)) .and. (p == (k+1)) ) then
      cm11(j,k) = cm11(j,k) + dble(j*k)*dble(rho(dr,dr))
    endif
  enddo ;   enddo
  cm11(j,k) = cm11(j,k)*((0.5d0*dble(da*db))/sqrt(dble(j*(j+1)*k*(k+1))))
enddo ;   enddo
 corrmat(1:da-1,1:db-1) = cm11 ;   deallocate( cm11 )
!------------------------------------------------- Diagonal-Symmetric (1,2)
allocate( cm12(1:(da-1),1:((db*(db-1))/2)) ) ;   cm12 = 0.d0
do j = 1, da-1  ! Loop for the rows
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1  ! Loops for the columns
    do m = 1, j ;   cm12(j,k) = cm12(j,k) + dble(rho((m-1)*db+q,(m-1)*db+p)) ;   enddo
    m = j+1 ;   cm12(j,k) = cm12(j,k) - dble(j)*dble(rho((m-1)*db+q,(m-1)*db+p))
    cm12(j,k) = cm12(j,k)*(dble(da*db)/sqrt(dble(2*j*(j+1))))
  enddo ;   enddo
enddo
 corrmat(1:da-1,db:ddb) = cm12 ;   deallocate( cm12 )
!------------------------------------------------- Diagonal-Asymmetric (1,3)
allocate( cm13(1:(da-1),1:((db*(db-1))/2)) ) ;   cm13 = 0.d0
do j = 1, da-1  ! Loop for the rows
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1  ! Loops for the columns
    do m = 1, j ;   cm13(j,k) = cm13(j,k) + aimag(rho((m-1)*db+q,(m-1)*db+p)) ;   enddo
    m = j+1 ;   cm13(j,k) = cm13(j,k) - dble(j)*aimag(rho((m-1)*db+q,(m-1)*db+p))
    cm13(j,k) = cm13(j,k)*((dble(da*db))/sqrt(dble(2*j*(j+1))))
  enddo ;   enddo
enddo
 corrmat(1:da-1,ddb+1:db**2-1) = cm13 ;   deallocate( cm13 )
!------------------------------------------------- Symmetric-Diagonal (2,1)
allocate( cm21(1:((da*(da-1))/2),1:(db-1)) ) ;   cm21 = 0.d0
do k = 1, db-1  ! Loop for the columns
  j = 0
  do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1  ! Loops for the rows
    do p = 1, k ;   cm21(j,k) = cm21(j,k) + dble(rho((n-1)*db+p,(m-1)*db+p)) ;   enddo
    p = k+1 ;   cm21(j,k) = cm21(j,k) - dble(k)*dble(rho((n-1)*db+p,(m-1)*db+p))
    cm21(j,k) = cm21(j,k)*((dble(da*db))/sqrt(dble(2*k*(k+1))))
  enddo ;   enddo
enddo
 corrmat(da:dda,1:db-1) = cm21 ;   deallocate( cm21 )
!------------------------------------------------- Symmetric-Symmetric (2,2)
allocate( cm22(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm22 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    cm22(j,k) = cm22(j,k) + dble(rho((n-1)*db+q,(m-1)*db+p)) + dble(rho((n-1)*db+p,(m-1)*db+q))
    cm22(j,k) = cm22(j,k)*(0.5d0*dble(da*db))
  enddo ;   enddo
enddo ;   enddo
 corrmat(da:dda,db:ddb) = cm22 ;   deallocate( cm22 )
!------------------------------------------------- Symmetric-Asymmetric (2,3)
allocate( cm23(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm23 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    cm23(j,k) = cm23(j,k) + aimag(rho((n-1)*db+q,(m-1)*db+p)) - aimag(rho((n-1)*db+p,(m-1)*db+q))
    cm23(j,k) = cm23(j,k)*(0.5d0*dble(da*db))
  enddo ;   enddo
enddo ;   enddo
 corrmat(da:dda,ddb+1:db**2-1) = cm23 ;   deallocate( cm23 )
!------------------------------------------------- Asymmetric-Diagonal (3,1)
allocate( cm31(1:((da*(da-1))/2),1:(db-1)) ) ;   cm31 = 0.d0
do k = 1, db-1  ! Loop for the columns
  j = 0
  do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1  ! Loops for the rows
    do p = 1, k ;   cm31(j,k) = cm31(j,k) + aimag(rho((n-1)*db+p,(m-1)*db+p)) ;   enddo
    p = k+1 ;   cm31(j,k) = cm31(j,k) - dble(k)*aimag(rho((n-1)*db+p,(m-1)*db+p))
    cm31(j,k) = cm31(j,k)*((dble(da*db))/sqrt(dble(2*k*(k+1))))
  enddo ;   enddo
enddo
 corrmat(dda+1:da**2-1,1:db-1) = cm31 ;   deallocate( cm31 )
!------------------------------------------------- Asymmetric-Symmetric (3,2)
allocate( cm32(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm32 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    cm32(j,k) = cm32(j,k) + aimag(rho((n-1)*db+q,(m-1)*db+p)) + aimag(rho((n-1)*db+p,(m-1)*db+q))
    cm32(j,k) = cm32(j,k)*(0.5d0*dble(da*db))
  enddo ;   enddo
enddo ;   enddo
 corrmat(dda+1:da**2-1,db:ddb) = cm32 ;   deallocate( cm32 )
!------------------------------------------------- Asymmetric-Asymmetric (3,3)
allocate( cm33(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm33 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    cm33(j,k) = cm33(j,k) + dble(rho((n-1)*db+p,(m-1)*db+q)) - dble(rho((n-1)*db+q,(m-1)*db+p))
    cm33(j,k) = cm33(j,k)*(0.5d0*dble(da*db))
  enddo ;   enddo
enddo ;   enddo
 corrmat(dda+1:da**2-1,ddb+1:db**2-1) = cm33 ;   deallocate( cm33 )

end
!###################################################################################################################################
!                                     BLOCH VECTOR & CORRELATION MATRIX (un-optimized)
!###################################################################################################################################
subroutine bloch_vector_gellmann_unopt(d, rho, bv)  ! Returns the BLOCH VECTOR using the generalized Gell Mann matrices as basis.
! UN-OPTIMIZED calculation
! Care must be taken regarding the choice of indexes for the components, as shown below
! Ref. J. Maziero, "Computing partial traces and reduced density matrices", arXiv:1601.07458.
implicit none
integer :: d  ! Dimension of the Bloch vector and of the density matrix
complex(8), intent(in) :: rho(1:d,1:d)  ! Density matrix
real(8), intent(out) :: bv(1:(d**2-1))  ! Bloch vector
integer :: j, k, l  ! Auxiliary variables for counters
complex(8), allocatable :: ggmm(:,:)  ! For the generalized Gell Mann matrices
real(8) :: trace_he  ! For the trace function (for Hermitian matices)

allocate( ggmm(1:d,1:d) )
 
! First d-1 components of the Bloch vector, corresponding to the digonal generators
do j = 1, d-1 ;   call gellmann(d, 1, j, j, ggmm) ;   bv(j) = (dble(d)/2.d0)*trace_he(d, matmul(rho,ggmm)) ;   enddo
! Components of the Bloch vector corresponding to the non-diagonal symmetric generators
j = d-1
do k = 1, d-1 ;   do l = k+1, d  ! 2nd group of generators
  j = j + 1 ;   call gellmann(d, 2, k, l, ggmm) ;   bv(j) = (dble(d)/2.d0)*trace_he(d, matmul(rho,ggmm))
enddo ;   enddo
! Components of the Bloch vector corresponding to the non-diagonal asymmetric generators
do k = 1, d-1 ;   do l = k+1, d  ! 3rd group of generators
  j = j + 1 ;   call gellmann(d, 3, k, l, ggmm) ;   bv(j) = (dble(d)/2.d0)*trace_he(d, matmul(rho,ggmm))
enddo ;   enddo

deallocate( ggmm )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine corrmat_gellmann_unopt(da, db, rho, corrmat)  ! Returns the CORRELATION MATRIX associated to a density operator, using Gell Mann basis
implicit none
integer :: da, db  ! Dimensions of the sub-systems
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite density operator
real(8) :: corrmat(1:da**2-1,1:db**2-1)  ! The correlation matrix
! Auxiliary matrices (sub-matrices to form the correlation matrix)
real(8), allocatable :: cm11(:,:), cm12(:,:), cm13(:,:), cm21(:,:), cm22(:,:), cm23(:,:), cm31(:,:), cm32(:,:), cm33(:,:)
integer :: dda, ddb, d  ! Auxiliary variables for the sub-matrices handling
integer :: j, k, m, n, p, q  ! Auxiliary variables for counters
complex(8), allocatable :: ggmma(:,:), ggmmb(:,:), kp(:,:)
real(8) :: trace_he  ! For the trace function (for Hermitian matices)

dda = ((da*(da+1))/2) - 1 ;   ddb = ((db*(db+1))/2) - 1 ;  d = da*db ;   allocate( ggmma(1:da,1:da), ggmmb(1:db,1:db), kp(1:d,1:d) )

!------------------------------------------------- Diagonal-Diagonal (1,1)
allocate( cm11(1:da-1,1:db-1) ) ;   cm11 = 0.d0
do j = 1, da-1 ;   do k = 1, db-1
  call gellmann(da, 1, j, j, ggmma) ; call gellmann(db, 1, k, k, ggmmb) ; call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
  cm11(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
enddo ;   enddo
 corrmat(1:da-1,1:db-1) = cm11 ;   deallocate( cm11 )
!------------------------------------------------- Diagonal-Symmetric (1,2)
allocate( cm12(1:(da-1),1:((db*(db-1))/2)) ) ;   cm12 = 0.d0
do j = 1, da-1  ! Loop for the rows
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1  ! Loops for the columns
    call gellmann(da, 1, j, j, ggmma); call gellmann(db, 2, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm12(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo
 corrmat(1:da-1,db:ddb) = cm12 ;   deallocate( cm12 )
!------------------------------------------------- Diagonal-Asymmetric (1,3)
allocate( cm13(1:(da-1),1:((db*(db-1))/2)) ) ;   cm13 = 0.d0
do j = 1, da-1  ! Loop for the rows
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1  ! Loops for the columns
    call gellmann(da, 1, j, j, ggmma); call gellmann(db, 3, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm13(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo
 corrmat(1:da-1,ddb+1:db**2-1) = cm13 ;   deallocate( cm13 )
!------------------------------------------------- Symmetric-Diagonal (2,1)
allocate( cm21(1:((da*(da-1))/2),1:(db-1)) ) ;   cm21 = 0.d0
do k = 1, db-1  ! Loop for the columns
  j = 0
  do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1  ! Loops for the rows
    call gellmann(da, 2, m, n, ggmma); call gellmann(db, 1, k, k, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm21(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo
 corrmat(da:dda,1:db-1) = cm21 ;   deallocate( cm21 )
!------------------------------------------------- Symmetric-Symmetric (2,2)
allocate( cm22(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm22 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    call gellmann(da, 2, m, n, ggmma); call gellmann(db, 2, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm22(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo ;   enddo
 corrmat(da:dda,db:ddb) = cm22 ;   deallocate( cm22 )
!------------------------------------------------- Symmetric-Asymmetric (2,3)
allocate( cm23(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm23 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    call gellmann(da, 2, m, n, ggmma); call gellmann(db, 3, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm23(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo ;   enddo
 corrmat(da:dda,ddb+1:db**2-1) = cm23 ;   deallocate( cm23 )
!------------------------------------------------- Asymmetric-Diagonal (3,1)
allocate( cm31(1:((da*(da-1))/2),1:(db-1)) ) ;   cm31 = 0.d0
do k = 1, db-1  ! Loop for the columns
  j = 0
  do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1  ! Loops for the rows
    call gellmann(da, 3, m, n, ggmma); call gellmann(db, 1, k, k, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm31(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo
 corrmat(dda+1:da**2-1,1:db-1) = cm31 ;   deallocate( cm31 )
!------------------------------------------------- Asymmetric-Symmetric (3,2)
allocate( cm32(1:((da*(da-1))/2),1:((db*(db-1))/2)) ) ;   cm32 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    call gellmann(da, 3, m, n, ggmma); call gellmann(db, 2, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm32(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo ;   enddo
 corrmat(dda+1:da**2-1,db:ddb) = cm32 ;   deallocate( cm32 )
!------------------------------------------------- Asymmetric-Asymmetric (3,3)
allocate( cm33(1:(da*(da-1))/2,1:(db*(db-1))/2) ) ;   cm33 = 0.d0
j = 0
do m = 1, da-1 ;   do n = m+1, da ;   j = j + 1
  k = 0
  do p = 1, db-1 ;   do q = p+1, db ;   k = k + 1
    call gellmann(da, 3, m, n, ggmma); call gellmann(db, 3, p, q, ggmmb); call kronecker_product_c(ggmma, da, da, ggmmb, db, db, kp)
    cm33(j,k) = (dble(d)/4.d0)*trace_he(d, matmul(kp,rho))
  enddo ;   enddo
enddo ;   enddo
 corrmat(dda+1:da**2-1,ddb+1:db**2-1) = cm33 ;   deallocate( cm33 )

end
!###################################################################################################################################