!###################################################################################################################################
!                                     VECTOR NORMS & INNER PRODUCTS
!###################################################################################################################################
complex(8) function inner_prod(d, psi, phi)  ! Returns the Euclidean inner product of two complex VECTORS psi and phi
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! Complex vectors whose inner product is to be computed
integer :: j  ! Auxiliary variable

inner_prod = (0.d0,0.d0) ;   do j = 1, d ;   inner_prod = inner_prod + conjg(psi(j))*phi(j) ;   enddo
! complexity: ~ O(16d**2), for large d

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function inner_prod_r(d, psi, phi)  ! Returns the Euclidean inner product of two real VECTORS psi and phi
implicit none
integer :: d  ! Dimension of the vectors
real(8) :: psi(1:d), phi(1:d)  ! Real vectors whose inner product is to be computed
integer :: j  ! Auxiliary variable

inner_prod_r = 0.d0 ;   do j = 1, d ;   inner_prod_r = inner_prod_r + psi(j)*phi(j) ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm(d, psi)  ! Returns the Euclidean norm of the complex VECTOR psi
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: psi(1:d)  ! Complex vector whose norm is to be computed
complex(8) :: inner_prod  ! Inner product function to be used

norm = sqrt(dble(inner_prod(d, psi, psi)))
! complexity: ~ O(16d**2), for large d

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_r(d, psi)  ! Returns the Euclidean norm of the real VECTOR psi
implicit none
integer :: d  ! Dimension of the vector
real(8) :: psi(1:d)  ! Real vector whose norm is to be computed
real(8) :: inner_prod_r  ! Inner product function to be used

norm_r = sqrt( inner_prod_r(d, psi, psi) )

end
!###################################################################################################################################
!                                                FIDELITIES
!###################################################################################################################################
real(8) function fidelity_pp(d, psi, phi)  ! Returns the fidelity between two pure states
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! Complex vectors whose fidelity is to be computed
complex(8) :: inner_prod, ip  ! Inner product function

ip = inner_prod(d, psi, phi) ;   fidelity_pp = (dreal(ip))**2.d0 + (dimag(ip))**2.d0

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function fidelity_pm(d, psi, rho)  ! Returns the fidelity between a pure and a mixed state
implicit none
integer :: d  ! Dimension of the vector and matrix
complex(8) :: psi(1:d)  ! Pure state vector
complex(8) :: rho(1:d,1:d)  ! Density matrix
complex(8) :: fid  ! Auxiliary variable
integer :: j, k  ! Auxiliary variables for counters

fid = 0.d0
do j = 1, d ;   do k = 1, d ;   fid = fid + conjg(psi(j))*rho(j,k)*psi(k) ;   enddo ;   enddo
fidelity_pm = dble(fid)

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function fidelity_mm(d, rho, zeta)  ! Returns the fidelity between two mixed states
implicit none
integer :: d  ! Dimension of the density matrices
complex(8) :: rho(1:d,1:d), zeta(1:d,1:d)  ! Density matrices
real(8) :: r(1:d), z(1:d)  ! Vectors for the eigenvalues of rho and zeta, respectively. 
complex(8) :: A(1:d,1:d)  ! Auxiliary matrix. We use r also for the eigenvalues if A=sqrt(rho)*zeta*sqrt(rho)
integer :: j, k, l  ! Auxiliary variables for counters
complex(8) :: op(1:d,1:d)  ! For the outer product
complex(8) ::  inner_prod  ! For the inner product function

call lapack_zheevd('V', d, rho, r) ;   call lapack_zheevd('V', d, zeta, z)

A = 0.d0
do j = 1, d ;   do k = 1, d ;   do l = 1, d
  call outer_product(d, rho(:,j), rho(:,l), op)
  A = A + sqrt(r(j)*r(l))*z(k)*inner_prod(d, rho(:,j), zeta(:,k))*inner_prod(d, zeta(:,k), rho(:,l))*op
enddo ;   enddo ;   enddo

call lapack_zheevd('N', d, A, r) ;   fidelity_mm = 0.d0 ;   do j = 1, d ;   fidelity_mm = fidelity_mm + sqrt(r(j)) ;   enddo
fidelity_mm = (fidelity_mm)**2.d0

end
!###################################################################################################################################
!                                                NORMS
!###################################################################################################################################
real(8) function norm_hs(m, n, A)  ! Returns the Hilbert-Schmidt norm (or 2-norm) of a complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n)  ! Matrix of which want to compute the 2-norm
integer :: j, k  ! Auxiliary variables for counters

norm_hs = 0.d0
do j = 1, m ;   do k = 1, n
  norm_hs = norm_hs + (dble(A(j,k)))**2.d0 + (aimag(A(j,k)))**2.d0
enddo ;   enddo
norm_hs = sqrt(norm_hs)

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_tr(d, A) ! Returns the trace norm (or 1-norm) of an hermitian matrix A
! ||A||_1 = \sum_j |a_j|, where a_j are the eigenvalues of A 
implicit none
integer :: d  ! Dimension of the matrix A
complex(8) :: A(1:d,1:d)  ! Matrix of which want to compute the 1-norm
real(8) :: egv(1:d)  ! Eigenvalues of A
integer :: j  ! Auxiliary variable for counters

call lapack_zheevd('N', d, A, egv) ;   norm_tr = 0.d0 ;   do j = 1, d ;   norm_tr = norm_tr + abs(egv(j)) ; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_tr_ge(d, A) ! Returns the trace norm (or 1-norm) of a general normal matrix A
! ||A||_1 = \sum_j |a_j|, where a_j are the eigenvalues of A 
implicit none
integer :: d  ! Dimension of the matrix A
complex(8) :: A(1:d,1:d)  ! Matrix of which want to compute the 1-norm
complex(8) :: egv(1:d)  ! Eigenvalues of A
integer :: j  ! Auxiliary variable for counters

call lapack_zgeev('N', d, A, egv) ;   norm_tr_ge = 0.d0 ;   do j = 1, d ;   norm_tr_ge = norm_tr_ge + abs(egv(j)) ; enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
real(8) function norm_l1(m, n, A)  ! Returns the l1-norm of a complex mxn complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n) ! The matrix
integer :: j, k

norm_l1 = 0.d0 ;   do j=1,m ;   do k=1,n ;   norm_l1 = norm_l1 + abs(A(j,k)) ;   enddo ;   enddo

end
!###################################################################################################################################
!                                                DISTANCES
!###################################################################################################################################
real(8) function d_hs(d, rho, zeta)  ! Returns the Hilbert-Schmidt (2-norm) distance between two density matrices
implicit none
integer :: d  ! Dimension of the density matrices
complex(8) :: rho(1:d,1:d), zeta(1:d,1:d)  ! Density matrices
real(8) :: norm_hs  ! Function for the Hilbert-Schmidt norm

d_hs = norm_hs(d, d, rho - zeta)

end
!###################################################################################################################################