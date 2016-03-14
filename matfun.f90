!###################################################################################################################################
!                                                    ARRAY DISPLAY
!###################################################################################################################################
subroutine array_display_rr(nr, nc, A)  ! Displays a real array on the screen
implicit none
integer :: nr, nc  ! No. of rows and columns of the matrix
real(8) :: A(1:nr,1:nc)  ! Complex matrix
real(8) :: Ar(1:nr,1:nc)
integer :: j, k  ! Auxiliary variables for counters
character(10) :: string

string = "(   F10.5)" ;   write(string(2:4),'(I3)') nc
Ar = real(A,4) ;   do j = 1, nr ;   print string, Ar(j,:) ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine array_display_cr(nr, nc, A)  ! Displays the real part of a complex array on the screen
implicit none
integer :: nr, nc  ! No. of rows and columns of the matrix
complex(8) :: A(1:nr,1:nc)  ! Complex matrix
real(8) :: Ar(1:nr,1:nc)
integer :: j, k  ! Auxiliary variables for counters
character(10) :: string

string = "(   F10.5)" ;   write(string(2:4),'(I3)') nc
Ar = real(A,4) ;   do j = 1, nr ;   print string, Ar(j,:) ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine array_display(nr, nc, A)  ! Displays the real and imaginary parts (separately) of a complex array on the screen
implicit none
integer :: nr, nc  ! No. of rows and columns of the matrix
complex(8) :: A(1:nr,1:nc)  ! Complex matrix
real(8) :: Ar(1:nr,1:nc), Ai(1:nr,1:nc)
integer :: j, k  ! Auxiliary variables for counters
character(10) :: string

string = "(   F10.5)" ;   write(string(2:4),'(I3)') nc
write(*,*) 'Real part'
Ar = real(A,4) ;   do j = 1, nr ;   print string, Ar(j,:) ;   enddo
write(*,*) 'Imaginary part'
Ai = aimag(A) ;   do j = 1, nr ;   print string, Ai(j,:) ;   enddo

end
!###################################################################################################################################
!                                                     Auxiliary Matrices
!###################################################################################################################################
subroutine ginibre(optg, m, n, G)  ! Returns a m x n complex matrix from the Ginibre ensemble
implicit none
integer :: m, n  ! No. of rows and columns of G
complex(8) :: G(1:m,1:n)  ! The Ginibre matrix
real(8), allocatable :: grn(:)  ! Vector of gaussianily distributed random numbers
integer :: j, k  ! Auxiliary variables for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate( grn(1:2*m) )
do j = 1, n
  call rng_gauss(optg, 2*m, grn) ;   forall( k = 1:m ) G(k,j) = grn(k) + (0.d0,1.d0)*grn(m+k)  ! Generates G column by column
enddo
deallocate( grn )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine identity_c(d, identity)  ! Returns the complex dxd identity matrix
implicit none
integer :: d  ! Dimension of the identity matrix
complex(8) :: identity(1:d,1:d)  ! Identity matrix
integer :: j, k  ! Auxiliary variable for counters

forall ( j = 1:d, k = 1:d, j /= k ) identity(j,k) = (0.d0,0.d0) ;   forall ( j = 1:d, k = 1:d, j == k ) identity(j,k) = (1.d0,0.d0)
  
end
!###################################################################################################################################
!                                                      Matrix functions
!###################################################################################################################################
subroutine rand_perm(optg, d, rperm)  ! Returns a random permutation of {1,2,...,d-1,d}
! Used e.g. in the normalization and trigonometric methods for rpvg
implicit none
integer :: d  ! Dimension of the random permutation vector
integer :: rperm(1:d)  ! Random permutation vector
integer :: j  ! Counter for do (auxiliary variable)
integer, allocatable :: counter(:)  ! Counter for the no. of times a component of rand_perm is randomly choosed (auxiliary variable)
integer :: ind  ! Identify the component of rand_perm choosed (auxiliary variable)
real(8) :: rn(1)  ! Random numbers
character(10), dimension(5) :: optg  ! Options for the generators

allocate( counter(1:d) )  ! Allocate memory for the counter for the componets of rand_perm
 counter = 0
 
do j = 1, d
  do
    call rng(optg, 1, rn)  ! Returns one random number using the method choosed via option op_rng 
    if ( rn(1) <= 0.d0 ) rn(1) = 1.d-10  ! These two lines are a precaution for the determination of ind (avoid rn = 0 and rn = 1) 
    if ( rn(1) >= 1.d0 ) rn(1) = 1.d0 - 1.d-10
    ind = aint( dble(d)*rn(1) + 1.d0 )  ! aint(x) returns the largest integer smaller than x (i.e., ind is in [1,d])
    counter(ind) = counter(ind) + 1  ! No. of times the value ind appeared in rperm
    if ( counter(ind) >= 2 ) cycle  ! cycle returns the pointer to the top of the do
    if ( counter(ind) == 1 ) then
      rperm(j) = ind ;   exit  ! exit the do going to the next j, i.e., the next component of rand_perm
    endif
  enddo
enddo
 
deallocate( counter )  ! Deallocate the memory used with the counter

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine projector(vec, d, proj)  ! Returns a PROJECTOR on the provided COMPLEX vector 
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: vec(1:d)  ! Vector we want the projector on
complex(8) :: proj(1:d,1:d)  ! Projector on vec
integer :: j, k  ! Auxiliary variables for counters

forall ( j=1:d, k=1:d, j <= k ) proj(j,k) = vec(j)*conjg(vec(k))  ! Elements in the diagonal and above
forall ( j=1:d, k=1:d, j > k ) proj(j,k) = conjg(proj(k,j))  ! Elements below the diagonal

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine projector_re(vec, d, proj)  ! Returns a PROJECTOR on the provided REAL vector 
implicit none
integer :: d  ! Dimension of the vector
real(8) :: vec(1:d)  ! Vector we want the projector on
real(8) :: proj(1:d,1:d)  ! Projector on vec
integer :: j, k  ! Auxiliary variables for counters

forall ( j=1:d, k=1:d, j <= k ) proj(j,k) = vec(j)*vec(k)  ! Elements in the diagonal and above
forall ( j=1:d, k=1:d, j > k ) proj(j,k) = proj(k,j)  ! Elements below the diagonal

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kron_prod_pauli_mat(ord_pm, n, kp_pauli_mat)  ! Returns the Kronecker's product of n Pauli matrices. 
! The especification of the order of the matrices is given as input within the vector vec_ord_pm, whose dimension is n
implicit none
integer, intent(in):: n
integer, intent(in):: ord_pm(1:n)
complex(kind=8), intent(out) :: kp_pauli_mat(1:2**n,1:2**n)
integer, parameter :: d2 = 2
complex(kind=8) :: M2(1:d2,1:d2)
complex(kind=8), allocatable :: M1(:,:), M1_kp_M2(:,:)
integer :: d1, i
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)

do i = 1, n - 1
       if (ord_pm(i+1) == 0) then ;   M2 = sigma_0 ;   else if (ord_pm(i+1) == 1) then ;   M2 = sigma_1
  else if (ord_pm(i+1) == 2) then ;   M2 = sigma_2 ;   else if (ord_pm(i+1) == 3) then ;   M2 = sigma_3
       end if
  d1 = 2**i
  if( i == 1 ) then               ! Used to initiate the sequence of tensor products
    allocate(M1(1:d1,1:d1))
         if (ord_pm(i) == 0) then ;   M1 = sigma_0 ;   else if (ord_pm(i) == 1) then ;   M1 = sigma_1
    else if (ord_pm(i) == 2) then ;   M1 = sigma_2 ;   else if (ord_pm(i) == 3) then ;   M1 = sigma_3
         end if
  end if
  allocate(M1_kp_M2(1:d1*d2,1:d1*d2)) ;   call kronecker_product_C(M1, d1, d1, M2, d2, d2, M1_kp_M2)  
  deallocate(M1) ;   allocate(M1(1:d1*d2,1:d1*d2)) ;   M1 = M1_kp_M2 ;   deallocate(M1_kp_M2)
end do

 kp_pauli_mat = M1 ;   deallocate(M1)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_id_bra(m, n, ket, kp)  ! Returns the tensor product of the mxm identity with the bra (adjoint) of the ket
implicit none
integer :: m  ! Dimension of the identity (given as input)
complex(8) :: id(1:m,1:m)  ! Identify
integer :: n  ! Dimension of the vector
complex(8) :: ket(1:n), bra(1,1:n)  ! Vector and its adjoint (we use a 2d array for bra for practical issues)
complex(8) :: kp(1:m,1:m*n)  ! The Kronecker product
integer :: j  ! Auxiliary variable for counters

forall(j=1:n) bra(1,j) = conjg(ket(j))
kp = 0.d0 ;   forall(j=1:m) kp(j,((j-1)*n+1):j*n) = bra(1,1:n)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_id_ket(m, n, ket, kp)  ! Returns the (optimized) tensor product of the mxm identity with a ket
implicit none
integer :: m  ! Dimension of the identity (given as input)
complex(8) :: id(1:m,1:m)  ! Identify
integer :: n  ! Dimension of the vector
complex(8) :: ket(1:n), keta(1,1:n)  ! Vector and its array version
complex(8) :: kp(1:m*n,1:m)  ! The Kronecker product
integer :: j  ! Auxiliary variable for counters

forall(j=1:n) keta(j,1) = ket(j)
kp = 0.d0 ;   forall(j=1:m) kp(((j-1)*n+1):j*n , j) = keta(1:n,1)

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine kronecker_product_c(M1, nr1, nc1, M2, nr2, nc2, M1_kp_M2)  ! Returns the tensor product of two general complex matrices
implicit none
integer :: nr1, nc1, nr2, nc2  ! Number of rows and columns of the two matrices
complex(8) :: M1(1:nr1,1:nc1), M2(1:nr2,1:nc2)  ! Matrices to take the tensor product of
complex(8) :: M1_kp_M2(1:nr1*nr2,1:nc1*nc2)  ! Matrix containing the tensor product of M1 and M2
integer :: i, j  ! Auxiliary variables for counters

M1_kp_M2 = (0.d0,0.d0) ;   forall ( i = 1:nr1 , j = 1:nc1 ) M1_kp_M2(nr2*(i-1)+1 : nr2*i , nc2*(j-1)+1 : nc2*j)  =  M1(i,j)*M2

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine gram_schmidt(m, n, A, B)  ! Returns, in the coloumns of B, an orthonormal basis obtained from linearly independent vectors 
! given as imput in the columns of A
! Ref: Golub, G. H., and Van Loan, C. F. (2013). Matrix Computations (4th ed.). Baltimore: The Johns Hopkins University Press.
implicit none
integer :: m, n  ! Dimensions of the matrices A and B
complex(8) :: A(1:m,1:n)  ! Matrix with the linearly independent vectors
complex(8) :: B(1:m,1:n)  ! Matrix with the orthonormal basis
real(8) :: norm  ! Vector norm function
complex(8) :: inner_prod  ! Function for the inner product of two complex vectors
integer :: j, k  ! Auxiliary variables for counters

B(:,1) = A(:,1)/norm(m, A(:,1))
do j = 2, n
  B(:,j) = A(:,j)
  do k = 1, j-1
    B(:,j) = B(:,j) - inner_prod(m, B(:,k), A(:,j))*B(:,k)
  enddo
  B(:,j) = B(:,j)/norm(m, B(:,j))
enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine gram_schmidt_modified(m, n, A, B)  ! Returns, in the coloumns of B, an orthonormal basis obtained from linearly 
! independent vectors given as imput in the columns of A
! Ref: Golub, G. H., and Van Loan, C. F. (2013). Matrix Computations (4th ed.). Baltimore: The Johns Hopkins University Press.
implicit none
integer :: m, n  ! Dimensions of the matrices A and B
complex(8) :: A(1:m,1:n)  ! Matrix with the linearly independent vectors
complex(8) :: B(1:m,1:n)  ! Matrix with the orthonormal basis
real(8) :: norm  ! Vector norm function
complex(8) :: inner_prod  ! Function for the inner product of two complex vectors
integer :: j, k  ! Auxiliary variables for counters

do j = 1, n
  B(:,j) = A(:,j)/norm(m, A(:,j))
  if ( j < n ) then ;   do k = j+1, n ;   A(:,k) = A(:,k) - inner_prod(m,B(:,j),A(:,k)) *B(:,j) ;   enddo ;   endif
enddo
! The complexity is (m=n): ~ O(8n**4+16n**3+n**2) ~ O(8n**4) for large n

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine adjoint(m, n, A, Ad)  ! Returns the adjoint (conjugate transpose) of a m by n complex matrix A
implicit none
integer :: m, n  ! Dimensions of the matrix A
complex(8) :: A(1:m,1:n)  ! Matrix we want to compute the adjoint of
complex(8) :: Ad(1:n,1:m)  ! Adjoint of the matrix A
integer :: j, k  ! Auxiliary variables for counters

forall( j = 1:m, k = 1:n ) Ad(k,j) = conjg(A(j,k))

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine outer_product(d, psi, phi, op)  ! Returns the outer product of two vectors
implicit none
integer :: d  ! Dimension of the vectors
complex(8) :: psi(1:d), phi(1:d)  ! The vectors
complex(8) :: op(1:d,1:d)  ! The outer product matrix
integer :: j, k  ! Auxiliary variables for counters

forall (j=1:d, k=1:d) op(j,k) = psi(j)*conjg(phi(k))

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine matmul_AAd(m, n, A, AAd)  ! Returns the product of a matrix A by its adjoint
implicit none
integer :: m, n  ! Dimensions of the matrix
complex(8) :: A(1:m,1:n)  ! The Matrix
complex(8) :: AAd(1:m,1:m)  ! Product of A by its adjoint
integer :: j, k, l  ! Auxiliary variables for counters

AAd = 0.d0
do j = 1, m ;   do k = 1, m
  do l = 1, n ;   AAd(j,k) = AAd(j,k) + A(j,l)*conjg(A(k,l)) ;   enddo
enddo ;   enddo

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine mm(m, n, o, A, B, AB) ! Returns the product of the matrices A and B. But only multiplies elements greater than 1.d-15,
! to avoid error propragation.
implicit none
integer :: m, n, o ! Dimensions of the matrices
complex(8) :: A(1:m,1:n), B(1:n,1:o), AB(1:m,1:o)  ! The Matrices
integer :: j, k, l  ! Auxiliary variables for counters

AB = (0.d0,0.d0)
do j = 1, m ;   do k = 1, o
  do l = 1, n ;   if ( (abs(A(j,l)) > 1.d-15) .and. (abs(B(l,k)) > 1.d-15) )  AB(j,k) = AB(j,k) + A(j,l)*B(l,k) ;   enddo
enddo ;   enddo

end
!###################################################################################################################################
!                                                     Scalar functions
!###################################################################################################################################
real(8) function trace_he(d, A)  ! Returns the trace of a Hermitian matrix A
implicit none
integer :: d  ! Dimension of the matrix
complex(8) :: A(1:d,1:d)  ! Matrix whose trace is computed
integer :: j  ! Auxiliary variable for counters

trace_he = 0.d0 ;   do j = 1, d ;   trace_he = trace_he + dble(A(j,j)) ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function stokes_parameter(rho, ord_pm, n) !  Returns the Stokes' parameter Tr(rho (sigma_j1 \otimes cdots \otimes sigma_jn),
! given a n-qubit density matrix and a sequence of n Pauli matrices
implicit none
integer :: n  ! No. of qubits
integer :: ord_pm(n)  ! Vector with the indexes for the order of the Pauli matrices
complex(8) :: rho(1:2**n,1:2**n)  ! Density matrix
complex(8), allocatable :: kp_pauli_mat(:,:) ! Matrix for the Kronecker product of nqb Pauli matrices 
integer :: j, k  ! Auxiliary variables for counters

allocate ( kp_pauli_mat(1:2**n,1:2**n) )

call kron_prod_pauli_mat(ord_pm, n, kp_pauli_mat)
stokes_parameter = 0.d0
do j = 1, 2**n ;   do k = 1, 2**n ;   stokes_parameter = stokes_parameter + dble(kp_pauli_mat(j,k)*rho(k,j)) ;   end do ;   end do
 
deallocate ( kp_pauli_mat )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine stokes_parameters_2qb(rho, ma, mb, corr)  ! Computes the Stokes parameters for a two-qubit density matrix
implicit none
complex(8) :: rho(1:4,1:4)  ! Density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations (magnetizations) of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations (magnetizations) of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
integer :: ord_pm(1:2)  ! For the order of the Pauli matrix, to be sent to the subroutine which computes the Stokes parameters
real(8) :: stokes_parameter  ! For the Stokes parameter function
integer :: j, k  ! Auxiliary variables for counters

ord_pm(2) = 0 ;   do j = 1, 3;   ord_pm(1) = j ;   ma(j) = stokes_parameter(rho, ord_pm, 2) ;   enddo
ord_pm(1) = 0 ;   do j = 1, 3;   ord_pm(2) = j ;   mb(j) = stokes_parameter(rho, ord_pm, 2) ;   enddo
do j = 1, 3 ;   ord_pm(1) = j ;   do k = 1, 3 ;    ord_pm(2) = k ;   corr(j,k) = stokes_parameter(rho, ord_pm, 2) ;   enddo ; enddo
     
end
!###################################################################################################################################