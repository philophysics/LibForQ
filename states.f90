!###################################################################################################################################
!                                                 Some popular Quantum STATES
!###################################################################################################################################
subroutine bell_basis(psi_p, psi_m, phi_p, phi_m)  ! Defines the Bell basis states
implicit none
complex(8) :: psi_p(4), psi_m(4), phi_p(4), phi_m(4) 

psi_p(1) = 0.d0 ;  psi_p(2) = 1.d0/sqrt(2.d0) ;   psi_p(3) = 1.d0/sqrt(2.d0) ;   psi_p(4) = 0.d0
psi_m(1) = 0.d0 ;  psi_m(2) = 1.d0/sqrt(2.d0) ;   psi_m(3) = -1.d0/sqrt(2.d0) ;   psi_m(4) = 0.d0
phi_p(1) = 1.d0/sqrt(2.d0) ;  phi_p(2) = 0.d0 ;   phi_p(3) = 0.d0 ;   phi_p(4) = 1.d0/sqrt(2.d0)
phi_m(1) = 1.d0/sqrt(2.d0) ;  phi_m(2) = 0.d0 ;   phi_m(3) = 0.d0 ;   phi_m(4) = -1.d0/sqrt(2.d0)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_bds(c1, c2, c3, rho)  ! Two-qubit Bell-diagonal state
implicit none
real(8) :: c1, c2, c3  ! Components of the correlation vector
complex(8) :: rho(1:4,1:4)  ! For the Bell-diagobal density matrix
complex(8), allocatable :: kp(:,:)  ! For the Kronecker product between elements of the Pauli group
complex(8), allocatable :: sigma_0(:,:), sigma_1(:,:), sigma_2(:,:), sigma_3(:,:)

allocate( kp(1:4,1:4), sigma_0(1:2,1:2), sigma_1(1:2,1:2), sigma_2(1:2,1:2), sigma_3(1:2,1:2) )
call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)

rho = 0.d0
call kronecker_product_c(sigma_0, 2, 2, sigma_0, 2, 2, kp) ;   rho = rho + kp  ! Identity term
call kronecker_product_c(sigma_1, 2, 2, sigma_1, 2, 2, kp) ;   rho = rho + c1*kp  ! sigma_x term
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, kp) ;   rho = rho + c2*kp  ! sigma_y term
call kronecker_product_c(sigma_3, 2, 2, sigma_3, 2, 2, kp) ;   rho = rho + c3*kp  ! sigma_z term
rho = (1.d0/4.d0)*rho

deallocate( kp, sigma_0, sigma_1, sigma_2, sigma_3 )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_werner(da, x, rho)  ! Returns the Werner's state
! Ref: S. J. Akhtarshenas, H. Mohammadi, S. Karimi, and Z. Azmi, Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
integer :: da ! Dimension of sub-system a (db = da and d = da*db = da^2)
complex(8) :: rho(1:da**2,1:da**2)  ! The density matrix (output)
real(8) :: x   ! For the "mixedness" parameter x (x \in [-1,1])
real(8) :: y  ! Auxiliary variable for computing the matrix elements
integer :: k, l  ! Auxiliary variables for counters

rho = 0.d0
! 'Identity' term
y = (dble(da)-x)/dble(da*(da**2-1))
do k = 1, da**2 ;   rho(k,k) = y ;   enddo
! Off-diagonal terms
y = (dble(da)*x - 1.d0)/dble(da*(da**2-1))
do k = 1, da ;   do l = 1, da
  if ( k == l ) then ;   rho((k-1)*da+l,(l-1)*da+k) = rho((k-1)*da+l,(l-1)*da+k) + y
  else if ( k /= l ) then ;   rho((k-1)*da+l,(l-1)*da+k) = y
  endif
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_isotropic(da, x, rho)  ! Returns the Isotropic state
! Ref: S. J. Akhtarshenas, H. Mohammadi, S. Karimi, and Z. Azmi, Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
integer :: da ! Dimension of sub-system a (db = da and d = da*db = da^2)
complex(8) :: rho(1:da**2,1:da**2)  ! The density matrix (output)
real(8) :: x   ! For the "mixedness" parameter x (x \in [-1,1])
real(8) :: y  ! Auxiliary variable for computing the matrix elements
integer :: k, l  ! Auxiliary variables for counters

rho = 0.d0
! 'Identity' term
y = (1.d0 - x)/dble(da**2-1)
do k = 1, da**2 ;   rho(k,k) = y ;   enddo
! Off-diagonal terms
y = (x*dble(da**2) - 1.d0)/dble(da**2-1)
do k = 1, da ;   do l = 1, da
  if ( k == l ) then ;   rho((k-1)*da+k,(l-1)*da+l) = rho((k-1)*da+k,(l-1)*da+l) + y/dble(da)
  else if ( k /= l ) then ;   rho((k-1)*da+k,(l-1)*da+l) = y/dble(da)
  endif
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_GHZ(nqb, GHZ)  ! Returns the ket for the GHZ state
implicit none
integer :: nqb  ! No. of qubits
integer :: d  ! Dimension
complex(8) :: GHZ(1:2**nqb)  ! Ket for the GHZ state (of a number nqb of qubits)
real(8) :: x  ! Auxiliary variable

x = 1.d0/dsqrt(2.d0) ;   d = 2**nqb ;   GHZ = 0.d0 ;   GHZ(1) = x ;   GHZ(d) = x

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine state_W(nqb, W)  ! Returns the ket for the W state
implicit none
integer :: nqb  ! No. of qubits
integer :: d  ! Dimension of the ket
complex(8) :: W(1:2**nqb)  ! Ket for the W state (of a number nqb of qubits)
real(8) :: x  ! Auxiliary variable
integer :: j  ! Auxiliary variable for counters
 
d = 2**nqb
if ( nqb == 1 ) then ;  x = 1.d0/sqrt(dble(d)) ;   W = x ;   return ;   endif  ! For one qubit
! For two or more qubits
x = 1.d0/sqrt(dble(nqb)) ;   W = 0.d0 ;   W(2) = x ;   W(2**(nqb-1)+1) = x 
if ( nqb > 2 ) then ;   do j = nqb-1, 2, -1 ;   W(2**(nqb-j)+1) = x ;   enddo ;   endif

end
!###################################################################################################################################