!###################################################################################################################################
!                                             Random density matrix generators - RDMG
!###################################################################################################################################
subroutine rdm_std(optg, d, rdm) ! Generates a random density matrix using the standard method: rdm = \sum_j p_j U|c_j><c_j|U^†
! Ref: Maziero, J. (2015). Random sampling of quantum states: A survey of methods. Braz. J. Phys. 45, 575.
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: ru(:,:)  ! Random unitary matrix
real(8), allocatable :: rpv(:)  ! Random probability vector
integer :: j, k, l  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate ( ru(1:d,1:d), rpv(1:d) ) ;   call rpvg(optg, d, rpv)
call rug(optg, d, ru)  ! Allocate memory for and get these random variables

rdm = (0.d0,0.d0)  ! Generates the rdm
do j = 1, d ;   do k = 1, d  ;  do l = 1, d ;   rdm(j,k) = rdm(j,k) + rpv(l)*ru(j,l)*conjg(ru(k,l)) ;   enddo ;   enddo ;   enddo

deallocate ( ru, rpv )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_ginibre(optg, d, rdm) ! Generates a random density matrix normalizing G*G^†, with G being a Ginibre matrix
! Ref:  \.{Z}yczkowski, K., and Sommers, H.-J. (2001). Induced measures in the space of mixed quantum states. 
!       J. Phys. A: Math. Gen. 34, 7111.
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: G(:,:), GGd(:,:)  ! For the Ginibre matrix and its product with its adjoint
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm function
character(10), dimension(5) :: optg  ! Options for the generators

allocate ( G(1:d,1:d), GGd(1:d,1:d) ) ;   call ginibre(optg, d, d, G) ;   call  matmul_AAd(d, d, G, GGd) 
rdm = GGd/((norm_hs(d, d, G))**2.d0)  ! Defines the density matrix
deallocate ( G, GGd )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_bures(optg, d, rdm) ! Generates a random density matrix normalizing (id+U)G*G^†(id+U^†), with G being a Ginibre matrix and U a random unitary
! Ref: Al Osipov, V., Sommers, H.-J., and \.{Z}yczkowski, K. (2010). Random Bures mixed states and the distribution of their purity. 
!      J. Phys. A: Math. Theor. 43, 055302.
implicit none
integer :: d  ! Dimension of the matrices
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: G(:,:)  ! For the Ginibre matrix
complex(8), allocatable :: U(:,:)  ! For the random unitary matrix
!complex(8), allocatable :: id(:,:)  ! For the indentity matrix
complex(8), allocatable :: A(:,:), AAd(:,:)  ! Auxiliary matrices
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm function
integer :: j  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate ( G(1:d,1:d), U(1:d,1:d), A(1:d,1:d), AAd(1:d,1:d) )
call ginibre(optg, d, d, G) ;   call rug(optg, d, U) ;   forall ( j = 1:d ) U(j,j) = U(j,j) + 1.d0  ! U -> id+U
A = matmul(U,G) ;   call  matmul_AAd(d, d, A, AAd)
rdm = AAd/((norm_hs(d, d, A))**2.d0)  ! Defines the density matrix
deallocate ( G, U, A, AAd )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_ptrace(optg, d, rdm) ! Generates a random density matrix via partial tracing over a random state vector
! Ref: Mej\'ia, J., Zapata, C., and Botero, A. (2015). The difference between two random mixed quantum states: Exact and asymptotic 
!      spectral analysis. arXiv:1511.07278.
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rdm(1:d,1:d)  ! Random density matrix 
complex(8), allocatable :: rsv(:)  ! For the random state vector
complex(8), allocatable :: proj(:,:)  ! For the projector
character(10), dimension(5) :: optg  ! Options for the generators

allocate ( rsv(1:d*d), proj(1:d*d,1:d*d) )
call rsvg(optg, d*d, rsv) ;   call projector(rsv, d*d, proj)
call partial_trace_a_he(proj, d, d, rdm)
deallocate ( rsv, proj )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rdm_bds(rbds) ! Generates a random Bell-diagonal state
implicit none
complex(8) :: rbds(4,4)  ! Random Bell-diagonal density matrix 
complex(8), allocatable :: psi_p(:), psi_m(:), phi_p(:), phi_m(:)  ! For the Bell states
real(8), allocatable :: rpv(:)  ! For the random probability vector
complex(8), allocatable :: proj(:,:)  ! For the projectors
character(10), dimension(5) :: optg  ! Options for the generators

allocate( psi_p(4), psi_m(4), phi_p(4), phi_m(4), rpv(4), proj(4,4) ) 
call bell_basis(psi_p, psi_m, phi_p, phi_m)
optg = 'std' ;   call rng_init(optg) ;   call rpvg(optg, 4, rpv)
rbds = 0.d0
call projector(psi_p, 4, proj) ;   rbds = rbds + rpv(1)*proj
call projector(psi_m, 4, proj) ;   rbds = rbds + rpv(2)*proj
call projector(phi_p, 4, proj) ;   rbds = rbds + rpv(3)*proj
call projector(phi_m, 4, proj) ;   rbds = rbds + rpv(4)*proj
deallocate ( psi_p, psi_m, phi_p, phi_m, rpv, proj )

end
!###################################################################################################################################
!                                                   Calling subroutines - RDMG
!###################################################################################################################################
subroutine rdmg(optg, d, rdm)  ! Calls the choosed random density matrix generator
implicit none
integer :: d  ! Dimension of the random density matrix
complex(8) :: rdm(1:d,1:d)  ! The random density matrix
character(10), dimension(5) :: optg  ! Options for the generators

     if ( optg(5) == "std" ) then ;   call rdm_std(optg, d, rdm)
else if ( optg(5) == "ginibre" ) then ;   call rdm_ginibre(optg, d, rdm)
else if ( optg(5) == "bures" ) then ;   call rdm_bures(optg, d, rdm)
else if ( optg(5) == "ptrace" ) then ;   call rdm_ptrace(optg, d, rdm)
     endif

end
!###################################################################################################################################
!                                                   Basic tests for density matrices
!###################################################################################################################################
subroutine rho_tests(d, rho)
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The random density matrix
complex(8), allocatable :: rhoa(:,:)  ! For the adjoint of rho
real(8), allocatable :: W(:)  ! For the eigenvalues of rho
real(8) :: trace_he  ! For the trace function
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm

allocate( rhoa(1:d,1:d), W(1:d) )

write(*,*) 'Verifying if the trace is equal to one'
write(*,*) 'Tr(rho) = ', trace_he(d, rho)  

write(*,*) 'Verifying Hermiticity'
call adjoint(d, d, rho, rhoa) ;   write(*,*) '||rho-rhoa||_2 = ', norm_hs(d, d, rho-rhoa) 

write(*,*) 'Verifying positivity'
call lapack_zheevd('N', d, rho, W) ;   write(*,*) 'Eigvals = ', W

deallocate( rhoa, W )

end
!###################################################################################################################################