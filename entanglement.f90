!###################################################################################################################################
!                                                       Entanglement
!###################################################################################################################################
real(8) function concurrence_2qb(rho)  ! Returns the entanglement measure concurrence, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(1:2**2,1:2**2)  ! Density matrix we want to compute the concurrence
complex(8) :: R(1:2**2,1:2**2), rho_tilde(1:2**2,1:2**2), s2_kp_s2(1:2**2,1:2**2)  ! Auxiliary matrices
complex(8) :: egv(1:2**2) ! Eigenvalues of R = rho*rho^tilde
real(8) :: egv_max  ! The greater eigenvalue of R
complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)

call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
call kronecker_product_c(sigma_2, 2, 2, sigma_2, 2, 2, s2_kp_s2) ;   rho_tilde = matmul( matmul(s2_kp_s2,conjg(rho)) , s2_kp_s2 )
R = matmul(rho,rho_tilde) ;   call lapack_zgeev('N', 4, R, egv)

egv_max = max( real(egv(1)), real(egv(2)), real(egv(3)), real(egv(4)) )
 concurrence_2qb = max( 0.d0, (2.d0*sqrt(egv_max)-sqrt(real(egv(1)))-sqrt(real(egv(2)))-sqrt(real(egv(3)))-sqrt(real(egv(4)))))
     
end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function EoF_2qb(rho)  ! Returns the entanglement of formation, for two-qubit states
! Ref: W. K. Wootters, Entanglement of Formation of an Arbitrary State of Two Qubits, Phys.Rev.Lett. 80, 2245 (1998).
implicit none
complex(8) :: rho(1:4,1:4)  ! Density matrix we want to compute the concurrence
real(8) :: concurrence_2qb  ! For the concurrence function
real(8) :: pv(1:2), shannon  ! Probability vector and Shannon's entropy

pv(1) = (1.d0 + sqrt(1.d0 - concurrence_2qb(rho)**2.d0))/2.d0 ;   pv(2) = 1.d0 - pv(1) ;   EoF_2qb = shannon(2, pv)

end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function negativity(da, db, ssys, rho)  ! Returns the entanglement negativity of a bipartite systems
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
implicit none
character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
real(8) :: norm_tr  ! For the trace norm function

if (ssys == 'a') then
  call partial_transpose_a(da, db, rho, rho_pt)
else if (ssys == 'b') then
  call partial_transpose_b(da, db, rho, rho_pt)
endif

negativity = 0.5d0*(norm_tr(da*db, rho_pt) - 1.d0)
     
end
!----------------------------------------------------------------------------------------------------------------------------------
real(8) function log_negativity(da, db, ssys, rho)  ! Returns the entanglement logaritmic negativity of a bipartite systems
! Ref: G. Vidal and R.F. Werner, A computable measure of entanglement, Phys. Rev. A 65, 032314 (2002).
implicit none
character(1) :: ssys  ! Determines in which sub-system the transposition is to be applied (subsys = 'a' or 'b')
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_pt(1:da*db,1:da*db)  ! Bipartite original and partial transposed states
!real(8) :: neg ! Negativity
real(8) :: norm_tr  ! For the trace norm function
real(8) :: log2  ! For the log base two

if (ssys == 'a') then
  call partial_transpose_a(da, db, rho, rho_pt)
else if (ssys == 'b') then
  call partial_transpose_b(da, db, rho, rho_pt)
endif

log_negativity = log2( norm_tr(da*db, rho_pt) )
     
end
!###################################################################################################################################