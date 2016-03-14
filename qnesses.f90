!###################################################################################################################################
!                                    Tests quantumnesses measures using the two-qubit Werner state
!###################################################################################################################################
subroutine werner_2qb()  ! Tests for quantumnesses quantifiers using a two-qubit Werner's state
! Gnuplot command to see the results
! plot 'werner.dat' u 1:2 t 'Conc', '' u 1:3 t 'EoF', '' u 1:4 t 'Neg', '' u 1:5 t 'logNeg', '' u 1:6 t '2*Disc_hs', 
!      '' u 1:7 t 'Disc_oz_bds', '' u 1:8 t 'max(0,CHSH-2)', '' u 1:9 t 'coherence_l1n', '' u 1:10 t 'coherence_re'
implicit none
complex(8) :: rho(1:4,1:4)  ! The 2-qubit Werner density matrix
complex(8) :: identity(1:4,1:4) ! The 4x4 identiy matrix
complex(8) :: proj(1:4,1:4)  ! For the projector on a Bell state 
real(8) :: w, dw  ! For the weight appearing in the Werner states
real(8) :: concurrence_2qb, EoF_2qb, negativity, log_negativity  ! Entanglement quantifiers
real(8) :: discord_hs_2qb, discord_oz_bds  ! Discord quantifiers
real(8) :: CHSH_2qb  ! Non-locality
real(8) :: coh_l1n, coh_re  ! Coherence quantifiers
complex(8) :: psi_p(4), psi_m(4), phi_p(4), phi_m(4) 
open(unit = 11, file = 'werner.dat', status = 'unknown')

call bell_basis(psi_p, psi_m, phi_p, phi_m)

call identity_c(4, identity) ;   call projector(psi_m, 4, proj)

dw = 1.d0/100.d0 ;   w = 0.d0 - dw  ! Sets the width of the step in w
do
  w = w + dw ;   if( w > 1.d0 ) exit
  rho = (1.d0-w)*(identity/4.d0) + w*proj
  write(11,*) w, concurrence_2qb(rho), EoF_2qb(rho), negativity(2,2,'a',rho), log_negativity(2,2,'a',rho), & 
                 2.d0*discord_hs_2qb('a',rho), discord_oz_bds(rho), & 
                 max(CHSH_2qb(rho)-2.d0,0.d0), coh_l1n(4, rho), coh_re(4, rho)
enddo

end
!###################################################################################################################################
!                                                      Non-Locality
!###################################################################################################################################
real(8) function CHSH_2qb(rho)  ! Returns the CHSH parameter (its maximazation)
! Ref: R. Horodecki, P. Horodecki e M. Horodecki, “Violating Bell inequality by mixed spin-1 states: Necessary and sufficient 
!      condition”, Phys. Lett. A 200, 3402 (1995).
implicit none
complex(8) :: rho(1:4,1:4)  ! The density matrix
real(8) :: ma(1:3)  ! Vector for the polarizations of qubit a
real(8) :: mb(1:3)  ! Vector for the polarizations of qubit b
real(8) :: corr(1:3,1:3)  ! Matrix for the correlations
real(8) :: mK(1:3,1:3)  ! Auxiliary matrix
real(8) :: W(1:3)  ! For the eigenvalues of K

call stokes_parameters_2qb(rho, ma, mb, corr) ;   mK = matmul(transpose(corr),corr) ;   call lapack_dsyevd('N', 3, mK, W)
 CHSH_2qb = 2.d0*sqrt( W(2) + W(3) ) 

end
!###################################################################################################################################