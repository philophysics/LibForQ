!###################################################################################################################################
!                                                        Coherence
!###################################################################################################################################
real(8) function coh_l1n(d, rho)  ! Returns the l1-norm quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014).
implicit none
integer :: d ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
real(8) :: norm_l1  ! For the l1-norm function

 coh_l1n = norm_l1(d, d, rho)

end
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function coh_re(d, rho)  ! Returns the relative entropy of quantum coherence
! Ref: T. Baumgratz, M. Cramer e M. B. Plenio, Quantifying coherence, PRL 113, 140401 (2014).
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The density matrix
real(8) :: shannon, neumann  ! For the Shannon and von Neumman entropy functions
real(8) :: pv(1:d)  ! Auxiliary probability vector
integer :: j  ! Auxiliary variable for counters

forall(j=1:d) pv(j) = dble(rho(j,j)) ;   coh_re = shannon(d, pv) - neumann(d, rho)

end
!###################################################################################################################################