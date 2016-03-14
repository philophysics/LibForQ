!###################################################################################################################################
!                                                     State vector generators - RSVG
!###################################################################################################################################
subroutine rsv_std(optg, d, rsv)  ! Generates a random state vector (pure quantum state) using the standard method
! Ref: Maziero, J. (2015). Generating pseudo-random discrete probability distributions. Braz. J. Phys. 45, 377.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! Random state vector
real(8), allocatable :: rpv(:)  ! Random probability vector
real(8), allocatable :: rn(:)  ! Vector of random numbers
integer :: j  ! Auxiliary variable for counters
real(8), parameter :: pi = 4.d0*atan(1.d0)
character(10), dimension(5) :: optg  ! Options for the generators

allocate( rpv(1:d), rn(1:d) ) ;  call rpvg(optg, d, rpv) ;   call rng(optg, d, rn)  ! Allocate memory for and get rpv and rn
forall( j = 1 : d ) rsv(j) = sqrt(rpv(j))*exp( (0.d0,1.d0)*(rn(j)*2.d0*pi) )  ! The random phases theta_j = rn*2*pi are in [0,2*pi] 
deallocate( rpv, rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rsv_gauss(optg, d, rsv) ! Generates a random state vector (pure quantum state) using the standard method with gaussianily distributed complex coefficients
! Ref: Maziero, J. (2015). Random sampling of quantum states: A survey of methods. Braz. J. Phys. 45, 575.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! (output) Random state vector
real(8), allocatable :: grn(:)  ! Random vector of gaussianily distributed random numbers
real(8) :: norm  ! Auxiliary variable for normalization
integer :: j  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate( grn(1:2*d) ) ;   call rng_gauss(optg, 2*d, grn) ;   forall ( j = 1:d ) rsv(j) = grn(j) + (0.d0,1.d0)*grn(j+d)
rsv = rsv/norm(d, rsv)  ! Normalize the vector

deallocate( grn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rsv_ru(optg, d, rsv) ! Generates a random quantum state vector using the first column of a random unitary matrix
! Ref: Zyczkowski, K. (1999). Volume of the set of separable states. II, Phys. Rev. A 60, 3496.
implicit none
integer :: d  ! (input) Dimension of the random state vector
complex(8) :: rsv(1:d)  ! (output) Random state vector
complex(8), allocatable :: ru(:,:)  ! Random unitary matrix
integer :: j  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate( ru(1:d,1:d) ) ;  call rug(optg, d, ru) ;   rsv(:) = ru(:,1) ;   deallocate( ru )

end
!###################################################################################################################################
!                                                        Calling subroutines - RSVS
!###################################################################################################################################
subroutine rsvg(optg, d, rsv)  ! Calls the choosed random state vector generator
implicit none
integer :: d  ! Dimension of the vector
complex(8) :: rsv(1:d)  ! The random state vector
character(10), dimension(5) :: optg  ! Options for the generators

     if ( optg(4) == "std" ) then ;   call rsv_std(optg, d, rsv)
else if ( optg(4) == "gauss" ) then ;   call rsv_gauss(optg, d, rsv)
else if ( optg(4) == "ru" ) then ;   call rsv_ru(optg, d, rsv)
     endif

end
!###################################################################################################################################