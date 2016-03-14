!###################################################################################################################################
!                                   Random number generators (uniform in [0,1], or integers)
!###################################################################################################################################
!                                      Merseene-Twister RNG (the module is in a separated file)
subroutine sub_rand_seed(seed) ! Uses the computer time to generate a random seed for the initialization of the random number generators
implicit none
character(8) :: date
character(10) :: time
character(5) :: zone
integer, dimension(8) :: values
integer(4) :: seed

call date_and_time(date,time,zone,values) ;   seed = 1.d3*values(8)
  
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine init_mt() ! Initializes the Merseene Twister pseudo-random number generator
use mt95
implicit none
integer(4) :: seed
real(8) :: rn

call sub_rand_seed(seed) ;   call genrand_init(seed) ;   call genrand_real1(rn)
 
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_mt(d, rn) ! Calls the Mersenne Twister pseudo-random number generator
use mt95
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1]

call genrand_real1(rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
!                                                    Gnu's RNG
subroutine init_gnu()  ! Initialization for the GNU RNG RANDOM_NUMBER()
IMPLICIT NONE
INTEGER, ALLOCATABLE :: seed(:)
INTEGER :: i, n, un, istat, dt(8), pid, t(2), s
INTEGER(8) :: count, tms
  
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
! First try if the OS provides a random number generator
OPEN(newunit=un, file="/dev/urandom", access="stream", form="unformatted", action="read", status="old", iostat=istat)
IF (istat == 0) THEN
  read(un) seed
  close(un)
ELSE
  ! Fallback to XOR:ing the current time and pid. The PID is useful in case one launches multiple instances of the same program in parallel.
  CALL SYSTEM_CLOCK(count)
  IF (count /= 0) THEN
    t = TRANSFER(count, t)
  ELSE
    CALL DATE_AND_TIME(values=dt)
    tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 + dt(6) * 60 * 1000 + dt(7) * 1000 + dt(8)
    t = TRANSFER(tms, t)
  ENDIF
  s = ieor(t(1), t(2))
  pid = getpid() + 1099279 ! Add a prime
  s = ieor(s, pid)
  IF (n >= 3) THEN
     seed(1) = t(1) + 36269
     seed(2) = t(2) + 72551
     seed(3) = pid
     IF (n > 3) THEN
        seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
     ENDIF
  ELSE
     seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
  ENDIF
ENDIF
CALL RANDOM_SEED(put=seed)
          
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_gnu(d, rn) ! Calls the Gnu "standard" random number generator. See https://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fNUMBER.html 
! The runtime-library implements George Marsaglia's KISS (Keep It Simple Stupid) random number generator.
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1)

call random_number(rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
!                                                    Netlib's RNG
subroutine init_netlib() ! Initializes Netlib's pseudo-random number generator zufall
implicit none
integer :: seed

call sub_rand_seed(seed) ;   call zufalli(seed)
 
end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_netlib(d, rn) ! Calls Netlib's pseudo-random number generator zufall
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d) ! Vector whose components are random numbers uniformly distributed in [0,1]

call zufall(d,rn)

end
!------------------------------------------------------------------------------------------------------------------------------------
!###################################################################################################################################
!                                           RN with different DISTRIBUTIONS and/or DOMAINS)
!###################################################################################################################################
subroutine rng_gauss(optg, d, grn) ! Returns a vector of random (real) numbers with Gaussian (standard normal) probability distributions
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8), allocatable :: rn(:)  ! Vector of uniformly distributed random numbers
real(8) :: grn(1:d)  ! Vector of Gaussianily distributed random numbers
real(8), parameter :: pi = 4.d0*atan(1.d0)
real(8) :: logterm, angle  ! Auxiliary variables
integer :: j  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate( rn(1:2*d) ) ;   call rng(optg, 2*d, rn)  ! Allocate memory for and get the uniformly distributed random numbers

do j = 1, 2*d, 2   ! Gets the Gaussianily distributed random numbers
  if ( rn(j+1) < 1.d-15 ) rn(j+1) = 1.d-15 ;   logterm = sqrt(-2.d0*log(rn(j+1))) ;   angle = 2.d0*pi*rn(j)
  grn((j+1)/2) = logterm*cos(angle)
enddo

deallocate( rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_exp(optg, d, ern) ! Returns a vector of random (real) numbers with exponential probability distributions
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8), allocatable :: rn(:)  ! Vector of uniformly distributed random numbers
real(8) :: ern(1:d)  ! Vector of exponentianilly distributed random numbers
integer :: j  ! Auxiliary variable for counters
character(10), dimension(5) :: optg  ! Options for the generators

allocate( rn(1:d) ) ;   call rng(optg, d, rn)  ! Allocate memory for and get the uniformly distributed random numbers
 
! Gets the exponentianily distributed random numbers 
do j = 1, d ;   if ( rn(j) < 1.d-15 ) rn(j+1) = 1.d-15 ;   ern(j) = -log(rn(j)) ;   enddo

deallocate( rn )

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng_unif(optg, d, a, b, rn) ! Returns a vector of random (real) numbers with each component uniformly distributed in the interval [a,b]
! Ref: Katzgraber, H. G. (2010). Random numbers in scientific computing: An introduction, arXiv:1005.4117.
implicit none
integer :: d  ! Dimension of the vectors of random numbers
real(8) :: rn(1:d)  ! Vector of random numbers
real(8) :: a, b  ! New lower and upper limits to the domain of the random numbers
character(10), dimension(5) :: optg  ! Options for the generators

call rng(optg, d, rn) ;   rn = a + (b-a)*rn  ! Gets rn in [0,1] and sets the new limits [a,b]

end
!###################################################################################################################################
!                                                  Calling subroutines - RNG
!###################################################################################################################################
subroutine rng_init(optg) ! Calls the initializer subroutine for the choosed random number generator
implicit none
integer :: seed  ! To the RNG initialization
character(10), dimension(5) :: optg  ! Options for the generators

     if ( optg(1) == "std" ) then ;   call init_mt() ! Initializes the Mersenne-Twister rng
else if ( optg(1) == "gnu" ) then ;   call init_gnu() ! Initializes the Gnu's rng
else if ( optg(1) == "netlib" ) then ;   call init_netlib() ! Initializes the Netlib's rng
     end if

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine rng(optg, d, rn) ! Calls the choosed random number generator
implicit none
integer :: d  ! Dimension of the vector of random numbers
real(8) :: rn(1:d)  ! Vector of random numbers
character(10), dimension(5) :: optg  ! Options for the generators

     if (optg(1) == "std") then ;   call rng_mt(d, rn)  ! Calls the Mersenne-Twister rng
else if (optg(1) == "gnu") then ;   call rng_gnu(d, rn)  ! ! Calls the Gnu's rng
else if (optg(1) == "netlib") then ;   call rng_netlib(d, rn)  ! Calls the Netlib's rng
     endif
  
end
!###################################################################################################################################