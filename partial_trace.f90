!###################################################################################################################################
! The next four subroutines are used, and needed, for computing the partial trace for general multi-partite systems
! (for the GENERAL case)
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace(rho, d, di, nss, ssys, dr, rhor)  ! Returns the partial trace for general multipartite systems
implicit none
integer :: nss  ! Number of sub-systems
integer :: di(1:nss)  ! Vector specifying the dimensions of the sub-systems
integer :: ssys(1:nss)  ! Vector specifying the sub-systems to be traced out
integer :: d ! Total dimension (is the product of the sub-systems dimensions)
complex(8) :: rho(1:d,1:d)  ! Total matrix (given as input)
integer :: dr  ! Dimension of the reduced matrix (is the product of the dimensions of the sub-systems we shall not trace out)
complex(8) :: rhor(1:dr,1:dr)  ! Reduced matrix
complex(8), allocatable :: mat1(:,:), mat2(:,:), rhored(:,:)  ! Auxiliary matrices
integer :: j, k, l  ! Auxiliary variables for counters
integer :: da, db, dc  ! Auxiliary variables for the dimensions

! For bipartite systems
if ( nss == 2 ) then
  if ( ssys(1) == 0 ) then ; call partial_trace_a(rho, di(1), di(2), rhor)
  else if ( ssys(2) == 0 ) then ; call partial_trace_b(rho, di(1), di(2), rhor) 
  endif
  rhored = rhor
! For multipartite systems
else if ( nss >= 3 ) then
  ! Left partial traces
    l = 0 ;   do ;   if ( ssys(l+1) == 1 ) exit ;   l = l + 1 ;   enddo  ! l defines up to which position we shall trace out
    if ( l == 0 ) then
      mat1 = rho  ! This matrix shall be used below if l = 0
    else if ( l > 0 ) then
      if ( l == 1 ) then ;   da = di(1) ;   else ;   da = product(di(1:l)) ;  endif ;   db = d/da
      allocate( mat1(1:db,1:db) ) ;   call partial_trace_a(rho, da, db, mat1)
      d = db  ! After this operation, the matrix left over is mat1, whose dimension is d = db  
      rhored = mat1    
    endif
  ! Right partial traces
    k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo  ! k defines up to which position we shall trace out
    if ( k == (nss+1) ) then
      mat2 = mat1  ! This matrix shall be used below if k = nss+1
      rhored = mat2
    else if ( k < (nss+1) ) then
      if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;  endif ;   da = d/db
      allocate( mat2(1:da,1:da) ) ;   call partial_trace_b(mat1, da, db, mat2) ;   deallocate( mat1 )
      d = da  ! After this operation, the matrix left over is mat2, whose dimension is d = da
      rhored = mat2
    endif
    if ( l == 0 .and. k == (nss+1) ) deallocate( mat1 )  ! To avoid allocating mat1 two times in these cases
  ! Inner partial traces
    if ( (k-l) > 3 ) then  ! If (k-l) <= 3 there is no need to take inner partial traces
    do j = (l+2), (k-2)
      if ( ssys(j) == 0 ) then
        db = di(j) ;   if ( j == (k-2) ) then ;   dc = di(k-1) ;   else ;   dc = product(di(j+1:k-1)) ;  endif ;   da = d/(db*dc)
        allocate( mat1(1:da*dc,1:da*dc) ) ;   call partial_trace_3(mat2, da, db, dc, mat1) ;   rhored = mat1
        d = da*dc ;   deallocate( mat2 ) ;   allocate( mat2(1:d,1:d) ) ;   mat2 = mat1 ;   deallocate( mat1 )
      endif
    enddo
    endif
endif

rhor = rhored

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_a(rho, da, db, rho_b)  ! Returns the left partial trace for a bi-partite matrix
implicit none
integer :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_b(1:db,1:db)  !  Reduced matrix
integer :: j, k, l  ! Auxiliary variable for counters

rho_b = (0.d0,0.d0)
do j = 1, db ; do k = 1, db ; do l = 1, da ;   rho_b(j,k) = rho_b(j,k) + rho((j-1)*da+l,(k-1)*da+l) ;   enddo ; enddo ; enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_b(rho, da, db, rho_a)  ! Returns the right partial trace for a bi-partite matrix
implicit none
integer :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_a(1:da,1:da)  !  Reduced matrix
integer :: j, k, l  ! Auxiliary variables for counters

rho_a = (0.d0,0.d0)
do j = 1, da ; do k = 1, da ; do l = 1, db ;   rho_a(j,k) = rho_a(j,k) + rho((j-1)*db+l,(k-1)*db+l) ;   enddo ; enddo ; enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_3(rho, da, db, dc, rho_ac)  ! Returns the inner partial trace for a three-partite matrix
implicit none
integer :: da, db, dc  ! Dimension of the three sub-systems (the dimension of the whole system is d = da*db*dc)
complex(8) :: rho(1:da*db*dc,1:da*db*dc)  ! Three-partite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_ac(1:da*dc,1:da*dc)  ! Bipartite reduced matrix
integer :: j, k, l, m, o  ! Auxiliary variables for counters
integer :: cj, ck, ccj, cck  ! Auxiliary variables for the dimensions

rho_ac = (0.d0,0.d0)
do j = 1, da ;   do l = 1, dc ;   cj = (j-1)*dc+l ;   ccj = (j-1)*db*dc+l
do m = 1, da ;   do o = 1, dc ;   ck = (m-1)*dc+o ;   cck = (m-1)*db*dc+o
  do k = 1, db ;   rho_ac(cj,ck) = rho_ac(cj,ck) + rho(ccj+(k-1)*dc,cck+(k-1)*dc) ;   enddo
enddo ;   enddo ;   enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
!###################################################################################################################################
! The next four subroutines are used, and needed, for computing the partial trace for general multi-partite systems 
! (for the HERMITIAN case)
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_he(rho, d, di, nss, ssys, dr, rhor)  ! Returns the partial trace for general multipartite systems
implicit none
integer :: nss  ! Number of sub-systems
integer :: di(1:nss)  ! Vector specifying the dimensions of the sub-systems
integer :: ssys(1:nss)  ! Vector specifying the sub-systems to be traced out
integer :: d ! Total dimension (is the product of the sub-systems dimensions)
complex(8) :: rho(1:d,1:d)  ! Total matrix (given as input)
integer :: dr  ! Dimension of the reduced matrix (is the product of the dimensions of the sub-systems we shall not trace out)
complex(8) :: rhor(1:dr,1:dr)  ! Reduced matrix
complex(8), allocatable :: mat1(:,:), mat2(:,:), rhored(:,:)  ! Auxiliary matrices
integer :: j, k, l  ! Auxiliary variables for counters
integer :: da, db, dc  ! Auxiliary variables for the dimensions

! For bipartite systems
if ( nss == 2 ) then
  if ( ssys(1) == 0 ) then ; call partial_trace_a_he(rho, di(1), di(2), rhor)
  else if ( ssys(2) == 0 ) then ; call partial_trace_b_he(rho, di(1), di(2), rhor) 
  endif
  rhored = rhor
! For multipartite systems
else if ( nss >= 3 ) then
  ! Left partial traces
    l = 0 ;   do ;   if ( ssys(l+1) == 1 ) exit ;   l = l + 1 ;   enddo  ! l defines up to which position we shall trace out
    if ( l == 0 ) then
      mat1 = rho  ! This matrix shall be used below if l = 0
    else if ( l > 0 ) then
      if ( l == 1 ) then ;   da = di(1) ;   else ;   da = product(di(1:l)) ;  endif ;   db = d/da
      allocate( mat1(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, mat1)
      !!!!! To test the PTr with Bloch parametrization, comment the last line and uncomment the next one
      !allocate( mat1(1:db,1:db) ) ;   call partial_trace_a_bloch(rho, da, db, mat1)
      d = db  ! After this operation, the matrix left over is mat1, whose dimension is d = db  
      rhored = mat1    
    endif
  ! Right partial traces
    k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo  ! k defines up to which position we shall trace out
    if ( k == (nss+1) ) then
      mat2 = mat1  ! This matrix shall be used below if k = nss+1
      rhored = mat2
    else if ( k < (nss+1) ) then
      if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;  endif ;   da = d/db
      allocate( mat2(1:da,1:da) ) ;   call partial_trace_b_he(mat1, da, db, mat2) ;   deallocate( mat1 )
      !!!!! To test the PTr with Bloch parametrization, comment the last line and uncomment the next one
      !allocate( mat2(1:da,1:da) ) ;   call partial_trace_b_bloch(mat1, da, db, mat2) ;   deallocate( mat1 )
      d = da  ! After this operation, the matrix left over is mat2, whose dimension is d = da
      rhored = mat2
    endif
    if ( l == 0 .and. k == (nss+1) ) deallocate( mat1 )  ! To avoid allocating mat1 two times in these cases
  ! Inner partial traces
    if ( (k-l) > 3 ) then  ! If (k-l) <= 3 there is no need to take inner partial traces
    do j = (l+2), (k-2)
      if ( ssys(j) == 0 ) then
        db = di(j) ;   if ( j == (k-2) ) then ;   dc = di(k-1) ;   else ;   dc = product(di(j+1:k-1)) ;  endif ;   da = d/(db*dc)
        allocate( mat1(1:da*dc,1:da*dc) ) ;   call partial_trace_3_he(mat2, da, db, dc, mat1) ;   rhored = mat1
        d = da*dc ;   deallocate( mat2 ) ;   allocate( mat2(1:d,1:d) ) ;   mat2 = mat1 ;   deallocate( mat1 )
      endif
    enddo
    endif
endif

rhor = rhored

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_a_he(rho, da, db, rho_b)  ! Returns the left partial trace (over a), for a bi-partite matrix (the HERMITIAN case)
implicit none
integer :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_b(1:db,1:db)  !  Reduced matrix
integer :: j, k, l  ! Auxiliary variable for counters

rho_b = 0.d0
do j = 1, db ;   do k = j, db
  do l = 1, da ;   rho_b(j,k) = rho_b(j,k) + rho((j-1)*da+l,(k-1)*da+l) ;   enddo
  if ( j /= k ) rho_b(k,j) = conjg(rho_b(j,k))
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_b_he(rho, da, db, rho_a)  ! Returns the right partial trace (over b), for a bi-partite matrix (the HERMITIAN case)
implicit none
integer :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_a(1:da,1:da)  !  Reduced matrix
integer :: j, k, l  ! Auxiliary variables for counters

rho_a = 0.d0
do j = 1, da ;   do k = j, da
  do l = 1, db ;   rho_a(j,k) = rho_a(j,k) + rho((j-1)*db+l,(k-1)*db+l) ;   enddo
  if ( j /= k ) rho_a(k,j) = conjg(rho_a(j,k))
enddo ;   enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine partial_trace_3_he(rho, da, db, dc, rho_ac)  ! Returns the inner partial trace for a three-partite matrix (for the HERMITIAN case)
implicit none
integer :: da, db, dc  ! Dimension of the three sub-systems (the dimension of the whole system is d = da*db*dc)
complex(8) :: rho(1:da*db*dc,1:da*db*dc)  ! Three-partite matrix (computational basis representation of the ragarded operator)
complex(8) :: rho_ac(1:da*dc,1:da*dc)  ! Bipartite reduced matrix
integer :: j, k, l, m, o  ! Auxiliary variables for counters
integer :: cj, ck, ccj, cck  ! Auxiliary variables for the dimensions

rho_ac = (0.d0,0.d0)
doj: do j = 1, da ;   dol: do l = 1, dc ;   cj = (j-1)*dc+l ;   ccj = (j-1)*db*dc+l
dom: do m = 1, da ;   doo: do o = 1, dc ;   ck = (m-1)*dc+o
  if ( cj > ck ) then ;   if ( o < dc ) then ;   cycle doo ;   else if ( o == dc ) then ;  cycle dom ;   endif ;   endif
  cck = (m-1)*db*dc+o
  dok: do k = 1, db ;   rho_ac(cj,ck) = rho_ac(cj,ck) + rho(ccj+(k-1)*dc,cck+(k-1)*dc) ;   enddo dok
  rho_ac(ck,cj) = conjg(rho_ac(cj,ck))
enddo doo ;   enddo dom ;   enddo dol ;   enddo doj

end
!-----------------------------------------------------------------------------------------------------------------------------------
!###################################################################################################################################