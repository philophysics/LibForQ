!###################################################################################################################################
!                                       PARTIAL TRANSPOSE
!###################################################################################################################################
subroutine partial_transpose_a(da, db, rho, rho_ta)  ! Returns its partial transpose with relation to system A
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_ta(1:da*db,1:da*db)  ! Bipartite original and transposed states
integer :: ja, jb, ka, kb  ! Auxiliary variable for counters

forall ( ja = 1:da, ka = 1:da, jb = 1:db, kb = 1:db )
  rho_ta((ka-1)*db+jb,(ja-1)*db+kb) = rho((ja-1)*db+jb,(ka-1)*db+kb)
end forall

end
!------------------------------------------------------------------------------------------------------------------------------------
subroutine partial_transpose_b(da, db, rho, rho_tb)  ! Returns the partial transpose with relation to system B
implicit none
integer :: da, db ! Dimensions of the subsystems
complex(8) :: rho(1:da*db,1:da*db), rho_tb(1:da*db,1:da*db)  ! Bipartite original and transposed states
integer :: ja, jb, ka, kb  ! Auxiliary variable for counters

forall ( ja = 1:da, ka = 1:da, jb = 1:db, kb = 1:db )
  rho_tb((ja-1)*db+kb,(ka-1)*db+jb) = rho((ja-1)*db+jb,(ka-1)*db+kb)
end forall

end
!###################################################################################################################################