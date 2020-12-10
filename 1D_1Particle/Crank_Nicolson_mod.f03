        Module CrankNicolsonMod
        Use sparseMod
!
        Contains
!
!       A module for solving Crank Nicolson time-propogation in 1D using
!       a 3 point and 5 point stencils.
!
        Subroutine CrankNicolson_3point(A,B,M,dx,dt)
!
!       Takes the number of grid points as well as the time and space
!       step size in the calculation and computes
!       the matrices A and B necessay for a crank nicolson time
!       propogation using a three point stencil.
!
          Implicit None
          Integer, Intent(In):: M
          Complex, Dimension(M,M), Intent(InOut):: A,B
          Complex:: a1, a2, b1, b2
          Complex:: CN
          Real:: hbar = 1.05457817e-34, me=9.109e-31, dx, dt
          Integer:: i
!
          CN = Complex(0,(dt*hbar)/(4*me*dx**2))
          a1 = 1+2*CN
          a2 = -CN
          Do i = 1,M
            A(i,i) = a1
            If (i<M) Then
              A(i+1,i) = a2
              A(i,i+1) = a2
            End If
          End Do
!
          b1 = 1-2*CN
          b2 = CN
          Do i = 1,M
            B(i,i) = b1
            If (i<M) Then
              B(i+1,i) = b2
              B(i,i+1) = b2
            End If
          End Do
!
        End Subroutine CrankNicolson_3point
!
        Subroutine CrankNicolson_5point(A,B,M,dx,dt)
!
!       does same as three point method but uses a five point stencil
!       instead.
!
          Implicit None
          Integer, Intent(In):: M
          Complex, Dimension(M,M), Intent(InOut):: A,B
          Complex:: a1, a2, a3, b1, b2 , b3
          Complex:: CN
          Real:: hbar = 1.05457817e-34, me=9.109e-31, dx, dt
          Integer:: i
!
          CN = Complex(0,(dt*hbar)/(24*me*dx**2))
          a1 = 1+30*CN
          a2 = -16*CN
          a3 = CN
          Do i = 1,M
            A(i,i) = a1
            If (i<M-1) Then
              A(i+1,i) = a2
              A(i,i+1) = a2
              A(i+2,i) = a3
              A(i,i+2) = a3
            ElseIf (i.eq.M-1) Then
              A(i+1,i) = a2
              A(i,i+1) = a2
            End If
          End Do
!
          b1 = 1-30*CN
          b2 = 16*CN
          b3 = -CN
          Do i = 1,M
            B(i,i) = b1
            If (i<M-1) Then
              B(i+1,i) = b2
              B(i,i+1) = b2
              B(i+2,i) = b3
              B(i,i+2) = b3
            ElseIf (i.eq.M-1) Then
              B(i+1,i) = b2
              B(i,i+1) = b2
            End If
          End Do
!
        End Subroutine CrankNicolson_5point
!
        Subroutine CN_Solve(psit2,psit1,A,B)
!
!       Takes the matrices A and B and neede for Crank-Nicolson time
!       propagation and solves for psit2 from psit1.
!
          Implicit None
          Complex, Dimension(:), Intent(In):: psit1
          Complex, Dimension(:), Intent(InOut):: psit2
          Complex, Dimension(:), Allocatable:: Bpsi
          Complex, Dimension(:,:), Intent(InOut):: A,B
        End Subroutine CN_Solve
!
        End Module CrankNicolsonMod
