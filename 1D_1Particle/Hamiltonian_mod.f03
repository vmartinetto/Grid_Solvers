      Module HamiltonianMod
      Use sparseMod
!
      Contains
!
!     A module for computing and solving the Hamiltonian using finite
!     difference methods, lapack, and sparse matrix objects from the
!     sparseMod file
!
      Subroutine Hamiltonian_5point(H,strt,endd,M)
!
!     Takes the number of points M, and the step size, dx, as well as the
!     start and end of the calculation region, strt and endd, and returns
!     a band matrix object from the sparse mod module.
!
        Implicit None
        Class(bandMatrix), Intent(InOut):: H
        Real, Intent(In):: strt,endd
        Integer, Intent(In):: M
        Real:: Potential, dx
        Real, Dimension(:,:), Allocatable:: band
        Real, Dimension(M):: v1, v2, main, x
        Integer:: i
!
        Allocate(band(M,5))
!
        dx = (endd-strt)/M
        v1 = 1/(24*dx**2)
        v1(1) = 0
        v1(2) = 0
        v2 = -2/(3*dx**2)
        v2(1) = 0
        x = (/(strt+i*dx , i=0,M-1)/)
        main = 0  
        Do i = 1,M
          main(i) = Potential(x(i)) + 15/(12*dx**2)
        End Do
        band(:,1) = v1 
        band(:,2) = v2
        band(:,3) = main
        band(:,4) = v2
        band(:,5) = v1
        Call H%set(band,(/M,M/),3)
!
        Return
      End Subroutine Hamiltonian_5point
!
      Subroutine Hamiltonian_3point(H,strt,endd,M)
!
!     Does the same as the 5 point stencil method but represents the hamiltonian
!     With a 3 point stencil .
!
        Implicit None
        Class(bandMatrix), Intent(InOut):: H
        Real, Intent(In):: strt, endd
        Integer, Intent(In):: M
        Integer:: i
        Real:: Potential, dx
        Real, Dimension(:,:), Allocatable:: band
        Real, Dimension(M):: off, main, x
!
        Allocate(band(M,3))
!
        dx = (endd-strt)/M
        off = -1/(2*dx**2)
        off(1) = 0
        x = (/(strt+i*dx , i=0,M-1)/)
        main = 0  
        Do i = 1,M
          main(i) = Potential(x(i)) + 1/(dx**2)
        End Do
        band(:,1) = off 
        band(:,2) = main
        band(:,3) = off
        Call H%set(band,(/M,M/),2)
        
!
      End Subroutine Hamiltonian_3point
!
      Subroutine Hsolve(H,eVecs,eVals)
!
!     takes the Hamiltonian as a band matrix object and converts it to 
!     upper packed form consitent with the LaPack documentation for band
!     matrices. It then solves the eigenvectors and eigenvalues and returns
!     them in eVecs and eVals
!
        Implicit None
        Class(bandMatrix), Intent(In):: H
        Integer:: info, N, Lwork, Liwork, k
        Integer, Dimension(:), Allocatable:: Iwork
        Real, Dimension(:,:), Allocatable:: A,Q,Z
        Real, Dimension(:), Allocatable:: D,E,work,W, Lworkvec
        Real, Dimension(:,:), Intent(InOut):: evecs
        Real, Dimension(:), Intent(InOut):: evals
!
!     Turn the band matrix into a upper packed form A for Lapack to use
!
        N = H%nDimDense(1)
        Allocate(A(H%maindiagonal,H%nDimSparse(1)))
        A = H%uppersymetricpack()
!
!     Reduce the band matrix into a tridiagonal matrix
!
        Allocate(Q(N,N), D(H%nDimSparse(1)),E(H%nDimSparse(1)-1), &
                 work(H%nDimSparse(1)),W(N),Z(N,N))
        Call SSBTRD('V', 'U', N, H%mainDiagonal-1, &
                    A, H%mainDiagonal, D, E, Q, N, &
                    work, info)
!
!     Solve the eigenvectors and eigenvalues of the reduced tridiagonal.
!     E need to be set to EE because SSTEGR calls for the off-diagonal to
!     be N elements but SSBTRD calss for N-1. The output Eigenvectors are
!     in Z and the output values in W. 
!
        k = int(log(float(N))/log(2.0))
        If (2**k < N) Then
          k = k+1
        End If
        Lwork = 1 + 3*N + 2*N*k + 4*N**2
        Liwork = 6 + 6*N + 5*N*k
        Allocate(Lworkvec(Lwork),Iwork(Liwork))
        Call SSTEDC('V', N, D, E, Q, N, Lworkvec, Lwork, &
                    Iwork, Liwork, info) 
!
!     The eigenvectors of the tri-diagonal stored in Z need to be converted 
!     back to the eigenvectors of the full matrix H. This is done by the multiplication
!     evecs = matmul(Q,Z). The eigenvalues are the same
!
        evecs = Q
        evals = D
!
      End Subroutine
!
      End Module HamiltonianMod
