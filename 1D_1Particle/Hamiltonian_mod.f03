Include 'sparseMod.f03'
      Module HamiltonianMod
      Use sparseMod
!
      Contains
!
!     A module for computing and solving the Hamiltonian using finite
!     difference methods, lapack, and sparse matrix objects from the
!     sparseMod file
!
      Subroutine Hamiltonian_5point(H,strt,endd,dx,M)
!
!     Takes the number of points. M, and the step size, dx, as well as the
!     start and end of the calculation region, strt and endd, and returns
!     a band matrix object from the sparse mod module.
!
        Implicit None
        Class(bandMatrix), Intent(InOut):: H
        Real, Intent(In):: strt,endd,dx
        Integer, Intent(In):: M
        Real:: Potential
        Real, Dimension(:,:), Allocatable:: band
        Real, Dimension(M):: v1, v2, main, x
        Integer:: i
!
        Allocate(band(M,5))
!
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
      Subroutine Hsolve(H,eVecs,eVals)
!
!     takes the Hamiltonian as a band matrix object and converts it to 
!     upper packed form consitent with the LaPack documentation for band
!     matrices. It then solves the eigenvectors and eigenvalues and returns
!     them in eVecs and eVals
!
        Implicit None
        Class(bandMatrix), Intent(In):: H
        Integer:: info
        Real, Dimension(:,:), Allocatable:: A,Q
        Real, Dimension(:), Allocatable:: D,E,work
        Real, Dimension(:,:), Intent(InOut):: eVecs
        Real, Dimension(:), Intent(InOut):: evals
!
!     Turn the band matrix into a upper packed form A for Lapack to use
!
        Allocate(A(H%maindiagonal,H%nDimSparse(1)))
        A = H%uppersymetricpack()
!
!     Reduce the band matrix into a tridiagonal matrix
!
        Allocate(Q(H%nDimDense(1),H%nDimDense(1)), D(H%nDimSparse(1)), &
                 E(H%nDimSparse(1)-1), work(H%nDimSparse(1)))
        Call SSBTRD('V', 'U', H%nDimDense(1), H%mainDiagonal-1, &
                    A, H%mainDiagonal, D, E, Q, H%nDimDense(1), &
                    work, info)
!
!     Solve the eigenvectors and eigenvalues of the reduced tridiagonal.
!     E need to be set to EE because SSTEGR calls for the off-diagonal to
!     be N elements but SSBTRD calss for N-1.
!
        Allocate(EE(H%nDimSparse(1))
        EE(1,)
        Call SSTEGR('V', 'A', H%nDimDense(1), D, EE)
!
      End Subroutine
!
      End Module HamiltonianMod
