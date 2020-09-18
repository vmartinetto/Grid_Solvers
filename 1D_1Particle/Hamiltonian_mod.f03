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
      End Module HamiltonianMod
