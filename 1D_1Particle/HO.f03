Include 'Hamiltonian_mod.f03'
      Program Harmonic
      Use sparseMod
      Use HamiltonianMod
!
!     All this file is for is providing the starting point, strt, and
!     the end of the calculation region, endd, as well as 
!     the size of each step, h, and defining the potential, as
!     a simple function Potential(x) to be used to make the hamiltonian with
!     the specified finite difference method.
!
      Implicit None
      Type(bandMatrix):: H
      Integer:: N, funit = 10, i, M = 3000
      Real:: dx
      Real, Dimension(:,:), Allocatable:: evecs
      Real, Dimension(:), Allocatable:: evals
!
      Call Hamiltonian_5point(H,-5.0,5.0,M)
!
      N = H%nDimDense(1)
      Allocate(evecs(N,N),evals(N))
      Call Hsolve(H,evecs,evals)
!
      write(*,*)
      write(*,*) 'GS eval:'
      Write(*,*) evals(1)
!
!     write evec to external file for plotting with gnuplot
!
      Open(unit = funit, file = 'HOGS.dat')
      Do i = 1, M
        Write(funit,*) i,evecs(i,1)
      End Do
      Close(unit = funit)
      End Program Harmonic
!
      Function Potential(x) Result(V)
!
!     change the definition of V to change the potential of the hamiltonian
!
        Implicit None
        Real, Intent(In):: x
        Real:: V, k 
!
        k = 3.3
        V = .5*k*x**2
!
        Return
      End Function Potential
