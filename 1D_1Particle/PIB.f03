Include 'SparseMod.f03'
Include 'Hamiltonian_mod.f03'
Include 'Crank_Nicolson_mod.f03'
      Program PIB
      Use sparseMod
      Use HamiltonianMod
      Use CrankNicolsonMod
!
!     All this file is for is providing the starting point, strt, and
!     the end of the calculation region, endd, as well as 
!     the size of each step, h, and defining the potential, as
!     a simple function Potential(x) to be used to make the hamiltonian with
!     the specified finite difference method.
!
      Implicit None
      Type(bandMatrix):: H5, H3
      Integer:: N, funit = 10, i, M = 1000
      Real:: dx
      Real, Dimension(:,:), Allocatable:: evecs,evec3
      Real, Dimension(:), Allocatable:: evals,evals3
!
      Call Hamiltonian_5point(H5,0.0,1.0,M)
!
      N = H5%nDimDense(1)
      Allocate(evecs(N,N),evals(N))
      Call Hsolve(H5,evecs,evals)
!
      write(*,*)
      write(*,*) 'GS eval:'
      Write(*,*) evals(1)
!
!     write evec to external file for plotting with gnuplot
!
      Open(unit = funit, file = 'PIBGS.dat')
      Do i = 1, M
        Write(funit,*) i,evecs(i,1)
      End Do
      Close(unit = funit)
      Open(unit = funit, file = 'PIBFES.dat')
      Do i = 1, M
        Write(funit,*) i,evecs(i,2)
      End Do
      Close(unit = funit)
!
      M = 3000
      Call Hamiltonian_3point(H3,0.0,1.0,M)
      N = H3%nDimDense(1)
      Allocate(evals3(N),evec3(N,N))
      Call Hsolve(H3,evecs,evals)
!
      write(*,*)
      write(*,*) 'GS eval:'
      Write(*,*) evals(1)
!
      End Program PIB
!
      Function Potential(x) Result(V)
!
!     change the definition of V to change the potential of the hamiltonian
!
        Implicit None
        Real, Intent(In):: x
        Real:: V
!
        V = 0
!
        Return
      End Function Potential
