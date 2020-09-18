Include 'Hamiltonian_mod.f03'
      Program PIB
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
!
      Call Hamiltonian_5point(H,0.0,2.0,1.0,3)
      Call H%print
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
        V = x
!
        Return
      End Function Potential
