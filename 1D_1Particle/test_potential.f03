      Program potentialTest
!
!     This is a test to test the scope of functions in Fortran
!
!
      Implicit None
      Interface
        Function V_vector(x) Result(y)
          Real, Dimension(:), Intent(In):: x
          Real:: V
          Real, Dimension(Size(x)):: y
          Integer:: i
        End Function
      End Interface

      Real, Dimension(3):: x,y
!
      x = (/1.0,2.0,3.0/)
      y = V_vector(x)
      Write(*,*) y
!
      End Program potentialTest
!
      Function V(x) Result(y)
!
!     Define potential here
!
        Implicit None
        Real, Intent(In):: x
        Real:: y
!
        y = 4*x
!
        Return
      End Function V
!
      Function V_vector(x) Result(y)
!
!     Uses previously designated V(x) to make and return a vector of the
!     potential across a range of x
!
        Implicit None
        Real, Dimension(:), Intent(In):: x
        Real:: V
        Real, Dimension(Size(x)):: y
        Integer:: i
!
        Do i = 1, Size(x)
          y(i) = V(x(i))
        End Do
!
        Return
      End Function V_vector
