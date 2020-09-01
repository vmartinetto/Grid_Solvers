      Module sparseMod
!
!     A module for doing sparse linear algebra operations
!
      Implicit None
!
!     define the sparse vector type
!
      Type, Public:: sparseVector
        Real, Dimension(:), Allocatable:: realVector
        Integer, Dimension(:), Allocatable:: indexVector
        Integer:: nDimDense, nDimSparse
        Contains
          Procedure, Public:: set => sparseVector_set
      End Type sparseVector
!
!     Procedure Interfaces
!

!
!     Subroutines
!
      Contains

      Subroutine sparseVector_set(mySV,realV,indexV,nDimDense)
!
!     This subroutine is used to set up a sparse vector object given
!     the vector of real values and their indices as well as the
!     diensions of the dense vector.
!
        Implicit None
        Class(sparseVector), Intent(InOut):: mySV
        Real, Dimension(:), Intent(In):: realV
        Integer, Dimension(:), Intent(In):: indexV
        Integer, Intent(In):: nDimDense
!
        mySV%nDimDense = nDimDense
        mySV%nDimSparse = Size(realV)
        Allocate(mySV%realVector(mySV%nDimSparse), &
          mySV%indexVector(mySV%nDimSparse))
        mySV%realVector = realV
        mySV%indexVector = indexV
!
        Return
      End Subroutine sparseVector_set
!
      End Module sparseMod
