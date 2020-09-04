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
          Procedure, Public:: print => sparseVector_print
          Procedure, Public:: sddot_product => sparse_dense_dot_product
          Procedure, Public:: ssdot_product => sparse_sparse_dot_product
      End Type sparseVector
!
!     Procedure Interfaces
!
      Interface dot_product
        module procedure sparse_dense_dot_product
        module procedure sparse_sparse_dot_product
      End Interface dot_product
!
      Interface outer_product
        Module Procedure outer_product
        Module Procedure dense_sparse_outer_product
      End Interface outer_product
!
!     Subroutines
!
      Contains
!
!---------------------- VECTOR OPERATIONS -------------------------------
!
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
      Subroutine sparseVector_print(mySV,iOut,header)
!
!     This function will pring the sparse vector contents
!
        Implicit None
        Class(sparseVector), Intent(In):: mySV
        Integer, Optional, Intent(In):: iOut
        Character(len=*), Optional, Intent(In):: header
        Integer:: myiOut, i
!
1000    Format(3x,'i=',I4,' v(i)=',f12.5)
!
        myiOut = 6
        If(Present(iOut)) myiOut = iOut
        If(Present(header)) Write(myiOut,'(A)') Trim(header)
        Do i = 1,mySV%nDimSparse
          Write(myiOut,1000) mySV%indexVector(i),mySV%realVector(i)
        End Do
!
        Return
      End Subroutine sparseVector_print
!
      Function sparse_dense_dot_product(mySV,dVec) Result(dotProduct)
!
!     This function does the dot product between a sparse vector object
!     and a FORTRAN intrinsic dense vector
!
        Implicit None
        Real:: dotProduct
        Class(sparseVector), Intent(In):: mySV
        Real, Dimension(:), Intent(In):: dVec
        Integer:: i
!
        dotProduct = 0
        Do i = 1,mySV%nDimsparse
          dotProduct = dotProduct + mySV%realVector(i)* &
            dVec(mySV%indexVector(i))
        End Do
!
        Return
      End Function sparse_dense_dot_product
!
      Function sparse_sparse_dot_product(mySV1,mySV2) Result(dotProduct)
!
!     This function takes the dot_product of two sparse vector objects
!     that have their indices sorted in order
!
        Implicit None
        Real:: dotProduct
        Integer:: i, j
        Class(sparseVector), Intent(In):: mySV1,mySV2
!
!     Iterates throught the indices of the sparse vectors checking if
!     they match, if they do the multiplication is added to the dot
!     product. If not, the index that is smaller is increased. continues
!     until one of the vectors has been iterated through.
!
        dotProduct = 0
        i = 1
        j = 1
        Do While ((i.le.mySV1%nDimSparse).and.(j.le.mySV1%nDimSparse))
          If (mySV1%indexVector(i).eq.mySV2%indexVector(j)) Then
            dotProduct = dotProduct + mySV1%realVector(i)* &
              mySV2%realVector(j)
            i = i + 1
            j = j + 1
          Else If (mySV1%indexVector(i).gt.mySV2%indexVector(j)) Then
            j = j + 1
          Else
            i = i + 1
          End If
        End Do
!
        Return
      End Function sparse_sparse_dot_product
!
      Function outer_product(vector1,vector2) Result(Matrix)
!
!     This function takes the outer product of two FORTRAN intrinsic
!     vectors one MX1 and the other 1XN and returns the MXN matrix that is produced
!
        Implicit None
        Real, Dimension(:), Intent(In):: vector1, vector2
        Real, Dimension(:,:), Allocatable:: Matrix
        Integer:: i
!
        Allocate(Matrix(Size(vector1),Size(vector2)))
        Do i = 1,Size(vector2)
          Matrix(:,i) = vector1*vector2(i)
        End Do
!
        Return
      End Function
!
      Function dense_sparse_outer_product(vector,mySV) Result(Matrix)
!
!     takes the outer product of a dense MX1 FORTRAN column vector and a
!     sparse 1XN row vector resulting in a MXN semi-sparse Matrix
!
        Implicit None
        Real, Dimension(:,:), Allocatable:: Matrix
        Real, Dimension(:), Intent(In):: vector
        Class(sparseVector), Intent(In):: mySV
        Integer:: i
!
        Allocate(Matrix(Size(vector),mySV%nDimDense))
        Matrix = 0
        Do i = 1,mySV%nDimSparse
          Matrix(:,mySV%indexVector(i)) = vector*mySV%realVector(i)
        End Do
!
        Return
      End Function dense_sparse_outer_product
!
!----------------------- MATRIX OPERATIONS -----------------------------
!
      Subroutine print_matrix_full_real(amat)
!
!     Prints the FORTRAN intrinsic matrix, amat, to the standard output.
!     the input matrix should be unpacked.
!
        Implicit None
        Real, Dimension(:,:), Intent(In):: amat
        Integer:: i, j, iFirst, Ilast
        Integer, Parameter:: iOut=6, nColumns = 5
!
1000    Format(1X,A)
2000    Format(5X,5(7X,I7))
2010    Format(1X,I7,5F14.6)
!
        Do iFirst = 1, Size(amat,2), nColumns
          iLast = Min(iFirst+nColumns-1,Size(amat,2))
          Write(iOut,2000) (i, i=iFirst,iLast) ! implied do?
          Do i = 1, Size(amat,1)
            Write(iOut,2010) i, (Amat(i,j), j=iFirst , iLast)
          End Do
        End Do
!
        Return
        End Subroutine print_matrix_full_real
!
      End Module sparseMod
