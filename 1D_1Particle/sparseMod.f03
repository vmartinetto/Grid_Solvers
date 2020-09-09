      Module sparseMod
!
!     A module for doing sparse linear algebra operations
!
      Implicit None
!-------------------------- Types --------------------------------------
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
!     Define Band matrix type
!
      Type, Public:: bandMatrix
        Integer, Dimension(2):: nDimDense, nDimSparse
        Real, Dimension(:,:), Allocatable:: bandDiags
        Integer:: mainDiagonal
        Contains
          Procedure, Public:: set => bandMatrix_set
          Procedure, Public:: set1D => bandMatrix_oneDiag_set
          Procedure, Public:: print => bandMatrix_print
      End Type bandMatrix
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
        Module Procedure sparse_dense_outer_product
        Module Procedure sparse_sparse_outer_product
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
      Function sparse_dense_outer_product(mySV,vector1) Result(Matrix)
!
!     Takes the outer product of MX1 sparse vector object and a 1XN
!     FORTRAN intrinsic vector returning a resulting MXN Matrix
!
        Implicit None
        Real, Dimension(:,:), Allocatable:: Matrix
        Real, Dimension(:), Intent(In):: vector1
        Class(sparseVector), Intent(In):: mySV
        Integer:: i, j
!
        Allocate(Matrix(mySV%nDimDense,Size(vector1)))
        Matrix = 0
        Do j = 1,Size(vector1)
          Do i = 1,mySV%nDimSparse
            Matrix(mySV%indexVector(i),j) = mySV%realVector(i)* &
              vector1(j)
          End Do
        End Do
!
        Return
      End Function sparse_dense_outer_product
!
      Function sparse_sparse_outer_product(mySV1,mySV2) Result(Matrix)
!
!     Takes the outer product of a MX1 and a 1XN sparse vector object
!     and returns the resulting matrix
!
        Implicit None
        Real, Dimension(:,:), Allocatable:: Matrix
        Class(sparseVector), Intent(In):: mySV1, mySV2
        Integer:: i, j
!
        Allocate(Matrix(mySV1%nDimDense,mySV2%nDimDense))
        Matrix = 0
        Do j = 1, mySV1%nDimSparse
          Do i = 1, mySV2%nDimSparse
            Matrix(mySV1%indexVector(i),mySV2%indexVector(j)) = &
              mySV1%realVector(i)*mySV2%realVector(j)
          End Do
        End Do
!
        Return
      End Function sparse_sparse_outer_product
!
!----------------------- MATRIX OPERATIONS -----------------------------
!
      Subroutine bandMatrix_oneDiag_set(myBM,diag,nDimDense)
!
!     This subroutine is used to set up a band matrix when the only band
!     in the matix is along the main diagonal, as such this suboutine
!     only works for square matrices.
!
        Implicit None
        Class(bandMatrix), Intent(InOut):: myBM
        Real, Dimension(:), Intent(In):: diag
        Integer, Dimension(2), Intent(In):: nDimDense
        Integer:: i
!
        myBM%mainDiagonal = 1
        myBM%nDimDense = nDimDense
        myBM%nDimSparse = (/Size(diag),1/)
        Allocate(myBM%bandDiags(myBM%nDimSparse(1),myBM%nDimSparse(2)))
        Do i = 1,Size(diag)
          myBM%bandDiags(i,1) = diag(i)
        End Do
!
        Return
      End Subroutine bandMatrix_oneDiag_set
!
      Subroutine bandMatrix_set(myBM,diags,nDimDense,mainpos)
!     
!     This subroutine is used to set up a band matrix object given the
!     diagonals of the band matrix, the position of the main diagonal in
!     the set of column vectors, and the dimensions of the dense matrix 
!     as a list of two values, (/rows,columns/). Currently only works for
!     square matices.
!
        Implicit None
        Integer, Intent(In):: mainpos
        Class(bandMatrix), Intent(InOut):: myBM
        Real, Dimension(:,:), Intent(In):: diags
        Integer, Dimension(2), Intent(In):: nDimDense
!
        myBM%mainDiagonal = mainpos
        myBM%nDimDense = nDimDense
        myBM%nDimSparse = Shape(diags)
        Allocate(myBM%bandDiags(myBM%nDimSparse(1),myBM%nDimSparse(2)))
        myBM%bandDiags = diags
!
        Return
      End Subroutine bandMatrix_set
!
      Subroutine bandMatrix_print(myBM,iOut,header)
!
!     Prints out the bands of the band matrix and specifies thir
!     starting position in the upper left hand corner of the matrix
!
        Implicit None
        Class(bandMatrix), Intent(In):: myBM
        Integer, Optional, Intent(In):: iOut
        Character(len=*), Optional, Intent(In):: header
        Integer:: myiOut, i
!
1000    Format(3x,'i=',I4,'j= ',I4,' mat(i,j)= ')
!
        myiOut = 6
        If(Present(iOut)) myiOut = iOut
        If(Present(header)) Write(myiOut,'(A)') Trim(header)
        Do i = 1,myBM%nDimSparse(2)
          If (i.lt.myBM%mainDiagonal) Then
            Write(myiOut,1000) myBM%mainDiagonal-(i-1),1
            Write(myiOut,*) myBM%bandDiags(:,i)
          Else If (i.eq.myBM%mainDiagonal) Then
            Write(myiOut,1000) 1,1
            Write(myiOut,*) myBM%bandDiags(:,i)
          Else
            Write(myiOut,1000) 1,myBM%mainDiagonal+(i-1)
            Write(myiOut,*) myBM%bandDiags(:,i)
          End If
        End Do
!
        Return
      End Subroutine bandMatrix_print
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
