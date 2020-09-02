include "sparseMod.f03"      
      Program sparseTest
      Use sparseMod
!
!     Test file to check the sparse linear algebra module
!
      Implicit None
      Type(sparseVector):: sparseVec1
      Integer:: i
      Integer:: nDimDense, nDimSparse
      Integer, Dimension(:), Allocatable:: indexVec1
      Real, Dimension(:), Allocatable:: realVec1
      Real, Dimension(:), Allocatable:: denseVec1
!
!     Input from the user
!
      Write(*,*) 'What is the dimension of the vector space'
      Read(*,*) nDimDense
      Write(*,*) 'How many non-zero elements are ther'
      Read(*,*) nDimSparse
      Allocate(realVec1(nDimSparse),indexVec1(nDimSparse))
      Write(*,*) 'Please enter the non-zero elements and their', & 
      ' indices in pairs'
      Do i = 1, nDimSparse
        Read(*,*) realVec1(i),indexVec1(i)
      End Do
!
!     Set up a sparse vector and test the print function
!
      Call sparseVec1%set(realVec1,indexVec1,nDimDense)
      Call sparseVec1%print(header='Vec1= ')
!
!     Test the sparse dense dot_product
!
      Allocate(denseVec1(nDimDense))
      Call Random_Number(denseVec1)
      Write(*,*)
      Write(*,*) 'This is the test dense vector'
      Do i = 1,nDimDense
        Write(*,*) denseVec1(i)
      End Do
      Write(*,*)
      Write(*,*) 'this sparse dense dot:', dot_product(sparseVec1,denseVec1)
!
      End Program sparseTest
