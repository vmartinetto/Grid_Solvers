include "sparseMod.f03"      
      Program sparseTest
      Use sparseMod
!
!     Test file to check the sparse linear algebra module
!
      Implicit None
      Type(sparseVector):: sparseVec1, sparseVec2
      Integer:: i
      Integer:: nDimDense, nDimSparse1, nDimSparse2
      Integer, Dimension(:), Allocatable:: indexVec1, indexVec2
      Real, Dimension(:), Allocatable:: realVec1, realVec2
      Real, Dimension(:), Allocatable:: denseVec1
!
!     Input from the user
!
      Write(*,*) 'What is the dimension of the vector space'
      Read(*,*) nDimDense
      Write(*,*) 'How many non-zero elements are ther'
      Read(*,*) nDimSparse1
      Allocate(realVec1(nDimSparse1),indexVec1(nDimSparse1))
      Write(*,*) 'Please enter the non-zero elements and their', & 
      ' indices in pairs'
      Do i = 1, nDimSparse1
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
      Write(*,*) 'sparse dense dot product:', dot_product(sparseVec1,denseVec1)
!
!     Test the sparse spare dot product
!
      Write(*,*) 'The second sparse vector has ', nDimDense, 'elements'
      Write(*,*) 'How many non-zero elements are ther'
      Read(*,*) nDimSparse2
      Allocate(realVec2(nDimSparse2),indexVec2(nDimSparse2))
      Write(*,*) 'Please enter the non-zero elements and their', & 
      ' indices in pairs'
      Do i = 1, nDimSparse2
        Read(*,*) realVec2(i),indexVec2(i)
      End Do
!
!     Set up a sparse vector and test the print function
!
      Call sparseVec1%set(realVec1,indexVec1,nDimDense)
      Call sparseVec1%print(header='Vec1= ')
!
      End Program sparseTest
