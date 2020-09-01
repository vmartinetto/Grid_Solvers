include "sparseMod.f03"      
      Program sparseTest
      Use sparseMod
!
!     Test file to check the sparse linear algebra module
!
      Implicit None
      Integer:: i
      Integer:: nDimDense, nDimSparse
      Integer, Dimension(:), Allocatable:: indexVec1
      Real, Dimension(:), Allocatable:: realVec1
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






      End Program sparseTest
