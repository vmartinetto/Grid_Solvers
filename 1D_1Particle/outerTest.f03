include "sparseMod.f03"
      Program outerTest
      Use sparseMod
!
!     Test file to check the cross product functions of the sparse
!     linear algebra module
!
      Implicit None
      Type(sparseVector):: sparseVec1
      Real, Dimension(3):: vector1, vector2
      Real, Dimension(3,3):: v1xv2, v2xv1, dsmat
!
!     set sparseVectors
!
      Call sparseVec1%set((/6.0,7.0/),(/1,2/),3)
!
!     set the vectors 1 and 2
!
      vector1 = (/1,2,3/)
      vector2 = (/2,3,4/)
!
!     do some dense dense outer products
!
      v1xv2 = outer_product(vector1,vector2)
      v2xv1 = outer_product(vector2,vector1)
!
!     print the dense dense outer product matrices
!
      Write(*,*) 'vector 1: ', vector1
      Write(*,*) 'vector 2: ', vector2
      Write(*,*) 'outer(vec1,vec2):'
      Call print_matrix_full_real(v1xv2)
      Write(*,*) 'outer(vec2,vec1):'
      Call print_matrix_full_real(v2xv1)
!
!     do the dense sparse outer product
!
      dsmat = outer_product(vector2,sparseVec1)
!
!     print the dense sparse matrix
!
      Write(*,*) 'Sparse Vector 1: '
      Call sparseVec1%print(header='sp1= ')
      Write(*,*) 'outer product of vec2 sparse1:'
      Call print_matrix_full_real(dsmat)
!
      End Program outerTest
