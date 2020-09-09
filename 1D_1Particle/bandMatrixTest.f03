Include "sparseMod.f03"
      Program bandTest
      Use sparseMod
!
!     Test file for the band Matrix object from the sparse linear
!     algebra modue 
!
      Implicit None
      Type(bandMatrix):: bandMat
!
!     Test the set subroutine
!
      Call bandMat%set1D((/1.0,2.0,3.0/),(/3,3/))
      Call bandMat%print(header='These are the band')
!
      End Program bandTest
