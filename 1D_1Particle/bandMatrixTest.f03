Include "sparseMod.f03"
      Program bandTest
      Use sparseMod
!
!     Test file for the band Matrix object from the sparse linear
!     algebra modue 
!
      Implicit None
      Real, Dimension(:), Allocatable:: band1, band2, band3
      Real, Dimension(:,:), Allocatable:: bands, packed
      Type(bandMatrix):: bandMat, denseBandMat
!
!     Test the set subroutine
!
      Call bandMat%set1D((/1.0,2.0,3.0/),(/3,3/))
      Call bandMat%print(header='These are the band')
!
!     Make a denser band matrix and print it
!
      Allocate(band1(4),band2(4),band3(4))
      band1 = (/0.0,1.0,2.0,3.0/)
      band2 = (/4.0,5.0,6.0,7.0/)
      band3 = (/0.0,8.0,9.0,10.0/)
      Allocate(bands(4,3))
      bands(:,1) = band1
      bands(:,2) = band2
      bands(:,3) = band3
!
!     Set and test on the denser band matrix
!
      Call denseBandMat%set(bands,(/4,4/),2)
      Call denseBandMat%print(header='These are the bands of the test 4X4')
!
!
!
      Allocate(packed(denseBandMat%mainDiagonal, &
        denseBandMat%nDimSparse(1)))
      packed = denseBandMat%uppersymetricpack()
      Call print_matrix_full_real(packed)
!
      End Program bandTest
