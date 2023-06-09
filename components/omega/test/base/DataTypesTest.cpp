//===-- Test driver for OMEGA data types -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA data types
///
/// This driver tests the definition of various file types for the OMEGA
/// model. In particular, it tests the length of various data types and
/// outputs a PASS if all of them are the expected length. It also tests
/// that a build with SINGLE_PRECISION converts the defaul real type to
/// single precision (4-byte) floating point.
///
//
//===-----------------------------------------------------------------------===/

#include <iostream>

#include "mpi.h"
#include "DataTypes.h"

int main(int argc, char *argv[])
{

   // initialize environments
   MPI_Init(&argc, &argv);
   yakl::init();

    // declare variables of each supported type
   OMEGA::I4 MyInt4 = 1;
   OMEGA::I8 MyInt8 = 2;
   OMEGA::R4 MyR4   = 3.0;
   OMEGA::R8 MyR8   = 4.0000000000001;
   OMEGA::Real MyReal = 5.000001;
   int SizeTmp = 0;

   // Check expected size (in bytes) for data types
   SizeTmp = sizeof(MyInt4); 
   if (SizeTmp == 4)
      std::cout << "Size of I4: PASS" << std::endl;
   else
      std::cout << "Size of I4: FAIL " << SizeTmp << std::endl;

   SizeTmp = sizeof(MyInt8); 
   if (SizeTmp == 8)
      std::cout << "Size of I8: PASS" << std::endl;
   else
      std::cout << "Size of I8: FAIL " << SizeTmp << std::endl;

   SizeTmp = sizeof(MyR4); 
   if (SizeTmp == 4)
      std::cout << "Size of R4: PASS" << std::endl;
   else
      std::cout << "Size of R4: FAIL " << SizeTmp << std::endl;

   SizeTmp = sizeof(MyR8); 
   if (SizeTmp == 8)
      std::cout << "Size of R8: PASS" << std::endl;
   else
      std::cout << "Size of R8: FAIL " << SizeTmp << std::endl;

   SizeTmp = sizeof(MyReal); 
#ifdef SINGLE_PRECISION
   if (SizeTmp == 4)
      std::cout << "Size of Real is 4: PASS" << std::endl;
   else
      std::cout << "Size of Real is 4: FAIL " << SizeTmp << std::endl;
#else
   if (SizeTmp == 8)
      std::cout << "Size of Real is 8: PASS" << std::endl;
   else
      std::cout << "Size of Real is 8: FAIL " << SizeTmp << std::endl;
#endif

   // Test creation of device arrays and copying to/from host
   // by initializing on the device, copying to host and comparing with
   // a reference host array.

   int NumCells = 100;
   int NumVertLvls = 100;
   int NumTracers = 4;
   int NumTimeLvls = 2;
   int NumExtra = 2;

   using yakl::c::parallel_for;
   using yakl::c::Bounds;

   // Test for 1DI4
   OMEGA::Array1DI4 TestArr1DI4("test",NumCells);
   OMEGA::ArrayHost1DI4 RefArr1DI4("ref",NumCells);

   for (int i=0; i < NumCells; ++i){
      RefArr1DI4(i) = i;
   }

   parallel_for(Bounds<1>(NumCells), YAKL_LAMBDA (int i) {
      TestArr1DI4(i) = i;
   });

   yakl::fence();
   auto TestHost1DI4 = TestArr1DI4.createHostCopy();

   int icount=0;
   for (int i=0; i < NumCells; ++i){
      if (TestHost1DI4(i) != RefArr1DI4(i)) ++icount;
   }
   TestHost1DI4.deallocate();
   RefArr1DI4.deallocate();
   TestArr1DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DI4 test: FAIL" << std::endl;

   // Test for 2DI4
   OMEGA::Array2DI4 TestArr2DI4("test",NumCells,NumVertLvls);
   OMEGA::ArrayHost2DI4 RefArr2DI4("ref",NumCells,NumVertLvls);

   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr2DI4(j,i) = i+j;
   }
   }

   parallel_for(Bounds<2>(NumCells,NumVertLvls), YAKL_LAMBDA (int j,int i) {
      TestArr2DI4(j,i) = i+j;
   });

   yakl::fence();
   auto TestHost2DI4 = TestArr2DI4.createHostCopy();

   icount=0;
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost2DI4(j,i) != RefArr2DI4(j,i)) ++icount;
   }
   }
   TestHost2DI4.deallocate();
   RefArr2DI4.deallocate();
   TestArr2DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DI4 test: FAIL" << std::endl;

   // Test for 3DI4
   OMEGA::Array3DI4 TestArr3DI4("test",NumTracers,NumCells,NumVertLvls);
   OMEGA::ArrayHost3DI4 RefArr3DI4("ref",NumTracers,NumCells,NumVertLvls);

   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr3DI4(k,j,i) = i+j+k;
   }
   }
   }

   parallel_for(Bounds<3>(NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int k, int j,int i) {
      TestArr3DI4(k,j,i) = i+j+k;
   });

   yakl::fence();
   auto TestHost3DI4 = TestArr3DI4.createHostCopy();

   icount=0;
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost3DI4(k,j,i) != RefArr3DI4(k,j,i)) ++icount;
   }
   }
   }
   TestHost3DI4.deallocate();
   RefArr3DI4.deallocate();
   TestArr3DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DI4 test: FAIL" << std::endl;

   // Test for 4DI4
   OMEGA::Array4DI4 TestArr4DI4("test",NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost4DI4 RefArr4DI4("ref",NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr4DI4(m,k,j,i) = i+j+k+m;
   }
   }
   }
   }

   parallel_for(Bounds<4>(NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int m, int k, int j,int i) {
      TestArr4DI4(m,k,j,i) = i+j+k+m;
   });

   yakl::fence();
   auto TestHost4DI4 = TestArr4DI4.createHostCopy();

   icount=0;
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost4DI4(m,k,j,i) != RefArr4DI4(m,k,j,i)) ++icount;
   }
   }
   }
   }
   TestHost4DI4.deallocate();
   RefArr4DI4.deallocate();
   TestArr4DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DI4 test: FAIL" << std::endl;

   // Test for 5DI4
   OMEGA::Array5DI4 TestArr5DI4("test",NumExtra,NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost5DI4 RefArr5DI4("ref",NumExtra,NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr5DI4(n,m,k,j,i) = i+j+k+m+n;
   }
   }
   }
   }
   }

   parallel_for(Bounds<5>(NumExtra,NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int n, int m, int k, int j,int i) {
      TestArr5DI4(n,m,k,j,i) = i+j+k+m+n;
   });

   yakl::fence();
   auto TestHost5DI4 = TestArr5DI4.createHostCopy();

   icount=0;
   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost5DI4(n,m,k,j,i) != RefArr5DI4(n,m,k,j,i)) ++icount;
   }
   }
   }
   }
   }
   TestHost5DI4.deallocate();
   RefArr5DI4.deallocate();
   TestArr5DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DI4 test: FAIL" << std::endl;

   // Test for 1DI8
   OMEGA::Array1DI8 TestArr1DI8("test",NumCells);
   OMEGA::ArrayHost1DI8 RefArr1DI8("ref",NumCells);

   for (int i=0; i < NumCells; ++i){
      RefArr1DI8(i) = i;
   }

   parallel_for(Bounds<1>(NumCells), YAKL_LAMBDA (int i) {
      TestArr1DI8(i) = i;
   });

   yakl::fence();
   auto TestHost1DI8 = TestArr1DI8.createHostCopy();

   icount=0;
   for (int i=0; i < NumCells; ++i){
      if (TestHost1DI8(i) != RefArr1DI8(i)) ++icount;
   }
   TestHost1DI8.deallocate();
   RefArr1DI8.deallocate();
   TestArr1DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DI8 test: FAIL" << std::endl;

   // Test for 2DI8
   OMEGA::Array2DI8 TestArr2DI8("test",NumCells,NumVertLvls);
   OMEGA::ArrayHost2DI8 RefArr2DI8("ref",NumCells,NumVertLvls);

   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr2DI8(j,i) = i+j;
   }
   }

   parallel_for(Bounds<2>(NumCells,NumVertLvls), YAKL_LAMBDA (int j,int i) {
      TestArr2DI8(j,i) = i+j;
   });

   yakl::fence();
   auto TestHost2DI8 = TestArr2DI8.createHostCopy();

   icount=0;
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost2DI8(j,i) != RefArr2DI8(j,i)) ++icount;
   }
   }
   TestHost2DI8.deallocate();
   RefArr2DI8.deallocate();
   TestArr2DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DI8 test: FAIL" << std::endl;

   // Test for 3DI8
   OMEGA::Array3DI8 TestArr3DI8("test",NumTracers,NumCells,NumVertLvls);
   OMEGA::ArrayHost3DI8 RefArr3DI8("ref",NumTracers,NumCells,NumVertLvls);

   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr3DI8(k,j,i) = i+j+k;
   }
   }
   }

   parallel_for(Bounds<3>(NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int k, int j,int i) {
      TestArr3DI8(k,j,i) = i+j+k;
   });

   yakl::fence();
   auto TestHost3DI8 = TestArr3DI8.createHostCopy();

   icount=0;
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost3DI8(k,j,i) != RefArr3DI8(k,j,i)) ++icount;
   }
   }
   }
   TestHost3DI8.deallocate();
   RefArr3DI8.deallocate();
   TestArr3DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DI8 test: FAIL" << std::endl;

   // Test for 4DI8
   OMEGA::Array4DI8 TestArr4DI8("test",NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost4DI8 RefArr4DI8("ref",NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr4DI8(m,k,j,i) = i+j+k+m;
   }
   }
   }
   }

   parallel_for(Bounds<4>(NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int m, int k, int j,int i) {
      TestArr4DI8(m,k,j,i) = i+j+k+m;
   });

   yakl::fence();
   auto TestHost4DI8 = TestArr4DI8.createHostCopy();

   icount=0;
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost4DI8(m,k,j,i) != RefArr4DI8(m,k,j,i)) ++icount;
   }
   }
   }
   }
   TestHost4DI8.deallocate();
   RefArr4DI8.deallocate();
   TestArr4DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DI8 test: FAIL" << std::endl;

   // Test for 5DI8
   OMEGA::Array5DI8 TestArr5DI8("test",NumExtra,NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost5DI8 RefArr5DI8("ref",NumExtra,NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr5DI8(n,m,k,j,i) = i+j+k+m+n;
   }
   }
   }
   }
   }

   parallel_for(Bounds<5>(NumExtra,NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int n, int m, int k, int j,int i) {
      TestArr5DI8(n,m,k,j,i) = i+j+k+m+n;
   });

   yakl::fence();
   auto TestHost5DI8 = TestArr5DI8.createHostCopy();

   icount=0;
   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost5DI8(n,m,k,j,i) != RefArr5DI8(n,m,k,j,i)) ++icount;
   }
   }
   }
   }
   }
   TestHost5DI8.deallocate();
   RefArr5DI8.deallocate();
   TestArr5DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DI8 test: FAIL" << std::endl;

   // Test for 1DR4
   OMEGA::Array1DR4 TestArr1DR4("test",NumCells);
   OMEGA::ArrayHost1DR4 RefArr1DR4("ref",NumCells);

   for (int i=0; i < NumCells; ++i){
      RefArr1DR4(i) = i;
   }

   parallel_for(Bounds<1>(NumCells), YAKL_LAMBDA (int i) {
      TestArr1DR4(i) = i;
   });

   yakl::fence();
   auto TestHost1DR4 = TestArr1DR4.createHostCopy();

   icount=0;
   for (int i=0; i < NumCells; ++i){
      if (TestHost1DR4(i) != RefArr1DR4(i)) ++icount;
   }
   TestHost1DR4.deallocate();
   RefArr1DR4.deallocate();
   TestArr1DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DR4 test: FAIL" << std::endl;

   // Test for 2DR4
   OMEGA::Array2DR4 TestArr2DR4("test",NumCells,NumVertLvls);
   OMEGA::ArrayHost2DR4 RefArr2DR4("ref",NumCells,NumVertLvls);

   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr2DR4(j,i) = i+j;
   }
   }

   parallel_for(Bounds<2>(NumCells,NumVertLvls), YAKL_LAMBDA (int j,int i) {
      TestArr2DR4(j,i) = i+j;
   });

   yakl::fence();
   auto TestHost2DR4 = TestArr2DR4.createHostCopy();

   icount=0;
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost2DR4(j,i) != RefArr2DR4(j,i)) ++icount;
   }
   }
   TestHost2DR4.deallocate();
   RefArr2DR4.deallocate();
   TestArr2DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DR4 test: FAIL" << std::endl;

   // Test for 3DR4
   OMEGA::Array3DR4 TestArr3DR4("test",NumTracers,NumCells,NumVertLvls);
   OMEGA::ArrayHost3DR4 RefArr3DR4("ref",NumTracers,NumCells,NumVertLvls);

   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr3DR4(k,j,i) = i+j+k;
   }
   }
   }

   parallel_for(Bounds<3>(NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int k, int j,int i) {
      TestArr3DR4(k,j,i) = i+j+k;
   });

   yakl::fence();
   auto TestHost3DR4 = TestArr3DR4.createHostCopy();

   icount=0;
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost3DR4(k,j,i) != RefArr3DR4(k,j,i)) ++icount;
   }
   }
   }
   TestHost3DR4.deallocate();
   RefArr3DR4.deallocate();
   TestArr3DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DR4 test: FAIL" << std::endl;

   // Test for 4DR4
   OMEGA::Array4DR4 TestArr4DR4("test",NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost4DR4 RefArr4DR4("ref",NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr4DR4(m,k,j,i) = i+j+k+m;
   }
   }
   }
   }

   parallel_for(Bounds<4>(NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int m, int k, int j,int i) {
      TestArr4DR4(m,k,j,i) = i+j+k+m;
   });

   yakl::fence();
   auto TestHost4DR4 = TestArr4DR4.createHostCopy();

   icount=0;
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost4DR4(m,k,j,i) != RefArr4DR4(m,k,j,i)) ++icount;
   }
   }
   }
   }
   TestHost4DR4.deallocate();
   RefArr4DR4.deallocate();
   TestArr4DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DR4 test: FAIL" << std::endl;

   // Test for 5DR4
   OMEGA::Array5DR4 TestArr5DR4("test",NumExtra,NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost5DR4 RefArr5DR4("ref",NumExtra,NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr5DR4(n,m,k,j,i) = i+j+k+m+n;
   }
   }
   }
   }
   }

   parallel_for(Bounds<5>(NumExtra,NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int n, int m, int k, int j,int i) {
      TestArr5DR4(n,m,k,j,i) = i+j+k+m+n;
   });

   yakl::fence();
   auto TestHost5DR4 = TestArr5DR4.createHostCopy();

   icount=0;
   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost5DR4(n,m,k,j,i) != RefArr5DR4(n,m,k,j,i)) ++icount;
   }
   }
   }
   }
   }
   TestHost5DR4.deallocate();
   RefArr5DR4.deallocate();
   TestArr5DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DR4 test: FAIL" << std::endl;

   // Test for 1DR8
   OMEGA::Array1DR8 TestArr1DR8("test",NumCells);
   OMEGA::ArrayHost1DR8 RefArr1DR8("ref",NumCells);

   for (int i=0; i < NumCells; ++i){
      RefArr1DR8(i) = i;
   }

   parallel_for(Bounds<1>(NumCells), YAKL_LAMBDA (int i) {
      TestArr1DR8(i) = i;
   });

   yakl::fence();
   auto TestHost1DR8 = TestArr1DR8.createHostCopy();

   icount=0;
   for (int i=0; i < NumCells; ++i){
      if (TestHost1DR8(i) != RefArr1DR8(i)) ++icount;
   }
   TestHost1DR8.deallocate();
   RefArr1DR8.deallocate();
   TestArr1DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DR8 test: FAIL" << std::endl;

   // Test for 2DR8
   OMEGA::Array2DR8 TestArr2DR8("test",NumCells,NumVertLvls);
   OMEGA::ArrayHost2DR8 RefArr2DR8("ref",NumCells,NumVertLvls);

   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr2DR8(j,i) = i+j;
   }
   }

   parallel_for(Bounds<2>(NumCells,NumVertLvls), YAKL_LAMBDA (int j,int i) {
      TestArr2DR8(j,i) = i+j;
   });

   yakl::fence();
   auto TestHost2DR8 = TestArr2DR8.createHostCopy();

   icount=0;
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost2DR8(j,i) != RefArr2DR8(j,i)) ++icount;
   }
   }
   TestHost2DR8.deallocate();
   RefArr2DR8.deallocate();
   TestArr2DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DR8 test: FAIL" << std::endl;

   // Test for 3DR8
   OMEGA::Array3DR8 TestArr3DR8("test",NumTracers,NumCells,NumVertLvls);
   OMEGA::ArrayHost3DR8 RefArr3DR8("ref",NumTracers,NumCells,NumVertLvls);

   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr3DR8(k,j,i) = i+j+k;
   }
   }
   }

   parallel_for(Bounds<3>(NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int k, int j,int i) {
      TestArr3DR8(k,j,i) = i+j+k;
   });

   yakl::fence();
   auto TestHost3DR8 = TestArr3DR8.createHostCopy();

   icount=0;
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost3DR8(k,j,i) != RefArr3DR8(k,j,i)) ++icount;
   }
   }
   }
   TestHost3DR8.deallocate();
   RefArr3DR8.deallocate();
   TestArr3DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DR8 test: FAIL" << std::endl;

   // Test for 4DR8
   OMEGA::Array4DR8 TestArr4DR8("test",NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost4DR8 RefArr4DR8("ref",NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr4DR8(m,k,j,i) = i+j+k+m;
   }
   }
   }
   }

   parallel_for(Bounds<4>(NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int m, int k, int j,int i) {
      TestArr4DR8(m,k,j,i) = i+j+k+m;
   });

   yakl::fence();
   auto TestHost4DR8 = TestArr4DR8.createHostCopy();

   icount=0;
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost4DR8(m,k,j,i) != RefArr4DR8(m,k,j,i)) ++icount;
   }
   }
   }
   }
   TestHost4DR8.deallocate();
   RefArr4DR8.deallocate();
   TestArr4DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DR8 test: FAIL" << std::endl;

   // Test for 5DR8
   OMEGA::Array5DR8 TestArr5DR8("test",NumExtra,NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost5DR8 RefArr5DR8("ref",NumExtra,NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr5DR8(n,m,k,j,i) = i+j+k+m+n;
   }
   }
   }
   }
   }

   parallel_for(Bounds<5>(NumExtra,NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int n, int m, int k, int j,int i) {
      TestArr5DR8(n,m,k,j,i) = i+j+k+m+n;
   });

   yakl::fence();
   auto TestHost5DR8 = TestArr5DR8.createHostCopy();

   icount=0;
   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost5DR8(n,m,k,j,i) != RefArr5DR8(n,m,k,j,i)) ++icount;
   }
   }
   }
   }
   }
   TestHost5DR8.deallocate();
   RefArr5DR8.deallocate();
   TestArr5DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DR8 test: FAIL" << std::endl;

   // Test for 1DReal
   OMEGA::Array1DReal TestArr1DReal("test",NumCells);
   OMEGA::ArrayHost1DReal RefArr1DReal("ref",NumCells);

   for (int i=0; i < NumCells; ++i){
      RefArr1DReal(i) = i;
   }

   parallel_for(Bounds<1>(NumCells), YAKL_LAMBDA (int i) {
      TestArr1DReal(i) = i;
   });

   yakl::fence();
   auto TestHost1DReal = TestArr1DReal.createHostCopy();

   icount=0;
   for (int i=0; i < NumCells; ++i){
      if (TestHost1DReal(i) != RefArr1DReal(i)) ++icount;
   }
   TestHost1DReal.deallocate();
   RefArr1DReal.deallocate();
   TestArr1DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DReal test: FAIL" << std::endl;

   // Test for 2DReal
   OMEGA::Array2DReal TestArr2DReal("test",NumCells,NumVertLvls);
   OMEGA::ArrayHost2DReal RefArr2DReal("ref",NumCells,NumVertLvls);

   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr2DReal(j,i) = i+j;
   }
   }

   parallel_for(Bounds<2>(NumCells,NumVertLvls), YAKL_LAMBDA (int j,int i) {
      TestArr2DReal(j,i) = i+j;
   });

   yakl::fence();
   auto TestHost2DReal = TestArr2DReal.createHostCopy();

   icount=0;
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost2DReal(j,i) != RefArr2DReal(j,i)) ++icount;
   }
   }
   TestHost2DReal.deallocate();
   RefArr2DReal.deallocate();
   TestArr2DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DReal test: FAIL" << std::endl;

   // Test for 3DReal
   OMEGA::Array3DReal TestArr3DReal("test",NumTracers,NumCells,NumVertLvls);
   OMEGA::ArrayHost3DReal RefArr3DReal("ref",NumTracers,NumCells,NumVertLvls);

   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr3DReal(k,j,i) = i+j+k;
   }
   }
   }

   parallel_for(Bounds<3>(NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int k, int j,int i) {
      TestArr3DReal(k,j,i) = i+j+k;
   });

   yakl::fence();
   auto TestHost3DReal = TestArr3DReal.createHostCopy();

   icount=0;
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost3DReal(k,j,i) != RefArr3DReal(k,j,i)) ++icount;
   }
   }
   }
   TestHost3DReal.deallocate();
   RefArr3DReal.deallocate();
   TestArr3DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DReal test: FAIL" << std::endl;

   // Test for 4DReal
   OMEGA::Array4DReal TestArr4DReal("test",NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost4DReal RefArr4DReal("ref",NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr4DReal(m,k,j,i) = i+j+k+m;
   }
   }
   }
   }

   parallel_for(Bounds<4>(NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int m, int k, int j,int i) {
      TestArr4DReal(m,k,j,i) = i+j+k+m;
   });

   yakl::fence();
   auto TestHost4DReal = TestArr4DReal.createHostCopy();

   icount=0;
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost4DReal(m,k,j,i) != RefArr4DReal(m,k,j,i)) ++icount;
   }
   }
   }
   }
   TestHost4DReal.deallocate();
   RefArr4DReal.deallocate();
   TestArr4DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DReal test: FAIL" << std::endl;

   // Test for 5DReal
   OMEGA::Array5DReal TestArr5DReal("test",NumExtra,NumTimeLvls,NumTracers,
                                       NumCells,NumVertLvls);
   OMEGA::ArrayHost5DReal RefArr5DReal("ref",NumExtra,NumTimeLvls,NumTracers,
                                         NumCells,NumVertLvls);

   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      RefArr5DReal(n,m,k,j,i) = i+j+k+m+n;
   }
   }
   }
   }
   }

   parallel_for(Bounds<5>(NumExtra,NumTimeLvls,NumTracers,NumCells,NumVertLvls), 
                YAKL_LAMBDA (int n, int m, int k, int j,int i) {
      TestArr5DReal(n,m,k,j,i) = i+j+k+m+n;
   });

   yakl::fence();
   auto TestHost5DReal = TestArr5DReal.createHostCopy();

   icount=0;
   for (int n=0; n < NumExtra; ++n){
   for (int m=0; m < NumTimeLvls; ++m){
   for (int k=0; k < NumTracers; ++k){
   for (int j=0; j < NumCells; ++j){
   for (int i=0; i < NumVertLvls; ++i){
      if (TestHost5DReal(n,m,k,j,i) != RefArr5DReal(n,m,k,j,i)) ++icount;
   }
   }
   }
   }
   }
   TestHost5DReal.deallocate();
   RefArr5DReal.deallocate();
   TestArr5DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DReal test: FAIL" << std::endl;

   // finalize environments
   //MPI_Status status;
   yakl::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
