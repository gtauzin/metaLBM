#include <mpi.h>
#include <iostream>
#include <string>

#include "metaLBM/Commons.h"
#include "metaLBM/Computation.cuh"
#include "metaLBM/DynamicArray.cuh"
#include "Functor.h"

using namespace lbm;


int main (int argc, char** argv) {
  constexpr int N = 10;
  MPI_Init(&argc, &argv);

  cudaSetDevice ( 0 );

  DynamicArray<double, Architecture::CPU> hArray(N, 2.5);


  for(int i = 0; i < N; ++i) std::cout << "Before copy: " << hArray[i] << std::endl;

  DynamicArray<double, Architecture::GPU> dArray(N, 3.5);

  //dArray.copyTo(hArray);

  for(int i = 0; i < N; ++i) std::cout << "After copy: " << hArray[i] << std::endl;

  Functor testFunctor(dArray, N);

  testFunctor.compute();

  // dim3 dimBlock(128, 1, 1);
  // dim3 dimGrid(1, 1, 1);
  // kernel_1D<<<dimBlock, dimGrid >>>(MathVector<unsigned int, 3>{{0}},
  //                                   MathVector<unsigned int, 3>{{10}},
  //                                   [&] DEVICE (const MathVector<unsigned int, 3>& iP){
  //                                     dArray[iP[0]] +=2;
  //                                   });

  CUDA_CALL( cudaDeviceSynchronize(); )


  dArray.copyTo(hArray);

  for(int i = 0; i < N; ++i) std::cout << "After kernel: " << hArray[i] << std::endl;


  MPI_Finalize();

  return EXIT_SUCCESS;
}
