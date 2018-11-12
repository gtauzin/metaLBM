#pragma once

#include <mpi.h>

#include "MPIInitializer.h"

namespace lbm {

/// Constructs a one to one map between process and GPU devices on a
/// multi-GPU node
struct CUDAInitializer {
  CUDAInitializer() {
    MPI_Comm localComm;
    MPI_Info info;

    LBM_MPI_CALL(MPI_Info_create(&info));
    LBM_MPI_CALL(MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED,
                                     MPIInit::rank[d::X], info,
                                     &localComm));

    int localRank;
    LBM_MPI_CALL(MPI_Comm_rank(localComm, &localRank));

    int numberDevices;
    LBM_CUDA_CALL(cudaGetDeviceCount(&numberDevices));
    LBM_CUDA_CALL(cudaSetDevice(localRank % numberDevices));

    LBM_MPI_CALL(MPI_Comm_free(&localComm));
    LBM_MPI_CALL(MPI_Info_free(&info));

  }

  /// Removes all device allocations
  ~CUDAInitializer() {
    LBM_CUDA_CALL(cudaDeviceReset());
  }
};

}  // end namespace lbm
