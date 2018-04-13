#pragma once

#include <mpi.h>

namespace lbm {

/// Constructs a one to one map between process and GPU devices on a
/// multi-GPU node
struct CUDAInitializer {
    CUDAInitializer () {
        MPI_Comm localComm;
        MPI_Info info;

        int globalRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &globalRank);

        MPI_Info_create(&info);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalRank,
                             info, &localComm);

        int localRank;
        MPI_Comm_rank(localComm, &localRank);

        int numberDevices;
        LBM_CUDA_CALL( cudaGetDeviceCount(&numberDevices); )
        LBM_CUDA_CALL( cudaSetDevice(localRank % numberDevices); )

        MPI_Comm_free(&localComm);
        MPI_Info_free(&info);
    }

    /// Removes all device allocations
    ~CUDAInitializer() {
        cudaDeviceReset();
    }
}; // end struct CUDAInitializer

} // end namespace lbm
