# Build Targets
if(SERIAL)
  add_subdirectory(serial)
endif()

if(MPI OR MPI_OMP)
  add_subdirectory(mpi)
endif()

if(CUDA)
  add_subdirectory(cuda)
endif()

if(MPI_CUDA)
  add_subdirectory(cuda)
endif()
