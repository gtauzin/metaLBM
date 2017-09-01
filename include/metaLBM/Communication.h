#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include <string>

#include "Commons.h"
#include "Options.h"
#include "DynamicArray.cuh"
#include "Lattice.h"
#include "Domain.h"
#include "Boundary.h"
#include "Computation.cuh"

namespace lbm {

  template<class T, LatticeType latticeType, AlgorithmType algorithmType,
           MemoryLayout memoryLayout,
           PartitionningType partitionningType, unsigned int Dimension>
  class Communication {};


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic, PartitionningType::Generic, 0>
    : public Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull> {
  protected:
    const MathVector<int, 3> rankMPI;
    const MathVector<int, 3> sizeMPI;
    const std::string processorName;

    T * RESTRICT haloComputedData;

    const unsigned int rightXRankMPI;
    const unsigned int leftXRankMPI;
    MPI_Status statusMPI[4];
    MPI_Request requestMPI[4];

    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>()
      , rankMPI(rankMPI_in)
      , sizeMPI(sizeMPI_in)
      , processorName(processorName_in)
      , rightXRankMPI((rankMPI_in[d::X] + sizeMPI_in[d::X] - 1) % sizeMPI_in[d::X])
      , leftXRankMPI((rankMPI_in[d::X] + 1) % sizeMPI_in[d::X])
      , statusMPI()
      , requestMPI()
      , haloComputedData(haloComputedData_in)
    {}

    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyY;
    using Boundary<T, BoundaryType::Periodic, AlgorithmType::Pull>::applyZ;

    void printInputs() {
      std::cout << "MPI #" << rankMPI[d::X] << " of " << sizeMPI[d::X]
                << " running on host " << processorName << std::endl
                << "Right MPI #" << rightXRankMPI
                << ", left MPI #" << leftXRankMPI << std::endl;
    }

    DEVICE HOST
    MathVector<int, 3> getRankMPI() {
      return rankMPI;
    }

    DEVICE HOST
    void sendGlobalToLocal(DynamicArray<T, Architecture::CPU>& globalArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           unsigned int numberComponents) {
      SCOREP_INSTRUMENT_ON("Communication<6>::sendGlobalToLocal")

      MPI_Scatter(globalArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                  localHostArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                  0, MPI_COMM_WORLD);

      localDeviceArray.copyFrom(localHostArray);
    }

    DEVICE HOST
    void sendLocalToGlobal(DynamicArray<T, Architecture::GPU>& localDeviceArray,
                           DynamicArray<T, Architecture::CPU>& localHostArray,
                           DynamicArray<T, Architecture::CPU>& globalArray,
                           unsigned int numberComponents) {
      SCOREP_INSTRUMENT_ON("Communication<6>::sendLocalToGlobal")

      localDeviceArray.copyTo(localHostArray);

      MPI_Gather(localHostArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                 globalArray.data(), numberComponents*lD::volume(), MPI_DOUBLE,
                 0, MPI_COMM_WORLD);
    }

    DEVICE HOST
    T reduce(const DynamicArray<T, Architecture::GPU>& localDeviceArray,
             DynamicArray<T, Architecture::CPU>& localHostArray) {

      localDeviceArray.copyTo(localHostArray);

      T localSum = (T) 0;
      for(unsigned int iZ = lD::start()[d::Z]; iZ < lD::end()[d::Z]; ++iZ) {
        for(unsigned int iY = lD::start()[d::Y]; iY < lD::end()[d::Y]; ++iY) {
          for(unsigned int iX = lD::start()[d::X]; iX < lD::end()[d::X]; ++iX) {
            localSum
              += localHostArray[lD::getIndex(MathVector<unsigned int, 3>({iX, iY, iZ}))];
          }
        }
      }


      MPI_Barrier(MPI_COMM_WORLD);

      T globalSum = (T) 0;
      MPI_Reduce(&localSum, &globalSum, 1, MPI_DOUBLE,
                 MPI_SUM, 0, MPI_COMM_WORLD);

      MPI_Barrier(MPI_COMM_WORLD);

      return globalSum;
    }

    DEVICE HOST
    void setHaloComputedData(T * RESTRICT haloComputedData_in) {
      haloComputedData = haloComputedData_in;
    }
  };


  template<class T, LatticeType latticeType>
    class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::SoA, PartitionningType::OneD, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic, 0> {
  protected:
    DEVICE HOST
    void sendAndReceiveHaloX(DynamicArray<T, Architecture::GPU>& haloDeviceArray,
                             DynamicArray<T, Architecture::CPU>& haloHostArray) {
      SCOREP_INSTRUMENT_ON("Communication<5, MemoryLayout::SoA>::sendAndReceiveHaloX")

      haloDeviceArray.copyTo(haloHostArray);

      UnrolledFor<0, L::dimQ>::Do([&] HOST DEVICE (unsigned int iQ) {
          sendToRightBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X]-1,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);
          receivedFromLeftBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({0,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);

          sendToLeftBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);

          receivedFromRightBeginX = hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X],
                  hMLD::start()[d::Y], hMLD::start()[d::Z]}), iQ);


          MPI_Irecv(haloHostArray.data()+receivedFromLeftBeginX, sizeStripeX,
                    MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestMPI[0]);

          MPI_Irecv(haloHostArray.data()+receivedFromRightBeginX, sizeStripeX,
                    MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestMPI[1]);

          MPI_Isend(haloHostArray.data()+sendToRightBeginX, sizeStripeX,
                    MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestMPI[2]);

          MPI_Isend(haloHostArray.data()+sendToLeftBeginX, sizeStripeX,
                    MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestMPI[3]);

          MPI_Waitall(4, requestMPI, statusMPI);

        });

        haloDeviceArray.copyFrom(haloHostArray);

    }

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic,
                        0>::haloComputedData;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic, 0>::applyZ;

  public:
    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic,
                      PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                     processorName_in, haloComputedData_in)
      , sizeStripeX(hMLD::volume()/hMLD::length()[d::X])
      , sendToRightBeginX(0)
      , receivedFromLeftBeginX(0)
      , sendToLeftBeginX(0)
      , receivedFromRightBeginX(0)

    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::setHaloComputedData;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestMPI;

    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::SoA, L::dimQ> hMLD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;
  };


  template<class T, LatticeType latticeType>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::AoS, PartitionningType::OneD, 0>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           MemoryLayout::Generic, PartitionningType::Generic, 0> {
  protected:
    DEVICE HOST
    void sendAndReceiveHaloX(DynamicArray<T, Architecture::GPU>& haloDeviceArray,
                             DynamicArray<T, Architecture::CPU>& haloHostArray) {
      SCOREP_INSTRUMENT_ON("Communication<5, MemoryLayout::AoS>::sendAndReceiveHaloX")

      haloDeviceArray.copyTo(haloHostArray);

      MPI_Irecv(haloHostArray.data()+receivedFromLeftBeginX, sizeStripeX,
                MPI_DOUBLE, leftXRankMPI, 17, MPI_COMM_WORLD, &requestMPI[0]);

      MPI_Irecv(haloHostArray.data()+receivedFromRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 23, MPI_COMM_WORLD, &requestMPI[1]);

      MPI_Isend(haloHostArray.data()+sendToRightBeginX, sizeStripeX,
                MPI_DOUBLE, rightXRankMPI, 17, MPI_COMM_WORLD, &requestMPI[2]);

      MPI_Isend(haloHostArray.data()+sendToLeftBeginX, sizeStripeX,
                MPI_DOUBLE, leftXRankMPI, 23, MPI_COMM_WORLD, &requestMPI[3]);

      MPI_Waitall(4, requestMPI, statusMPI);

      haloDeviceArray.copyFrom(haloHostArray);

    }
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic,
                        0>::haloComputedData;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic, 0>::applyY;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic, PartitionningType::Generic, 0>::applyZ;

  public:
    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      MemoryLayout::Generic,
                      PartitionningType::Generic, 0>(rankMPI_in, sizeMPI_in,
                                                     processorName_in, haloComputedData_in)
      , sizeStripeX(L::dimQ*hMLD::volume()*L::halo()[d::X]/hMLD::length()[d::X])
      , sendToRightBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X]-1,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , receivedFromLeftBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({0,
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , sendToLeftBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
      , receivedFromRightBeginX(hMLD::getIndex(MathVector<unsigned int, 3>({L::halo()[d::X]+lD::length()[d::X],
                hMLD::start()[d::Y], hMLD::start()[d::Z]}), 0))
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::setHaloComputedData;

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::rightXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::leftXRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::statusMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        MemoryLayout::Generic,
                        PartitionningType::Generic, 0>::requestMPI;

    typedef Domain<DomainType::Halo, PartitionningType::Generic,
                   MemoryLayout::AoS, L::dimQ> hMLD;

    unsigned int sizeStripeX;
    unsigned int sendToRightBeginX;
    unsigned int receivedFromLeftBeginX;
    unsigned int sendToLeftBeginX;
    unsigned int receivedFromRightBeginX;
  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 1>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::OneD, 0> {
  public:
    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 0>(rankMPI_in, sizeMPI_in,
                                                                processorName_in,
                                                                haloComputedData_in)
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::setHaloComputedData;

    DEVICE HOST
    inline void periodic(DynamicArray<T, Architecture::GPU>& haloDeviceArray,
                         DynamicArray<T, Architecture::CPU>& haloHostArray) {
      SCOREP_INSTRUMENT_ON("Communication<6>::periodic")

        sendAndReceiveHaloX(haloDeviceArray, haloHostArray);
    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 0>::sendAndReceiveHaloX;
  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 2>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::OneD, 1> {
  public:
    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 1>(rankMPI_in, sizeMPI_in,
                                                                processorName_in,
                                                                haloComputedData_in)
      , endY(MathVector<unsigned int, 3>({hD::end()[d::X],
              hD::start()[d::Y]+1, hD::end()[d::Z]}))
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::setHaloComputedData;

    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP) {
      applyY(haloComputedData, iP);
    }

    DEVICE HOST
    inline void periodic(DynamicArray<T, Architecture::GPU>& haloDeviceArray,
                         DynamicArray<T, Architecture::CPU>& haloHostArray) {
      SCOREP_INSTRUMENT_ON("Communication<6>::periodic")

        //record event (compute_done) in default stream (0)
        //launch local periodic boundary kernel in default stream (0)
      Computation_::Do(hD::start(), endY, *this);

      //wait for compute to be done -> synchronize event compute_done

      Communication<T, latticeType, AlgorithmType::Pull,
                    memoryLayout, PartitionningType::OneD, 1>::periodic(haloDeviceArray,
                                                                        haloHostArray);

    }

  protected:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::haloComputedData;

  private:
    const MathVector<unsigned int, 3> endY;

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 1>::applyY;

  };


  template<class T, LatticeType latticeType,
           MemoryLayout memoryLayout>
  class Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 3>
    : public Communication<T, latticeType, AlgorithmType::Pull,
                           memoryLayout, PartitionningType::OneD, 2> {
  public:
    DEVICE HOST
    Communication(const MathVector<int, 3>& rankMPI_in,
                  const MathVector<int, 3>& sizeMPI_in,
                  const std::string& processorName_in,
                  T * RESTRICT haloComputedData_in)
      : Communication<T, latticeType, AlgorithmType::Pull,
                      memoryLayout, PartitionningType::OneD, 2>(rankMPI_in, sizeMPI_in,
                                                                processorName_in,
                                                                haloComputedData_in)
      , endZ(MathVector<unsigned int, 3>({hD::end()[d::X],
              hD::end()[d::Y], hD::start()[d::Z]+1}))
    {}

    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::printInputs;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::getRankMPI;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::sendGlobalToLocal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::sendLocalToGlobal;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::reduce;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::setHaloComputedData;

    DEVICE HOST
    inline void operator()(const MathVector<unsigned int, 3>& iP) {
      applyZ(haloComputedData, iP);
    }

    DEVICE HOST
    inline void periodic(DynamicArray<T, Architecture::GPU>& haloDeviceArray,
                         DynamicArray<T, Architecture::CPU>& haloHostArray) {
      SCOREP_INSTRUMENT_ON("Communication<6>::periodic")

      Computation_::Do(hD::start(), endZ, *this);

      Communication<T, latticeType, AlgorithmType::Pull,
                    memoryLayout, PartitionningType::OneD, 2>::periodic(haloDeviceArray,
                                                                        haloHostArray);

    }

  private:
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::haloComputedData;
    using Communication<T, latticeType, AlgorithmType::Pull,
                        memoryLayout, PartitionningType::OneD, 2>::applyZ;
    const MathVector<unsigned int, 3> endZ;
  };


  typedef Communication<dataT, latticeT, algorithmT,
                        memoryL, partitionningT, L::dimD> Communication_;
}

#endif // COMMUNICATION_H
