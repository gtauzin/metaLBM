#pragma once

#include <string>

#include "Commons.h"
#include "Domain.h"
#include "Field.h"
#include "MathVector.h"
#include "Options.h"
#include "Writer.h"

namespace lbm {

template <class T, InputOutput inputOutput>
class FieldReader {};

template <class T>
class FieldReader<T, InputOutput::HDF5>
    : public FieldWriter<T, InputOutput::HDF5> {
 private:
  using Base = FieldWriter<T, InputOutput::HDF5>;

 public:
  FieldReader(const std::string& filePrefix_in,
              const MathVector<int, 3>& rankMPI_in,
              const std::string& name_in = "field")
      : Base(filePrefix_in, rankMPI_in, name_in) {}

  inline void openFile(const unsigned int iteration) {
    std::string fileName = Base::getFileName(iteration);
    open(fileName);
  }

  inline void closeFile() { Base::statusHDF5 = H5Fclose(Base::fileHDF5); }

  template <unsigned int NumberComponents, Architecture architecture>
  void readField(Field<T, NumberComponents, architecture, true>& field) {
    {LBM_INSTRUMENT_ON(
        "Reader<T, InputOutput::HDF5, "
        "readerFileFromat>::readField<NumberComponents>",
        3)}

    std::string fieldName = field.fieldName;
    Base::propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);

    for (auto iC = 0; iC < NumberComponents; ++iC) {
      if (NumberComponents > 1) {
        fieldName = field.fieldName + dName[iC];
      }
      Base::fileSpaceHDF5 = H5Screate_simple(
          L::dimD,
          Project<hsize_t, unsigned int, L::dimD>::Do(gSD::pLength()).data(),
          NULL);

      Base::dataSetHDF5 = H5Dcreate2(Base::fileHDF5, (fieldName).c_str(),
                                     H5T_NATIVE_DOUBLE, Base::fileSpaceHDF5,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5);

      Base::dataSpaceHDF5 = H5Screate_simple(
          L::dimD,
          Project<hsize_t, unsigned int, L::dimD>::Do(lSD::pLength()).data(),
          NULL);

      Base::fileSpaceHDF5 = H5Dget_space(Base::dataSetHDF5);

      H5Sselect_hyperslab(
          Base::fileSpaceHDF5, H5S_SELECT_SET,
          Project<hsize_t, unsigned int, L::dimD>::Do(
              gSD::pOffset(Base::rankMPI))
              .data(),
          NULL,
          Project<hsize_t, unsigned int, L::dimD>::Do(lSD::pLength()).data(),
          NULL);

      H5Pset_dxpl_mpio(Base::propertyListHDF5, H5FD_MPIO_COLLECTIVE);

      Base::statusHDF5 = H5Dread(
          Base::dataSetHDF5, H5T_NATIVE_DOUBLE, Base::dataSpaceHDF5,
          Base::fileSpaceHDF5, Base::propertyListHDF5, field.getLocalData(iC));

      Base::statusHDF5 = H5Dclose(Base::dataSetHDF5);
      Base::statusHDF5 = H5Sclose(Base::dataSpaceHDF5);
      Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5);
    }

    Base::statusHDF5 = H5Pclose(Base::propertyListHDF5);
  }

  template <unsigned int NumberComponents, Architecture architecture>
  void readField(Field<T, NumberComponents, architecture, false>& field) {}

  inline void open(const std::string& fileName) {
    Base::propertyListHDF5 = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(Base::propertyListHDF5, MPI_COMM_WORLD, MPI_INFO_NULL);

    Base::fileHDF5 =
        H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, Base::propertyListHDF5);
    H5Pclose(Base::propertyListHDF5);

    if (!Base::fileHDF5) {
      std::cout << "Could not open file " << fileName << std::endl;
    }
  }
};

template <class T, InputOutput inputOutput>
class DistributionReader {};

template <class T>
class DistributionReader<T, InputOutput::HDF5>
    : public FieldReader<T, InputOutput::HDF5> {
 private:
  using Base = FieldReader<T, InputOutput::HDF5>;

 public:
  DistributionReader(const std::string& filePrefix_in,
                     const MathVector<int, 3>& rankMPI_in)
      : Base(filePrefix_in, rankMPI_in, "distribution") {}

  template <Architecture architecture>
  void readDistribution(Distribution<T, architecture>& distribution) {
    {LBM_INSTRUMENT_ON(
        "Reader<T, InputOutput::HDF5, "
        "readerFileFromat>::readField<NumberComponents>",
        3)}

    Base::propertyListHDF5 = H5Pcreate(H5P_DATASET_XFER);

    for (auto iC = 0; iC < L::dimQ; ++iC) {
      /* Base::fileSpaceHDF5 = H5Screate_simple(L::dimD, */
      /*                                        Project<hsize_t, */
      /*                                        unsigned int,
       * L::dimD>::Do(gSD::pLength()).data(), */
      /*                                        NULL); */

      Base::dataSetHDF5 =
          H5Dopen1(Base::fileHDF5,
                   (distribution.fieldName + std::to_string(iC)).c_str());

      /* Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5); */

      Base::dataSpaceHDF5 = H5Screate_simple(
          L::dimD,
          Project<hsize_t, unsigned int, L::dimD>::Do(lSD::pLength()).data(),
          NULL);

      Base::fileSpaceHDF5 = H5Dget_space(Base::dataSetHDF5);

      H5Sselect_hyperslab(
          Base::dataSpaceHDF5, H5S_SELECT_SET,
          Project<hsize_t, unsigned int, L::dimD>::Do(
              gSD::pOffset(Base::rankMPI))
              .data(),
          NULL,
          Project<hsize_t, unsigned int, L::dimD>::Do(lSD::pLength()).data(),
          NULL);

      H5Pset_dxpl_mpio(Base::propertyListHDF5, H5FD_MPIO_COLLECTIVE);

      Base::statusHDF5 =
          H5Dread(Base::dataSetHDF5, H5T_NATIVE_DOUBLE, Base::dataSpaceHDF5,
                  Base::fileSpaceHDF5, Base::propertyListHDF5,
                  distribution.getLocalData(iC));

      Base::statusHDF5 = H5Dclose(Base::dataSetHDF5);
      Base::statusHDF5 = H5Sclose(Base::dataSpaceHDF5);
      Base::statusHDF5 = H5Sclose(Base::fileSpaceHDF5);
    }

    Base::statusHDF5 = H5Pclose(Base::propertyListHDF5);
  }

  using Base::closeFile;
  using Base::openFile;
};

typedef FieldReader<dataT, InputOutput::HDF5> FieldReader_;
typedef DistributionReader<dataT, InputOutput::HDF5> DistributionReader_;

}  // namespace lbm
