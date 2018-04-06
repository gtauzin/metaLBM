#pragma once

#include <string>
#include <fstream>
#include <ostream>

#include "Commons.h"
#include "Options.h"
#include "Writer.h"
#include "Field.h"
#include "Initialize.h"

namespace lbm {

   /**
   * FieldList containing all fields.
   *
   * @tparam T datatype.
   * @tparam Architecture on which the code is executed.
   */

  template <class T, Architecture architecture>
  class FieldList {
  public:
    Field<T, 1, architecture, true> density;
    Field<T, L::dimD, architecture, true> velocity;
    Field<T, L::dimD, architecture, writeForce> force;
    Field<T, 1, architecture, writeAlpha> alpha;
    Field<T, 2*L::dimD-3, architecture, writeVorticity> vorticity;
    FieldWriter_& fieldWriter;
    Curl<double, Architecture::CPU, PartitionningType::OneD,
         L::dimD, L::dimD> curlVelocity;

    FieldList(const MathVector<int, 3>& rankMPI_in, const unsigned int numberElements_in,
              FieldWriter_& fieldWriter_in)
      : density(initLocalDensity<T, architecture>(numberElements_in, rankMPI_in))
      , velocity(initLocalVelocity<T, architecture>(numberElements_in))
      , force(initLocalForce<T, architecture>(numberElements_in, rankMPI_in))
      , alpha(initLocalAlpha<T, architecture>(numberElements_in))
      , vorticity("vorticity", numberElements_in, 0)
      , curlVelocity(velocity.getLocalData(), vorticity.getLocalData(),
                     numberElements_in, Cast<unsigned int,
                     ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                     gFD::offset(rankMPI_in))
      , fieldWriter(fieldWriter_in)
    {}

    inline void writeFields() {
      if(writeVorticity) curlVelocity.executeSpace();

      fieldWriter.writeField(density);
      fieldWriter.writeField(velocity);
      fieldWriter.writeField(alpha);
      fieldWriter.writeField(force);
      fieldWriter.writeField(vorticity);
    }

  };


} // namespace lbm
