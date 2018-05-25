#pragma once

#include <fstream>
#include <ostream>
#include <string>

#include "Commons.h"
#include "Field.h"
#include "Initialize.h"
#include "Options.h"
#include "Writer.h"

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
  Field<T, 2 * L::dimD - 3, architecture, writeVorticity> vorticity;
  FieldWriter_& fieldWriter;
  Curl<double, Architecture::CPU, PartitionningType::OneD, L::dimD, L::dimD>
      curlVelocity;

  FieldList(FieldWriter_& fieldWriter_in,
            const Stream<architecture>& stream_in)
    : density(initLocalDensity<T, architecture>(stream_in))
    , velocity(initLocalVelocity<T, architecture>(stream_in))
    , force(initLocalForce<T, architecture>(stream_in))
    , alpha(initLocalAlpha<T, architecture>(stream_in))
    , vorticity("vorticity")
    , curlVelocity(velocity.getLocalData(FFTWInit::numberElements),
                   vorticity.getLocalData(FFTWInit::numberElements),
                   Cast<unsigned int, ptrdiff_t, 3>::Do(gSD::sLength()).data(),
                   gFD::offset(MPIInit::rank))
    , fieldWriter(fieldWriter_in)
  {}

  inline void writeFields() {
    fieldWriter.writeField(density);
    fieldWriter.writeField(velocity);
    fieldWriter.writeField(alpha);
    fieldWriter.writeField(force);
    fieldWriter.writeField(vorticity);
  }
};

}  // namespace lbm
