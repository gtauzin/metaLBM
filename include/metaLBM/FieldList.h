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
  Field<T, 1, architecture, writeT> T2;
  Field<T, 1, architecture, writeT> T3;
  Field<T, 1, architecture, writeT> T4;
  Field<T, 1, architecture, writeT> T2_approx;
  Field<T, 1, architecture, writeT> T3_approx;
  Field<T, 1, architecture, writeT> T4_approx;
  Field<T, 2 * L::dimD - 3, architecture, writeVorticity> vorticity;
  FieldWriter_& fieldWriter;


  FieldList(FieldWriter_& fieldWriter_in,
            const Stream<architecture>& stream_in)
    : density(initDensity<T, architecture>(stream_in))
    , velocity(initVelocity<T, architecture>(stream_in))
    , force(initForce<T, architecture>(stream_in))
    , alpha(initAlpha<T, architecture>(stream_in))
    , T2("T2")
    , T3("T3")
    , T4("T4")
    , T2_approx("T2_approx")
    , T3_approx("T3_approx")
    , T4_approx("T4_approx")
    , vorticity("vorticity")
    , fieldWriter(fieldWriter_in)
  {}

  inline void writeFields() {
    fieldWriter.writeField(density);
    fieldWriter.writeField(velocity);
    if(writeAlpha) fieldWriter.writeField(alpha);
    if(writeT) {
      fieldWriter.writeField(T2);
      fieldWriter.writeField(T3);
      fieldWriter.writeField(T4);

      fieldWriter.writeField(T2_approx);
      fieldWriter.writeField(T3_approx);
      fieldWriter.writeField(T4_approx);
    }
    if(writeForce) fieldWriter.writeField(force);
    if(writeVorticity) fieldWriter.writeField(vorticity);
  }
};

}  // namespace lbm
