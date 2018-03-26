#ifndef FIELDLIST_H
#define FIELDLIST_H

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

    FieldList(const MathVector<int, 3>& rankMPI_in, FieldWriter_& fieldWriter_in)
      : density(initLocalDensity<T, architecture>(rankMPI_in))
      , velocity(initLocalVelocity<T, architecture>())
      , force(initLocalForce<T, architecture>(rankMPI_in))
      , alpha(initLocalAlpha<T, architecture>())
      , vorticity("vorticity", 0)
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


} // namespace lbm

#endif // FIELDLIST_H
