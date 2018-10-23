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
  Field<T, L::dimD, architecture, true> strainDiagonal;
  Field<T, 2 * L::dimD - 3, architecture, true> strainSymmetric;
  Field<T, L::dimD, architecture, writeForce> force;
  Field<T, 1, architecture, writeAlpha> alpha;
  Field<T, 1, architecture, writeAlpha> numberIterations;
  Field<T, 1, architecture, writeKinetics> T2;
  Field<T, 1, architecture, writeKinetics> T3;
  Field<T, 1, architecture, writeKinetics> T4;
  Field<T, 1, architecture, writeKinetics> T2_approx;
  Field<T, 1, architecture, writeKinetics> T3_approx;
  Field<T, 1, architecture, writeKinetics> T4_approx;
  Field<T, 3, architecture, writeKinetics> fNeq_5_6_8;
  Field<T, 3, architecture, writeKinetics> f1_5_6_8;
  Field<T, L::dimD, architecture, writeKinetics> piNeqDiagonal;
  Field<T, 2 * L::dimD - 3, architecture, writeKinetics> piNeqSymmetric;
  Field<T, L::dimD, architecture, writeKinetics> pi1Diagonal;
  Field<T, 2 * L::dimD - 3, architecture, writeKinetics> pi1Symmetric;
  Field<T, 1, architecture, writeKinetics> squaredQContractedPiNeq;
  Field<T, 1, architecture, writeKinetics> cubedQContractedPiNeq;
  Field<T, 1, architecture, writeKinetics> squaredQContractedPi1;
  Field<T, 1, architecture, writeKinetics> cubedQContractedPi1;
  Field<T, 2 * L::dimD - 3, architecture, writeVorticity> vorticity;
  FieldWriter_& fieldWriter;


  FieldList(FieldWriter_& fieldWriter_in,
            const Stream<architecture>& stream_in)
    : density(initDensity<T, architecture>(stream_in))
    , velocity(initVelocity<T, architecture>(stream_in))
    , strainDiagonal("strainDiagonal")
    , strainSymmetric("strainSymmetric")
    , force(initForce<T, architecture>(stream_in))
    , alpha(initAlpha<T, architecture>(stream_in))
    , numberIterations("numberIterations")
    , T2("T2")
    , T3("T3")
    , T4("T4")
    , T2_approx("T2_approx")
    , T3_approx("T3_approx")
    , T4_approx("T4_approx")
    , fNeq_5_6_8("fNeq_5_6_8")
    , f1_5_6_8("f1_5_6_8")
    , piNeqDiagonal("piNeqDiagonal")
    , piNeqSymmetric("piNeqSymmetric")
    , pi1Diagonal("pi1Diagonal")
    , pi1Symmetric("pi1Symmetric")
    , squaredQContractedPiNeq("squaredQContractedPiNeq")
    , cubedQContractedPiNeq("cubedQContractedPiNeq")
    , squaredQContractedPi1("squaredQContractedPi1")
    , cubedQContractedPi1("cubedQContractedPi1")
    , vorticity("vorticity")
    , fieldWriter(fieldWriter_in)
  {}

  inline void writeFields() {
    fieldWriter.writeField(density);
    fieldWriter.writeField(velocity);
    fieldWriter.writeField(strainDiagonal);
    fieldWriter.writeField(strainSymmetric);

    if(writeAlpha) {
      fieldWriter.writeField(alpha);
      fieldWriter.writeField(numberIterations);
    }
    if(writeKinetics) {
      fieldWriter.writeField(T2);
      fieldWriter.writeField(T3);
      fieldWriter.writeField(T4);

      fieldWriter.writeField(T2_approx);
      fieldWriter.writeField(T3_approx);
      fieldWriter.writeField(T4_approx);

      fieldWriter.writeField(fNeq_5_6_8);
      fieldWriter.writeField(f1_5_6_8);

      fieldWriter.writeField(piNeqDiagonal);
      fieldWriter.writeField(piNeqSymmetric);
      fieldWriter.writeField(pi1Diagonal);
      fieldWriter.writeField(pi1Symmetric);

      fieldWriter.writeField(squaredQContractedPiNeq);
      fieldWriter.writeField(cubedQContractedPiNeq);
      fieldWriter.writeField(squaredQContractedPi1);
      fieldWriter.writeField(cubedQContractedPi1);
    }
    if(writeForce) fieldWriter.writeField(force);
    if(writeVorticity) fieldWriter.writeField(vorticity);
  }
};

}  // namespace lbm
