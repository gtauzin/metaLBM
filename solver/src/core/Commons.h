#ifndef COMMONS_H
#define COMMONS_H

#ifdef __CUDACC__
#define _HOST_ __host__
#define _DEVICE_ __device__
#else
#define _HOST_
#define _DEVICE_
#endif

#define CACHE_LINE 64

namespace lbm {

  // halo updateIndex
  inline int idxL(const int iX, const int iY, const int iZ) {
    return sizeZ_h*(sizeY_h*iX + iY) + iZ;
  }

  inline int idx_inL(const int iX, const int iY, const int iZ) {
    return sizeZ_h*(sizeY_h*(iX-L::hX) + (iY-L::hY)) + (iZ-L::hZ);
  }

  inline int idx_outL(const int iX, const int iY, const int iZ) {
    return sizeZ_h*(sizeY_h*(iX+L::hX) + (iY+L::hY)) + (iZ+L::hZ);
  }

  inline int idx_inF(const int iX, const int iY, const int iZ) {
    return L::lZ_l * (L::lY_l * (iX-L::hX) + (iY-L::hY)) + (iZ-L::hZ);
  }

  // local updateIndex
  inline int idx_lF(const int iX, const int iY, const int iZ) {
    return L::lZ_l * (L::lY_l*iX + iY) + iZ;
  }

  // it means that Global must inherit from Local
  // global updateIndex OneD partitionning
  inline int idx_gF(const int iX, const int iY, const int iZ) {
    const int iX_l = iX - iX/L::lX_l * L::lX_l;
    return iX/L::lX_l * volume_l() + idx_lF(iX_l, iY, iZ);
  }

#if defined(_AOS)

  // idxL = halo updateIndex
  inline int idxPop(const int iX, const int iY, const int iZ,
                                  const int iQ) {
    return idxL(iX, iY, iZ) * L::dimQ + iQ;
  }

  // idx_gF = global updateIndex
  inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    return idx_gF(iX, iY, iZ) * L::dimQ + iQ;
  }

#elif defined(_SOA)
  // idxL = halo updateIndex
  inline int idxPop(const int iX, const int iY, const int iZ,
                                    const int iQ) {
    return iQ * sizeX_h * sizeY_h * sizeZ_h + idxL(iX, iY, iZ);
  }

  // idx_gF = global updateIndex
  inline int idxPop_gF(const int iX, const int iY, const int iZ,
                                     const int iQ) {
    const int iX_l = iX - iX/L::lX_l * L::lX_l;
    return iX/L::lX_l * L::dimQ * volume_l() + idxPop_lF(iX_l, iY, iZ, iQ);
  }

}

#endif // COMMONS_H
