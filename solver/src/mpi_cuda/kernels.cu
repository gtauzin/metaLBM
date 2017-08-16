extern "C" {
#include "kernels.h"
#include "utils_gpu.h"
}

#include "commons.h"

namespace lbm {

   __constant__ double _celerityX[9] = {0.0,     -1.0,  -1.0,   -1.0,   0.0,    1.0,   1.0,    1.0,   0.0};
   __constant__ double _celerityY[9] = {0.0,      1.0,   0.0,   -1.0,  -1.0,   -1.0,   0.0,    1.0,   1.0};


  __global__ void propagateDevice(double * nxt,
                                  double * prv,
                                  const  int startX,
                                  const  int   endX,
                                  const  int startY,
                                  const  int   endY) {
    int iX = threadIdx.x+blockIdx.x*blockDim.x;
    int iY = threadIdx.y+blockIdx.y*blockDim.y;

    if(iX >= startX && iX < endX && iY >= startY && iY < endY) {
      for(int i = 0; i < 9; ++i){
        int idx_lattice = idxL(iX, iY);

        int iX_next = iX + _celerityX[i];
        int iY_next = iY + _celerityY[i];

        if(iX == haloX && _celerityX[i] == -1) iX_next = haloX+lengthX_g-1;
        if(iX == haloX+lengthX_g-1 && _celerityX[i] == 1) iX_next = haloX;
        if(iY == haloY && _celerityY[i] == -1) iY_next = haloY+lengthY_g-1;
        if(iY == haloY+lengthY_g-1 && _celerityY[i] == 1) iY_next = haloY;
        int idx_lattice_next = idxL(iX_next, iY_next);

        nxt[idxPop(idx_lattice_next, i)] = prv[idxPop(idx_lattice, i)];
      }

    }

    //  __syncthreads();
  }

  extern "C"
  void dummyPropagateDevice(double * nxt,
                            double * prv,
                            const  int startX,
                            const  int   endX,
                            const  int startY,
                            const  int   endY) {

    //  fprintf(stdout, "%d, %d, %d\n", size, size/NX, size/NX/NY);
    int sizeX = startX+endX;
    int sizeY = startY+endY;
    int size = 9*sizeX*sizeY;

    double * d_nxt;
    double * d_prv;
    cudaErrorCheck(cudaMalloc((void **) &d_nxt, size));
    cudaErrorCheck(cudaMalloc((void **) &d_prv, size));

    cudaErrorCheck(cudaMemcpy(d_prv, prv, size, cudaMemcpyHostToDevice));
    cudaErrorCheck(cudaMemcpy(d_nxt, nxt, size, cudaMemcpyHostToDevice));

    dim3 block(BLOCK_SIZE, BLOCK_SIZE, 1);
    dim3 grid((sizeX+block.x-1)/block.x, (sizeY+block.y-1)/block.y, 1);

    propagateDevice<<<grid, block>>>(d_nxt, d_prv, startX, endX, startY, endY);

    cudaErrorCheck(cudaMemcpy(prv, d_prv, size, cudaMemcpyDeviceToHost));
    cudaErrorCheck(cudaMemcpy(nxt, d_nxt, size, cudaMemcpyDeviceToHost));

    cudaErrorCheck(cudaFree(d_nxt));
    cudaErrorCheck(cudaFree(d_prv));
  }

}
