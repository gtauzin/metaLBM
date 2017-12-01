#include <cstdio>

constexpr int N = 10;

__constant__ float a[N];

__global__ void kernel(float *out)
{
    if (threadIdx.x < N)
        out[threadIdx.x] = a[threadIdx.x];

}

int main()
{
    const size_t size = size_t(N) * sizeof(float);
    const float aVals[N]={ 1., 2., 3., 4., 5., 6., 7., 8., 9., 10. };
    float aHost[N];

    cudaMemcpyToSymbol("a", &avals[0], sz, size_t(0),cudaMemcpyHostToDevice);

    float *aDevice;
    cudaMalloc((void **)&aDevice, sz);

    kernel<<<dim3(1),dim3(16)>>>(aDevice);

    cudaMemcpy(&aHost[0],aDevice,size,cudaMemcpyDeviceToHost);

    for(int i = 0; i < N; i++) {
        printf("%d %f\n", i, aHost[i]);
    }
}
