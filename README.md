To compile:

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE:STRING=Release -DSERIAL=ON -DNPROCS:INT=1 -DNTHREADS:INT=1 -DDATA_STRUCT:STRING=-D_SOA -DLOG:STRING= .. &&
make lbm