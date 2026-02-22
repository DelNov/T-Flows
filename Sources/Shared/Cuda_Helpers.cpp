#include <cstdio>
#if T_FLOWS_GPU == 1
#include <cuda_runtime.h>
#endif

/*-----------------------------------------------------------------------------+
|                                                                              |
|  GPU is defined                                                              |
|                                                                              |
+-----------------------------------------------------------------------------*/
#if T_FLOWS_GPU == 1

/*============================================================================*/
extern "C" void cuda_alloc_double_(void **     ptr,
                                   const int & N) {
/*-----------------------------------------------------------------------------+
|   Allocates device memory for N doubles                                      |
+-----------------------------------------------------------------------------*/

  cudaError_t err = cudaMalloc(ptr, N * sizeof(double));
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMalloc failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*============================================================================*/
extern "C" void cuda_alloc_copyin_double_(void **        ptr,
                                          const double * val_host,
                                          const int    & N) {
/*-----------------------------------------------------------------------------+
|   Allocates device memory and copies doubles from host to device             |
+-----------------------------------------------------------------------------*/

  cudaError_t err = cudaMalloc(ptr, N * sizeof(double));
  if (err != cudaSuccess) {
    fprintf(stderr, "cudaMalloc (double) failed: %s\n",
            cudaGetErrorString(err));
    return;
  }

  err = cudaMemcpy(*ptr, val_host, N * sizeof(double), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy H->D (double) failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*============================================================================*/
extern "C" void cuda_alloc_int_(void **     ptr,
                                const int & N) {
/*-----------------------------------------------------------------------------+
|   Allocates device memory for N integers                                     |
+-----------------------------------------------------------------------------*/

  cudaError_t err = cudaMalloc(ptr, N * sizeof(int));
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaMalloc failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*============================================================================*/
extern "C" void cuda_alloc_copyin_int_(void **     ptr,
                                       const int * val_host,
                                       const int & N) {
/*-----------------------------------------------------------------------------+
|   Allocates device memory and copies integers from host to device            |
+-----------------------------------------------------------------------------*/

  cudaError_t err = cudaMalloc(ptr, N * sizeof(int));
  if (err != cudaSuccess) {
    fprintf(stderr, "cudaMalloc (int) failed: %s\n",
            cudaGetErrorString(err));
    return;
  }

  err = cudaMemcpy(*ptr, val_host, N * sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy H->D (int) failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*============================================================================*/
extern "C" void cuda_copyout_double_(double *      val_host,
                                     const void ** ptr_dev,
                                     const int &   N) {
/*-----------------------------------------------------------------------------+
|   Copies N doubles from device memory to host array                          |
+-----------------------------------------------------------------------------*/

  cudaError_t err =
    cudaMemcpy(val_host, *ptr_dev, N * sizeof(double), cudaMemcpyDeviceToHost);
  if (err != cudaSuccess) {
    fprintf(stderr, "cudaMemcpy D->H (double) failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*============================================================================*/
extern "C" void cuda_free_(void * ptr) {
/*-----------------------------------------------------------------------------+
|   Frees previously allocated device memory                                   |
+-----------------------------------------------------------------------------*/

  cudaError_t err = cudaFree(ptr);
  if(err != cudaSuccess) {
    fprintf(stderr, "cudaFree failed: %s\n",
            cudaGetErrorString(err));
  }
}

/*-----------------------------------------------------------------------------+
|                                                                              |
|  GPU is not defined                                                          |
|                                                                              |
+-----------------------------------------------------------------------------*/
#else

extern "C" void cuda_alloc_double_(void ** ptr, const int & N) {}

extern "C" void cuda_alloc_copyin_double_(void **        ptr,
                                          const double * val_host,
                                          const int    & N) {}

extern "C" void cuda_alloc_int_(void ** ptr, const int & N) {}

extern "C" void cuda_alloc_copyin_int_(void **     ptr,
                                       const int * val_host,
                                       const int & N) {}

extern "C" void cuda_copyout_double_(double *      val_host,
                                     const void ** ptr_dev,
                                     const int &   N) {}

extern "C" void cuda_free_(void * ptr) {}

#endif
