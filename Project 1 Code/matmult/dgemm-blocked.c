/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm (int n, double* A, double* B, double* C)
{
    // TODO: Implement the blocking optimization
    for (int i = 0; i < n-5; i+=2) {
        for (int j = 0; j < n-5; j+=2) {
            double c00 = C[i+j*n];
            double c01 = C[i+(j+1)*n];
            double c10 = C[(i+1)+j*n];
            double c11 = C[(i+1)+(j+1)*n];
            for (int k = 0; k < n; ++k) {
                double a0 = A[i+k*n];
                double a1 = A[(i+1)+k*n];
                double b0 = B[k+j*n];
                double b1 = B[k+(j+1)*n];
                c00 += a0 * b0;
                c01 += a0 * b1;
                c10 += a1 * b0;
                c11 += a1 * b1;
            }
            C[i+j*n] = c00;
            C[i+(j+1)*n] = c01;
            C[(i+1)+j*n] = c10;
            C[(i+1)+(j+1)*n] = c11;
        }
    }
    for (int i = n-5; i < n; ++i) {
        for (int j = n-5; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                C[i+j*n] += A[i+k*n] * B[k+j*n];
            }
        }
    }
}