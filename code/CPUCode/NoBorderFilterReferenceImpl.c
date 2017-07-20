// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "NoBorderFilterReferenceImpl.h"

#include "omp.h"

void noBorderFilterCPUParallel(input_type src[], output_type dst[], size_t dimx, size_t dimy, output_type Kernel[], size_t Kdim)
{
    output_type sum = 0;
    int start = (int) (Kdim) / 2;
	
    #pragma omp parallel for private(sum)
    for (size_t x = start; x < (dimx - start); x++) {
        for (size_t y = start; y < (dimy - start); y++) {
            sum = 0;
            for (int k = -start; k < start; k++) {
                for (int j = -start; j < start; j++) {
                    sum = sum + Kernel[(int)(k + start) * Kdim + (int) (j + start)] * (output_type) src[(int)(x - k) * dimy + y - j];
                }
            }
            dst[x * dimy + y] = sum;
        }
    }
}


void noBorderFilterCPUReferenceImpl(input_type src[], output_type dst[], size_t dimx, size_t dimy, output_type Kernel[], size_t Kdim)
{
    output_type sum = 0;
    int start = (int) (Kdim) / 2;

    for (size_t y = start; y < (dimy - start); y++) {
        for (size_t x = start; x < (dimx - start); x++) {
            sum = 0;
            for (int k = -start; k < start; k++) {
                for (int j = -start; j < start; j++) {
                    sum = sum + Kernel[(int)(k + start) * Kdim + (int) (j + start)] * (output_type) src[(int)(x - k) * dimy + y - j];
                }
            }
            dst[x * dimy + y] = sum;
        }
    }
}
