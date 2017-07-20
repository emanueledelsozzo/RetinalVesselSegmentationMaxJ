// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#ifndef NO_BORDER_FILTER_REFERENCE_IMPL
#define VESSEL_REFERENCE_IMPL

#include "NoBorderFilterTypes.h"

void noBorderFilterCPUParallel(input_type src[], output_type dst[], size_t dimx, size_t dimy, output_type Kernel[], size_t Kdim);

void noBorderFilterCPUReferenceImpl(input_type src[], output_type dst[], size_t dimx, size_t dimy, output_type Kernel[], size_t Kdim);

#endif
