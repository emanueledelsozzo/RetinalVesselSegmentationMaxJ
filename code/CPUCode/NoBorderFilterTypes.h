// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#ifndef NO_BORDER_FILTER_TYPES_H
#define NO_BORDER_FILTER_TYPES_H

#include <stdint.h>

#include "NoBorderFilter.h"

// Hack to synchronize CPUCode types with the ones defined in RunRules.settings
// defining a MaxFileConstant in the Manager supports only numeric values,
// not generic preprocessor "#define"s.

#ifndef NoBorderFilter_input_type
    #error NoBorderFilter_input_type undefined
#endif

#if NoBorderFilter_input_type == 0
    #define IPL_INPUT_IMG_DEPTH IPL_DEPTH_32S
    typedef int32_t input_type;
#elif NoBorderFilter_input_type == 1
    #define IPL_INPUT_IMG_DEPTH IPL_DEPTH_32F
    typedef float input_type;
#elif NoBorderFilter_input_type == 2
    #define IPL_INPUT_IMG_DEPTH IPL_DEPTH_16U
    typedef uint16_t input_type;
#elif NoBorderFilter_input_type == 3
    #define IPL_INPUT_IMG_DEPTH IPL_DEPTH_8U
    typedef uint8_t input_type;
#else
    #error Unrecognized input type
#endif

#ifndef NoBorderFilter_output_type
    #error NoBorderFilter_output_type undefined
#endif

#if NoBorderFilter_output_type == 0
    #define IPL_OUTPUT_IMG_DEPTH IPL_DEPTH_32S
    typedef int32_t output_type;
#elif NoBorderFilter_output_type == 1
    #define IPL_OUTPUT_IMG_DEPTH IPL_DEPTH_32F
    typedef float output_type;
#elif NoBorderFilter_output_type == 2
    #define IPL_OUTPUT_IMG_DEPTH IPL_DEPTH_16U
    typedef uint16_t output_type;
#elif NoBorderFilter_output_type == 3
    #define IPL_OUTPUT_IMG_DEPTH IPL_DEPTH_8U
    typedef uint8_t output_type;
#else
    #error Unrecognized output type
#endif

#endif /* NO_BORDER_FILTER_TYPES_H */
