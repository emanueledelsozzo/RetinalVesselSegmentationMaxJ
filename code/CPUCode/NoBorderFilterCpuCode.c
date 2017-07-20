// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

#include "omp.h"
#include "opencv/highgui.h"
#include "opencv/cv.h"

#include "NoBorderFilterTypes.h"

#define BENCHMARK
#define DEBUG

#include "benchmark_utils.h"
#include "NoBorderFilterReferenceImpl.h"

#define MAX_PRINT_MATRIX 16

#define IMPL_SIMPLE


//int kernel_zero[] = {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0};
//int kernel_line[] = {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0};

/* The kernel must be even */
static size_t kernelSize = NoBorderFilter_kernel_size;

static output_type kernel[16][16] = {
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 4, 3, 2, 1, -2, -5, -6, -5, -2, 1 ,2, 3, 4, 0, 0},
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0,  0,  0,  0,  0,  0, 0, 0, 0, 0, 0, 0}
    };


void print_matrix_in(input_type* src, size_t width, size_t height, FILE* f)
{
    for(size_t i = 0; i < MIN(MAX_PRINT_MATRIX, height); i++) {
        for(size_t j = 0; j < MIN(MAX_PRINT_MATRIX, width); j++) {
            fprintf(f, IPL_INPUT_IMG_DEPTH == IPL_DEPTH_32F ? "%3f" : "%3d ", src[i*width +j]);
        }
        fprintf(f, "\n");
    }
}


void print_matrix_out(output_type* src, size_t width, size_t height, FILE* f)
{
    for(size_t i = 0; i < MIN(MAX_PRINT_MATRIX, height); i++) {
	    for(size_t j = 0; j < MIN(MAX_PRINT_MATRIX, width); j++) {
            fprintf(f, IPL_OUTPUT_IMG_DEPTH == IPL_DEPTH_32F ? "%3f" : "%3d ", src[i*width +j]);
        }
        fprintf(f, "\n");
    }
}



void set_dfe_kernel_coefficients(max_engine_t* dfe) {
    int64_t* tmp = malloc(NoBorderFilter_kernel_size * NoBorderFilter_kernel_size * sizeof(int64_t));
    int k = 0;
    for(int j=0; j < NoBorderFilter_kernel_size; j++) {
        for(int i = 0; i < NoBorderFilter_kernel_size; i++) {
           tmp[k++] = (int64_t) kernel[i][j];
        }
    }

    NoBorderFilter_setCoefficients_actions_t coef_acn;
    coef_acn.inmem_NoBorderFilterKernel_coefficients_0 = (uint64_t*) &tmp[0];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_1 = (uint64_t*) &tmp[1 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_2 = (uint64_t*) &tmp[2 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_3 = (uint64_t*) &tmp[3 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_4 = (uint64_t*) &tmp[4 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_5 = (uint64_t*) &tmp[5 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_6 = (uint64_t*) &tmp[6 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_7 = (uint64_t*) &tmp[7 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_8 = (uint64_t*) &tmp[8 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_9 = (uint64_t*) &tmp[9 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_10 = (uint64_t*) &tmp[10 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_11 = (uint64_t*) &tmp[11 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_12 = (uint64_t*) &tmp[12 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_13 = (uint64_t*) &tmp[13 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_14 = (uint64_t*) &tmp[14 * NoBorderFilter_kernel_size];
    coef_acn.inmem_NoBorderFilterKernel_coefficients_15 = (uint64_t*) &tmp[15 * NoBorderFilter_kernel_size];

    NoBorderFilter_setCoefficients_run(dfe, &coef_acn);

    free(tmp);
}

void noBorderFilterDFERun(input_type* src, output_type* dst, int dimx, int dimx_aligned, int dimy, int bufsize, int kernel_height)
{

#ifdef BENCHMARK
    struct timeval loadStart, loadEnd, writeLMemStart, writeLMemEnd, kernelsStart, kernelsEnd, readLMemStart, readLMemEnd;
#endif

    // invoke dfe
    max_file_t* maxfile = NoBorderFilter_init();
    BENCH_GETTIME(&loadStart);
    max_engine_t* dfe = max_load(maxfile, "local:*");
    BENCH_GETTIME(&loadEnd);

    //int initial_offset = dimx_aligned * (kernel_height / 2 - 1) * sizeof(input_type);
    int initial_offset = dimx_aligned * (kernel_height / 2) * sizeof(input_type);


    set_dfe_kernel_coefficients(dfe);


    BENCH_GETTIME(&writeLMemStart);
    NoBorderFilter_writeLMem_actions_t write_acn;
    write_acn.param_size = bufsize;
    write_acn.param_start = initial_offset;
    write_acn.instream_fromcpu = src;
#ifdef DEBUG
    fprintf(stderr, "Loading data to LMem with parameters: size = %d, start = %d\n",
            write_acn.param_size,
            write_acn.param_start);
#endif
    NoBorderFilter_writeLMem_run(dfe, &write_acn);
    BENCH_GETTIME(&writeLMemEnd);

    // Strong assumption ~> kernel_height is divisible by the number of streams per run
    BENCH_GETTIME(&kernelsStart);
#ifdef DEBUG
    fprintf(stderr, "Using %d lanes\n", NoBorderFilter_number_of_lanes);
#endif
    int i;
    for(i = 0; i < kernel_height / NoBorderFilter_number_of_lanes; i++) {
        NoBorderFilter_oneKernel_actions_t kern_acn;
        kern_acn.param_size_x = dimx;
        kern_acn.param_size_x_aligned = dimx_aligned;
        kern_acn.param_size_y = dimy;
        kern_acn.param_start_in = i * NoBorderFilter_number_of_lanes * dimx_aligned * sizeof(input_type);
        kern_acn.param_start_out = (i % 2 ? 1 : 2) * bufsize * sizeof(output_type) + initial_offset;
        kern_acn.param_start_feedback = (i % 2 ? 2 : 1) * bufsize * sizeof(output_type) + initial_offset;
        kern_acn.param_run = i * NoBorderFilter_number_of_lanes;
#ifdef DEBUG
        fprintf(stderr, "Running kernel with parameters: start_in=%d, start_out=%d, start_feedback = %d, size_x = %d, size_y = %d, run = %d\n",
            kern_acn.param_start_in,
            kern_acn.param_start_out,
            kern_acn.param_start_feedback,
            kern_acn.param_size_x,
            kern_acn.param_size_y,
            kern_acn.param_run);
#endif
        NoBorderFilter_oneKernel_run(dfe, &kern_acn);
    }
    BENCH_GETTIME(&kernelsEnd);

    BENCH_GETTIME(&readLMemStart);
    NoBorderFilter_readLMem_actions_t read_acn;
    read_acn.param_size = bufsize;
    read_acn.param_start = (i % 2 ? 2 : 1) * bufsize * sizeof(output_type) + initial_offset;
    read_acn.outstream_tocpu = dst;
    NoBorderFilter_readLMem_run(dfe, &read_acn);
    BENCH_GETTIME(&readLMemEnd);

    max_unload(dfe);

#ifdef BENCHMARK
    print_duration(loadStart, loadEnd, "DFELoad");
    print_duration(writeLMemStart, writeLMemEnd, "DFEWriteLMem");
    print_duration(kernelsStart, kernelsEnd, "DFEComputation");
    print_duration(readLMemStart, readLMemEnd, "DFEReadLMem");
#endif

}

size_t align_to_burst(size_t size, size_t burst){
	size_t alignment = size % burst;
	if (alignment == 0){
		return size;
	}else{
		return size + (burst - alignment);
	}
}

int run_NoBorderFilter_dfe_wrapper(input_type* imgbuf, output_type** output, size_t width, size_t height)
{
    // base size ~> 288
    const size_t burst_size = 384;
    const size_t aligned_width = align_to_burst(width * sizeof(input_type), burst_size) / sizeof(input_type);
    const size_t aligned_full = aligned_width * height;

#ifdef BENCHMARK
    add_duration((float) width, "Width");
    add_duration((float) height, "Height");
#endif

#ifdef DEBUG
    fprintf(stderr, "DFE - width %ld, height %ld, aligned width %ld, aligned size %ld\n",
            width, height, aligned_width, aligned_full);
#endif

    input_type* input = malloc(aligned_full * sizeof(input_type));
    output_type* dfeout = malloc(aligned_full * sizeof(output_type));
    if(!input || !dfeout) {
        fprintf(stderr, "Malloc cannot allocate memory\n");
        return -1;
    }
    memset(input,  0, aligned_full * sizeof(input_type));
    memset(dfeout,  0, aligned_full * sizeof(output_type));

    // Align each row to the LMem burst

    #pragma omp parallel for
    for(size_t i = 0; i < height; i++) {
        memcpy(&(input[i * aligned_width]), &(imgbuf[i * width]), width * sizeof(input_type));
    }

    noBorderFilterDFERun(input, dfeout, width, aligned_width, height, aligned_full, kernelSize);

    free(input);

    // Do the very same thing to the output. Can do in place, actually.
    // memmove works even if src and dst overlap; memcpy is more efficient but does not work in this case!
    // this loop is not parallel due to the same reason
    for(size_t i = 0; i < height; i++) {
        memmove(&dfeout[i * width], &dfeout[i * aligned_width], width * sizeof(output_type));
    }

    *output = dfeout;

    return 0;
}

int run_NoBorderFilter(input_type* imgbuf, size_t width, size_t height, output_type** dfeout_ptr, output_type** cpuout_ptr)
{

#ifdef BENCHMARK
    struct timeval cpuStart, cpuEnd, dfeStart, dfeEnd, ompStart, ompEnd;
#endif

    output_type* cpuout = malloc(width * height * sizeof(output_type));
    output_type* ompout = malloc(width * height * sizeof(output_type));
    output_type* dfeout = NULL;

    if (!cpuout || !ompout) {
        fprintf(stderr, "Malloc cannot allocate memory\n");
        exit(1);
    }

    memset(cpuout, 0, width * height * sizeof(output_type));
    memset(ompout, 0, width * height * sizeof(output_type));

    BENCH_GETTIME(&dfeStart);
    run_NoBorderFilter_dfe_wrapper(imgbuf, &dfeout, width, height);
    BENCH_GETTIME(&dfeEnd);

    BENCH_GETTIME(&ompStart);
    noBorderFilterCPUParallel(imgbuf, ompout, height, width, (output_type*) kernel, kernelSize);
    BENCH_GETTIME(&ompEnd);

    BENCH_GETTIME(&cpuStart);
    noBorderFilterCPUReferenceImpl(imgbuf, cpuout, height, width, (output_type*) kernel, kernelSize);
    BENCH_GETTIME(&cpuEnd);

    int nwrong = 0, ntot = 0;
    for(uint32_t i=0; i < height; i++) {
        for(uint32_t j=0; j < width; j++) {
            if (dfeout[i * width + j] != cpuout[i * width + j] || cpuout[i * width + j] != ompout[i * width + j]) {
                nwrong++;
            }
            ntot++;
        }
    }

#ifdef DEBUG
    fprintf(stderr, "Input: \n");
    print_matrix_in(imgbuf, width, height, stderr);

    fprintf(stderr, "Output: \n");
    print_matrix_out(cpuout, width, height, stderr);

    fprintf(stderr, "DFE Output: \n");
    print_matrix_out(dfeout, width, height, stderr);

    fprintf(stderr, "Nwrong: %d/%d\n", nwrong, ntot);
#endif

#ifdef BENCHMARK
    print_duration(cpuStart, cpuEnd, "cpu");
    print_duration(ompStart, ompEnd, "omp");
    print_duration(dfeStart, dfeEnd, "dfe");
#endif

    *dfeout_ptr = dfeout;
    *cpuout_ptr = cpuout;
    free(ompout);

    return 0;
}

void print_help(char* progname)
{
    fprintf(stdout, "Usage: \n");
    fprintf(stdout, " * Random input: %s --[random|sequential] --width <width> --height <height>\n", progname);
    fprintf(stdout, " * Image input: %s --file filename.jpg\n", progname);
}

int generate_imgbuf(input_type** imgbuf, size_t width, size_t height, int isSequential)
{
    input_type* input = malloc(sizeof(input_type) * width * height);
    if(!input)
        return -1;

    int k = 0;
    for(size_t i = 0; i < height; ++i) {
        for(size_t j = 0; j < width; ++j) {
            input[i * width + j] = isSequential ? k++ : ( random() % 100 );
        }
    }

    *imgbuf = input;

    return 0;
}

int load_from_file(const char* filename, IplImage** imgbuf)
{
    IplImage* source = cvLoadImage(filename, CV_LOAD_IMAGE_ANYDEPTH | CV_LOAD_IMAGE_COLOR);
    int mpixels = source->height * source->width;
    float lin_scale = sqrt((float)(mpixels) / 0.8e6);
    fprintf(stderr, "lin_scale = %f\n", lin_scale);
    if(lin_scale < 0.9 || lin_scale > 1.1) {
        int new_width = (int)(source->width / lin_scale);
        int new_height = (int)(source->height / lin_scale);
        IplImage* newimg = cvCreateImage(cvSize(new_width, new_height), source->depth, source->nChannels);
        cvResize(source, newimg, CV_INTER_LINEAR);
        cvReleaseImage(&source);
        source = newimg;
    }

    int crop_size = abs(source->width - source->height)/2;
    fprintf(stderr, "source->width: %d, source->height: %d, crop_size: %d\n", source->width, source->height, crop_size);
    cvSetImageROI(source, cvRect(crop_size, 0, source->width - 2 * crop_size, source->height));

    *imgbuf = source;

    return 1;
}

int main(int argc, char** argv)
{
    static struct option long_options[] = {
        { "random",      no_argument,       NULL, 'r' },
        { "sequential",  no_argument,       NULL, 's' },
        { "file",        required_argument, NULL, 'f' },
        { "width",       required_argument, NULL, 'w' },
        { "height",      required_argument, NULL, 'h' },
        { "causes",      no_argument,       NULL, 'c' },
        { NULL,          0,                 NULL, 0 }
    };
    int option_index = 0;

    int ch;
    enum {MODE_UNSPEC, MODE_RANDOM, MODE_FILE, MODE_SEQUENTIAL} mode = MODE_UNSPEC;
    size_t width = 0, height = 0;
    int opt_print_causes = 0;
    char* filename = NULL;
    input_type *imgbuf;
    output_type *cpuout, *dfeout;

    while((ch = getopt_long(argc, argv, "rsf:w:h:c", long_options, &option_index)) != -1) {
        switch(ch) {
            case 'r':
                mode = MODE_RANDOM;
                break;
            case 's':
                mode = MODE_SEQUENTIAL;
                break;
            case 'f':
                mode = MODE_FILE;
                filename = optarg;
                break;
            case 'w':
                width = atoi(optarg);
                break;
            case 'h':
                height = atoi(optarg);
                break;
            case 'c':
                opt_print_causes = 1;
                break;
            default:
                fprintf(stderr, "Unknown option\n");
                print_help(argv[0]);
                exit(1);
        }
    }

    IplImage* img32 = NULL;

    if(mode == MODE_RANDOM || mode == MODE_SEQUENTIAL) {
        if(generate_imgbuf(&imgbuf, width, height, mode == MODE_SEQUENTIAL) < 0) {
            fprintf(stderr, "Error generating image buffer\n");
            exit(1);
        }
    } else if(mode == MODE_FILE) {
        IplImage* source;
        if(load_from_file(filename, &source) < 0) {
            fprintf(stderr, "Error loading image from file\n");
            exit(1);
        }
        IplImage* red = cvCreateImage(cvGetSize(source), source->depth, 1);
        IplImage* green = cvCreateImage(cvGetSize(source), source->depth, 1);
        IplImage* blue = cvCreateImage(cvGetSize(source), source->depth, 1);

        img32 = cvCreateImage(cvGetSize(source),IPL_INPUT_IMG_DEPTH, 1);

        cvSplit(source, blue, green, red, NULL);

        cvReleaseImage(&source);
        cvReleaseImage(&blue);
        cvReleaseImage(&red);

        cvSaveImage("original.png", green, NULL);

        cvConvertScale(green, img32, 1, 0);


        width = green->width;
        height = green->height;

        cvReleaseImage(&green);

        cvSaveImage("rescaled.png", img32, NULL);

        cvGetRawData(img32, (unsigned char**) &imgbuf, NULL, NULL);
    } else {
        print_help(argv[0]);
        exit(1);
    }

    int ret = run_NoBorderFilter(imgbuf, width, height, &dfeout, &cpuout);

    // process the image and output the file, finally
    if (mode == MODE_FILE) {
        IplImage* dfeoutImg = cvCreateImage(cvGetSize(img32), IPL_OUTPUT_IMG_DEPTH, 1);
        output_type* outbuf;
        cvGetRawData(dfeoutImg, (unsigned char**) &outbuf, NULL, NULL);

        memcpy(outbuf, dfeout, sizeof(output_type) * width * height);
        cvSaveImage("dfeout.png", dfeoutImg, NULL);

        memcpy(outbuf, cpuout, sizeof(output_type) * width * height);
        cvSaveImage("cpuout.png", dfeoutImg, NULL);

        cvReleaseImage(&dfeoutImg);
        //cvReleaseImage(&img32);
    }

    //free(cpuout);
    //free(dfeout);

#ifdef BENCHMARK
    if(opt_print_causes) {
        print_causes();
    }
    print_durations();
#endif

    return ret;
}
