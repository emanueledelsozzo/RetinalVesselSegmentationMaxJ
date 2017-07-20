// Authors:
// Emanuele Del Sozzo (emanuele.delsozzo@polimi.it), Marcello Pogliani (marcello.pogliani@polimi.it)

#include <stdio.h>

#include "benchmark_utils.h"
#include <stdlib.h>

typedef struct duration {
    double val;
    char* why;
    struct duration* next;
} duration_t;

duration_t* durations = NULL;

void print_duration(struct timeval start, struct timeval end, char* why)
{
    double duration = ((end.tv_sec-start.tv_sec)*1000000
            + end.tv_usec - start.tv_usec)/1000.0;
    add_duration(duration, why);
}

void add_duration(double val, char* why)
{
    if(durations) {
        duration_t* dur = durations;
        while(dur->next)
            dur = dur->next;
        dur->next = malloc(sizeof(duration_t));
        dur->next->next = NULL;
        dur->next->val = val;
        dur->next->why = why;
    } else {
        durations = malloc(sizeof(duration_t));
        durations->next = NULL;
        durations->val = val;
        durations->why = why;
    }
}

void print_durations()
{
    duration_t* dur = durations;
    while(dur) {
        fprintf(stdout, "%lf", dur->val);
        dur = dur->next;
        if(dur) {
            fprintf(stdout, ",");
        }
    }
    fprintf(stdout, "\n");
}

void print_causes()
{
    duration_t* dur = durations;
    while(dur) {
        fprintf(stdout, "%s", dur->why);
        dur = dur->next;
        if(dur) {
            fprintf(stdout, ",");
        }
    }
    fprintf(stdout, "\n");
}

