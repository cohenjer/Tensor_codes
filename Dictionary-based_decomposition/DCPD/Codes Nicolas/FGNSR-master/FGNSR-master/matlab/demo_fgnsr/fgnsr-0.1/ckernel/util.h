#ifndef UTIL_H_INCDEF
#define UTIL_H_INCDEF

/****************************************************************/
/* Data structures */
/****************************************************************/

/* Adapt for portability and functionality via chains of #define ... */
struct datum {
    double tag_cpu;
    double  tag_wc;
};
typedef struct datum DATUM;

/****************************************************************/
/* Memory management */
/****************************************************************/

void *FGNSRcalloc(size_t count, size_t size);
void *FGNSRmalloc(size_t size);
void FGNSRfree(void **ptr);

/****************************************************************/
/* Time measurment */
/****************************************************************/

/* Obtain portable time stamp */
void FGNSRtimestamp(DATUM *stamp);

/* Wallclock time in seconds */
double FGNSRwctime(DATUM start, DATUM end);

/* CPU time in seconds */
double FGNSRcputime(DATUM start, DATUM end);

#endif /* UTIL_H_INCDEF */
