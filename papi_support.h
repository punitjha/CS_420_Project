#ifndef _PAPI_SUPPORT_H_
#define _PAPI_SUPPORT_H_

#ifndef NOPAPI
#include <papi.h>
#endif
#include "papi_helpers.h"

/** Initialize PAPI */
void papi_init();

/** Mark a PAPI event to be recorded */
void papi_prepare_counter(int eventcode);

/** Start recording the PAPI counters */
void papi_start();

/** Stop recording the PAPI counters */
void papi_stop_and_report();

#endif
