#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>
#include "papi_support.h"

double start_time;

#if NOPAPI
void papi_init() { }
void papi_hw_details() {}
void papi_prepare_counter(int event_code) { }
#else

int eventset = PAPI_NULL;
int ins_counter_active = 0;

/* Initializes PAPI system. */
void papi_init()
{
  /* Init the PAPI library. */
	int retval = PAPI_library_init(PAPI_VER_CURRENT);
	if(retval != PAPI_VER_CURRENT)
    handle_error("PAPI library init error.", retval);

  /* Initialize the eventset. */
  check_error(PAPI_create_eventset(&eventset),
    "PAPI eventset initialization error");

  printf("Intializing PAPI.\n");
}

/* The following are functions to initialize the various counters. */

void papi_prepare_counter(int event_code)
{
  check_error(PAPI_add_event(eventset, event_code),
      "PAPI add events error.");

  char eventCodeStr[PAPI_MAX_STR_LEN];
  PAPI_event_code_to_name(event_code, eventCodeStr);

  printf("Prepared counter: '%s'.\n", eventCodeStr);
}

#endif

/* Starts the counters. */
void papi_start()
{
#ifndef NOPAPI
  printf("Starting PAPI counters.\n");

  check_error(PAPI_start(eventset),
      "PAPI start counters error.");
#endif

  start_time = current_time();
}

/* Stops the counters and prints the results. */
void papi_stop_and_report()
{
  /* Print stats. */
  printf("Execution time: %1.7f\n", current_time() - start_time);

#ifndef NOPAPI
  /* Will hold counter values. */
  long_long values[2];

  /* Place counter values into the above values arrays. */
  check_error(PAPI_stop(eventset, values),
      "PAPI stop counters error.\n");
#endif

#ifndef NOPAPI
  printf("Event 1: %lld\nEvent 2: %lld\n",
    values[0],
    values[1]);

  /* Reset the counters. */
  printf("Resetting PAPI.\n");
  check_error(PAPI_reset(eventset),
      "Resetting counters failed.");
#endif
}
