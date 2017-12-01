#include <stdlib.h>
#include <stdio.h>
#ifndef NOPAPI
#include <papi.h>
#endif
#include <string.h>
#include <errno.h>
#include <sys/time.h>
#include "papi_helpers.h"

#ifdef NOPAPI
char* papi_error_to_string(int err_code) { return 0; }
void handle_error(char* msg, int err_code) { }
void check_error(int err_code, char* msg) { }
#else

/* Returns a string description of the given error code. */
char* papi_error_to_string(int error_code)
{
  switch (error_code)
  {
    case PAPI_EINVAL:
      return "One or more of the arguments is invalid.";
    case PAPI_ENOMEM:
      return "Insufficient memory to complete the operation.";
    case PAPI_ENOEVST:
      return "The event set specified does not exist.";
    case PAPI_EISRUN:
      return "The event set is currently counting events.";
    case PAPI_ECNFLCT:
      return "The underlying counter hardware can not count this event and other events" \
             " in the event set simultaneously.";
    case PAPI_ENOEVNT:
      return "The PAPI preset is not available on the underlying hardware.";
    case PAPI_EBUG:
      return "Internal error, please send mail to the developers.";
    case PAPI_ESYS:
    {
      char* msg_buff = malloc(100 * sizeof(char));

      sprintf(msg_buff,
          "A system or C library call failed inside PAPI, errno is %d\n",
          errno);

      return msg_buff;
    }
    default:
    {
      char* msg_buff = malloc(100 * sizeof(char));

      sprintf(msg_buff, "Return code: %d", error_code);

      return msg_buff;
    }
  }
}

/* Prints an error message and exits. */
void handle_error(char* msg, int error_code)
{
  char* papi_msg = papi_error_to_string(error_code);
  fprintf(stderr, "%s\nMessage: %s\n",
    msg, papi_msg);
  free(papi_msg);
  exit(1);
}

/* Meant to wrap calls to PAPI libraries. */
void check_error(int code, char* msg)
{
  if (code != PAPI_OK)
    handle_error(msg, code);
}

#endif

/* Returns the current time. */
double current_time()
{
	struct timeval tv;
	int retval = gettimeofday(&tv, NULL);
	if(retval < 0)
	{
		fprintf(stderr, "Timer gettimeofday error.\n");
		exit(1);
	}
	return tv.tv_sec * 1.0 + tv.tv_usec * 1.0e-6;
}

