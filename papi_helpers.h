#ifndef _PAPI_HELPERS_H
#define _PAPI_HELPERS_H

/* Returns a string description of the given error code. */
char* papi_error_to_string(int error_code);

/* Prints an error message and exits. */
void handle_error(char* msg, int error_code);

/* Meant to wrap calls to PAPI libraries. */
void check_error(int code, char* msg);

/* Returns the current time. */
double current_time();

#endif

