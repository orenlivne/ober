#include <stdio.h>

long smallest_index_long(long val, long *data, long length, long start, long step) {
    long i;
	for (i = start; (i >= 0) && (i < length); i += step) {
		if (data[i] == val)
			return(i);
	}	
    return(-1L);
}

long smallest_index_byte(long val, char *data, long length, long start, long step) {
    long i;
	for (i = start; (i >= 0) && (i < length); i += step) {
		if (data[i] == val)
			return(i);
	}
    return(-1L);
}
