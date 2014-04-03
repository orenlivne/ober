#include <stdio.h>
/* Apply a median filter of size filter_size to each row of the m x n matrix data. Fill outside with 0's. */
long median_filter(char *data, long m, long n, long filter_size) {
    long i, j;
    char d = *data;
    for (i = 0; i < m; i++) {
    	for (j = 0; j < n; j++) {
    		printf("%ld %ld %d\n", i, j, d);
    		d++;
    	}
    }
    return 0;
}
