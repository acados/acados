
#include "common_header.h"
#include <stdio.h>

/** A simple helper function. */
void printMatrix(	const char* name,
					const real_t* mat,
					uint nRows,
					uint nCols
					)
{
    uint r, c;
    printf("%s: \n", name);
    for (r = 0; r < nRows; ++r)
    {
        for(c = 0; c < nCols; ++c)
        	printf("\t%f", mat[r * nCols + c]);
        printf("\n");
    }
}
