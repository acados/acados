#include <stdio.h>
#include "acados/utils/types.h"
#include "examples/c/acados_gnuplot/acados_gnuplot.h"

#define _GNU_SOURCE

extern FILE *popen(const char *command, const char *type);

void acados_gnuplot(real_t **data, int_t n_data, real_t T,
    int_t N, char **labels, int_t layout_x, int_t layout_y) {

    double t_grid[N];
    for (int_t i = 0; i < N; i++) t_grid[i] = i*T/N;

    FILE *gnuplotPipe, *tempDataFile;
    char *temp_file = NULL;

    double x, y;
    int i;
    gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe) {
        fprintf(gnuplotPipe, "set multiplot layout %i,%i\n", layout_x, layout_y);
        for (int_t k = 0; k < n_data; k++) {
            temp_file = labels[k];
            // Plot x1
            tempDataFile = fopen(temp_file, "w");
            fprintf(gnuplotPipe, "set grid ytics\n");
            fprintf(gnuplotPipe, "set grid xtics\n");
            fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
            fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", temp_file);
            fflush(gnuplotPipe);
            for (i=0; i < N; i++) {
                x = t_grid[i];
                y = data[k][i];
                fprintf(tempDataFile, "%lf %lf\n", x, y);
            }
            fclose(tempDataFile);
        }

    printf("Press enter to continue...");
    getchar();
    remove(temp_file);

    fprintf(gnuplotPipe, "exit gnuplot\n");
    } else {
        printf("gnuplot not found...");
    }
}
