#include <stdio.h>
#include <plotting.h>
void plot_states_controls(real_t *w, real_t T, int_t NN, int_t NX, int_t NU, FILE *gnuplotPipe) {
      double t_grid[NN];
      for (int_t i = 0; i < NN; i++) t_grid[i] = i*T;

      FILE *tempDataFile;
      char *x1_temp_file;
      char *x2_temp_file;
      char *x3_temp_file;
      char *x4_temp_file;
      char *x5_temp_file;
      char *x6_temp_file;
      char *x7_temp_file;
      char *x8_temp_file;
      char *x9_temp_file;
      char *x10_temp_file;
      char *x11_temp_file;
      char *x12_temp_file;

      char *u1_temp_file;
      char *u2_temp_file;
      char *u3_temp_file;
      char *u4_temp_file;

      double x, y;
      int i;
      x1_temp_file = "q1";
      if (gnuplotPipe) {
          fprintf(gnuplotPipe, "set multiplot layout 4,4\n");

          // Plot q1
          tempDataFile = fopen(x1_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x1_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot q2
          x2_temp_file = "q2";
          tempDataFile = fopen(x2_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x2_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+1];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot q3
          x3_temp_file = "q3";
          tempDataFile = fopen(x3_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x3_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+2];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot q4
          x4_temp_file = "q4";
          tempDataFile = fopen(x4_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x4_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+3];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);




          // Plot w1
          x5_temp_file = "w1";
          tempDataFile = fopen(x5_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x5_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+4];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot w2
          x6_temp_file = "w2";
          tempDataFile = fopen(x6_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x6_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+5];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot w3
          x7_temp_file = "w3";
          tempDataFile = fopen(x7_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x7_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+6];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          x8_temp_file = "empty";
          tempDataFile = fopen(x8_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x8_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = 0;
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u1
          x9_temp_file = "u1";
          tempDataFile = fopen(x9_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x9_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+7];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u2
          x10_temp_file = "u2";
          tempDataFile = fopen(x10_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x10_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+8];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u3
          x11_temp_file = "u3";
          tempDataFile = fopen(x11_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x11_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+9];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u4
          x12_temp_file = "u4";
          tempDataFile = fopen(x12_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with lines lt rgb \"blue\"\n", x12_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+10];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);




          // Plot u1r
          u1_temp_file = "u1r";
          tempDataFile = fopen(u1_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u1_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+11];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u2r
          u2_temp_file = "u2r";
          tempDataFile = fopen(u2_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u2_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+12];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u3r
          u3_temp_file = "u3r";
          tempDataFile = fopen(u3_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u3_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+13];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);

          // Plot u2r
          u4_temp_file = "u4r";
          tempDataFile = fopen(u4_temp_file, "w");
          fprintf(gnuplotPipe, "set grid ytics\n");
          fprintf(gnuplotPipe, "set grid xtics\n");
          fprintf(gnuplotPipe, "set xlabel \"%s\"\n", "time [s]");
          fprintf(gnuplotPipe, "plot \"%s\" with steps lt rgb \"red\" \n", u4_temp_file);
          fflush(gnuplotPipe);
          for (i=0; i < NN; i++) {
              x = t_grid[i];
              y = w[i*(NX+NU)+14];
              fprintf(tempDataFile, "%lf %lf\n", x, y);
          }
          fclose(tempDataFile);



      } else {
          printf("gnuplot not found...");
      }
}
