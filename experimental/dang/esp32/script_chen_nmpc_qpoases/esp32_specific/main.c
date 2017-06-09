#include "freertos/FreeRTOS.h"
#include "esp_wifi.h"
#include "esp_system.h"
#include "esp_event.h"
#include "esp_event_loop.h"
#include "nvs_flash.h"
#include "driver/gpio.h"

/* Begin acados code */
#include <stdio.h>
#include <stdlib.h>

#include "aux_d.h"

#include "ocp_qp_condensing_qpoases.h"
#include "sim_erk_integrator.h"
#include "timing.h"
#include "types.h"
#include "chen_model.h"

#define NN 13
#define NX 2
#define NU 1

#ifdef DEBUG
static void print_states_controls(real_t *w, int_t N) {
    printf("node\tx\t\t\t\tu\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e\t%+e\n", i, w[i*(NX+NU)], w[i*(NX+NU)+1], w[i*(NX+NU)+2]);
    }
    printf("%4d\t%+e %+e\n", N, w[N*(NX+NU)], w[N*(NX+NU)+1]);
}
#endif  // DEBUG

static void shift_states(real_t *w, real_t *x_end, int_t N) {
    for (int_t i = 0; i < N; i++) {
        for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
    }
    for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = x_end[j];
}

static void shift_controls(real_t *w, real_t *u_end, int_t N) {
    for (int_t i = 0; i < N-1; i++) {
        for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
    }
    for (int_t j = 0; j < NU; j++) w[(N-1)*(NX+NU)+NX+j] = u_end[j];
}

/* End acados code */

esp_err_t event_handler(void *ctx, system_event_t *event)
{
    return ESP_OK;
}

void main_memory_task(void *pv)
{
  // function to test memory size
  // result on board Nano32 running ESP32:
  //  More stack size == less heap size for malloc(). Their sum is about 180 kB.

  unsigned long occupy_heap_size=1024;

  // void *dump_heap_size = malloc(occupy_heap_size); // for debug 115111: already big
  while(1) {
    int *dump_heap_size = malloc(occupy_heap_size);
    if(dump_heap_size == NULL) {
      printf("Pointer %p, failed to allocate %lu bytes.\n", dump_heap_size, occupy_heap_size);
      break;
    }
    else {
    printf("Pointer %p, Allocated %lu bytes.\n", dump_heap_size, occupy_heap_size);
    free(dump_heap_size);
    occupy_heap_size += 1024;
    }
  }
}

void main_task(void *pv)
{
    int loopnumber = 0; // for debug

    /* Begin acados code */
    // Problem data
    int_t   N                   = NN;
    real_t  x0[NX]              = {0.5, 0};
    real_t  w[NN*(NX+NU)+NX]    = {0};  // States and controls stacked
    real_t  Q[NX*NX]            = {0};
    real_t  R[NU*NU]            = {0};
    real_t  xref[NX]            = {0};
    real_t  uref[NX]            = {0};
    int_t   max_sqp_iters       = 1;
    int_t   max_iters           = 10;
    real_t  x_end[NX]           = {0};
    real_t  u_end[NU]           = {0};

    for (int_t i = 0; i < NX; i++) w[i] = x0[i];
    for (int_t i = 0; i < NX; i++) Q[i*(NX+1)] = 1.0;
    for (int_t i = 0; i < NU; i++) R[i*(NU+1)] = 0.05;

    // Integrator structs
    real_t T = 0.1;
    sim_in  sim_in;
    sim_out sim_out;
    sim_in.nSteps = 10;
    sim_in.step = T/sim_in.nSteps;
    sim_in.VDE_forw = &VDE_fun;
    sim_in.nx = NX;
    sim_in.nu = NU;

    sim_in.sens_forw = true;
    sim_in.sens_adj = false;
    sim_in.sens_hess = false;
    sim_in.nsens_forw = NX+NU;

    sim_in.x = malloc(sizeof(*sim_in.x) * (NX));
    sim_in.u = malloc(sizeof(*sim_in.u) * (NU));
    sim_in.S_forw = malloc(sizeof(*sim_in.S_forw) * (NX*(NX+NU)));
    for (int_t i = 0; i < NX*(NX+NU); i++) sim_in.S_forw[i] = 0.0;
    for (int_t i = 0; i < NX; i++) sim_in.S_forw[i*(NX+1)] = 1.0;

    sim_out.xn = malloc(sizeof(*sim_out.xn) * (NX));
    sim_out.S_forw = malloc(sizeof(*sim_out.S_forw) * (NX*(NX+NU)));

    sim_info erk_info;
    sim_out.info = &erk_info;

    sim_RK_opts rk_opts;
    sim_erk_create_arguments(&rk_opts, 4);
    void *erk_work;
    int_t work_space_size = sim_erk_calculate_workspace_size(&sim_in, &rk_opts);
    erk_work = (void *) malloc(work_space_size);

    int_t nx[NN+1] = {0};
    int_t nu[NN] = {0};
    int_t nb[NN+1] = {0};
    int_t nc[NN+1] = {0};
    for (int_t i = 0; i < N; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[N] = NX;

    real_t *pA[N];
    real_t *pB[N];
    real_t *pb[N];
    real_t *pQ[N+1];
    real_t *pS[N];
    real_t *pR[N];
    real_t *pq[N+1];
    real_t *pr[N];
    real_t *px[N+1];
    real_t *pu[N];
    real_t *ppi[N];
    real_t *plam[N+1];
    real_t *px0[1];
    d_zeros(&px0[0], nx[0], 1);
    for (int_t i = 0; i < N; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        d_zeros(&ppi[i], nx[i], 1);
        d_zeros(&plam[i], nb[i]+nc[i], 1);
    }
    d_zeros(&pq[N], nx[N], 1);
    d_zeros(&px[N], nx[N], 1);
    d_zeros(&plam[N], nb[N]+nc[N], 1);

    /* End acados code */

    gpio_set_direction(GPIO_NUM_16, GPIO_MODE_OUTPUT);
    int level = 0;
    while (true) {
        gpio_set_level(GPIO_NUM_16, level);
        level = !level;
        vTaskDelay(100 / portTICK_PERIOD_MS);
        printf("\n\n New loop %d\n\n",++loopnumber);

/* Begin acados code */
// Allocate OCP QP variables
ocp_qp_in qp_in;
qp_in.N = N;
ocp_qp_out qp_out;
ocp_qp_condensing_qpoases_args args;
real_t *work = NULL;
qp_in.nx = nx;
qp_in.nu = nu;
qp_in.nb = nb;
qp_in.nc = nc;
for (int_t i = 0; i < N; i++) {
    pQ[i] = Q;
    pR[i] = R;
}
pQ[N] = Q;
qp_in.Q = (const real_t **) pQ;
qp_in.S = (const real_t **) pS;
qp_in.R = (const real_t **) pR;
qp_in.q = (const real_t **) pq;
qp_in.r = (const real_t **) pr;
qp_in.A = (const real_t **) pA;
qp_in.B = (const real_t **) pB;
qp_in.b = (const real_t **) pb;
qp_in.lb = (const real_t **) px0;
qp_out.x = px;
qp_out.u = pu;
qp_out.pi = ppi;
qp_out.lam = plam;
/* End acados code */

printf("Free heap size before initializing qpoases: %d\n",esp_get_free_heap_size()); // for debug

/* Begin acados code */
acado_timer timer;
real_t total_time = 0;
acado_tic(&timer);
for (int_t iter = 0; iter < max_iters; iter++) {
    // printf("\n------ ITERATION %d ------\n", iter);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iters; sqp_iter++) {
        for (int_t i = 0; i < N; i++) {
            // Pass state and control to integrator
            for (int_t j = 0; j < NX; j++) sim_in.x[j] = w[i*(NX+NU)+j];
            for (int_t j = 0; j < NU; j++) sim_in.u[j] = w[i*(NX+NU)+NX+j];
            sim_erk(&sim_in, &sim_out, &rk_opts, 0, erk_work);
            // Construct QP matrices
            for (int_t j = 0; j < NX; j++) {
                pq[i][j] = Q[j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
            }
            for (int_t j = 0; j < NU; j++) {
                pr[i][j] = R[j*(NX+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
            }
            for (int_t j = 0; j < NX; j++) {
                pb[i][j] = sim_out.xn[j] - w[(i+1)*(NX+NU)+j];
                for (int_t k = 0; k < NX; k++) pA[i][j*NX+k] = sim_out.S_forw[j*(NX)+k];
            }
            for (int_t j = 0; j < NU; j++)
                for (int_t k = 0; k < NX; k++) pB[i][j*NX+k] = sim_out.S_forw[NX*NX + NX*j+k];
        }
        for (int_t j = 0; j < NX; j++) {
            px0[0][j] = (x0[j]-w[j]);
        }
        for (int_t j = 0; j < NX; j++) {
            pq[N][j] = Q[j*(NX+1)]*(w[N*(NX+NU)+j]-xref[j]);
        }
        int status = ocp_qp_condensing_qpoases(&qp_in, &qp_out, &args, NULL, work);
        if (status) {
            printf("qpOASES returned error status %d\n", status);
            return -1;
        }
        for (int_t i = 0; i < N; i++) {
            for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] += qp_out.x[i][j];
            for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] += qp_out.u[i][j];
        }
        for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] += qp_out.x[N][j];
    }
    for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
    shift_states(w, x_end, N);
    shift_controls(w, u_end, N);
}
#ifdef DEBUG
print_states_controls(&w[0], N);
#endif  // DEBUG
total_time = acado_toc(&timer);  // in seconds
printf("Average of %.3f ms per iteration.\n", 1e3*total_time/max_iters);
/* End acados code */

printf("Free heap size: %d\n",esp_get_free_heap_size()); // for debug

    }
}

void app_main(void)
{
    nvs_flash_init();
    // tcpip_adapter_init();
    // ESP_ERROR_CHECK( esp_event_loop_init(event_handler, NULL) );
    // wifi_init_config_t cfg = WIFI_INIT_CONFIG_DEFAULT();
    // ESP_ERROR_CHECK( esp_wifi_init(&cfg) );
    // ESP_ERROR_CHECK( esp_wifi_set_storage(WIFI_STORAGE_RAM) );
    // ESP_ERROR_CHECK( esp_wifi_set_mode(WIFI_MODE_STA) );
    // wifi_config_t sta_config = {
    //     .sta = {
    //         .ssid = "access_point_name",
    //         .password = "password",
    //         .bssid_set = false
    //     }
    // };
    // ESP_ERROR_CHECK( esp_wifi_set_config(WIFI_IF_STA, &sta_config) );
    // ESP_ERROR_CHECK( esp_wifi_start() );
    // ESP_ERROR_CHECK( esp_wifi_connect() );
    xTaskCreate(main_task, "main", 10*1024, NULL, 5, NULL);
    // xTaskCreate(main_memory_task, "main_memory", 10*1024, NULL, 5, NULL);
}
