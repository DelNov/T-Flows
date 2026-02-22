#include <cstdio>
#include <fstream>
#include <iterator>
#include <string>
#if T_FLOWS_AMGX == 1
#include <amgx_c.h>
#include <cuda_runtime.h>
#include "Amgx_Shared.h"
#endif

/*-----------------------------------------------------------------------------+
|                                                                              |
|  AMGX is defined                                                             |
|                                                                              |
+-----------------------------------------------------------------------------*/
#if T_FLOWS_AMGX == 1

/*------------------------------------------------------------------+
|  Global variables.  Now I know this is a bad pracitice, but they  |
|  are global just inside the C++ space, which is rather limitied.  |
+------------------------------------------------------------------*/

/* Plain C/C++ pointers */
int    * amgx_a_row;
int    * amgx_a_col;
double * amgx_a_val;
double * amgx_x;
double * amgx_b;

/* AMGX handles (in order in which they are defined) */
AMGX_config_handle    amgx_conf_hand;
AMGX_resources_handle amgx_rsrc_hand;
AMGX_matrix_handle    amgx_a_hand;
AMGX_vector_handle    amgx_x_hand;
AMGX_vector_handle    amgx_b_hand;
AMGX_solver_handle    amgx_solv_hand;

/* AMGX precision
   dDDI means: (d)evice
               (D)ouble precision matrix  (alternative: (F)loat, 32 bit)
               (D)ouble precision vector  (alternative: (F)loat, 32 bit)
               (I)nteger indices          (alternative: (L)ong,  64 bit)  */
const AMGX_Mode amgx_precision = AMGX_mode_dDDI;

/*============================================================================*/
static void check_amgx(AMGX_RC rc, const char * where) {
/*----------------------------------------------------------------------------*/

  if(rc != AMGX_RC_OK) {
    fprintf(stderr, "AMGX error in %s: rc = %d\n", where, rc);
    exit(EXIT_FAILURE);
  }
}

/*============================================================================*/
static void amgx_silent_print(const char* msg, int length) {
/*-----------------------------------------------------------------------------+
|  Suppress message from AMGX                                                  |
+-----------------------------------------------------------------------------*/

  (void)msg; (void)length;  // swallow everything
}

/*============================================================================*/
extern "C" int amgx_solver_(
    const int & verbose,        // AMGX_QUIET or AMGX_LOUD
    const int & call_count,     // AMGX_FIRST_CALL, AMGX_FINAL_CALL
    const int & matrix_state,   // AMGX_MATRIX_SAME or AMGX_MATRIX_CHANGED
    const int & N,              // number of unknowns
    const int & nnz,            // number of non-zeros
    void **     a_row_dev_ptr,  // size: N+1
    void **     a_col_dev_ptr,  // size: nnz
    void **     a_val_dev_ptr,  // size: nnz
    const int & S,              // start of amgx_x
    void **     x_dev_ptr,      // size: N
    void **     b_dev_ptr,      // size: N
    int &       iters,          // performed iterations
    double &    res) {          // achieved residual
/*-----------------------------------------------------------------------------+
|  amgx_solver - Minimal Fortran - C++ bridge to NVIDIA AMGX
|
|  This routine is an interface layer between T-Flows (Fortran/OpenACC) and
|  NVIDIA AMGX (C API). The linear system (CSR matrix A and vectors x, b) is
|  assembled and kept on the GPU on the Fortran side; this function only
|  creates AMGX objects, uploads/updates data as needed, runs the solve, and
|  returns iteration count and achieved residual.
|
|  AMGX API objects used here form a small hierarchy and must be created in a
|  specific order:
|
|    1) Configuration   (AMGX_config_handle)
|    2) Resources       (AMGX_resources_handle)  depends on configuration
|    3) Matrix A        (AMGX_matrix_handle)     depends on resources
|    4) Vectors x and b (AMGX_vector_handle)     depend on resources
|    5) Solver          (AMGX_solver_handle)     depends on (resources+config)
|                                                and is set up for a
|                                                particular matrix A
|  The routine itself has four logical phases:
|
|    A) Initialization / creation (call_count == AMGX_FIRST_CALL)
|       - Connect device pointers received from Fortran to C/C++ pointers.
|       - Convert CSR indexing from Fortran’s 1-based to AMGX’s 0-based
|         (done once; the CSR arrays are mutated in-place on the device).
|       - AMGX_initialize(), read amgx.json, create config/resources/matrix/
|         vectors, and create+setup the solver.
|       - Optionally suppress AMGX internal printing via a print callback.
|
|    B) Data upload/update (every call)
|       - Upload current A, x, and b into AMGX objects (GPU-to-GPU).
|         (In this version we upload A every time for simplicity.)
|
|    C) Solve (every call)
|       - AMGX_solver_solve()
|       - Query solve status, number of iterations, and residual norm.
|       - Download solution back into the caller’s device vector x
|         (GPU-to-GPU via AMGX_vector_download).
|
|    D) Finalization / destruction (call_count == AMGX_FINAL_CALL)
|       - Destroy AMGX solver, vectors, matrix, resources, config
|       - AMGX_finalize()
|
|  Notes:
|    - Mode is AMGX_mode_dDDI: device memory, double values (matrix & vector),
|      32-bit integer indices.
|    - The CSR index shift is performed only once because the CSR arrays live
|      on the GPU and persist across calls.
+-----------------------------------------------------------------------------*/

  /*-------------------------------------+
  |                                      |
  |  Just a welcome message, marvelling  |
  |     the success of this coupling     |
  |                                      |
  +-------------------------------------*/
  if(call_count == AMGX_FIRST_CALL) {
    if(verbose > 0) {
      printf("Hello from Calling_Amgx!\n\n");
      printf("Everything works as it should.  Fortan side has everything\n");
      printf("on the device and passes it all to C++ for solving it.\n");
      printf("\nThe AMGX functions I am using are:\n");
      printf("- AMGX_initialize\n");
      printf("- AMGX_get_api_version\n");
      printf("- AMGX_config_create\n");
      printf("- AMGX_resources_create_simple\n");
      printf("- AMGX_matrix_create\n");
      printf("- AMGX_vector_create\n");
      printf("- AMGX_matrix_upload_all\n");
      printf("- AMGX_matrix_replace_coefficients\n");
      printf("- AMGX_vector_upload\n");
      printf("- AMGX_solver_create\n");
      printf("- AMGX_solver_setup\n");
      printf("- AMGX_solver_solve\n");
      printf("- AMGX_solver_get_status\n");
      printf("- AMGX_solver_get_iterations_number\n");
      printf("- AMGX_solver_calculate_residual_norm\n");
      printf("- AMGX_vector_download\n");
      printf("- AMGX_solver_destroy\n");
      printf("- AMGX_matrix_destroy\n");
      printf("- AMGX_vector_destroy\n");
      printf("- AMGX_resources_destroy\n");
      printf("- AMGX_config_destroy\n");
      printf("- AMGX_finalize\n\n");
      printf("and these AMGX's data types:\n");
      printf("- AMGX_config_handle\n");
      printf("- AMGX_resources_handle\n");
      printf("- AMGX_matrix_handle\n");
      printf("- AMGX_vector_handle\n");
      printf("- AMGX_solver_handle\n\n");
    }
  }

  /*-------------------------------------------------+
  |                                                  |
  |  Connect to the linear system of equations from  |
  |  the Fortran side with standard C/C++ pointers.  |
  |                                                  |
  +-------------------------------------------------*/
  if(call_count == AMGX_FIRST_CALL) {

    amgx_a_row = static_cast <int *>    (* a_row_dev_ptr);
    amgx_a_col = static_cast <int *>    (* a_col_dev_ptr);
    amgx_a_val = static_cast <double *> (* a_val_dev_ptr);
    amgx_x     = static_cast <double *> (* x_dev_ptr) + S;
    amgx_b     = static_cast <double *> (* b_dev_ptr);
  }

  /*-----------------------------------------------------------+
  |                                                            |
  |  Shift the indices to match the C/C++ numbering from zero  |
  |                                                            |
  +-----------------------------------------------------------*/
  if(call_count   == AMGX_FIRST_CALL ||
     matrix_state == AMGX_MATRIX_CHANGED) {

    #pragma acc parallel loop deviceptr(amgx_a_row)
    for(int c = 0; c <= N; c++) {
      amgx_a_row[c]--;
    }
    #pragma acc parallel loop deviceptr(amgx_a_col)
    for(int i = 0; i < nnz; i++) {
      amgx_a_col[i]--;
    }
  }

  /*----------------------+
  |                       |
  |  AMGX Initialization  |
  |                       |
  +----------------------*/
  if(call_count == AMGX_FIRST_CALL) {

    check_amgx(AMGX_initialize(), "AMGX_initialize");

    /* Query API version (non-essential, but cute) */
    if(verbose > 0) {
      int api_major = 0;
      int api_minor = 0;
      check_amgx(AMGX_get_api_version(&api_major, &api_minor),
                "AMGX_get_api_version");
      printf("AMGX API version: %d.%d\n", api_major, api_minor);
    }

    /* Control the output from AMGX (non-essential, but very useful) */
    if(verbose == 0) {
      check_amgx(AMGX_register_print_callback(amgx_silent_print),
                "AMGX_register_print_callback");
    }

    /*------------------------+
    |  1. AMGX Configuration  |
    +------------------------*/
    amgx_conf_hand = nullptr;

    /* Read .json file ... */
    std::ifstream in("amgx.json");
    std::string cfg_json((std::istreambuf_iterator<char>(in)),
                          std::istreambuf_iterator<char>());
    if(!in) {
      std::fprintf(stderr, "Could not open amgx.json\n");
      std::exit(EXIT_FAILURE);
    }

    /* ... and create config */
    check_amgx(AMGX_config_create(&amgx_conf_hand, cfg_json.c_str()),
              "AMGX_config_create");

    /*--------------------+
    |  2. AMGX Resources  |  (depend on amgx_conf_hand)
    +--------------------*/
    amgx_rsrc_hand = nullptr;

    // Simple resources: single GPU, device 0
    check_amgx(AMGX_resources_create_simple(&amgx_rsrc_hand, amgx_conf_hand),
              "AMGX_resources_create_simple");

    /*---------------------------------+
    |  3. Create matrix on the device  |  (depends on: amgx_rsrc_hand)
    +---------------------------------*/
    amgx_a_hand = nullptr;

    check_amgx(AMGX_matrix_create(&amgx_a_hand, amgx_rsrc_hand, amgx_precision),
              "AMGX_matrix_create");

    /*----------------------------------+
    |  4. Create vectors on the device  |  (depends on: amgx_rsrc_hand)
    +----------------------------------*/
    amgx_x_hand = nullptr;
    amgx_b_hand = nullptr;

    check_amgx(AMGX_vector_create(&amgx_x_hand, amgx_rsrc_hand, amgx_precision),
              "AMGX_vector_create(x)");
    check_amgx(AMGX_vector_create(&amgx_b_hand, amgx_rsrc_hand, amgx_precision),
               "AMGX_vector_create(b)");

  }  /* call_count == AMGX_FIRST_CALL */

  /*------------------------------------------+
  |                                           |
  |  Upload matrix and vectors to the device  |  (before creating solver)
  |                                           |
  +------------------------------------------*/
  if(call_count == AMGX_FIRST_CALL) {

    check_amgx(AMGX_matrix_upload_all(
               amgx_a_hand, N, nnz, 1, 1,
               amgx_a_row, amgx_a_col, amgx_a_val, nullptr),
              "AMGX_matrix_upload_all");
    if(verbose > 0) printf("Matrix uploaded to AMGX.\n");

  /*--------------------------------------------+
  |                                             |
  |  Replace matrix coefficients on the device  |
  |                                             |
  +--------------------------------------------*/
  } else {

    if(matrix_state == AMGX_MATRIX_CHANGED) {
      check_amgx(AMGX_matrix_replace_coefficients(
                 amgx_a_hand, N, nnz, amgx_a_val, nullptr),
                "AMGX_matrix_replace_coefficients");
      if(verbose > 0) printf("Matrix replaced in AMGX.\n");

      check_amgx(AMGX_solver_setup(amgx_solv_hand, amgx_a_hand),
                "AMGX_solver_setup");
    }
  }

  /* "Upload" host vectors to AMGX' workspace */
  check_amgx(AMGX_vector_upload(amgx_x_hand, N, 1, amgx_x),
            "AMGX_vector_upload(x)");
  check_amgx(AMGX_vector_upload(amgx_b_hand, N, 1, amgx_b),
            "AMGX_vector_upload(b)");

  if(verbose > 0) printf("Vectors x and b uploaded to AMGX.\n");

  if(call_count == AMGX_FIRST_CALL) {
    /*-----------------------+
    |  5. Create the solver  |  (when matrix and vectors are on the device)
    +-----------------------*/

    amgx_solv_hand = nullptr;
    check_amgx(AMGX_solver_create(&amgx_solv_hand,
                                   amgx_rsrc_hand,
                                   amgx_precision,
                                   amgx_conf_hand),
              "AMGX_solver_create");
    check_amgx(AMGX_solver_setup(amgx_solv_hand, amgx_a_hand),
              "AMGX_solver_setup");

  }  /* call_count == AMGX_FIRST_CALL */

  /*-----------------+
  |                  |
  |  Run the solver  |
  |                  |
  +-----------------*/
  check_amgx(AMGX_solver_solve(amgx_solv_hand, amgx_b_hand, amgx_x_hand),
            "AMGX_solver_solve");

  AMGX_SOLVE_STATUS status;
  check_amgx(AMGX_solver_get_status(amgx_solv_hand, &status),
            "AMGX_solver_get_status");

  if(verbose > 0) printf("Solve status: %d (0=SUCCESS)\n", status);

  /* Fetch the number of performed iterations */
  iters = 0;
  check_amgx(AMGX_solver_get_iterations_number(amgx_solv_hand, &iters),
            "AMGX_solver_get_iterations_number");

  if(verbose > 0) printf("Number of iterations: %d\n", iters);

  /* Fetch the number of performed iterations */
  res = 0.0;
  check_amgx(AMGX_solver_calculate_residual_norm(
             amgx_solv_hand, amgx_a_hand,amgx_b_hand,amgx_x_hand, &res),
            "AMGX_solver_calculate_residual_norm");

  /* "Download" solution back to caller's memory space */
  check_amgx(AMGX_vector_download(amgx_x_hand, amgx_x),
            "AMGX_vector_download");

  /*------------------------------------------------------------+
  |                                                             |
  |  Shift the indices to match the Fortran numbering from one  |
  |                                                             |
  +------------------------------------------------------------*/
  if(call_count   == AMGX_FIRST_CALL ||
     matrix_state == AMGX_MATRIX_CHANGED) {

    #pragma acc parallel loop deviceptr(amgx_a_row)
    for(int c = 0; c <= N; c++) {
      amgx_a_row[c]++;
    }
    #pragma acc parallel loop deviceptr(amgx_a_col)
    for(int i = 0; i < nnz; i++) {
      amgx_a_col[i]++;
    }
  }

  /*--------------------+
  |                     |
  |  AMGX Finalization  |
  |                     |
  +--------------------*/

  /* Destroy AMGX objects (in order opposite in which they were defined) */
  if(call_count == AMGX_FINAL_CALL) {
    check_amgx(AMGX_solver_destroy   (amgx_solv_hand), "AMGX_solver_destroy");
    check_amgx(AMGX_vector_destroy   (amgx_b_hand),    "AMGX_vector_destroy");
    check_amgx(AMGX_vector_destroy   (amgx_x_hand),    "AMGX_vector_destroy");
    check_amgx(AMGX_matrix_destroy   (amgx_a_hand),    "AMGX_matrix_destroy");
    check_amgx(AMGX_resources_destroy(amgx_rsrc_hand), "AMGX_resources_destroy");
    check_amgx(AMGX_config_destroy   (amgx_conf_hand), "AMGX_config_destroy");

    check_amgx(AMGX_finalize(), "AMGX_finalize");
  }

  if(verbose > 0) printf("\nGoodbye from Calling_Amgx.\n");

  return 0;
}

/*-----------------------------------------------------------------------------+
|                                                                              |
|  AMGX is not defined                                                         |
|                                                                              |
+-----------------------------------------------------------------------------*/
#else

static void check_amgx(const int rc, const char * where) {}

extern "C" int amgx_solver_(const int & verbose,        // 0 or 1
                            const int & N,
                            const int & nnz,
                            void **     a_row_dev_ptr,  // size: N+1
                            void **     a_col_dev_ptr,  // size: nnz
                            void **     a_val_dev_ptr,  // size: nnz
                            const int & S,              // start of x
                            void **     x_dev_ptr,      // size: N
                            void **     b_dev_ptr,      // size: N
                            int &       iters,          // performed iterations
                            double &    res) {          // achieved residual
  check_amgx(0, " ");
  return 0;
}

#endif
