/*------------------------------------------------------------------------------------------------*/
/*  DFTB+: general package for performing fast atomistic simulations                              */
/*  Copyright (C) 2006 - 2024  DFTB+ developers group                                             */
/*                                                                                                */
/*  See the LICENSE file for terms of usage and distribution.                                     */
/*------------------------------------------------------------------------------------------------*/

#include <assert.h>
#include <complex.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "dftbplus.h"
#include "testhelpers.h"

// Dimension of the dense matrices
int basis_size;

// Number of spin channels
int n_spin;

// Number of k-points
int n_kpts;

// Is this a real or complex set of matrices
_Bool is_hs_real;

// Has the density matrix been read yet for each k/spin
_Bool *is_dm_stored;

double *dm_real = NULL;
double complex *dm_cmplx = NULL;
double *s_real = NULL;
double complex *s_cmplx = NULL;
// Not storing these variables in this code:
//double *h_real = NULL;
//double complex *h_cmplx = NULL;


/*-----------------------------------------*/
/*  Callback functions to attach to DFTB+  */
/*-----------------------------------------*/

/*
Access memory for hamiltonian inside DFTB+ and modify H in a simple
way by adding the overlap times the diagonal of the density
matrix. For a self-consistent calculation, this modification prevent
the total energy from being calculated, but the band energy is
available.
Note: we're not storing H in this routine, so just working with the
pointer to the current iK and iS part of H.
*/
void h_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data) {

  printf("Hamiltonian matrix\n");

  // Offset due to multiple k/s elements
  int offset = ((iS - 1) * n_kpts + iK - 1 ) * basis_size * basis_size;

  /* Note: relying on S being available before H in callback, but the
     DM is available only after the iteration, so will be used in the
     next iteration if SCC */

  if (is_hs_real != false)  {

    double *h_local = blacs_data;

    // Print the H matrix in DFTB+
    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
        // h_local[i * basis_size + j] = hamiltonian_real[i][j];  // write to matrix in DFTB+
        // h_real[i*basis_size + j] = h_local[i*basis_size + j]; // read from matrix in by DFTB+
        double h = (i < j) ? h_local[i*basis_size + j] : h_local[j*basis_size + i];
        printf("%9.6f ", h);
      }
      printf("\n");
    }

    // Modify relevant triangle of H
    if (is_dm_stored[(iS-1) * n_kpts + iK -1] != false) {
      for (int i = 0; i < basis_size; ++i) {
        for (int j = i; j < basis_size; ++j) {
          // Modify H by adding a small amount of diagonal density
          // matrix weighted overlap to it:
          h_local[i*basis_size + j] = 0.95 * h_local[i*basis_size + j]
            + 0.05 * s_real[i*basis_size + j + offset] * 0.5
            * ( dm_real[i*basis_size + i + offset] + dm_real[j*basis_size + j + offset] );
        }
      }
    }

  } else {

    double complex *h_local = blacs_data;

    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
        // h_local[i * basis_size + j] = hamiltonian_real[i][j];  // write to matrix in DFTB+
        // h_cmplx[i*basis_size + j] = h_local[i*basis_size + j]; // read from matrix in DFTB+

        // As only one triangle of matrix is populated inside DFTB+:
        double complex h = (i < j) ? h_local[i*basis_size + j] : conj(h_local[j*basis_size + i]);
        printf("(%9.6f,%9.6f) ", creal(h), cimag(h));
      }
      printf("\n");
    }

    // Modify relevant triangle of H
    if (is_dm_stored[(iS-1) * n_kpts + iK -1] != false) {
      for (int i = 0; i < basis_size; ++i) {
        for (int j = i; j < basis_size; ++j) {
          // Modify H by adding a small amount of diagonal density
          // matrix weighted overlap to it:
          h_local[i*basis_size + j] = 0.95 * h_local[i*basis_size + j]
            + 0.05 * s_cmplx[i*basis_size + j + offset] * 0.5
            * ( dm_cmplx[i*basis_size + i + offset] + dm_cmplx[j*basis_size + j + offset] );
        }
      }
    }

  }

}


/* Access memory for overlap inside DFTB+ and take a local copy. Note:
 DFTB+ will trigger this routine once for each spin and k-combination
 at each new geometry, but not neccessarily inside self-consistent
 cycles */
void s_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data) {

  int s_offset = ((iS - 1) * n_kpts + iK - 1 ) * basis_size * basis_size;

  printf("Overlap matrix at ik=%d is=%d\n", iK, iS);
  if (is_hs_real != false)  {

    double *s_local = blacs_data;

    for (int i = 0; i < basis_size; ++i) {
      for (int j = i; j < basis_size; ++j) {
	// relevant triangle from matrix in DFTB+
	s_real[i*basis_size + j + s_offset] = s_local[i*basis_size + j];
      }
    }
    // print both triangles
    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
        double s = (i < j) ? s_real[i*basis_size + j + s_offset] : s_real[j*basis_size + i + s_offset];
        printf("%9.6f ", s);
      }
      printf("\n");
    }

  } else {

    double complex *s_local = blacs_data;

    for (int i = 0; i < basis_size; ++i) {
      for (int j = i; j < basis_size; ++j) {
	// relevant triangle from matrix in DFTB+
	s_cmplx[i*basis_size + j + s_offset] = s_local[i*basis_size + j];
      }
    }
    // print both triangles
    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
	// Elements from other triangle by conjugation:
        double complex s = (i < j) ? s_cmplx[i*basis_size + j + s_offset] : conj(s_cmplx[j*basis_size + i + s_offset]);
        printf("(%9.6f,%9.6f) ", creal(s), cimag(s));
      }
      printf("\n");
    }

  }

}


/* Access memory for density matrix inside DFTB+ and take a local
 copy. Note: DFTB+ will trigger this routine once for each spin and
 k-combination after solving H matrix */
void dm_callback(void *aux_ptr, int iK, int iS, int *blacs_descr,
                void *blacs_data) {

  int dm_offset = ((iS - 1) * n_kpts + iK - 1 ) * basis_size * basis_size;

  printf("Density matrix at ik=%d is=%d\n", iK, iS);
  if (is_hs_real != false)  {

    double *dm_local = blacs_data;

    for (int i = 0; i < basis_size; ++i) {
      for (int j = i; j < basis_size; ++j) {
	// relevant triangle from matrix in DFTB+
	dm_real[i*basis_size + j + dm_offset] = dm_local[i*basis_size + j];
      }
    }
    // print both triangles
    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
        double s = (i < j) ? dm_real[i*basis_size + j + dm_offset] : dm_real[j*basis_size + i + dm_offset];
        printf("%9.6f ", s);
      }
      printf("\n");
    }

  } else {

    double complex *dm_local = blacs_data;

    for (int i = 0; i < basis_size; ++i) {
      for (int j = i; j < basis_size; ++j) {
	// relevant triangle from matrix in DFTB+
	dm_cmplx[i*basis_size + j + dm_offset] = dm_local[i*basis_size + j];
      }
    }
    for (int i = 0; i < basis_size; ++i) {
      for (int j = 0; j < basis_size; ++j) {
        // Elements from other triangle by conjugation:
        double complex s = (i < j) ? dm_cmplx[i*basis_size + j + dm_offset] : conj(dm_cmplx[j*basis_size + i + dm_offset]);
        printf("(%9.6f,%9.6f) ", creal(s), cimag(s));
      }
      printf("\n");
    }

  }

  /* Now have the density matrix for this k-point and spin
     combination: */
  is_dm_stored[(iS-1) * n_kpts + iK -1]  = true;

}


/*--------------------------------------------------------*/
/*  Main code to invoke DFTB+ and access matrix elements  */
/*--------------------------------------------------------*/
int main() {
  FILE *atf = fopen("autotest.tag", "w+");

  DftbPlus calculator;
  DftbPlusInput input;

  int major, minor, patch;
  dftbp_api(&major, &minor, &patch);
  printf("API version %d.%d.%d\n", major, minor, patch);

  _Bool instsafe = dftbp_is_instance_safe();
  printf(instsafe ? "API is instance safe\n" : "API is NOT instance safe\n");

  /* Initialize DFTB+ input tree from input in external file */
  dftbp_init(&calculator, NULL);
  dftbp_get_input_from_file(&calculator, "dftb_in.hsd", &input);
  dftbp_process_input(&calculator, &input);

  is_hs_real = dftbp_is_hs_real(&calculator);
  if (is_hs_real != false) { // Note, can't guarantee all compilers make true == 1, but false == 0 is safe
    printf("Real calculation\n");
    // workaround for Intel(R) oneAPI 2024.0, as otherwise get value
    // outside bool {0,1} for following print:
    printf("is_hs_real : %d\n", is_hs_real);
    fprintf(atf, "is_hs_real       :integer:0:\n%d\n", 1);
  } else {
    printf("Complex calculation\n");
    fprintf(atf, "is_hs_real       :integer:0:\n%d\n", 0);
  }

  basis_size = dftbp_get_basis_size(&calculator);
  printf("Matrix size %d x %d\n", basis_size, basis_size);
  fprintf(atf, "basis_size       :integer:0:\n%d\n", basis_size);

  n_spin = dftbp_get_nr_spin(&calculator);
  printf("Spin channels %d\n", n_spin);
  fprintf(atf, "n_spin       :integer:0:\n%d\n", n_spin);

  n_kpts = dftbp_nr_kpoints(&calculator);
  printf("K-points %d\n", n_kpts);
  fprintf(atf, "n_kpts       :integer:0:\n%d\n", n_kpts);

  is_dm_stored = malloc(n_spin * n_kpts * sizeof(bool));
  for (int i = 0; i < n_spin * n_kpts; ++i) {
    is_dm_stored[i] = false;
  }

  if (is_hs_real != false) {
    dm_real = calloc( basis_size * basis_size * n_spin * n_kpts, sizeof(double)); // storage for multiple spin or k
    // h_real = calloc( basis_size * basis_size, sizeof(double)); // not going to store this
    s_real = calloc( basis_size * basis_size * n_spin * n_kpts, sizeof(double));
  } else {
    dm_cmplx = calloc( basis_size * basis_size * n_spin * n_kpts, sizeof(double complex)); // storage for multiple spin or k
    // h_cmplx = calloc( basis_size * basis_size, sizeof(double complex)); // not going to store this
    s_cmplx = calloc( basis_size * basis_size * n_spin * n_kpts, sizeof(double complex));
  }

  dftbp_register_s_callback(&calculator, s_callback, 0);
  dftbp_register_h_callback(&calculator, h_callback, 0);
  dftbp_register_dm_callback(&calculator, dm_callback, 0);

  /* Request energy, which will invoke s_callback, h_callback and then
   dm_callback. If DFTB+ is running SCC, this will be internally
   iterative, but the returned energy will not be correct for the
   (modified) hamiltonian */
  double mermin_energy;
  dftbp_get_energy(&calculator, &mermin_energy);
  printf("mermin_energy : %9.6f\n", mermin_energy);
  fprintf(atf, "mermin_energy       :real:0:\n%f\n", mermin_energy);

  /*  Clean up */
  dftbp_final(&calculator);

  if (is_hs_real != false) {
    free(dm_real);
    //free(h_real);
    free(s_real);
  } else {
    free(dm_cmplx);
    //free(h_cmplx);
    free(s_cmplx);
  }

  fclose(atf);
  return 0;

}
