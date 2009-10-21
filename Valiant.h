/*
 * Valiant.h
 *
 *  Created on: Sep 29, 2009
 *      Author: yye00
 */

#ifndef VALIANT_H_
#define VALIANT_H_

static char help[] = "Valiant, the Ensemble Kalman filter and Ensemble Optimization Library/executable.\n\n";

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "petscda.h"
#include "petscksp.h"
#include "petscdmmg.h"

#define VALIANT_DEBUG             1
#define VALIANT_USE_X_VIEWER      0
#define VALIANT_USE_ASCII_VIEWER  0
#define VALIANT_USE_MATLAB_VIEWER 1

#if defined(PETSC_HAVE_MATLAB_ENGINE)
#undef VALIANT_USE_MATLAB_VIEWER
#define VALIANT_USE_MATLAB_VIEWER 1
#else
#undef VALIANT_USE_MATLAB_VIEWER
#define VALIANT_USE_MATLAB_VIEWER 0
#endif

/* Some Math Stuff */
#define PI 3.14159265358979323846264338327
#define EPSILON 1e-30
#define  ABS(a) ((a) < 0.0 ? -(a) : (a))

typedef struct {
  Mat Le, HLe, HLeT;
  Mat HLeHLeT, HLeHLeTCD, HLeHLeTCDInv;
  Mat CD;
  Mat Ke;

  /* Vector pointers */
  Vec ** EnsembleVecs;         /* Double pointer pointing to individual vectors */
  Vec * EnsembleObservations;  /* Pointer to observations of each ensemble */
  Vec * EnsembleColumn;        /* Pointer to scattered-gathered big vector: [Phi, LnK11, LnK22, LnK33..]^T */

  /* Ensemble observation error vectors */
  /* This is epsilon j, only need one   */
  Vec EnsembleObsError;

  /* Observations */
  Vec Observations;

  /* Mean vectors */
  Vec Mean, ObsMean;

  /* Parameters */
  PetscInt NumberOfVecsPerEnsemble;
  PetscInt NumberOfEnsembles;

  /* File locations */
  char *SimPathPrefix;
  char *ObsPathPrefix;

} PEnKFRun;

/* Create and Destroy functions */
extern PetscErrorCode ValiantPEnKFCreate(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantPEnKFDestroy(PEnKFRun *MyPEnKF);

/* Load vectors from Defiant */
extern PetscErrorCode ValiantDefiant2PhLoadVecs(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantDefiant3PhLoadVecs(PEnKFRun *MyPEnKF);

/* Load ADCIRC vector, just to freak out people */
extern PetscErrorCode ValiantDefiantADCIRCLoadVecs(PEnKFRun *MyPEnKF);

/* Write vectors back to Defiant */
extern PetscErrorCode ValiantDefiant2PhWriteVecs(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantDefiant3PhWriteVecs(PEnKFRun *MyPEnKF);

/* Scatter and Gather functions */
extern PetscErrorCode ValiantPEnKFScatterForward(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantPEnKFScatterReverse(PEnKFRun *MyPEnKF);

/* Perform the actual data assimilation */
extern PetscErrorCode ValiantPEnKFAssimilate(PEnKFRun *MyPEnKF);

/* Calculate the measurement error covariance matrix */
extern PetscErrorCode ValiantPEnKFComputeCD(PEnKFRun *MyPEnKF);

/* Perturb the forecast measurement errors with this function */
extern PetscErrorCode ValiantPEnKFPerturb(PEnKFRun *MyPEnKF);

/* Calculate the CD matrix with small % random number */
extern PetscErrorCode ValiantPEnKFComputeCDRandomPercentage(PEnKFRun *MyPEnKF);

/* Inverse normal perturbation function. Crude but effective */
extern PetscErrorCode ValiantPEnKFPerturbNIG(PEnKFRun *MyPEnKF);

/* Aggregate functions */
extern PetscErrorCode ValiantPEnKF2Ph(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantPEnKF3Ph(PEnKFRun *MyPEnKF);

#endif /* VALIANT_H_ */
