/*
 * Valiant.h
 *
 *  Created on: Sep 29, 2009
 *      Author: yye00
 */

#ifndef VALIANT_H_
#define VALIANT_H_

static char help[] = "Valiant, the Ensemble Kalman filter and Ensemble Optimization Library/executable.\n\n";

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
  Mat Le, HLe, HLeT, HLeHLeT;
  Mat EMatrix, EMatrixT, CD;
  Mat HLeHLeTCDInv;
  Mat TempI, TempMat;
  Mat Ke;

  /* Vector pointers */
  Vec ** EnsembleVecs;
  Vec * EnsembleObservations;
  Vec * EnsembleColumn;

  /* Ensemble observation error vectors */
  Vec *EnsembleObsError;

  /* Observations */
  Vec Observations;

  /* Mean vectors */
  Vec Mean, ObsMean;

  /* LU factorization variables */
  MatFactorInfo  info;
  IS perm,iperm;

  /* Randcom context */
  PetscRandom rctx;

  /* Error Code */
  PetscErrorCode ierr;

  /* Parameters */
  PetscInt NumberOfVecsPerEnsemble;
  PetscInt NumberOfEnsembles;
  PetscInt SizeOfStateVector;
  PetscInt SizeOfObservations;
  PetscInt NumberOfEntries;

  /* File locations */
  char *FileDirectory;

} PEnKFRun;

/* Create and Destroy functions */
extern PetscErrorCode ValiantPEnKFCreate(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantPEnKFDestroy(PEnKFRun *MyPEnKF);

/* Load vectors from Defiant */
extern PetscErrorCode ValiantDefiant2PhLoadVecs(PEnKFRun *MyPEnKF);
extern PetscErrorCode ValiantDefiant3PhLoadVecs(PEnKFRun *MyPEnKF);

#endif /* VALIANT_H_ */
