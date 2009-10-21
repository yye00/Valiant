/*
 * Valiant.c
 *
 *  Created on: Sep 29, 2009
 *      Author: yye00
 */

#include "Valiant.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PEnKFRun SimPEnKF;

  /* Initialize PETSc */
  PetscInitialize(&argc, &argv, (char *) 0, help);
  char SimPrefix[] = "/home/yye00/workspace/Defiant/Simulation_";
  char ObsPrefix[] = "/home/yye00/workspace/Defiant/Observations/";

  SimPEnKF.NumberOfEnsembles = 6;
  SimPEnKF.SimPathPrefix = SimPrefix;
  SimPEnKF.ObsPathPrefix = ObsPrefix;

  /* 2Phase takes six state vecs in an ensemble */
  SimPEnKF.NumberOfVecsPerEnsemble = 6;

  ierr = ValiantPEnKF2Ph(&SimPEnKF);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);

  return 0;
}
