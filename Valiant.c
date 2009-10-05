/*
 * Valiant.c
 *
 *  Created on: Sep 29, 2009
 *      Author: yye00
 */

#include "Valiant.h"

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFCreate"
extern PetscErrorCode ValiantPEnKFCreate(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i, j;

  PetscFunctionBegin;

  /* allocate the pointers to the ensemble vector pointers */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles *sizeof (Vec *),&MyPEnKF->EnsembleVecs);CHKERRQ(ierr);
  /* allocate the pointers to the ensemble observations */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles * sizeof (Vec),&MyPEnKF->EnsembleObservations);CHKERRQ(ierr);
  /* Create the ensemble observation error vectors */
  ierr = PetscMalloc (MyPEnKF->NumberOfEnsembles *sizeof (Vec),&MyPEnKF->EnsembleObsError);CHKERRQ(ierr);
  /* allocate the vectors for the full ensemble vecs */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles *sizeof (Vec),&MyPEnKF->EnsembleColumn);CHKERRQ(ierr);

  /* Now for every ensemble allocate the vecs */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
    ierr = PetscMalloc (NMyPEnKF->umberOfVecsPerEnsemble *sizeof (Vec),&MyPEnKF->EnsembleVecs[i]); CHKERRQ(ierr);

  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* for every ensemble iterate over all the state vectors */
    for(j = 0; j< MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleVecs[i][j]));CHKERRQ(ierr);
    }
    /* create the Observation vectors */
    ierr = VecCreate (PETSC_COMM_WORLD,&(MyPEnKF->EnsembleObservations[i]));CHKERRQ(ierr);
    /* Create the ensemble observation error vectors */
    ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleObsError[i]));CHKERRQ(ierr);
    /* create the Full Ensemble vectors */
    ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleColumn[i]));CHKERRQ(ierr);
  }

  /* Create the  Real world observations vectors */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Observations);CHKERRQ(ierr);
  /* Setup the Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Mean);CHKERRQ(ierr);
  /* Setup the Observations Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->ObsMean);CHKERRQ(ierr);

  /* Allocate memory for the scatters */
  ierr = PetscMalloc (NumberOfVecsPerEnsemble *sizeof (IS), &isVecs);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFDestroy"
extern PetscErrorCode ValiantPEnKFDestroy(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Clean up Matrices */
   ierr = MatDestroy(MyPEnKF->Le);
   ierr = MatDestroy(MyPEnKF->HLe);
   ierr = MatDestroy(MyPEnKF->HLeT);
   ierr = MatDestroy(MyPEnKF->HLeHLeT);
   ierr = MatDestroy(MyPEnKF->EMatrix);
   ierr = MatDestroy(MyPEnKF->EMatrixT);
   ierr = MatDestroy(MyPEnKF->CD);
   ierr = MatDestroy(MyPEnKF->HLeHLeTCDInv);
   ierr = MatDestroy(MyPEnKF->TempI);
   ierr = MatDestroy(MyPEnKF->TempMat);
   ierr = MatDestroy(MyPEnKF->Ke);

   /* Clean up Vectors */
   for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
   {
     for(j = 0; j< MyPEnKF->NumberOfVecsPerEnsemble; j++)
     {
       ierr = VecDestroy (MyPEnKF->EnsembleVecs[i][j]);CHKERRQ(ierr);
     }
     ierr = PetscFree (MyPEnKF->EnsembleVecs[i]);CHKERRQ(ierr);
     ierr = VecDestroy (MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);
     ierr = VecDestroy (MyPEnKF->EnsembleObsError[i]);CHKERRQ(ierr);
     ierr = VecDestroy (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);
   }

   ierr = VecDestroy(MyPEnKF->Observations);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->Mean);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->ObsMean);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->TempVector);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->VecOnes);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
