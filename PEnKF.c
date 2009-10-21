/*
 * PEnKF.c
 *
 *  Created on: Sep 29, 2009
 *      Author: yye00
 */

#include "Valiant.h"

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKF2Ph"
extern PetscErrorCode ValiantPEnKF2Ph(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ValiantPEnKFCreate(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantDefiant2PhLoadVecs(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFScatterForward(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFAssimilate(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFScatterReverse(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantDefiant2PhWriteVecs(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFDestroy(MyPEnKF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKF3Ph"
extern PetscErrorCode ValiantPEnKF3Ph(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = ValiantPEnKFCreate(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantDefiant3PhLoadVecs(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFScatterForward(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFAssimilate(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFScatterReverse(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantDefiant3PhWriteVecs(MyPEnKF);CHKERRQ(ierr);
  ierr = ValiantPEnKFDestroy(MyPEnKF);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFCreate"
extern PetscErrorCode ValiantPEnKFCreate(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i, j;

  PetscFunctionBegin;

  /* allocate the pointers to the ensemble vector pointers */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles *sizeof (Vec *),&MyPEnKF->EnsembleVecs);CHKERRQ(ierr);CHKMEMQ;
  /* allocate the pointers to the ensemble observations */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles * sizeof (Vec),&MyPEnKF->EnsembleObservations);CHKERRQ(ierr);CHKMEMQ;
  /* Create the ensemble observation error vectors */
  ierr = PetscMalloc (MyPEnKF->NumberOfEnsembles *sizeof (Vec),&MyPEnKF->EnsembleObsError);CHKERRQ(ierr);CHKMEMQ;
  /* allocate the vectors for the full ensemble vecs */
  ierr = PetscMalloc(MyPEnKF->NumberOfEnsembles *sizeof (Vec),&MyPEnKF->EnsembleColumn);CHKERRQ(ierr);CHKMEMQ;

  /* Now for every ensemble allocate the vecs */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
    ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (Vec),&(MyPEnKF->EnsembleVecs[i])); CHKERRQ(ierr);CHKMEMQ;

  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* for every ensemble iterate over all the state vectors */
    for(j = 0; j< MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleVecs[i][j]));CHKERRQ(ierr);CHKMEMQ;
      ierr = VecSetType(MyPEnKF->EnsembleVecs[i][j], VECMPI);CHKERRQ(ierr);CHKMEMQ;
    }
    /* create the Observation vectors */
    ierr = VecCreate (PETSC_COMM_WORLD,&(MyPEnKF->EnsembleObservations[i]));CHKERRQ(ierr);CHKMEMQ;
    ierr = VecSetType(MyPEnKF->EnsembleObservations[i], VECMPI);CHKERRQ(ierr);CHKMEMQ;
    /* create the Full Ensemble vectors */
    ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleColumn[i]));CHKERRQ(ierr);CHKMEMQ;
    ierr = VecSetType(MyPEnKF->EnsembleColumn[i], VECMPI);CHKERRQ(ierr);CHKMEMQ;
  }

  /* Create the  Real world observations vectors */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Observations);CHKERRQ(ierr);CHKMEMQ;
  /* Setup the Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Mean);CHKERRQ(ierr);CHKMEMQ;
  /* Setup the Observations Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->ObsMean);CHKERRQ(ierr);CHKMEMQ;

  /* Now for some matrices */
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->Le);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLe);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeT);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeT);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeTCDInv);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->CD);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->Ke);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFDestroy"
extern PetscErrorCode ValiantPEnKFDestroy(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i, j;

  PetscFunctionBegin;
  /* Clean up Matrices */
   ierr = MatDestroy(MyPEnKF->Le);
   ierr = MatDestroy(MyPEnKF->HLe);
   ierr = MatDestroy(MyPEnKF->HLeT);
   ierr = MatDestroy(MyPEnKF->HLeHLeT);
   ierr = MatDestroy(MyPEnKF->HLeHLeTCD);
   ierr = MatDestroy(MyPEnKF->HLeHLeTCDInv);
   ierr = MatDestroy(MyPEnKF->CD);
   ierr = MatDestroy(MyPEnKF->Ke);

   /* Clean up Vectors */
   for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
   {
     for(j = 0; j< MyPEnKF->NumberOfVecsPerEnsemble; j++)
     {
       ierr = VecDestroy (MyPEnKF->EnsembleVecs[i][j]);CHKERRQ(ierr);CHKMEMQ;
     }
     ierr = PetscFree (MyPEnKF->EnsembleVecs[i]);CHKERRQ(ierr);CHKMEMQ;
     ierr = VecDestroy (MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);CHKMEMQ;
     ierr = VecDestroy (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);CHKMEMQ;
   }

   ierr = VecDestroy(MyPEnKF->Observations);CHKERRQ(ierr);CHKMEMQ;
   ierr = VecDestroy(MyPEnKF->Mean);CHKERRQ(ierr);CHKMEMQ;
   ierr = VecDestroy(MyPEnKF->ObsMean);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantDefiant2PhLoadVecs"
extern PetscErrorCode ValiantDefiant2PhLoadVecs(PEnKFRun *MyPEnKF)
{
  /* Keep in mind for 3 phases: NumberOfVecsPerEnsemble is 6 */
  PetscErrorCode ierr;
  PetscInt i;
  PetscViewer viewer;

  char PathBuffer[256];
  char FileNameBuffer[256];
  struct stat st;
  PetscFunctionBegin;

  /* For each ensemble */
  for(i=0;i<MyPEnKF->NumberOfEnsembles;i++) {
    /* go to the ensemble simulation directory */
    sprintf(PathBuffer, "%s%05d", MyPEnKF->SimPathPrefix, i);
    /* if we can get into it start loading */
    if(stat(PathBuffer,&st)==0) {
      /* Load the porosity */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Phi.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][0]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][1]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][2]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][3]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][4]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the Saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][5]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the Observation vector */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Observations.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

    } else
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Simulation directory: %s does NOT exist", PathBuffer);CHKERRQ(ierr);CHKMEMQ;
  }

  sprintf(FileNameBuffer, "%s%s", MyPEnKF->ObsPathPrefix, "/FieldObservations.Valiant");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->Observations);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ValiantDefiant3PhLoadVecs"
extern PetscErrorCode ValiantDefiant3PhLoadVecs(PEnKFRun *MyPEnKF)
{
  /* Keep in mind for 3 phases: NumberOfVecsPerEnsemble is 9 */
  PetscErrorCode ierr;
  PetscInt i;
  PetscViewer viewer;

  char PathBuffer[256];
  char FileNameBuffer[256];
  struct stat st;
  PetscFunctionBegin;

  /* For each ensemble */
  for(i=0;i<MyPEnKF->NumberOfEnsembles;i++) {
    /* go to the ensemble simulation directory */
    sprintf(PathBuffer, "%s%05d", MyPEnKF->SimPathPrefix, i);
    /* if we can get into it start loading */
    if(stat(PathBuffer,&st)==0) {
      /* Load the porosity */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Phi.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][0]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][1]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][2]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][3]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the water pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][4]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the oil pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Po.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][5]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the gas pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pg.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][6]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the water saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][7]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the oil saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/So.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][8]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the Observation vector */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Observations.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

    } else
      ierr = PetscPrintf(PETSC_COMM_WORLD, "Simulation directory: %s does NOT exist", PathBuffer);CHKERRQ(ierr);CHKMEMQ;
  }

  sprintf(FileNameBuffer, "%s%s", MyPEnKF->ObsPathPrefix, "/FieldObservations.Valiant");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->Observations);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantDefiant2PhWriteVecs"
extern PetscErrorCode ValiantDefiant2PhWriteVecs(PEnKFRun *MyPEnKF)
{
  /* Keep in mind for 3 phases: NumberOfVecsPerEnsemble is 6 */
  PetscErrorCode ierr;
  PetscInt i;
  PetscViewer viewer;

  char PathBuffer[256];
  char FileNameBuffer[256];
  struct stat st;
  PetscFunctionBegin;

  /* For each ensemble */
  for(i=0;i<MyPEnKF->NumberOfEnsembles;i++) {
    /* go to the ensemble simulation directory */
    sprintf(PathBuffer, "%s%05d", MyPEnKF->SimPathPrefix, i);
    /* if we can get into it start loading */
    if(stat(PathBuffer,&st)==0) {
      /* Load the porosity */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Phi.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][0],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][1],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][2],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][3],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][4],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the Saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][5],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantDefiant3PhWriteVecs"
extern PetscErrorCode ValiantDefiant3PhWriteVecs(PEnKFRun *MyPEnKF)
{
  /* Keep in mind for 3 phases: NumberOfVecsPerEnsemble is 10 */
  PetscErrorCode ierr;
  PetscInt i;
  PetscViewer viewer;

  char PathBuffer[256];
  char FileNameBuffer[256];
  struct stat st;
  PetscFunctionBegin;

  /* For each ensemble */
  for(i=0;i<MyPEnKF->NumberOfEnsembles;i++) {
    /* go to the ensemble simulation directory */
    sprintf(PathBuffer, "%s%05d", MyPEnKF->SimPathPrefix, i);
    /* if we can get into it start loading */
    if(stat(PathBuffer,&st)==0) {
      /* Load the porosity */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Phi.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][0],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][1],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][2],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][3],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the water pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][4],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the oil pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Po.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][5],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the gas pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pg.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][6],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the water saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][7],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;

      /* Load the oil saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/So.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView(MyPEnKF->EnsembleVecs[i][8],viewer);CHKERRQ(ierr);CHKMEMQ;
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);CHKMEMQ;
    }
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFScatterForward"
extern PetscErrorCode ValiantPEnKFScatterForward(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i,j;
  PetscInt ColumnVecGlobalSize;
  PetscInt StateVecGlobalSize, StateVecLocalSize, StateVecLower, StateVecHigher;
  PetscInt ObsVecGlobalSize,   ObsVecLocalSize,   ObsVecLower,     ObsVecHigher;
  PetscInt FullVecGlobalSize, FullVecLocalSize,   FullVecLower,   FullVecHigher;
  PetscInt *indices, *RowIndices, *FullRowIndices;
  /* Vector scatters */
  IS *isVecs, isObs, intoObs;
  VecScatter VecToEnsembleVec, ObsToEnsembleVec;

  PetscFunctionBegin;

  /* Retrieve global, local sizes and extents of the state vectors */
  ierr = VecAssemblyBegin(MyPEnKF->EnsembleVecs[0][0]);
  ierr = VecAssemblyEnd(MyPEnKF->EnsembleVecs[0][0]);
  ierr = VecGetSize (MyPEnKF->EnsembleVecs[0][0], &StateVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleVecs[0][0], &StateVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleVecs[0][0], &StateVecLower, &StateVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Allocate memory for the scatters */
  ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (IS), &isVecs);CHKERRQ(ierr);CHKMEMQ;

  /* For every vector, create the scatter into the ensemble vector */
  for( i =0; i < MyPEnKF->NumberOfVecsPerEnsemble; i++)
  {
    ierr = PetscMalloc (StateVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
    for (j = 0; j < StateVecLocalSize; j++)
    {
      indices[j] = i*StateVecGlobalSize + StateVecLower + j;
    }
    ierr = ISCreateGeneral (PETSC_COMM_WORLD, StateVecLocalSize, indices,&isVecs[i]);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisVecs is:\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = ISView(isVecs[i], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
  }

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize + ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&isObs);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisObs\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = ISView (isObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&intoObs);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nintoObs\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = ISView (intoObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Set the size of the column vector, assuming we already loaded the state vectors and observations */
  ColumnVecGlobalSize = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize+ObsVecGlobalSize;
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    ierr = VecSetSizes (MyPEnKF->EnsembleColumn[i], PETSC_DECIDE, ColumnVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecSetFromOptions (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);CHKMEMQ;
  }

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);CHKMEMQ;
  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);CHKMEMQ;
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &RowIndices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    RowIndices[i] = ObsVecLower + i;

  /* For all ensembles, perform scatter */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    for(j = 0; j < MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecScatterCreate (MyPEnKF->EnsembleVecs[i][j], isVecs[0],MyPEnKF->EnsembleColumn[i],isVecs[j], &VecToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nVecToEnsembleVec\n");CHKERRQ(ierr);CHKMEMQ;
      ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
      ierr = VecScatterBegin (VecToEnsembleVec, MyPEnKF->EnsembleVecs[i][j],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecScatterEnd (VecToEnsembleVec, MyPEnKF->EnsembleVecs[i][j],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);CHKMEMQ;

      /* Be nice and cleanup before leaving */
      ierr = VecScatterDestroy (VecToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
    }

    ierr = VecScatterCreate (MyPEnKF->EnsembleObservations[i], intoObs,MyPEnKF->EnsembleColumn[i],isObs, &ObsToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nObsToEnsembleVec\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
    ierr = VecScatterBegin (ObsToEnsembleVec, MyPEnKF->EnsembleObservations[i],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecScatterEnd (ObsToEnsembleVec, MyPEnKF->EnsembleObservations[i],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);CHKMEMQ;
    /* Be nice and cleanup before leaving */
    ierr = VecScatterDestroy (ObsToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "Ensemble vector after scatter\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = VecView (MyPEnKF->EnsembleColumn[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
  }
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFScatterReverse"
extern PetscErrorCode ValiantPEnKFScatterReverse(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i,j;
  PetscInt ColumnVecGlobalSize;
  PetscInt StateVecGlobalSize, StateVecLocalSize, StateVecLower, StateVecHigher;
  PetscInt ObsVecGlobalSize,   ObsVecLocalSize,   ObsVecLower,     ObsVecHigher;
  PetscInt FullVecGlobalSize, FullVecLocalSize,   FullVecLower,   FullVecHigher;
  PetscInt *indices, *RowIndices, *FullRowIndices;
  /* Vector scatters */
  IS *isVecs, isObs, intoObs;
  VecScatter VecToEnsembleVec, ObsToEnsembleVec;

  PetscFunctionBegin;

  /* Retrieve global, local sizes and extents of the state vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleVecs[0][0], &StateVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleVecs[0][0], &StateVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleVecs[0][0], &StateVecLower, &StateVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Allocate memory for the scatters */
  ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (IS), &isVecs);CHKERRQ(ierr);CHKMEMQ;

  /* For every vector, create the scatter into the ensemble vector */
  for( i =0; i < MyPEnKF->NumberOfVecsPerEnsemble; i++)
  {
    ierr = PetscMalloc (StateVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
    for (j = 0; j < StateVecLocalSize; j++)
    {
      indices[j] = i*StateVecGlobalSize + StateVecLower + j;
    }
    ierr = ISCreateGeneral (PETSC_COMM_WORLD, StateVecLocalSize, indices,&isVecs[i]);CHKERRQ(ierr);CHKMEMQ;
    ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisVecs is:\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = ISView(isVecs[i], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
  }

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize + ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&isObs);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisObs\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = ISView (isObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&intoObs);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscFree (indices);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nintoObs\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = ISView (intoObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);CHKMEMQ;
  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);CHKMEMQ;
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &RowIndices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    RowIndices[i] = ObsVecLower + i;

  /* Scatter in reverse mode */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    for(j = 0; j < MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecScatterCreate (MyPEnKF->EnsembleVecs[i][j], isVecs[0],MyPEnKF->EnsembleColumn[i],isVecs[j], &VecToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
#ifdef PENKF_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nEnsembleVec To Vec\n");CHKERRQ(ierr);CHKMEMQ;
      ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
      ierr = VecScatterBegin (VecToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleVecs[i][j],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);CHKMEMQ;
      ierr = VecScatterEnd (VecToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleVecs[i][j],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);CHKMEMQ;

      /* Be nice and cleanup before leaving */
      ierr = VecScatterDestroy (VecToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
#ifdef PENKF_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nPost Data Assimilation Vectors\n");CHKERRQ(ierr);CHKMEMQ;
      ierr = VecView (MyPEnKF->EnsembleVecs[i][j], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
    }
    ierr = VecScatterCreate (MyPEnKF->EnsembleObservations[i], intoObs,MyPEnKF->EnsembleColumn[i],isObs, &ObsToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
#ifdef PENKF_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nEnsembleVec To VecObservation\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
    ierr = VecScatterBegin (ObsToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleObservations[i],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecScatterEnd (ObsToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleObservations[i],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);CHKMEMQ;
    /* Be nice and cleanup before leaving */
    ierr = VecScatterDestroy (ObsToEnsembleVec);CHKERRQ(ierr);CHKMEMQ;
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFAssimilate"
extern PetscErrorCode ValiantPEnKFAssimilate(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscInt i, size;
  PetscInt FullVecLocalSize, FullVecGlobalSize;
  PetscInt FullVecLower, FullVecHigher;
  PetscInt ObsVecLocalSize, ObsVecGlobalSize;
  PetscInt ObsVecLower, ObsVecHigher;
  PetscInt *FullRowIndices, *ObsRowIndices;

  PetscScalar TempScalar;
  PetscScalar * avec;

  /* LU factorization variables */
  MatFactorInfo  info;
  IS perm,iperm;

  Vec VecOnes, TempVector;
  Mat TempI, TempMat;

  PetscFunctionBegin;
  /* Ge the the size, make things easier */
  ierr = MPI_Comm_size (PETSC_COMM_WORLD, &size);CHKERRQ(ierr);CHKMEMQ;

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);CHKMEMQ;

  /* Create the Le matrix */
  ierr = MatSetSizes (MyPEnKF->Le, FullVecLocalSize, PETSC_DECIDE, FullVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetType (MyPEnKF->Le, MATDENSE);CHKERRQ(ierr);CHKMEMQ;

  /* Create the HLe matrix */
  ierr = MatSetSizes (MyPEnKF->HLe, ObsVecLocalSize, PETSC_DECIDE, ObsVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetType (MyPEnKF->HLe, MATAIJ);CHKERRQ(ierr);CHKMEMQ;

  /* Setup the mean vectors */
  ierr = VecSetSizes (MyPEnKF->Mean,PETSC_DECIDE,FullVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetFromOptions (MyPEnKF->Mean);CHKERRQ(ierr);CHKMEMQ;
  /* Setup the observations mean vector */
  ierr = VecSetSizes (MyPEnKF->ObsMean,PETSC_DECIDE,ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetFromOptions (MyPEnKF->ObsMean);CHKERRQ(ierr);CHKMEMQ;

  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);CHKMEMQ;
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &ObsRowIndices);CHKERRQ(ierr);CHKMEMQ;
  for (i = 0; i < ObsVecLocalSize; i++)
    ObsRowIndices[i] = ObsVecLower + i;

  /* Compute the mean vectors, the Le, HLe and EMatrix matrices */
  for( i =0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* Get values and put them in the full ensemble mean vector */
    ierr = VecGetArray (MyPEnKF->EnsembleColumn[i], &avec);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecSetValues (MyPEnKF->Mean, FullVecLocalSize, FullRowIndices, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;

    /* Get values and put them in the observation mean vector */
    ierr = VecGetArray (MyPEnKF->EnsembleObservations[i], &avec);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecSetValues (MyPEnKF->ObsMean, ObsVecLocalSize, ObsRowIndices, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;

    /* now copy from the ensemble vectors into the Le matrix */
    ierr = VecGetArray (MyPEnKF->EnsembleColumn[i], &avec);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatSetValuesBlocked (MyPEnKF->Le, FullVecLocalSize, FullRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;

    /* now copy from the observation vectors into the HLe matrix */
    ierr = VecGetArray (MyPEnKF->EnsembleObservations[i], &avec);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatSetValuesBlocked (MyPEnKF->HLe, ObsVecLocalSize, ObsRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;
  }

  /* Perform assembly for mean vectors */
  ierr = VecAssemblyBegin(MyPEnKF->Mean);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MyPEnKF->Mean);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyBegin(MyPEnKF->ObsMean);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecAssemblyEnd(MyPEnKF->ObsMean);CHKERRQ(ierr);CHKMEMQ;

  /* Now we need to scale by -1.0/Ne  to subtract the actual mean */
  TempScalar = -1.0/MyPEnKF->NumberOfEnsembles;
  ierr = VecScale(MyPEnKF->Mean, TempScalar); CHKERRQ(ierr);CHKMEMQ;
  ierr = VecScale(MyPEnKF->ObsMean, TempScalar); CHKERRQ(ierr);CHKMEMQ;

  /* Now finalize the Le and HLe matrices */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* we are going to work with column i */
    ierr = VecGetArray (MyPEnKF->Mean, &avec); CHKERRQ(ierr);CHKMEMQ;
    ierr = MatSetValues (MyPEnKF->Le, FullVecLocalSize, FullRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecGetArray (MyPEnKF->ObsMean, &avec); CHKERRQ(ierr);CHKMEMQ;
    ierr = MatSetValues (MyPEnKF->HLe, ObsVecLocalSize, ObsRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);CHKMEMQ;
  }

  /* set everything back for the mean and observation mean */
  TempScalar = -1.0;
  ierr = VecScale(MyPEnKF->Mean, TempScalar); CHKERRQ(ierr);CHKMEMQ;
  ierr = VecScale(MyPEnKF->ObsMean, TempScalar); CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The mean vector is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = VecView (MyPEnKF->Mean,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The observation mean vector is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = VecView (MyPEnKF->ObsMean,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyEnd (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;

  /* Now Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyEnd (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;

  /* Scale Le and HLe by 1/sqrt(Ne-1) */
  TempScalar = 1.0/sqrt(MyPEnKF->NumberOfEnsembles-1);
  ierr = MatScale(MyPEnKF->Le, TempScalar);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatScale(MyPEnKF->HLe,TempScalar);CHKERRQ(ierr);CHKMEMQ;

  /* Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyEnd (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;

  /* Now Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyEnd (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;

  /* Before we transpose we make everything solvable */
#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLe matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLe,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Now we get the transpose matrix HLeT */
  ierr = MatTranspose(MyPEnKF->HLe, MAT_INITIAL_MATRIX,&MyPEnKF->HLeT);
  /* Now we have everything we need, let's find HLe*HLeT */
  ierr = MatMatMult(MyPEnKF->HLe,MyPEnKF->HLeT,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->HLeHLeT);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLEHLET matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLeHLeT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Now we compute the CD matrix */
  ierr = ValiantPEnKFComputeCD(MyPEnKF);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The CD matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->CD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  /* Sum HLeHLeT and CD, store in HLeHLeTCD */
  TempScalar = 1.0;
  ierr = MatDuplicate(MyPEnKF->HLeHLeT, MAT_COPY_VALUES, &MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAXPY(MyPEnKF->HLeHLeTCD, TempScalar, MyPEnKF->CD, DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "Le matrix looks like this:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->Le, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "HLe matrix looks like this:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLe, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "HLeT matrix looks like this:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLeT, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLeHLeT matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLeHLeT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The CD matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->CD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLeHLeTCD matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLeHLeTCD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif


  /* Now we invert HLeHLeTCD */
  /* It is very expensive to compute the inverse of a matrix and very
   * rarely needed in practice. We highly recommend avoiding  algorithms
   * that need it. The inverse of a matrix (dense or sparse) is essentially
   * always dense, so begin by creating a dense matrix B and fill it with
   * the identity matrix (ones along the diagonal), also create a dense
   * matrix X of the same size that will hold the solution. Then factor
   * the matrix you wish to invert with MatLUFactor() or MatCholeskyFactor(),
   * call the result A. Then call MatMatSolve(A,B,X) to compute the inverse into X */

  /* Get the factorization for the matrix */
  ierr = MatConvert(MyPEnKF->HLeHLeTCD,MATDENSE,MAT_REUSE_MATRIX,&MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetOption(MyPEnKF->HLeHLeTCD,MAT_ROW_ORIENTED,PETSC_FALSE);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetFromOptions(MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);CHKMEMQ;

  ierr = MatAssemblyBegin(MyPEnKF->HLeHLeTCD,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(MyPEnKF->HLeHLeTCD,MAT_FINAL_ASSEMBLY);

  /* copy the non-zero structure for the identity matrix */
  ierr = MatDuplicate(MyPEnKF->HLeHLeTCD, MAT_DO_NOT_COPY_VALUES,&TempI);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecCreate (PETSC_COMM_WORLD, &VecOnes);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetSizes (VecOnes, PETSC_DECIDE,ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecSetFromOptions (VecOnes);CHKERRQ(ierr);CHKMEMQ;
  TempScalar = 1.0;
  ierr = VecSet(VecOnes, TempScalar);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatZeroEntries(TempI);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatDiagonalSet(TempI, VecOnes,INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The TempI matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (TempI,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  ierr = MatGetOrdering(MyPEnKF->HLeHLeTCD,MATORDERING_NATURAL,&perm,&iperm);CHKERRQ(ierr);CHKMEMQ;
  /* Now for the factorization */
  if (size == 1) {
    ierr = MatGetFactor(MyPEnKF->HLeHLeTCD,MAT_SOLVER_PETSC,MAT_FACTOR_LU,&TempMat);CHKERRQ(ierr);CHKMEMQ;
  } else {
    ierr = MatGetFactor(MyPEnKF->HLeHLeTCD,MAT_SOLVER_PLAPACK,MAT_FACTOR_LU,&TempMat);CHKERRQ(ierr);CHKMEMQ;
  }
  ierr = MatLUFactorSymbolic(TempMat,MyPEnKF->HLeHLeTCD,perm,iperm,&info);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatLUFactorNumeric(TempMat,MyPEnKF->HLeHLeTCD,&info);CHKERRQ(ierr);CHKMEMQ;

  /* free the perm and iperm */
  ierr = ISDestroy (perm);CHKERRQ(ierr);CHKMEMQ;
  ierr = ISDestroy (iperm);CHKERRQ(ierr);CHKMEMQ;

  ierr = MatSetSizes(MyPEnKF->HLeHLeTCDInv, PETSC_DECIDE, PETSC_DECIDE, ObsVecGlobalSize, ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetType(MyPEnKF->HLeHLeTCDInv, MATDENSE);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetFromOptions(MyPEnKF->HLeHLeTCDInv);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyBegin(MyPEnKF->HLeHLeTCDInv, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatAssemblyEnd(MyPEnKF->HLeHLeTCDInv, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);CHKMEMQ;

  ierr = MatMatSolve(TempMat, TempI ,MyPEnKF->HLeHLeTCDInv);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatDestroy(TempMat);CHKERRQ(ierr);CHKMEMQ;

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLeHLeTCDInv matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView (MyPEnKF->HLeHLeTCDInv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  if (size == 1) {
    ierr = MatConvert(MyPEnKF->HLeT,MATSEQDENSE,MAT_REUSE_MATRIX,&MyPEnKF->HLeT);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatConvert(MyPEnKF->Le,MATSEQDENSE,MAT_REUSE_MATRIX,&MyPEnKF->Le);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatMatMult(MyPEnKF->HLeT,MyPEnKF->HLeHLeTCDInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &TempMat);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatMatMult(MyPEnKF->Le,TempMat,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->Ke);CHKERRQ(ierr);CHKMEMQ;
  } else {
    ierr = MatConvert(MyPEnKF->HLeT,MATMPIAIJ,MAT_REUSE_MATRIX,&MyPEnKF->HLeT);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatConvert(MyPEnKF->Le,MATMPIAIJ,MAT_REUSE_MATRIX,&MyPEnKF->Le);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatMatMult(MyPEnKF->HLeT,MyPEnKF->HLeHLeTCDInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &TempMat);CHKERRQ(ierr);CHKMEMQ;
    ierr = MatMatMult(MyPEnKF->Le,TempMat,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->Ke);CHKERRQ(ierr);CHKMEMQ;
  }

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The Inverse matrix is:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView(MyPEnKF->HLeHLeTCDInv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscPrintf (PETSC_COMM_WORLD, "Ke is now:\n");CHKERRQ(ierr);CHKMEMQ;
  ierr = MatView(MyPEnKF->Ke,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif

  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* find dj-Hyj */
    ierr = VecDuplicate(MyPEnKF->Observations, &TempVector);CHKERRQ(ierr);CHKMEMQ;
    ierr = VecCopy(MyPEnKF->Observations,TempVector);CHKERRQ(ierr);CHKMEMQ;

    /* Call the observational data perturbation function */
    /* This function should be adding random errors from */
    /* The same distributions as the measurement errors  */
    ierr = ValiantPEnKFPerturb(MyPEnKF);CHKERRQ(ierr);CHKMEMQ;

    /* Add the observation error: epsilon j */
    TempScalar = 1.0;
    //ierr = VecAXPY(TempVector,TempScalar,MyPEnKF->EnsembleObsError);CHKERRQ(ierr);CHKMEMQ;

    /* Subtract the forecast observations */
    TempScalar = -1.0;
    //ierr = VecAXPY(TempVector,TempScalar,MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);CHKMEMQ;

    ierr = MatMultAdd(MyPEnKF->Ke,TempVector,MyPEnKF->EnsembleColumn[i],MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);CHKMEMQ;
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "The updated ensemble vector is:\n");CHKERRQ(ierr);CHKMEMQ;
    ierr = VecView(MyPEnKF->EnsembleColumn[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);CHKMEMQ;
#endif
  }


  /* do some cleanup */
  ierr = PetscFree(FullRowIndices); CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscFree(ObsRowIndices); CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFComputeCD"
extern PetscErrorCode ValiantPEnKFComputeCD(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = ValiantPEnKFComputeCDRandomPercentage(MyPEnKF);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFPerturb"
extern PetscErrorCode ValiantPEnKFPerturb(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = ValiantPEnKFPerturbNIG(MyPEnKF);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFComputeCDRandomPercentage"
extern PetscErrorCode ValiantPEnKFComputeCDRandomPercentage(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  Vec TempVec, RandomVec;
  PetscScalar Percentage;
  PetscRandom    rctx;
  PetscInt ObsVecLocalSize, ObsVecGlobalSize;

  PetscFunctionBegin;
  /* create the temporary and random vector */
  ierr = VecDuplicate(MyPEnKF->Observations, &TempVec);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecDuplicate(MyPEnKF->Observations, &RandomVec);CHKERRQ(ierr);CHKMEMQ;

  /* copy the Observations to the temporary vector */
  ierr = VecCopy(MyPEnKF->Observations, TempVec);CHKERRQ(ierr);CHKMEMQ;
  Percentage = 0.05;
  /* scale it to a certain percentage */
  ierr = VecScale(TempVec, Percentage);CHKERRQ(ierr);CHKMEMQ;
  /* Create the random context */
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscRandomSetFromOptions(rctx);CHKERRQ(ierr);CHKMEMQ;
  /* Set the random vector */
  ierr = VecSetRandom(RandomVec, rctx);CHKERRQ(ierr);CHKMEMQ;
  ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);CHKMEMQ;
  /* multiply the %-tage with a random number */
  ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);CHKMEMQ;
  /* Compute the diagonal for CD */
  ierr = VecPointwiseMult(TempVec,TempVec,TempVec);CHKERRQ(ierr);CHKMEMQ;

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);CHKMEMQ;

  /* Setup CD, with sizes and everything */
  ierr = MatSetSizes (MyPEnKF->CD, ObsVecLocalSize, PETSC_DECIDE, ObsVecGlobalSize, ObsVecGlobalSize);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatSetType (MyPEnKF->CD, MATAIJ);CHKERRQ(ierr);CHKMEMQ;
  ierr = MatZeroEntries(MyPEnKF->CD);CHKERRQ(ierr);CHKMEMQ;

  /* zero the entries */
  ierr = MatZeroEntries(MyPEnKF->CD);CHKERRQ(ierr);CHKMEMQ;
  /* fill up with the diagonal */
  ierr = MatDiagonalSet(MyPEnKF->CD, TempVec,INSERT_VALUES);CHKERRQ(ierr);CHKMEMQ;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFPerturbNIG"
extern PetscErrorCode ValiantPEnKFPerturbNIG(PEnKFRun *MyPEnKF)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = 0;
  PetscFunctionReturn(0);
}
