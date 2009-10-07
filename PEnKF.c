/*
 * PEnKF.c
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
    ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (Vec),&MyPEnKF->EnsembleVecs[i]); CHKERRQ(ierr);

  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* for every ensemble iterate over all the state vectors */
    for(j = 0; j< MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleVecs[i][j]));CHKERRQ(ierr);
    }
    /* create the Observation vectors */
    ierr = VecCreate (PETSC_COMM_WORLD,&(MyPEnKF->EnsembleObservations[i]));CHKERRQ(ierr);
    /* create the Full Ensemble vectors */
    ierr = VecCreate (PETSC_COMM_WORLD, &(MyPEnKF->EnsembleColumn[i]));CHKERRQ(ierr);
  }

  /* Create the  Real world observations vectors */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Observations);CHKERRQ(ierr);
  /* Setup the Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->Mean);CHKERRQ(ierr);
  /* Setup the Observations Mean Vector */
  ierr = VecCreate (PETSC_COMM_WORLD, &MyPEnKF->ObsMean);CHKERRQ(ierr);

  /* Now for some matrices */
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->Le);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLe);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeT);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeT);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->HLeHLeTCDInv);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->CD);CHKERRQ(ierr);
  ierr = MatCreate (PETSC_COMM_WORLD, &MyPEnKF->Ke);CHKERRQ(ierr);

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
       ierr = VecDestroy (MyPEnKF->EnsembleVecs[i][j]);CHKERRQ(ierr);
     }
     ierr = PetscFree (MyPEnKF->EnsembleVecs[i]);CHKERRQ(ierr);
     ierr = VecDestroy (MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);
     ierr = VecDestroy (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);
   }

   ierr = VecDestroy(MyPEnKF->Observations);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->Mean);CHKERRQ(ierr);
   ierr = VecDestroy(MyPEnKF->ObsMean);CHKERRQ(ierr);

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
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][0]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][1]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][2]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][3]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][4]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the Saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][5]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the Observation vector */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Observations.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

    }
  }

  sprintf(FileNameBuffer, "%s%s", MyPEnKF->ObsPathPrefix, "/FieldObservations.Valiant");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->Observations);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

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
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][0]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][1]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][2]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][3]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the water pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][4]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the oil pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Po.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][5]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the gas pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pg.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][6]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the water saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][7]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the oil saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/So.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleVecs[i][8]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the Observation vector */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Observations.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
      ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

    }
  }

  sprintf(FileNameBuffer, "%s%s", MyPEnKF->ObsPathPrefix, "/FieldObservations.Valiant");
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = VecLoad(viewer,PETSC_NULL,&MyPEnKF->Observations);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

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
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][0],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][1],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][2],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][3],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][4],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the Saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][5],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
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
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][0],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the LnK 11 22 33 */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK11.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][1],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK22.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][2],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/LnK33.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][3],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the water pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][4],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the oil pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Po.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][5],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the gas pressure */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Pg.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][6],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the water saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/Sw.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][7],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

      /* Load the oil saturation */
      sprintf(FileNameBuffer, "%s%s", PathBuffer, "/So.Valiant");
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,FileNameBuffer,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
      ierr = VecView(MyPEnKF->EnsembleVecs[i][8],viewer);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);
    }
  }

  PetscFunctionReturn(0);
}


#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFGather"
extern PetscErrorCode ValiantPEnKFGather(PEnKFRun *MyPEnKF)
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
  ierr = VecGetSize (MyPEnKF->EnsembleVecs[0][0], &StateVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleVecs[0][0], &StateVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleVecs[0][0], &StateVecLower, &StateVecHigher);CHKERRQ(ierr);

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);

  /* Allocate memory for the scatters */
  ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (IS), &isVecs);CHKERRQ(ierr);

  /* For every vector, create the scatter into the ensemble vector */
  for( i =0; i < MyPEnKF->NumberOfVecsPerEnsemble; i++)
  {
    ierr = PetscMalloc (StateVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
    for (j = 0; j < StateVecLocalSize; j++)
    {
      indices[j] = i*StateVecGlobalSize + StateVecLower + j;
    }
    ierr = ISCreateGeneral (PETSC_COMM_WORLD, StateVecLocalSize, indices,&isVecs[i]);CHKERRQ(ierr);
    ierr = PetscFree (indices);CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisVecs is:\n");CHKERRQ(ierr);
    ierr = ISView(isVecs[i], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  }

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize + ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&isObs);CHKERRQ(ierr);
  ierr = PetscFree (indices);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisObs\n");CHKERRQ(ierr);
  ierr = ISView (isObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&intoObs);CHKERRQ(ierr);
  ierr = PetscFree (indices);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nintoObs\n");CHKERRQ(ierr);
  ierr = ISView (intoObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Set the size of the column vector, assuming we already loaded the state vectors and observations */
  ColumnVecGlobalSize = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize+ObsVecGlobalSize;
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    ierr = VecSetSizes (MyPEnKF->EnsembleColumn[i], PETSC_DECIDE, ColumnVecGlobalSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);
  }

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);
  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &RowIndices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    RowIndices[i] = ObsVecLower + i;

  /* Scatter in reverse mode */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    for(j = 0; j < MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecScatterCreate (MyPEnKF->EnsembleVecs[i][j], isVecs[0],MyPEnKF->EnsembleColumn[i],isVecs[j], &VecToEnsembleVec);CHKERRQ(ierr);
#ifdef PENKF_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nEnsembleVec To Vec\n");CHKERRQ(ierr);
      ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
      ierr = VecScatterBegin (VecToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleVecs[i][j],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);
      ierr = VecScatterEnd (VecToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleVecs[i][j],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);

      /* Be nice and cleanup before leaving */
      ierr = VecScatterDestroy (VecToEnsembleVec);CHKERRQ(ierr);
#ifdef PENKF_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nPost Data Assimilation Vectors\n");CHKERRQ(ierr);
      ierr = VecView (MyPEnKF->EnsembleVecs[i][j], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
    }
    ierr = VecScatterCreate (MyPEnKF->EnsembleObservations[i], intoObs,MyPEnKF->EnsembleColumn[i],isObs, &ObsToEnsembleVec);CHKERRQ(ierr);
#ifdef PENKF_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nEnsembleVec To VecObservation\n");CHKERRQ(ierr);
    ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
    ierr = VecScatterBegin (ObsToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleObservations[i],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);
    ierr = VecScatterEnd (ObsToEnsembleVec, MyPEnKF->EnsembleColumn[i], MyPEnKF->EnsembleObservations[i],INSERT_VALUES, SCATTER_REVERSE);CHKERRQ(ierr);
    /* Be nice and cleanup before leaving */
    ierr = VecScatterDestroy (ObsToEnsembleVec);CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "ValiantPEnKFScatter"
extern PetscErrorCode ValiantPEnKFScatter(PEnKFRun *MyPEnKF)
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
  ierr = VecGetSize (MyPEnKF->EnsembleVecs[0][0], &StateVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleVecs[0][0], &StateVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleVecs[0][0], &StateVecLower, &StateVecHigher);CHKERRQ(ierr);

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);

  /* Allocate memory for the scatters */
  ierr = PetscMalloc (MyPEnKF->NumberOfVecsPerEnsemble *sizeof (IS), &isVecs);CHKERRQ(ierr);

  /* For every vector, create the scatter into the ensemble vector */
  for( i =0; i < MyPEnKF->NumberOfVecsPerEnsemble; i++)
  {
    ierr = PetscMalloc (StateVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
    for (j = 0; j < StateVecLocalSize; j++)
    {
      indices[j] = i*StateVecGlobalSize + StateVecLower + j;
    }
    ierr = ISCreateGeneral (PETSC_COMM_WORLD, StateVecLocalSize, indices,&isVecs[i]);CHKERRQ(ierr);
    ierr = PetscFree (indices);CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisVecs is:\n");CHKERRQ(ierr);
    ierr = ISView(isVecs[i], PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  }

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize + ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&isObs);CHKERRQ(ierr);
  ierr = PetscFree (indices);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nisObs\n");CHKERRQ(ierr);
  ierr = ISView (isObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Set up the scatter array for the Observation vectors */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &indices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    indices[i] = ObsVecLower + i;

  ierr = ISCreateGeneral (PETSC_COMM_WORLD, ObsVecLocalSize, indices,&intoObs);CHKERRQ(ierr);
  ierr = PetscFree (indices);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "\nintoObs\n");CHKERRQ(ierr);
  ierr = ISView (intoObs, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Set the size of the column vector, assuming we already loaded the state vectors and observations */
  ColumnVecGlobalSize = MyPEnKF->NumberOfVecsPerEnsemble*StateVecGlobalSize+ObsVecGlobalSize;
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    ierr = VecSetSizes (MyPEnKF->EnsembleColumn[i], PETSC_DECIDE, ColumnVecGlobalSize);CHKERRQ(ierr);
    ierr = VecSetFromOptions (MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);
  }

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);
  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &RowIndices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    RowIndices[i] = ObsVecLower + i;

  /* For all ensembles, perform scatter */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    for(j = 0; j < MyPEnKF->NumberOfVecsPerEnsemble; j++)
    {
      ierr = VecScatterCreate (MyPEnKF->EnsembleVecs[i][j], isVecs[0],MyPEnKF->EnsembleColumn[i],isVecs[j], &VecToEnsembleVec);CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
      ierr = PetscPrintf (PETSC_COMM_WORLD, "\nVecToEnsembleVec\n");CHKERRQ(ierr);
      ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
      ierr = VecScatterBegin (VecToEnsembleVec, MyPEnKF->EnsembleVecs[i][j],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
      ierr = VecScatterEnd (VecToEnsembleVec, MyPEnKF->EnsembleVecs[i][j],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);

      /* Be nice and cleanup before leaving */
      ierr = VecScatterDestroy (VecToEnsembleVec);CHKERRQ(ierr);
    }

    ierr = VecScatterCreate (MyPEnKF->EnsembleObservations[i], intoObs,MyPEnKF->EnsembleColumn[i],isObs, &ObsToEnsembleVec);CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "\nObsToEnsembleVec\n");CHKERRQ(ierr);
    ierr = VecScatterView (VecToEnsembleVec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
    ierr = VecScatterBegin (ObsToEnsembleVec, MyPEnKF->EnsembleObservations[i],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    ierr = VecScatterEnd (ObsToEnsembleVec, MyPEnKF->EnsembleObservations[i],MyPEnKF->EnsembleColumn[i],INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
    /* Be nice and cleanup before leaving */
    ierr = VecScatterDestroy (ObsToEnsembleVec);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "Ensemble vector after scatter\n");CHKERRQ(ierr);
    ierr = VecView (MyPEnKF->EnsembleColumn[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
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
  ierr = MPI_Comm_size (PETSC_COMM_WORLD, &size);CHKERRQ(ierr);

  /* Retrive global size, local size lower and higher extent of the full vector */
  ierr = VecGetSize (MyPEnKF->EnsembleColumn[0], &FullVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleColumn[0], &FullVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleColumn[0], &FullVecLower, &FullVecHigher);CHKERRQ(ierr);

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);
  ierr = VecGetOwnershipRange (MyPEnKF->EnsembleObservations[0], &ObsVecLower, &ObsVecHigher);CHKERRQ(ierr);

  /* Create the Le matrix */
  ierr = MatSetSizes (MyPEnKF->Le, FullVecLocalSize, PETSC_DECIDE, FullVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);
  ierr = MatSetType (MyPEnKF->Le, MATDENSE);CHKERRQ(ierr);

  /* Create the HLe matrix */
  ierr = MatSetSizes (MyPEnKF->HLe, ObsVecLocalSize, PETSC_DECIDE, ObsVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);
  ierr = MatSetType (MyPEnKF->HLe, MATDENSE);CHKERRQ(ierr);

  /* Setup the mean vectors */
  ierr = VecSetSizes (MyPEnKF->Mean,PETSC_DECIDE,FullVecGlobalSize);CHKERRQ(ierr);
  ierr = VecSetFromOptions (MyPEnKF->Mean);CHKERRQ(ierr);
  /* Setup the observations mean vector */
  ierr = VecSetSizes (MyPEnKF->ObsMean,PETSC_DECIDE,ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecSetFromOptions (MyPEnKF->ObsMean);CHKERRQ(ierr);

  /* set the indices for the local vecsetvalues */
  ierr = PetscMalloc (FullVecLocalSize*sizeof (PetscInt), &FullRowIndices);CHKERRQ(ierr);
  for(i = 0; i < FullVecLocalSize; i++)
    FullRowIndices[i] = FullVecLower + i;

  /* Allocate memory for all the number of rows we will have in the HLe matrix */
  /* The number of rows is the local size of the observations vector */
  ierr = PetscMalloc (ObsVecLocalSize * sizeof (PetscInt), &ObsRowIndices);CHKERRQ(ierr);
  for (i = 0; i < ObsVecLocalSize; i++)
    ObsRowIndices[i] = ObsVecLower + i;

  /* Now we have good sizes we set sizes of matrices */
  /* Setup the Le matrix */
  ierr = MatSetSizes (MyPEnKF->Le, FullVecLocalSize, PETSC_DECIDE, FullVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);
  ierr = MatSetType (MyPEnKF->Le, MATDENSE);CHKERRQ(ierr);

  /* Setup the HLe matrix */
  ierr = MatSetSizes (MyPEnKF->HLe, ObsVecLocalSize, PETSC_DECIDE, ObsVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);
  ierr = MatSetType (MyPEnKF->HLe, MATDENSE);CHKERRQ(ierr);

  /* Compute the mean vectors, the Le, HLe and EMatrix matrices */
  for( i =0; i < MyPEnKF->NumberOfVecsPerEnsemble; i++)
  {
    /* Get values and put them in the full ensemble mean vector */
    ierr = VecGetArray (MyPEnKF->EnsembleColumn[i], &avec);CHKERRQ(ierr);
    ierr = VecSetValues (MyPEnKF->Mean, FullVecLocalSize, FullRowIndices, avec, ADD_VALUES);CHKERRQ(ierr);

    /* Get values and put them in the observation mean vector */
    ierr = VecGetArray (MyPEnKF->EnsembleObservations[i], &avec);CHKERRQ(ierr);
    ierr = VecSetValues (MyPEnKF->ObsMean, ObsVecLocalSize, ObsRowIndices, avec, ADD_VALUES);CHKERRQ(ierr);

    /* now copy from the ensemble vectors into the Le matrix */
    ierr = VecGetArray (MyPEnKF->EnsembleColumn[i], &avec);CHKERRQ(ierr);
    ierr = MatSetValuesBlocked (MyPEnKF->Le, FullVecLocalSize, FullRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);

    /* now copy from the observation vectors into the HLe matrix */
    ierr = VecGetArray (MyPEnKF->EnsembleObservations[i], &avec);CHKERRQ(ierr);
    ierr = MatSetValuesBlocked (MyPEnKF->HLe, ObsVecLocalSize, ObsRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);
  }

  /* Perform assembly for mean vectors */
  ierr = VecAssemblyBegin(MyPEnKF->Mean);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MyPEnKF->Mean);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(MyPEnKF->ObsMean);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(MyPEnKF->ObsMean);CHKERRQ(ierr);

  /* Now we need to scale by -1.0/Ne  to subtract the actual mean */
  TempScalar = -1.0/MyPEnKF->NumberOfEnsembles;
  ierr = VecScale(MyPEnKF->Mean, TempScalar); CHKERRQ(ierr);
  ierr = VecScale(MyPEnKF->ObsMean, TempScalar); CHKERRQ(ierr);

  /* Now finalize the Le and HLe matrices */
  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* we are going to work with column i */
    ierr = VecGetArray (MyPEnKF->Mean, &avec); CHKERRQ(ierr);
    ierr = MatSetValues (MyPEnKF->Le, FullVecLocalSize, FullRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);
    ierr = VecGetArray (MyPEnKF->ObsMean, &avec); CHKERRQ(ierr);
    ierr = MatSetValues (MyPEnKF->HLe, ObsVecLocalSize, ObsRowIndices,1,&i, avec, ADD_VALUES);CHKERRQ(ierr);
  }

  /* set everything back for the mean and observation mean */
  TempScalar = -1.0;
  ierr = VecScale(MyPEnKF->Mean, TempScalar); CHKERRQ(ierr);
  ierr = VecScale(MyPEnKF->ObsMean, TempScalar); CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The mean vector is:\n");CHKERRQ(ierr);
  ierr = VecView (MyPEnKF->Mean,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The observation mean vector is:\n");CHKERRQ(ierr);
  ierr = VecView (MyPEnKF->ObsMean,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Now Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Scale Le and HLe by 1/sqrt(Ne-1) */
  TempScalar = 1.0/sqrt(MyPEnKF->NumberOfEnsembles-1);
  ierr = MatScale(MyPEnKF->Le, TempScalar);CHKERRQ(ierr);
  ierr = MatScale(MyPEnKF->HLe,TempScalar);CHKERRQ(ierr);

  /* Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd (MyPEnKF->Le, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Now Assemble the Matrix, about to start our operations */
  ierr = MatAssemblyBegin (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd (MyPEnKF->HLe, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  /* Before we transpose we make everything solvable */
#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLe matrix is:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->HLe,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  /* Now we get the transpose matrix HLeT */
  ierr = MatTranspose(MyPEnKF->HLe, MAT_INITIAL_MATRIX,&MyPEnKF->HLeT);
  /* Now we have everything we need, let's find HLe*HLeT */
  ierr = MatMatMult(MyPEnKF->HLe,MyPEnKF->HLeT,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->HLeHLeT);CHKERRQ(ierr);

  /* Now we compute the CD matrix */
  ierr = ValiantPEnKFComputeCD(MyPEnKF);CHKERRQ(ierr);

  /* Sum HLeHLeT and CD, store in HLeHLeTCD */
  TempScalar = 1.0;
  ierr = MatDuplicate(MyPEnKF->HLeHLeT, MAT_COPY_VALUES, &MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);
  ierr = MatAYPX(MyPEnKF->HLeHLeTCD, TempScalar, MyPEnKF->CD, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "Le matrix looks like this:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->Le, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "HLe matrix looks like this:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->HLe, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "HLeT matrix looks like this:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->HLeT, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLeHLeT matrix is:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->HLeHLeT,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The CD matrix is:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->CD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The HLeHLeTCD matrix is:\n");CHKERRQ(ierr);
  ierr = MatView (MyPEnKF->HLeHLeTCD,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
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
  ierr = MatConvert(MyPEnKF->HLeHLeTCD,MATDENSE,MAT_REUSE_MATRIX,&MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);
  ierr = MatSetOption(MyPEnKF->HLeHLeTCD,MAT_ROW_ORIENTED,PETSC_FALSE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(MyPEnKF->HLeHLeTCD);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(MyPEnKF->HLeHLeTCD,MAT_FINAL_ASSEMBLY);
  ierr = MatAssemblyEnd(MyPEnKF->HLeHLeTCD,MAT_FINAL_ASSEMBLY);

  /* copy the non-zero structure for the identity matrix */
  ierr = MatDuplicate(MyPEnKF->HLeHLeTCD, MAT_DO_NOT_COPY_VALUES,&TempI);CHKERRQ(ierr);
  ierr = VecCreate (PETSC_COMM_WORLD, &VecOnes);CHKERRQ(ierr);
  ierr = VecSetSizes (VecOnes, PETSC_DECIDE,ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecSetFromOptions (VecOnes);CHKERRQ(ierr);
  TempScalar = 1.0;
  ierr = VecSet(VecOnes, TempScalar);CHKERRQ(ierr);
  ierr = MatZeroEntries(TempI);CHKERRQ(ierr);
  ierr = MatDiagonalSet(TempI, VecOnes,INSERT_VALUES);CHKERRQ(ierr);

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The TempI matrix is:\n");CHKERRQ(ierr);
  ierr = MatView (TempI,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  ierr = MatGetOrdering(MyPEnKF->HLeHLeTCD,MATORDERING_NATURAL,&perm,&iperm);CHKERRQ(ierr);
  /* Now for the factorization */
  if (size == 1) {
    ierr = MatGetFactor(MyPEnKF->HLeHLeTCD,MAT_SOLVER_PETSC,MAT_FACTOR_LU,&TempMat);CHKERRQ(ierr);
  } else {
    ierr = MatGetFactor(MyPEnKF->HLeHLeTCD,MAT_SOLVER_PLAPACK,MAT_FACTOR_LU,&TempMat);CHKERRQ(ierr);
  }
  ierr = MatLUFactorSymbolic(TempMat,MyPEnKF->HLeHLeTCD,perm,iperm,&info);CHKERRQ(ierr);
  ierr = MatLUFactorNumeric(TempMat,MyPEnKF->HLeHLeTCD,&info);CHKERRQ(ierr);

  /* free the perm and iperm */
  ierr = ISDestroy (perm);CHKERRQ(ierr);
  ierr = ISDestroy (iperm);CHKERRQ(ierr);

  ierr = MatMatSolve(TempMat, TempI ,MyPEnKF->HLeHLeTCDInv);CHKERRQ(ierr);
  ierr = MatDestroy(TempMat);CHKERRQ(ierr);

  if (size == 1) {
    ierr = MatConvert(MyPEnKF->HLeT,MATSEQDENSE,MAT_REUSE_MATRIX,&MyPEnKF->HLeT);CHKERRQ(ierr);
    ierr = MatConvert(MyPEnKF->Le,MATSEQDENSE,MAT_REUSE_MATRIX,&MyPEnKF->Le);CHKERRQ(ierr);
    ierr = MatMatMult(MyPEnKF->HLeT,MyPEnKF->HLeHLeTCDInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &TempMat);CHKERRQ(ierr);
    ierr = MatMatMult(MyPEnKF->Le,TempMat,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->Ke);CHKERRQ(ierr);
  } else {
    ierr = MatMatMult(MyPEnKF->HLeT,MyPEnKF->HLeHLeTCDInv,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &TempMat);CHKERRQ(ierr);
    ierr = MatMatMult(MyPEnKF->Le,TempMat,MAT_INITIAL_MATRIX,PETSC_DEFAULT, &MyPEnKF->Ke);CHKERRQ(ierr);
  }

#ifdef VALIANT_DEBUG
  ierr = PetscPrintf (PETSC_COMM_WORLD, "The Inverse matrix is:\n");CHKERRQ(ierr);
  ierr = MatView(MyPEnKF->HLeHLeTCDInv,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf (PETSC_COMM_WORLD, "Ke is now:\n");CHKERRQ(ierr);
  ierr = MatView(MyPEnKF->Ke,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif

  for (i = 0; i < MyPEnKF->NumberOfEnsembles; i++)
  {
    /* find dj-Hyj */
    ierr = VecDuplicate(MyPEnKF->Observations, &TempVector);CHKERRQ(ierr);
    ierr = VecCopy(MyPEnKF->Observations,TempVector);CHKERRQ(ierr);

    /* Call the observational data perturbation function */
    /* This function should be adding random errors from */
    /* The same distributions as the measurement errors  */
    ierr = ValiantPEnKFPerturb(MyPEnKF);CHKERRQ(ierr);

    /* Add the observation error: epsilon j */
    TempScalar = 1.0;
    ierr = VecAXPY(TempVector,TempScalar,MyPEnKF->EnsembleObsError);CHKERRQ(ierr);

    /* Subtract the forecast observations */
    TempScalar = -1.0;
    ierr = VecAXPY(TempVector,TempScalar,MyPEnKF->EnsembleObservations[i]);CHKERRQ(ierr);

    ierr = MatMultAdd(MyPEnKF->Ke,TempVector,MyPEnKF->EnsembleColumn[i],MyPEnKF->EnsembleColumn[i]);CHKERRQ(ierr);
#ifdef VALIANT_DEBUG
    ierr = PetscPrintf (PETSC_COMM_WORLD, "The updated ensemble vector is:\n");CHKERRQ(ierr);
    ierr = VecView(MyPEnKF->EnsembleColumn[i],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  }


  /* do some cleanup */
  ierr = PetscFree(FullRowIndices); CHKERRQ(ierr);
  ierr = PetscFree(ObsRowIndices); CHKERRQ(ierr);

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
  ierr = VecDuplicate(MyPEnKF->Observations, &TempVec);CHKERRQ(ierr);
  ierr = VecDuplicate(MyPEnKF->Observations, &RandomVec);CHKERRQ(ierr);

  /* copy the Observations to the temporary vector */
  ierr = VecCopy(MyPEnKF->Observations, TempVec);CHKERRQ(ierr);
  Percentage = 0.05;
  /* scale it to a certain percentage */
  ierr = VecScale(TempVec, Percentage);CHKERRQ(ierr);
  /* Create the random context */
  ierr = PetscRandomCreate(PETSC_COMM_WORLD,&rctx);
  /* Set the random vector */
  ierr = VecSetRandom(RandomVec, rctx);CHKERRQ(ierr);
  ierr = PetscRandomDestroy(rctx);CHKERRQ(ierr);
  /* multiply the %-tage with a random number */
  ierr = VecPointwiseMult(TempVec,TempVec,RandomVec);CHKERRQ(ierr);
  /* Compute the diagonal for CD */
  ierr = VecPointwiseMult(TempVec,TempVec,TempVec);CHKERRQ(ierr);

  /* Retrieve global, local sizes and extents of the observation vectors */
  ierr = VecGetSize (MyPEnKF->EnsembleObservations[0], &ObsVecGlobalSize);CHKERRQ(ierr);
  ierr = VecGetLocalSize (MyPEnKF->EnsembleObservations[0], &ObsVecLocalSize);CHKERRQ(ierr);

  /* Setup CD, with sizes and everything */
  ierr = MatSetSizes (MyPEnKF->CD, ObsVecLocalSize, PETSC_DECIDE, ObsVecGlobalSize,MyPEnKF->NumberOfEnsembles);CHKERRQ(ierr);
  ierr = MatSetType (MyPEnKF->CD, MATAIJ);CHKERRQ(ierr);
  ierr = MatZeroEntries(MyPEnKF->CD);CHKERRQ(ierr);

  /* zero the entries */
  ierr = MatZeroEntries(MyPEnKF->CD);CHKERRQ(ierr);
  /* fill up with the diagonal */
  ierr = MatDiagonalSet(MyPEnKF->CD, TempVec,INSERT_VALUES);CHKERRQ(ierr);

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
