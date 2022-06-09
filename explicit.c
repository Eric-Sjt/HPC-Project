static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <math.h>
int main(int argc, char **args)
{
    Vec uold,unew;
    PetscMPIInt size;
    PetscInt n=100,i,ix[n+2];
    PetscScalar dx=1.0/n,dt=0.01,y[100],x,ex;

    for(i=0;i<n+2;i++)ix[i]=i;

    /*初始化uold向量并全部赋值为0*/
    PetscCall(PetscInitialize(&argc,&args,(char*)0,help));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCall(VecCreate(PETSC_COMM_WORLD,&uold));
    PetscCall(VecSetSizes(uold,PETSC_DECIDE,n+2));
    PetscCall(VecSetFromOptions(uold));
    PetscCall(VecSet(uold,0));
    PetscCall(VecAssemblyBegin(uold));
    PetscCall(VecAssemblyEnd(uold));

    /*初始化unew向量并全部赋值为0*/
    PetscCall(VecDuplicate(uold,&unew));
    PetscCall(VecAssemblyBegin(unew));
    PetscCall(VecAssemblyEnd(unew));

    x=1.0/(2.0*n);
    for(i=1;i<n+1;i++)
    {
        ex=exp(x);
        PetscCall(VecSetValue(uold,i,ex,INSERT_VALUES));
	x=x+dx;
    }

    PetscCall(VecGetValues(uold,n,ix,y));
    for(i=0;i<n;i++)PetscPrintf(PETSC_COMM_WORLD,"%lf\n",y[i]);
    PetscCall(VecDestroy(&uold));
    PetscCall(VecDestroy(&unew));
}
