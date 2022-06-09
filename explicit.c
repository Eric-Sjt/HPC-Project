static char help[] = "Solves a tridiagonal linear system with KSP.\n\n";
#include <petscksp.h>
#include <stdio.h>
#include <math.h>

double fun(double x)
{
    double pi=acos(-1);
    double l=1.0;
    double result;
    result=sin(l*pi*x);
    return result;
}

double fun2(double x)
{
    double pi=acos(-1);
    double result;
    result=sin(pi*x)/pi/pi;
    return result;
}

int main(int argc, char **args)
{
    Vec                 uold,unew,uexact,f;
    Mat                 A;
    PetscMPIInt         size;
    PetscInt            n=100,i,ix[n];
    PetscScalar         dx=1.0/n,dt=0.00001,y[100],x,ex,t_start=0.0,t_end=2.0;
    PetscScalar         CFL=dt/(dx*dx);

	printf("%lf",CFL);

    for(i=0;i<n+2;i++)ix[i]=i;

    PetscCall(PetscInitialize(&argc,&args,(char*)0,help));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"dt = %f\n",dt));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,"CFL = %f\n",CFL));

    /*初始化uold向量并全部赋值为0，uold表示当前时刻状态*/
    PetscCall(VecCreate(PETSC_COMM_WORLD,&uold));
    PetscCall(VecSetSizes(uold,PETSC_DECIDE,n));
    PetscCall(VecSetFromOptions(uold));
    PetscCall(VecSet(uold,0));

    /*初始化unew向量并全部赋值为0，unew表示下一时刻状态*/
    PetscCall(VecDuplicate(uold,&unew));

    /*初始化uexact向量并全部赋值为0，exact表示解析解*/
    PetscCall(VecDuplicate(uold,&uexact));
    PetscCall(VecAssemblyBegin(uexact));
    PetscCall(VecAssemblyEnd(uexact));

    /*初始化f向量并全部赋值为0*/
    PetscCall(VecDuplicate(uold,&f));

    /*初始化三对角矩阵A*/
    PetscCall(MatCreate(PETSC_COMM_WORLD,&A));
    PetscCall(MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n));
    PetscCall(MatSetFromOptions(A));
    PetscCall(MatSetUp(A));

    /*对向量赋值*/
    x=1.0/(2.0*n);
    for(i=0;i<n;i++)
    {
        PetscCall(VecSetValue(uold,i,exp(x),INSERT_VALUES));
        PetscCall(VecSetValue(f,i,fun(x)*dt,INSERT_VALUES));
        PetscCall(VecSetValue(uexact,i,fun2(x),INSERT_VALUES));
        x+=dx;
    }
    PetscCall(VecAssemblyBegin(uold));
    PetscCall(VecAssemblyEnd(uold));

    PetscCall(VecAssemblyBegin(f));
    PetscCall(VecAssemblyEnd(f));

    PetscCall(VecAssemblyBegin(uexact));
    PetscCall(VecAssemblyEnd(uexact));

    /*对三对角矩阵A赋值*/
    PetscCall(MatSetValue(A,0,0,1-2*CFL,INSERT_VALUES));
    PetscCall(MatSetValue(A,0,1,CFL,INSERT_VALUES));
    PetscCall(MatSetValue(A,n-1,n-2,CFL,INSERT_VALUES));
    PetscCall(MatSetValue(A,n-1,n-1,1-2*CFL,INSERT_VALUES));
    for(i=1;i<n-1;i++)
    {
        PetscCall(MatSetValue(A,i,i-1,CFL,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i,1-2*CFL,INSERT_VALUES));
        PetscCall(MatSetValue(A,i,i+1,CFL,INSERT_VALUES));
    }

    PetscCall(MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY));


    PetscCall(VecAssemblyBegin(unew));
    PetscCall(VecAssemblyEnd(unew));

    for(i=0;i<t_end/dt;i++)
    {
        MatMultAdd(A,uold,f,unew);
        VecCopy(unew,uold);
    }

    VecView(unew,PETSC_VIEWER_STDOUT_WORLD);
    VecView(uexact,PETSC_VIEWER_STDOUT_WORLD);

    //PetscCall(VecGetValues(uold,n,ix,y));
    //PetscPrintf(PETSC_COMM_WORLD,"%d\n",i);
    //for(i=1;i<n+1;i++)PetscPrintf(PETSC_COMM_WORLD,"%lf\n",y[i]);
    PetscCall(VecDestroy(&uold));
    PetscCall(VecDestroy(&unew));
    PetscCall(VecDestroy(&uexact));
    PetscCall(VecDestroy(&f));
    PetscCall(MatDestroy(&A));
}
