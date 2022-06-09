static char help[] = "Solves a tridiagonal linear system.\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>
#include <math.h>

#define pi acos(-1)    /*定义pi（3.1415926...）的值*/

int main(int argc,char **args)
{
  Vec            x, z, b, tem;    /*设置所需向量*/
  Mat            A;    /*设置所需矩阵*/
  KSP            ksp;    /*设置求解方法*/
  PC             pc;    /*设置求解参数*/
  PetscErrorCode ierr;    /*检查错误信息*/
  PetscInt       i, ii, col[3], rstart, rend, nlocal, rank, iter = 0;
    /*其中i,ii是矩阵和向量的角标，col是三对角矩阵参数的位置，rstart和rend均为设置矩阵时需要的参数，
    nlocal和rank为程序并行化所需参数,iter是迭代次数*/
  PetscBool      restart = PETSC_FALSE;    /*增加重读标志，默认为False*/
  PetscInt       n = 128, start = 0, end, index;    /*这是将区域分成n块，start是起始边界，end是终止边界,index仅在读取存储基础数据时使用*/
  PetscReal      dx, dt = 0.00003, t = 0.0;    /*dx是空间步长，dt是时间步长，t是已经走过的时间*/
  PetscReal      p = 1.0, c = 1.0, k = 1.0;    /*设置初始的条件参数*/
  PetscReal      te = k/p/c, alpha, u0 = 0.0;      /*通过dt和dx求解alpha，方便后续计算，u0是初始条件*/
  PetscScalar    zero = 0.0, value[3], data[3];    /*zero是所有对象的默认值，value是设置三对角矩阵的参数，data存储基础数据：dx,dt,t*/
  PetscViewer    h5;    /*创建输出*/

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;    /*初始化Petsc*/
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,NULL,NULL);CHKERRQ(ierr);    /*开始读取选项参数*/
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);    /*从命令行读取n的值（若有）*/
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);    /*从命令行读取dt的值（若有）*/
  ierr = PetscOptionsGetBool(NULL,NULL,"-restart",&restart,NULL);CHKERRQ(ierr);    /*从命令行读取是否重启（若有）*/
  ierr = PetscOptionsEnd();CHKERRQ(ierr);    /*读取选项参数结束*/
  
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);    /*设置并行MPI参数*/
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);    /*将n的值打印出来，方便阅读输出文件时参考*/

  dx = 1.0/n;    /*计算出每一小格的长度*/
  alpha = te*dt*n*n;    /*计算出CFL的值（隐式格式中，CFL的值无限制）*/
  end = n;    /*更新end的值*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr);    /*将dx的值打印出来，方便阅读输出文件时参考*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dt = %f\n",dt);CHKERRQ(ierr);    /*将dt的值打印出来，方便阅读输出文件时参考*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"alpha = %f\n",alpha);CHKERRQ(ierr);    /*将alpha的值打印出来，方便阅读输出文件时参考*/
  ierr = PetscPrintf(PETSC_COMM_WORLD,"restart = %d\n",restart);CHKERRQ(ierr);    /*将restart的值打印出来，方便阅读输出文件时参考*/

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);    /*创建一个并行空间*/
  ierr = VecCreate(PETSC_COMM_WORLD,&tem);CHKERRQ(ierr);    /*创建临时向量*/
  ierr = VecSetSizes(x,PETSC_DECIDE,n+1);CHKERRQ(ierr);    /*创建一个长度n+1的矩阵*/
  ierr = VecSetSizes(tem, 3, PETSC_DECIDE);CHKERRQ(ierr);    /*创建长度为3的临时向量*/
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);    /*从选项数据库中配置向量*/
  ierr = VecSetFromOptions(tem);CHKERRQ(ierr);    /*获得参数*/
  ierr = VecDuplicate(x,&z);CHKERRQ(ierr);    /*将x的格式赋给z*/
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);    /*将x的格式赋给b*/

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);    /*设置并行x的起始终止点*/
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);    /*设置并行x的位置点*/

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);    /*在并行空间创建一个矩阵*/
  ierr = MatSetSizes(A,nlocal,nlocal,n+1,n+1);CHKERRQ(ierr);    /*设置矩阵的行数和列数*/
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);    /*从选项数据库中配置矩阵*/
  ierr = MatSetUp(A);CHKERRQ(ierr);    /*开始建立矩阵*/

  if (!rstart){    /*若rstart为0时，即为首行*/
    rstart = 1;    /*将rstart设为1*/
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1+2.0*alpha; value[1] = -alpha;    /*设置要用到的参数*/
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的第一行*/
  }
  
  if (rend == n+1){    /*最后一行*/
    rend = n;    /*将rend设为n*/
    i    = n; col[0] = n-1; col[1] = n; value[0] = -alpha; value[1] = 1+2.0*alpha;    /*设置要用到的参数*/
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的最后一行*/
  }

  value[0] = -alpha; value[1] = 1+2.0*alpha; value[2] = -alpha;    /*设置三对角矩阵除首尾两行外的其余行的三个值*/
  for (i=rstart; i<rend; i++){    /*除首尾两行外的行*/
    col[0] = i-1; col[1] = i; col[2] = i+1;    /*设置要用到的参数*/
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);    /*设置三对角矩阵的三对角值*/
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);    /*通知其余并行块将矩阵统一*/
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);    /*结束通知*/
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    /*打印矩阵，检查是否出错*/
  
  /*判断是否需要重读还是新建矩阵*/
  if(restart){    /*如果restart为True，表示重读，则从文件开始读入*/
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"implicit.h5", FILE_MODE_READ, &h5);CHKERRQ(ierr);    /*创建输入文件*/
      ierr = PetscObjectSetName((PetscObject) z, "implicit-vector");CHKERRQ(ierr);    /*将z输入的名字命名为explicit-vector*/
      ierr = PetscObjectSetName((PetscObject) tem, "implicit-necess-data");CHKERRQ(ierr);    /*将临时向量tem输入的名字命名为explicit-necess-data*/
      ierr = VecLoad(tem, h5);CHKERRQ(ierr);    /*将读入的数据加载到向量tem中*/
      ierr = VecLoad(z, h5);CHKERRQ(ierr);    /*将读入的数据加载到向量z中*/
      ierr = PetscViewerDestroy(&h5);CHKERRQ(ierr);    /*关闭输入*/
      index=0;    /*将索引初始化*/
      ierr = VecGetValues(tem,1,&index,&dx);CHKERRQ(ierr);    /*将第一个值赋给dx*/
      index=index+1;    /*索引移向下一位*/
      ierr = VecGetValues(tem,1,&index,&dt);CHKERRQ(ierr);    /*将第二个值赋给dt*/
      index=index+1;    /*索引移向下一位*/
      ierr = VecGetValues(tem,1,&index,&t);CHKERRQ(ierr);    /*将第三个值赋给t*/
      index= 0;    /*索引复位*/
  }
  else{    /*如果restart为False，表示新的开始，则开始进行向量初始化，构建向量*/
    ierr = VecSet(z,zero);CHKERRQ(ierr);    /*设置初始向量b*/
    if(rank == 0){    /*开始设置初始条件*/
      for(ii = 1; ii < n; ii++){    /*除首尾两个点外的其余点*/
        u0 = exp(ii*dx);    /*根据当前位置来获取初始值*/
        ierr = VecSetValues(z, 1, &ii, &u0, INSERT_VALUES);CHKERRQ(ierr);    /*将向量的对应位置的值进行修改*/
      }
    }
    ierr = VecAssemblyBegin(z);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
    ierr = VecAssemblyEnd(z);CHKERRQ(ierr);    /*结束通知*/
  }

  ierr = VecSet(b,zero);CHKERRQ(ierr);    /*设置初始向量u*/
  if(rank == 0){    /*开始设置初始条件*/
    for(ii = 1; ii < n; ii++){    /*除首尾两个点外的其余点*/
      u0 = dt*sin(ii*dx*pi);    /*根据当前位置来获取传热值*/
      ierr = VecSetValues(b, 1, &ii, &u0, INSERT_VALUES);CHKERRQ(ierr);    /*将向量的对应位置的值进行修改*/
    }
  }
  
  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);    /*结束通知*/
  
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);    /*创建ksp解空间*/
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);    /*设置方程左侧的系数*/
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);    /*设置矩阵求解的相关系数*/
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);    /*设置pc的默认参数*/
  ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);    /*设置各种误差值*/
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);    /*从选项数据库中配置ksp解空间*/

  while(PetscAbsReal(t)<2.0){    /*计算0-2时间内的传播*/

    t += dt;    /*时间向前走*/

    ierr = VecAXPY(z,1.0,b);CHKERRQ(ierr);    /*设置方程右边项的值*/
    ierr = KSPSolve(ksp,z,x);CHKERRQ(ierr);    /*求解方程*/

    ierr = VecSetValues(x, 1, &start, &zero, INSERT_VALUES);CHKERRQ(ierr);    /*设置边界条件*/
    ierr = VecSetValues(x, 1, &end, &zero, INSERT_VALUES);CHKERRQ(ierr);    /*设置边界条件*/
    ierr = VecAssemblyBegin(x);CHKERRQ(ierr);    /*统一向量更新*/
    ierr = VecAssemblyEnd(x);CHKERRQ(ierr);    /*结束更新*/

    ierr = VecCopy(x,z);CHKERRQ(ierr);    /*将x的值赋给z*/

    iter += 1;    /*记录迭代次数*/
     if((iter%10)==0){    /*如果迭代次数为10的倍数，即每迭代十次*/

       data[0] = dx; data[1] = dt; data[2] = t;    /*将值赋给数组*/
       ierr = VecSet(tem,zero);CHKERRQ(ierr);    /*初始化矩阵*/
       for(index=0;index<3;index++){    /*循环遍历数组，并将值赋给向量*/
        u0 = data[index];    /*将数组的值赋给u0*/
        ierr = VecSetValues(tem,1,&index,&u0,INSERT_VALUES);CHKERRQ(ierr);    /*将矩阵赋值给向量*/
       }
       ierr = VecAssemblyBegin(tem);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
       ierr = VecAssemblyEnd(tem);CHKERRQ(ierr);    /*结束通知*/

       ierr = PetscViewerCreate(PETSC_COMM_WORLD,&h5);CHKERRQ(ierr);    /*创建输出指针*/
       ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"implicit.h5", FILE_MODE_WRITE, &h5);CHKERRQ(ierr);    /*创建输出文件*/
       ierr = PetscObjectSetName((PetscObject) z, "implicit-vector");CHKERRQ(ierr);    /*将z输出的名字命名为explicit-vector*/
       ierr = PetscObjectSetName((PetscObject) tem, "implicit-necess-data");CHKERRQ(ierr);    /*将tem输出的名字命名为implicit-necess-data*/
       ierr = VecView(tem, h5);CHKERRQ(ierr);    /*tem输出到文件*/
       ierr = VecView(z, h5);CHKERRQ(ierr);    /*z输出到文件*/
       ierr = PetscViewerDestroy(&h5);CHKERRQ(ierr);    /*关闭输出*/
     }
  }

  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);    /*打印向量，获得结束时隐式方法的值*/
 
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);    /*关闭ksp*/
  ierr = VecDestroy(&tem);CHKERRQ(ierr);    /*关闭临时向量*/
  ierr = VecDestroy(&x);CHKERRQ(ierr);    /*关闭向量x*/
  ierr = VecDestroy(&b);CHKERRQ(ierr);    /*关闭向量u*/
  ierr = VecDestroy(&z);CHKERRQ(ierr);    /*关闭向量b*/
  ierr = MatDestroy(&A);CHKERRQ(ierr);    /*关闭矩阵A*/

  ierr = PetscFinalize();    /*结束并行*/
  return ierr;    /*程序结束*/
}

// EOF
