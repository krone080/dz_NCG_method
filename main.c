#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void usr_input(double ***A, unsigned* psize, double** b);

void res_output(double** A,unsigned size,double* b,double* x);

void NCG_method(double** A, unsigned size, double* b, double* x1, unsigned n);

double scal_prod(const double* a,const double* b,unsigned size);

int main (int argc,char* argv[])
 {
 double **A,*b;
 unsigned size;

 usr_input(&A,&size,&b);

 double* x=(double*) malloc(sizeof(double)*size);

 NCG_method(A,size,b,x,size);

 res_output(A,size,b,x);
 for(unsigned i=0;i<size;++i)
  free(A[i]);
 free(A);
 free(b);
 free(x);
 }

void usr_input(double*** A, unsigned *psize, double** b)
 {
 printf("Input size of matrix A: ");
 scanf("%i",psize);

 *A=(double**) malloc(sizeof(double*)*(*psize));
 for(unsigned i=0;i<*psize;++i)
  (*A)[i]=(double*) malloc(sizeof(double)*(*psize));
 *b=(double*) malloc(sizeof(double)*(*psize));

 for(unsigned i=0;i<*psize;++i)
  {
  printf("Fill matrix string A[%i]: ",i);
  fflush(stdout);
  for(unsigned j=0;j<*psize;++j)
   scanf("%lf",(*A)[i]+j);
  }

// printf("*(A[i]+j)=%0.3lf\n",*((*A)[1]+2));
// printf("A[1][2]=%0.3lf\n",(*A)[1][2]);

 printf("Fill vector b: ");
  for(unsigned j=0;j<*psize;++j)
   scanf("%lf",*b+j);
 }

void res_output(double** A,unsigned size,double* b,double* x)
 {
 printf("\nThe solution of linear system:\n");
 for(unsigned i=0;i<size;++i)
  {
  for(unsigned j=0;j<size;++j)
   printf("%.3lf\t",A[i][j]);
  printf("| %.3lf\n",b[i]);
  }

 printf("is:\n");
 for(unsigned i=0;i<size;++i)
  printf("%.3lf\n",x[i]);
 }

void NCG_method(double** A, unsigned size, double* b, double* x1, unsigned n)
 {
 double s;
 double *x0,*d0,*d1,*g0,*g1;
 x0=(double*) malloc(sizeof(double)*size);
 d0=(double*) malloc(sizeof(double)*size);
 d1=(double*) malloc(sizeof(double)*size);
 g0=(double*) malloc(sizeof(double)*size);
 g1=(double*) malloc(sizeof(double)*size);

 memset(x0,0,size*sizeof(double));
 memset(x1,0,size*sizeof(double));
 memset(d0,0,size*sizeof(double));
 memset(d1,0,size*sizeof(double));
 memcpy(g0,b,size*sizeof(double));
// printf("g0[0]=%.3lf; g0[1]=%.3lf;\n",g0[0],g0[1]);

 for(unsigned k=0;k<n;++k)
  {
  //Шаг 1
  double sum;
  for(unsigned i=0;i<size;++i)
   {
   sum=0;
   for(unsigned j=0;j<size;++j)
    sum+=A[i][j]*x0[j];
   g1[i]=sum-b[i];
   }

  //Шаг 2
  double coef=(scal_prod(g1,g1,size)/scal_prod(g0,g0,size));
  for(unsigned i=0;i<size;++i)
   d1[i]=-g1[i]+coef*d0[i];

  //Шаг 3. Массив d0 используем как временный вектор,
  //всё равно он нам больше не потребуется
  double divider=0.;
  memset(d0,0,size*sizeof(double));
  for(unsigned i=0;i<size;++i)
   {
   for(unsigned j=0;j<size;++j)
    d0[i]+=d1[j]*A[j][i];
   divider+=d0[i]*d1[i];
   }
  //На интуите ошибка в примере на этом месте: потерян минус
  s=-scal_prod(d1,g1,size)/divider;

  //Шаг 4
  for(unsigned i=0;i<size;++i)
   x1[i]=x0[i]+s*d1[i];

  //Сдвиги
  memcpy(x0,x1,size*sizeof(double));
  memcpy(d0,d1,size*sizeof(double));
  memcpy(g0,g1,size*sizeof(double));
  }
 //Освобождение памяти
 free(x0);
 free(d0);
 free(d1);
 free(g0);
 free(g1);
 }

double scal_prod(const double *a, const double *b, unsigned size)
 {
 double ret=0;
 for(unsigned i=0;i<size;++i)
  ret+=a[i]*b[i];
 return ret;
 }
