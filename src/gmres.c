/**************************************************************
 $Id: gmres.c,v 1.1 1998/06/20 06:01:59 ichiki Exp $
 GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 vol7(1986) pp.856-869
 implemented by Kengo Ichiki @ Caltech

 $Log: gmres.c,v $
 Revision 1.1  1998/06/20 06:01:59  ichiki
 Initial revision


**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrutil.h"

void mygmres(unsigned long n, double a[], double f[], double x[],
	     int m,
	     int itol, double tol,
	     int itmax, int *iter, double *err);
void back_sub(unsigned long m, double r[], double g[], double y[]);
double norm(unsigned long n, double x[]);
double inner(unsigned long n, double x[], double y[]);
void atimes(unsigned long n, double a[], double x[], double r[], int itrnsp);
void asolve(unsigned long n, double a[], double b[], double x[], int itrnsp);

/* q[] は不要、g[] をその場で計算すればよい */
/* r[] は対角部分とその下１列が H[] (h_i,j) と違うだけ */
/* しかも下１列は回転で 0 になるので、対角部分だけ保存しておけばよい */
void mygmres(unsigned long n, double a[], double f[], double x[],
	     int m,
	     int itol, double tol,
	     int itmax, int *iter, double *err)
{
  /* solve linear system A.x = f */
  /* n: dimension of this sysmtem */
  /* m: # of iteration at once */
  double norm(unsigned long n, double x[]);
  double inner(unsigned long n, double x[], double y[]);

  unsigned long
    i,j,k;
  double
    beta,
    hv,
    hjoj,
    rr,hh,
    cj,sj,
    g0,
    *tmp,*v,*h,*r,*g;

  tmp = dvector(0,n-1);
  v   = dvector(0,n*(m+1)-1);
  h   = dvector(0,m*m-1);
  r   = dvector(0,m-1);
  g   = dvector(0,m);

  /* 1. start: */
  /* compute r0 */
  for(i=0;i<n;i++){
    /* tmp = A.x0 */
    tmp[i]= 0.0;
    for(k=0;k<n;k++){
      tmp[i]+= a[i*n+k]*x[k];
    }
    tmp[i] = f[i]-tmp[i]; /* r0 */
  }
  /* compute v1 */
  beta = norm(n,tmp);
  for(i=0;i<n;i++){
    v[0*n+i] = tmp[i]/beta;
  }
  /* main loop */
  while( *iter <= itmax ){
    ++(*iter);
    /* 2. iterate: */
    for(j=0;j<m;j++){
      /* tmp = A.vj */
      for(i=0;i<n;i++){
	tmp[i] = 0.0;
	for(k=0;k<n;k++){
	  tmp[i]+= a[i*n+k]*v[j*n+k];
	}
      }
      /* h_i,j (i=1,...,j) */
      for(i=0;i<=j;i++){
	h[i*m+j] = inner(n,tmp,&v[i*n]);
      }
      /* vv_j+1 */
      for(k=0;k<n;k++){
	hv = 0.0;
	for(i=0;i<=j;i++){
	  hv += h[i*m+j]*v[i*n+k];
	}
	tmp[k] = tmp[k]-hv; /* vv_j+1 */
      }
      /* h_j+1,j */
      hjoj = norm(n,tmp); /* h_j+1,j */
      /* v_j+1 */
      for(k=0;k<n;k++){
	v[(j+1)*n+k] = tmp[k]/hjoj;
      }
      /* rotate */
      rr = h[j*m+j];
      hh = hjoj;
      hv = sqrt(rr*rr+hh*hh); /* temporary variable */
      cj =  rr/hv;
      sj = -hh/hv;
      r[j] = hv; /* diagonal element of R */
      if(j == 0){
	g[0] = cj*beta;
	g[1] = sj*beta;
      }else{
	g0 = g[j];
	g[j  ] = cj*g0;
	g[j+1] = sj*g0;
      }
    }
    /* 3. form the approximate solution */
    /* solve y_k */
    /* restore R matrix */
    for(i=0;i<m;i++){
      h[i*m+i] = r[i];
    }
    back_sub(m,h,g,r); /* use r[] as y_k */
    /* x_m */
    for(i=0;i<n;i++){
      for(k=0;k<m;k++){
	x[i] += v[k*n+i]*r[k];
      }
    }
    /* 4. restart */
    *err = g[m]; /* residual */
    /* if satisfied, */
    if( *err <= tol ) break;
    /* else */
    /* compute r_m */
    for(i=0;i<n;i++){
      /* tmp = A.x_m */
      tmp[i]= 0.0;
      for(k=0;k<n;k++){
	tmp[i]+= a[i*n+k]*x[k];
      }
      tmp[i] = f[i]-tmp[i]; /* r0 */
    }
    /* compute v1 */
    beta = norm(n,tmp);
    for(i=0;i<n;i++){
      v[0*n+i] = tmp[i]/beta;
    }
  }

  free_dvector(tmp,0,n-1);
  free_dvector(v  ,0,n*m-1);
  free_dvector(h  ,0,m*m-1);
  free_dvector(r  ,0,m*m-1);
  free_dvector(g  ,0,m);
}


void back_sub(unsigned long m, double r[], double g[], double y[])
{
  int i,j;

  for(j=m-1;j>=0;j--){
    y[j] = 0.0;
    for(i=j+1;i<m;i++){
      y[j]-=r[i*m+j]*y[i];
    }
    y[j]+=g[j];
    y[j]=y[j]/r[j*m+j];
  }
}

double norm(unsigned long n, double x[])
{
  unsigned long i;
  double ans;

  ans = 0.0;
  for (i=0;i<n;i++) ans += x[i]*x[i];
  return sqrt(ans);
}

double inner(unsigned long n, double x[], double y[])
{
  unsigned long i;
  double ans;

  ans = 0.0;
  for (i=0;i<n;i++) ans += x[i]*y[i];
  return ans;
}

void atimes(unsigned long n, double a[], double x[], double r[], int itrnsp)
{
  unsigned long i,j;

  if(itrnsp) {
    for (i=0;i<n;i++){
      r[i] = 0.0;
      for (j=0;j<n;j++){
	r[i] += a[i*n+j]*x[j];
      }
    }
  }else{
    for (i=0;i<n;i++){
      r[i] = 0.0;
      for (j=0;j<n;j++){
	r[i] += a[j*n+i]*x[j];
      }
    }
  }
}

void asolve(unsigned long n, double a[], double b[], double x[], int itrnsp)
{
  unsigned long i;

  for (i=0;i<n;i++){
    x[i] = a[i*n+i]*b[i];
  }
}
