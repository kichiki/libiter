/**************************************************************
 $Id: gmres.c,v 1.3 1998/06/20 22:31:02 ichiki Exp $
 GMRES(m) : Y.Saad & M.H.Schultz, SIAM J.Sci.Stat.Comput.
 vol7(1986) pp.856-869
 implemented by Kengo Ichiki @ Caltech

 $Log: gmres.c,v $
 Revision 1.3  1998/06/20 22:31:02  ichiki
 now debuging!

 Revision 1.2  1998/06/20 19:18:45  ichiki
 make portable subroutine 'myatimes', which multiply vector to matrix.

 Revision 1.1  1998/06/20 06:01:59  ichiki
 Initial revision

**************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "nrutil.h"
#include "mygmres.h"

/* q[] は不要、g[] をその場で計算すればよい */
/* r[] は対角部分とその下１列が H[] (h_i,j) と違うだけ */
/* しかも下１列は回転で 0 になるので、対角部分だけ保存しておけばよい */
void mygmres(unsigned long n, double f[], double x[],
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
    hv,
    rr,hh,
    r1,r2,
    g0,
    *tmp,*v,*h,*g,*c,*s;

  tmp = dvector(0,n-1);
  v   = dvector(0,n*(m+1)-1);
  h   = dvector(0,m*m-1);
  g   = dvector(0,m);
  c   = dvector(0,m-1);
  s   = dvector(0,m-1);

  /* 1. start: */
  /* compute r0 */
  /* tmp = A.x0 */
  myatimes(n,x,tmp);
  for(i=0;i<n;i++){
    tmp[i] = f[i]-tmp[i]; /* r0 */
  }
  /* compute v1 */
  g[0] = norm(n,tmp); /* beta */
  for(i=0;i<n;i++){
    v[0*n+i] = tmp[i]/g[0];
  }
  /* main loop */
  while( *iter <= itmax ){
    ++(*iter);
    /* 2. iterate: */
    for(j=0;j<m;j++){
      /* tmp = A.vj */
      myatimes(n,&v[j*n],tmp);
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
      hh = norm(n,tmp); /* h_j+1,j */
      /* v_j+1 */
      for(k=0;k<n;k++){
	v[(j+1)*n+k] = tmp[k]/hh;
      }
      /* rotate */
      for(i=0;i<j;i++){
	r1 = h[ i   *m+j];
	r2 = h[(i+1)*m+j];
	h[ i   *m+j] = c[i]*r1-s[i]*r2;
	h[(i+1)*m+j] = s[i]*r1+c[i]*r2;
      }
      rr = h[j*m+j];
      hv = sqrt(rr*rr+hh*hh); /* temporary variable */
      c[j] =  rr/hv;
      s[j] = -hh/hv;
      h[j*m+j] = hv; /* resultant (after rotated) element */

      g0 = g[j];
      g[j  ] = c[j]*g0;
      g[j+1] = s[j]*g0;
    }
    /* 3. form the approximate solution */
    /* solve y_k */
    back_sub(m,h,g,c); /* use c[] as y_k */
    /* x_m */
    for(i=0;i<n;i++){
      for(k=0;k<m;k++){
	x[i] += v[k*n+i]*c[k];
      }
    }
    /* 4. restart */
    *err = g[m]; /* residual */
    /* if satisfied, */
    if( *err <= tol ) break;
    /* else */
    /* compute r_m */
    /* tmp = A.x_m */
    myatimes(n,x,tmp);
    for(i=0;i<n;i++){
      tmp[i] = f[i]-tmp[i]; /* r0 */
    }
    /* compute v1 */
    g[0] = norm(n,tmp);
    for(i=0;i<n;i++){
      v[0*n+i] = tmp[i]/g[0];
    }
  }

  free_dvector(tmp,0,n-1);
  free_dvector(v  ,0,n*m-1);
  free_dvector(h  ,0,m*m-1);
  free_dvector(g  ,0,m);
  free_dvector(c  ,0,m-1);
  free_dvector(s  ,0,m-1);
}


void back_sub(unsigned long m, double r[], double g[], double y[])
{
  unsigned long i,j,jj;

  /*for(j=m-1;j>=0;j--){*/
  /* above for-loop fail, because j is unsigned!! */
  for(jj=0;jj<m;jj++){
    j = m-1-jj;
    y[j] = 0.0;
    for(i=j+1;i<m;i++){
      y[j]-=r[j*m+i]*y[i];
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

void myatimes(unsigned long n, double x[], double y[])
{
  extern double ITER_A[];
  unsigned long i,j;

  for(i=0;i<n;i++){
    y[i]= 0.0;
    for(j=0;j<n;j++){
      y[i] += ITER_A[i*n+j]*x[j];
    }
  }
}
