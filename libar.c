#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// compile to make shared library
// gcc -fPIC -c asym_rotor.c asym_rotor.o
// gcc -shared -lgsl -lgslcblas -Wl,-soname,libar.so.1 -o libar.so.1 asym_rotor.o

// parameters of the derivative
typedef struct {
  double F, mu1, mu2, mu3, I1, I2, I3; 
} params;

// evaluate the right hand side of the state variable
//
// as written in code this routine should be completely opaque
// for a clearer explanation and justification see
// for example: Ph. Dugourd et. al. Chem. Phys. Lett. 423, (2006) pp. 13-16

int 
dqt(double t, const double qt[], double dqt[], void *par) {
  params *st = (params *)par; 

  double F = st->F; // field in space fixed axes 
  double mu1 = st->mu1, mu2 = st->mu2, mu3 = st->mu3; 
  double I1 = st->I1, I2 = st->I2, I3 = st->I3;
  
  double q0 = qt[0], q1 = qt[1], q2 = qt[2], q3 = qt[3];
  double w1 = qt[4], w2 = qt[5], w3 = qt[6];
  
  double dq0 = 0.5 * (-q1*w1 - q2*w2 - q3*w3);
  double dq1 = 0.5 * (q0*w1 - q3*w2 + q2*w3);
  double dq2 = 0.5 * (q3*w1 + q0*w2 - q1*w3);
  double dq3 = 0.5 * (-q2*w1 + q1*w2 + q0*w3);

  // now compute the electric field vector in the body fixed axes
  // the orientation quaternion is converted to a rotation matrix Q
  // and the transpose of this matrix is multiplied by the field
  // vector in the space fixed axes

  double f1 = F * 2.0 * (q1*q3 + q0*q2);
  double f2 = F * 2.0 * (q2*q3 - q0*q1);
  double f3 = F * (q0*q0 - q1*q1 - q2*q2 - q3*q3);
  
  double dw1 = ((I2-I3)*w2*w3 - (mu3*f2 - mu2*f3))/I1;
  double dw2 = ((I3-I1)*w3*w1 - (mu1*f3 - mu3*f1))/I2;
  double dw3 = ((I1-I2)*w1*w2 - (mu2*f1 - mu1*f2))/I3;
  
  dqt[0] = dq0; dqt[1] = dq1; dqt[2] = dq2; dqt[3] = dq3;
  dqt[4] = dw1; dqt[5] = dw2; dqt[6] = dw3;
  return GSL_SUCCESS;
}

// numerically integrate the motion of an asymmetric rotor
// 
int asym_rotor(const double t0, const double tf, 
	       const double qt0[], 
	       double *qts, 
	       long max_tstep,
	       const double Fmax,
	       const double mu1, const double mu2, const double mu3,
	       const double I1, const double I2, const double I3) {

  // allocate a stepper that uses Runge-Kutta Cash-Karp (4,5) step size
  // this should not require a jacobean matrix
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rkck;
  
  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc(T, 7);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new(1e-6, 0.0);
  gsl_odeiv_evolve * e 
     = gsl_odeiv_evolve_alloc(7);

  // initial parameters
  params pi = { Fmax, mu1, mu2, mu3, I1, I2, I3 };
  
  gsl_odeiv_system sys = { dqt, NULL, 7, &pi };
  
  double t = t0;
  double h = 1e-6;
  
  double qt[7] = { qt0[0], qt0[1], qt0[2], qt0[3], qt0[4], qt0[5], qt0[6] };

  long tstep = 0;

  while(t < tf) {
    
    if(tstep >= max_tstep)
      break;

    int status = gsl_odeiv_evolve_apply(e, c, s,
					&sys,
					&t, tf,
					&h, qt);

    if(status != GSL_SUCCESS)
      break;

    qts[tstep*8 + 0] = t;
    qts[tstep*8 + 1] = qt[0];
    qts[tstep*8 + 2] = qt[1];
    qts[tstep*8 + 3] = qt[2];
    qts[tstep*8 + 4] = qt[3];
    qts[tstep*8 + 5] = qt[4];
    qts[tstep*8 + 6] = qt[5];
    qts[tstep*8 + 7] = qt[6];

    // printf("tstep %ld: t = %f\n", tstep, t);

    // print out the state
    //printf("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
    //   t, qt[0], qt[1], qt[2], qt[3], qt[4], qt[5], qt[6]);
    
    tstep++;
  }

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return tstep;
}

/* int main() { */
/*   double qt0[7] = {1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0}; */
/*   asym_rotor(0.0, 100.0, */
/* 	     qt0, */
/* 	     1.0, */
/* 	     0.0, 0.5, 1.0, */
/* 	     1.0, 2.0, 1.5); */
/*   printf("done!\n"); */
/*   return 0; */
/* } */

