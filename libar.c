#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// compile to make shared library
// gcc -fPIC -c asym_rotor.c asym_rotor.o
// gcc -shared -lgsl -lgslcblas -Wl,-soname,libar.so.1 -o libar.so.1 asym_rotor.o

// parameters of the derivative
typedef struct {
  double t0, t1, F, mu1, mu2, mu3, I1, I2, I3; 
} params;

// evaluate the right hand side of the state variable
//
// as written in code this routine should be completely opaque
// for a clearer explanation and justification see
// for example: Ph. Dugourd et. al. Chem. Phys. Lett. 423, (2006) pp. 13-16

int 
dqt(double t, const double qt[], double dqt[], void *par) {
  params *st = (params *)par; 

  double t0 = st->t0;
  double t1 = st->t1;
  
  double F = 0;
  double mu1 = st->mu1, mu2 = st->mu2, mu3 = st->mu3; 
  double I1 = st->I1, I2 = st->I2, I3 = st->I3;

  double q0 = qt[0], q1 = qt[1], q2 = qt[2], q3 = qt[3];
  double w1 = qt[4], w2 = qt[5], w3 = qt[6];

  // normalize q
  //double qn = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
  //q0 /= qn; q1 /= qn; q2 /= qn; q3 /= qn;

  //printf("qn = %f\n", qn);

  // we have a problem with dugourd
  //double dq0 = 0.5 * (-q1*w1 - q2*w2 - q3*w3);
  //double dq1 = 0.5 * (q0*w1 - q3*w2 + q2*w3);
  //double dq2 = 0.5 * (q3*w1 + q0*w2 - q1*w3);
  //double dq3 = 0.5 * (-q2*w1 + q1*w2 + q0*w3);
  
  double dq0 = 0.5 * (-q1*w1 - q2*w2 - q3*w3);
  double dq1 = 0.5 * (q0*w1 + q3*w2 - q2*w3);
  double dq2 = 0.5 * (-q3*w1 + q0*w2 + q1*w3);
  double dq3 = 0.5 * (q2*w1 - q1*w2 + q0*w3);

  // now compute the electric field vector in the body fixed axes
  // the orientation quaternion is converted to a rotation matrix Q
  // and the transpose of this matrix is multiplied by the field
  // vector in the space fixed axes

  double f1 = 0;
  double f2 = 0;
  double f3 = 0;

  // slowly increase the field to F between t0 and t1
  if(t < t0) {
    F = 0.0;
  } else if(t < t1) {
    F = st->F * (t - t0)/(t1 - t0);
  } else {
    F = st->F;
  }

  f1 = F * 2.0 * (q1*q3 + q0*q2);
  f2 = F * 2.0 * (q2*q3 - q0*q1);
  f3 = F * (q0*q0 - q1*q1 - q2*q2 + q3*q3);

  double dw1 = (I2-I3)*w2*w3/I1 - (mu3*f2 - mu2*f3)/I1;
  double dw2 = (I3-I1)*w3*w1/I2 - (mu1*f3 - mu3*f1)/I2;
  double dw3 = (I1-I2)*w1*w2/I3 - (mu2*f1 - mu1*f2)/I3;
  
  dqt[0] = dq0; dqt[1] = dq1; dqt[2] = dq2; dqt[3] = dq3;
  dqt[4] = dw1; dqt[5] = dw2; dqt[6] = dw3;
  
  return GSL_SUCCESS;
}

int jac(double t, const double qt[], double *dfdq, double dfdt[], void *par) {
    params *st = (params *)par; 

    double t0 = st->t0;
    double t1 = st->t1;
    double tau = t1 - t0;
        
    double F = 0;
    double mu1 = st->mu1, mu2 = st->mu2, mu3 = st->mu3; 
    double I1 = st->I1, I2 = st->I2, I3 = st->I3;

    double q0 = qt[0], q1 = qt[1], q2 = qt[2], q3 = qt[3];
    double w1 = qt[4], w2 = qt[5], w3 = qt[6];

    // slowly increase the field to F between t0 and t1
    if(t < t0) {
        F = 0.0;
    } else if(t < t1) {
        F = st->F * (t - t0)/tau;
    } else {
        F = st->F;
    }
    
    // it is necessary to use the jacobean - we wish to understand the stability 
    // of the different trajectories so this will be necessary anyway..
    
    printf(".1.");

    gsl_matrix_view dfdq_mat = gsl_matrix_view_array(dfdq, 7, 7);
    gsl_matrix *m = &dfdq_mat.matrix;

    gsl_matrix_set(m, 0, 0, 0);
    gsl_matrix_set(m, 0, 1, w1);
    gsl_matrix_set(m, 0, 2, w2);
    gsl_matrix_set(m, 0, 3, w3);
    gsl_matrix_set(m, 0, 4, q1);
    gsl_matrix_set(m, 0, 5, q2);
    gsl_matrix_set(m, 0, 6, q3);
    
    gsl_matrix_set(m, 1, 0, w1);
    gsl_matrix_set(m, 1, 1, 0);
    gsl_matrix_set(m, 1, 2, -w3);
    gsl_matrix_set(m, 1, 3, w2);
    gsl_matrix_set(m, 1, 4, q0);
    gsl_matrix_set(m, 1, 5, q3);
    gsl_matrix_set(m, 1, 6, -q2);
    
    gsl_matrix_set(m, 2, 0, w2);
    gsl_matrix_set(m, 2, 1, w3);
    gsl_matrix_set(m, 2, 2, 0);
    gsl_matrix_set(m, 2, 3, -w1);
    gsl_matrix_set(m, 2, 4, -q3);
    gsl_matrix_set(m, 2, 5, q0);
    gsl_matrix_set(m, 2, 6, q1);

    gsl_matrix_set(m, 3, 0, w3);
    gsl_matrix_set(m, 3, 1, -w2);
    gsl_matrix_set(m, 3, 2, w1);
    gsl_matrix_set(m, 3, 3, 0);
    gsl_matrix_set(m, 3, 4, q2);
    gsl_matrix_set(m, 3, 5, -q1);
    gsl_matrix_set(m, 3, 6, q0);

    gsl_matrix_set(m, 4, 0, mu3*2*F/I1 *  q1 + mu2*2*F/I1 *  q0);
    gsl_matrix_set(m, 4, 1, mu3*2*F/I1 *  q0 + mu2*2*F/I1 * -q1);
    gsl_matrix_set(m, 4, 2, mu3*2*F/I1 * -q3 + mu2*2*F/I1 * -q2);
    gsl_matrix_set(m, 4, 3, mu3*2*F/I1 * -q2 + mu2*2*F/I1 *  q3);
    gsl_matrix_set(m, 4, 4, 0);
    gsl_matrix_set(m, 4, 5, w3*(I2-I3)/I2);
    gsl_matrix_set(m, 4, 6, w2*(I2-I3)/I2);

    gsl_matrix_set(m, 5, 0, mu1*2*F/I2 * -q0 + mu3*2*F/I2 *  q2);
    gsl_matrix_set(m, 5, 1, mu1*2*F/I2 *  q1 + mu3*2*F/I2 *  q3);
    gsl_matrix_set(m, 5, 2, mu1*2*F/I2 *  q2 + mu3*2*F/I2 *  q0);
    gsl_matrix_set(m, 5, 3, mu1*2*F/I2 * -q3 + mu3*2*F/I2 *  q1);
    gsl_matrix_set(m, 5, 4, w3*(I3-I1)/I2);
    gsl_matrix_set(m, 5, 5, 0);
    gsl_matrix_set(m, 5, 6, w1*(I3-I1)/I2);

    gsl_matrix_set(m, 6, 0, mu2*2*F/I3 * -q2 + mu1*2*F/I3 * -q1);
    gsl_matrix_set(m, 6, 1, mu2*2*F/I3 * -q3 + mu1*2*F/I3 * -q0);
    gsl_matrix_set(m, 6, 2, mu2*2*F/I3 * -q0 + mu1*2*F/I3 *  q3);
    gsl_matrix_set(m, 6, 3, mu2*2*F/I3 * -q1 + mu1*2*F/I3 *  q2);
    gsl_matrix_set(m, 6, 4, w2*(I1-I2)/I3);
    gsl_matrix_set(m, 6, 5, w1*(I1-I2)/I3);
    gsl_matrix_set(m, 6, 6, 0);

    dfdt[0] = 0.0; dfdt[1] = 0.0; dfdt[2] = 0.0; dfdt[3] = 0.0;
    dfdt[4] = 0.0; dfdt[5] = 0.0; dfdt[6] = 0.0;
    
    if(t < t1 && t > t0) {
        F = st->F * (t - t0)/tau;
        dfdt[4] = -mu3*2*F/I1/tau * (q2*q3 - q0*q1) +
            mu2*F/I1/tau * (q0*q0 - q1*q1 - q2*q2 + q3*q3);
        dfdt[5] = -mu1*F/I2/tau * (q0*q0 - q1*q1 - q2*q2 + q3*q3) +
            mu3*2*F/I2/tau * (q1*q3 + q0*q2);
        dfdt[6] = -mu2*2*F/I3/tau * (q1*q3 + q0*q2) +
            mu1*2*F/I3/tau * (q2*q3 - q0*q1);
    }
   
    return GSL_SUCCESS;
}


// numerically integrate the motion of an asymmetric rotor
// 
int asym_rotor(const double t0, const double t1, const double t2, 
	       const double qt0[], 
	       double *qts, 
	       long max_tstep,
	       const double Fmax,
	       const double mu1, const double mu2, const double mu3,
	       const double I1, const double I2, const double I3) {

  // allocate a stepper that uses Runge-Kutta Cash-Karp (4,5) step size
  // this should not require a jacobean matrix
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk4imp;
  
  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc(T, 7);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new(1e-6, 1e-6);

  gsl_odeiv_evolve * e 
     = gsl_odeiv_evolve_alloc(7);

  // initial parameters
  params pi = { t0, t1, Fmax, mu1, mu2, mu3, I1, I2, I3 };
  
  gsl_odeiv_system sys = { dqt, jac, 7, &pi };
  
  double t = 0.0;
  double h = 1e-8;
  
  double qt[7] = { qt0[0], qt0[1], qt0[2], qt0[3], qt0[4], qt0[5], qt0[6] };

  long tstep = 0;

  printf("Error control erl/abs = 1e-6\n");


  while(t < t2) {
    
    if(tstep >= max_tstep)
      break;

    int status = gsl_odeiv_evolve_apply(e, c, s,
					&sys,
					&t, t2,
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

// numerically integrate the motion of an asymmetric rotor
// and calculate the projected dipole moment after the field is turned on
// this routine is much faster if all we want to know is the time averaged
// field projection - once this works we can apply a convergence test
int asym_rotor_muz(const double t0, const double t1, const double t2, 
	       const double qt0[], 
	       double *muz_bar, 
	       long max_tstep,
	       const double Fmax,
	       const double mu1, const double mu2, const double mu3,
	       const double I1, const double I2, const double I3) {

  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk4imp;
  
  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc(T, 7);
  gsl_odeiv_control * c 
    = gsl_odeiv_control_y_new(1e-6, 1e-6);

  gsl_odeiv_evolve * e 
     = gsl_odeiv_evolve_alloc(7);

  // initial parameters
  params pi = { t0, t1, Fmax, mu1, mu2, mu3, I1, I2, I3 };
  
  gsl_odeiv_system sys = { dqt, jac, 7, &pi };
  
  double t = 0.0;
  double h = 1e-8;
  
  double qt[7] = { qt0[0], qt0[1], qt0[2], qt0[3], qt0[4], qt0[5], qt0[6] };

  long tstep = 0;

  printf("Error control erl/abs = 1e-6\n");

  // init and field turnon
  while(t < t1) {
    if(tstep >= max_tstep)
      break;

    int status = gsl_odeiv_evolve_apply(e, c, s,
					&sys,
					&t, t2,
					&h, qt);

    if(status != GSL_SUCCESS)
      break;
  }

  double muz_avg = 0;
  double n = 0;

  // now take a time average of muz until t2
  while(t < t2) {

    if(tstep >= max_tstep)
      break;

    int status = gsl_odeiv_evolve_apply(e, c, s,
					&sys,
					&t, t2,
					&h, qt);

    if(status != GSL_SUCCESS)
      break;
 
    double muz = mu1 * 2 * (qt[1]*qt[3] + qt[0]*qt[2]) +
                 mu2 * 2 * (qt[2]*qt[3] - qt[0]*qt[1]) + 
                 mu3 * (qt[0]*qt[0] - qt[1]*qt[1] - qt[2]*qt[2] + qt[3]*qt[3]);

    // compute a running average
    muz_avg = (muz + n*muz_avg)/(n+1);

    // printf("tstep %ld: t = %f\n", tstep, t);

    // print out the state
    //printf("%.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
    //   t, qt[0], qt[1], qt[2], qt[3], qt[4], qt[5], qt[6]);
    
    tstep++;
  }

  *muz_bar = muz_avg;  

  gsl_odeiv_evolve_free(e);
  gsl_odeiv_control_free(c);
  gsl_odeiv_step_free(s);

  return tstep;
}


int test_byref(double *out) {
    *out = 5.0;
    return 1;
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

