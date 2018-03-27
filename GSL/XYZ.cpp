//gcc -lgsl -lgslcblas prog.c -o prog

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int func (double t, const double xyz[], double f[], void *params) {

  double w = *(double *)params;

  //Transformation to first order.
  //Equation Ux = x'
  f[0] = xyz[1];
  //Equation Uy = y'
  f[2] = xyz[3];
  //Equation dUx = 2 * w * y' - 3 * w * w * x;
  f[1] = 2 * w * xyz[3] + 3 * (w * w) * xyz[0];
  //Equation dUy = -2 * w * x'
  f[3] = -2 * w * xyz[1];
  //Equation Q = z'
  f[4] = xyz[5];
  //Equation dQ = - z * w * w;
  f[5] = -xyz[4] * w * w;

  return GSL_SUCCESS;
}

int jac (double t, const double z[], double *dfdy, double dfdt[], void *params) {

  //Vetor de elementos da derivada.
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  dfdt[2] = 0.0;
  dfdt[4] = 0.0;
  dfdt[5] = 0.0;

  return GSL_SUCCESS;
}

int main (void) {

  //Tipo de dado que define um sistema de EDO generico com parametros arbitrarios
  double w = 0.15;

  gsl_odeiv2_system sys = {func, jac, 6, &w};

  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 0.5, 1e-6, 10.0);

  int i;
  double t = 0.0, t1 = 10.0;
  double xyz[6] = { 1.0 , 1.0, 1.0, 1.0, 1.0, 1.0};

  for (i = 0; i <= 10; i++) {
      double ti = i * t1 / 20.0;

      //evolui o sistema norteador d iniciando a evolução em t e terminando em t1
      int status = gsl_odeiv2_driver_apply (d, &t, ti, xyz);

      if (status != GSL_SUCCESS) {
    	  printf ("error, return value=%d\n", status);
    	  break;
    	}

      printf("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, xyz[3], xyz[2], xyz[1], xyz[2], xyz[0], xyz[4] );
  }

  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z" );

  gsl_odeiv2_driver_free (d);
  return 0;
}
