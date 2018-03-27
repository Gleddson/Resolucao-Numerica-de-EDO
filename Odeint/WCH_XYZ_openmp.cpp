#include <iostream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <omp.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>


using namespace boost::numeric::odeint;
using namespace std;

typedef boost::array< double , 6 > state_type;
//typedef runge_kutta4< state_type > rk4;

//Paraller runge_kutta4
typedef runge_kutta4<
                  state_type , double ,
                  state_type , double ,
                  openmp_range_algebra
                > rk4;

double w = 0.15;
double x = 1.0, dx = 1.0, y = 1.0, dy = 1.0;

const size_t N = 0;

void ode( const state_type &XY , state_type &dUdt , double t )
{

  #pragma omp parallel for schedule(runtime)
  for (size_t aux = 0; aux <= N; aux++) {

      //Transformation to first order.
      //Equation Ux = x'
      dUdt[0] = XY[1];
      //Equation Uy = y'
      dUdt[2] = XY[3];
      //Equation dUx = 2 * w * y' - 3 * w * w * x;
      dUdt[1] = 2 * w * XY[3] + 3 * (w * w) * XY[0];
      //Equation dUy = -2 * w * x'
      dUdt[3] = -2 * w * XY[1];
      //Equation Q = z'
      dUdt[4] = XY[5];
      //Equation dQ = - z * w * w;
      dUdt[5] = -XY[4] * w * w;
    }
}

void observer( const state_type &XY , const double t ) {

    double xt = (dx / w) * sin(w * t) - ((2*dy/w) + 3*x) * cos( w * t ) + ((2*dy / w) + 4 * x);
    //double xt = (dx / w) * sin(w * t) - (3*x +(2*dy)/w)*cos(w*t) + (2/w)*(2*w*x + dy);
    //double yt =  ((2 * dx) / w) * cos(w*t) - ((4*dy)/w + 6*x)*sin(w*t) + (y - ((2 * dx)/w)) - (3 * y + 6 * w * x ) * t;
    double yt =  ((2 * dx) / w) * cos(w*t) + 2 * (3*x + (2*dy)/w)*sin(w*t) - (2*dx)/w - 3*(2*w*x + dy)* t + y;

    double zt = cos(w * t) + (1 / w) * sin(w * t);
    //printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z", "Y(t)", "X(t)", "Z(t)" );
    printf("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, XY[3], XY[2], XY[1], XY[2], XY[0], XY[4], yt, xt, zt );
}

//2 = y'
//0 = x
//4 = z

int main()
{
    using namespace std;
    using namespace boost::numeric::odeint;

    state_type XY = { 1.0 , 1.0, 1.0, 1.0, 1.0, 1.0};

    //omp_set_num_threads(10);
    int chunk_size = omp_get_max_threads();
    omp_set_schedule( omp_sched_static , chunk_size );

    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z", "Y(t)", "X(t)", "Z(t)" );
    //integrate_adaptive( rk4(), ode , XY , 0.0 , 5.0 , 0.001, observer );

    //integrate_n_steps( rk4() , ode , XY , start value , precision , number of intertions, observer);
    integrate_n_steps( rk4() , ode , XY , 0.0 , 1.0 , 1000, observer);
    printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z", "Y(t)", "X(t)", "Z(t)" );
}
