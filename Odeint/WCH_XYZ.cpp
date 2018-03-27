#include <iostream>
#include <vector>
#include <time.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <omp.h>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>


using namespace boost::numeric::odeint;
using namespace std;

typedef std::vector< double > state_type;
//typedef boost::array< double , 6 > state_type;

//typedef runge_kutta4< state_type > rk4;

//Paraller runge_kutta4
typedef runge_kutta4<
                  state_type , double ,
                  state_type , double ,
                  openmp_range_algebra
                > rk4;

double w = 0.15;
double x = 1.0, dx = 1.0, y = 1.0, dy = 1.0;

//const size_t N = 0;

struct ode
{
    void operator()( const state_type &XY , state_type &dUdt , double t )
    {

      const size_t N = XY.size();

      //#pragma omp parallel for schedule(runtime)
      //for(size_t i = 1 ; i < N - 1 ; ++i)
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
};

void observer( const state_type &XY , const double t ) {

    double xt = (dx / w) * sin(w * t) - ((2*dy/w) + 3*x) * cos( w * t ) + ((2*dy / w) + 4 * x);
    //double xt = (dx / w) * sin(w * t) - (3*x +(2*dy)/w)*cos(w*t) + (2/w)*(2*w*x + dy);
    //double yt =  ((2 * dx) / w) * cos(w*t) - ((4*dy)/w + 6*x)*sin(w*t) + (y - ((2 * dx)/w)) - (3 * y + 6 * w * x ) * t;
    double yt =  ((2 * dx) / w) * cos(w*t) + 2 * (3*x + (2*dy)/w)*sin(w*t) - (2*dx)/w - 3*(2*w*x + dy)* t + y;

    double zt = cos(w * t) + (1 / w) * sin(w * t);
    //printf("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", t, XY[3], XY[2], XY[1], XY[2], XY[0], XY[4], yt, xt, zt );

    if (omp_get_num_threads() > 1) {
      printf("num_threads in ode: %i\n", omp_get_num_threads() );
    }
}

//2 = y'
//0 = x
//4 = z

int main(int argc, char** argv)
{
    using namespace std;
    using namespace boost::numeric::odeint;

    clock_t tempoInicial, tempoFinal;
    double tempoGasto;

    if (argc < 1) {
        printf("Usage : WCH_XYZ <Number of threads>\n");
        exit(EXIT_FAILURE);
    }

    int number_threads = atoi(argv[1]);

    state_type XY(6);

    XY[0] = 1.0;
    XY[1] = 1.0;
    XY[2] = 1.0;
    XY[3] = 1.0;
    XY[4] = 1.0;
    XY[5] = 1.0;

    omp_set_num_threads(number_threads);
    int chunk_size = 1000 / omp_get_max_threads();
    omp_set_schedule( omp_sched_static , chunk_size );

    //printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z", "Y(t)", "X(t)", "Z(t)" );
    tempoInicial = clock();
    //integrate_adaptive( rk4(), ode , XY , 0.0 , 5.0 , 0.001, observer );
    //integrate_n_steps( rk4() , ode , XY , start value , precision , number of intertions, observer);
    //runge_kutta4<state_type>()
    integrate_n_steps( rk4() , ode() , XY , 0.0 , 0.00001 , 2000, observer);
    //integrate_n_steps( runge_kutta4<state_type>() , ode , XY , 0.0 , 0.00001 , 2000, observer);

    tempoFinal = clock();

    tempoGasto = (( tempoFinal - tempoInicial ));
    tempoGasto = (((float)tempoGasto)/CLOCKS_PER_SEC);

    //printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "t", "Y'", "X'", "Z'", "Y", "X", "Z", "Y(t)", "X(t)", "Z(t)" );

    printf("\nTempo em segundos: %.5fs\n\n", tempoGasto);
    return 0;
}
