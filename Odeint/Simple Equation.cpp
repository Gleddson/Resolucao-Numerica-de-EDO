#include <iostream>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

void rhs( const double y , double &dxdt , const double x )
{
    dxdt = y - 2*x/y;
}

void write_cout( const double &y , const double x )
{
    cout << x << '\t' << y << endl;
}

// state_type = double
typedef runge_kutta4< double > stepper_type;

int main()
{
    double y = 1.0;
    double x = 0.0;

    cout << 'x' << '\t' << 'y' << '\n';
    integrate_adaptive( stepper_type() , rhs , y , x , 4.0 , 0.25 , write_cout);
}
