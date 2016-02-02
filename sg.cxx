#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda, const double dx, const double dt, const int Nx, const double xmin);
void step(cmplx* const v, const double dt, const double dx, const double omega, const int N, const double xmin);
void writeToFile(const cmplx* const v, const string s, const double dx, const int Nx, const double xmin, const double alpha, const double lambda, const double omega, const double t);
//-----------------------------------
int main(){

	const int Nx = 300;
	const double xmin = -40.0;
	const double xmax = 40.0;
	const double Tend = 10.0*M_PI;

	const double dx = (xmax-xmin)/(Nx-1);
	const double dt = dx/500.0;
 	double t = 0;

	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);

	const double lambda = 10.0;
	const double omega = 0.2;
	const double alpha = sqrt(omega);

  stringstream strm;

	cmplx* psi0 = new cmplx[Nx];

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
			step(psi0,dt,dx,omega,Nx,xmin);
			t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
	delete[] psi0;
	return 0;
}
//-----------------------------------
void step(cmplx* const v, const double dt, const double dx, const double omega, const int N, const double xmin)
{
  cmplx* d=new cmplx[N];
  cmplx* d1=new cmplx[N];
  cmplx* u=new cmplx[N];
  cmplx* l=new cmplx[N];
  cmplx a  = cmplx(0.0,-dt/(4.0*dx*dx));
  cmplx a1 = cmplx(0.0,dt/(4.0*dx*dx));
  cmplx temp1;
  cmplx temp2;
  const double k=omega*omega;
  double x=xmin;

  for(int i=0;i<N-1;i++){
	d1[i] = cmplx(1.0,-(dt/(2.0*dx*dx)+dt*k*x*x/4.0));
	d[i]  = cmplx(1.0, (dt/(2.0*dx*dx)+dt*k*x*x/4.0));
	x+=dx;
	u[i] = a;
	l[i] = a;
  }
  v[0]=d1[0]*v[0]+a1*v[1];
  temp1=v[0];
  for(int i=1;i<N-1;i++){
	temp2=v[i];
	v[i] = a1*temp1+d1[i]*v[i]+a1*v[i+1];
	temp1=temp2;
  }
  v[N-1] = a1*temp1+d1[N-1]*v[N-1];

  for(int i=1;i<N;i++){
    d[i] = d[i] - u[i-1]*l[i]/d[i-1];
    v[i] = v[i] - v[i-1]*l[i]/d[i-1];
    l[i] = 0;
  }
  v[N-1] = v[N-1]/d[N-1];
  for(int i=N-2;i>=0;i--) v[i] = (v[i]-u[i]*v[i+1])/d[i];

  delete[] d;
  delete[] d1;
  delete[] u;
  delete[] l;
}
//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx, const int Nx, const double xmin, const double alpha, const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , -(h2 + h3)  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda, const double dx, const double dt, const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2.0)/2.0 );
	}
}
