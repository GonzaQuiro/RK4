/*					/
/					/
/ CODE CREATE BY DR GONZALO QUIROGA     /
/					/
/				       */
#include <iostream>
#include <math.h>
#include <fstream>
#include <iomanip>
#include <string> //def string
#include <sstream> //convert to string
#include <stdlib.h>  

#define dim 2
void rhs(double l, double m, double I, double mu,int n, double x[], double dx[]);
int sgn(double v);

using namespace std;

int main()
{

	double x[dim], dx[dim], k1[dim], k2[dim], k3[dim], k4[dim];
	double t, l, h, m, mu, I;
	double x_aux[dim];
	int n;

	int count, loop, Max_loop;

	cout << "Initial position: "; 
	cin >> x[0];
	cout << "Initial velocity: "; 
	cin >> x[1];
	cout << "Integer Pow: ";
	cin >> n;

	string number; 
	stringstream ss;
	ss << n;
	number = ss.str();
	
	string filename = "data"+number+".dat";
	string filename2= "energy"+number+".dat";
	
	ofstream output;
	ofstream output2;
	output.open(filename.c_str());
	output2.open(filename2.c_str());
	
	cout << "Final time: "; 
	cin >> t;
	cout << "Pendulum length: "; 
	cin >> l;
	cout << "Pendulum mass: ";
	cin >> m;
	cout << "Damping Coefficient: ";
	cin >> mu;
	//Integration step	
	h=0.001;
	//Inertia
	I=m*l*l;	
	
	Max_loop = t / h;
	loop = 0;
	while (loop <= Max_loop) //Heart of RK4
	{
		output << fixed << setprecision(4) << loop*h << " "<< fixed << setprecision(16)<< x[0]<< " "<< x[1]<< endl;		
		output2 << fixed << setprecision(4) << loop*h << " "<< fixed << setprecision(16)<< 0.5*I*x[1]*x[1]+m*9.8*l*(1-cos(x[0]))<< endl;
		
		rhs(l,m,I,mu,n, x, dx);

		for (count = 0; count <dim; count++)
		k1[count] = h*dx[count];

		for (count = 0; count <dim; count++)
		x_aux[count] = x[count] + 0.5*k1[count];

		rhs(l,m,I,mu,n, x_aux, dx);

		for (count = 0; count <dim; count++) 
		k2[count] = h*dx[count];

		for (count = 0; count <dim; count++)
		x_aux[count] = x[count] + 0.5*k2[count];

		rhs(l,m,I,mu,n, x_aux, dx);

		for (count = 0; count <dim; count++)
 		k3[count] = h*dx[count];

		for (count = 0; count <dim; count++)
		x_aux[count] = x[count] + k3[count];

		rhs(l,m,I,mu,n, x_aux, dx);

		for (count = 0; count <dim; count++)
		k4[count] = h*dx[count];

		for (count = 0; count <dim; count++)
		x[count] = x[count] + (k1[count] + 2 * k2[count]+ 2 * k3[count] + k4[count]) / 6;
		
		loop++;
	}
	output.close();
	output2.close();
	return 0;
}

void rhs(double l, double m, double I, double mu, int n, double x[], double dx[]) //
{		
	dx[0] = x[1];
	dx[1] = -sgn(x[1])*mu*pow(x[1]*x[1],n/2.0)-((m*9.8*l)/I)*sin(x[0]);
}

int sgn(double v) {
  return (v < 0) ? -1 : ((v > 0) ? 1 : 0);
}
