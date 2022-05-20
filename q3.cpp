#include <cstdio>
# include <iostream>
# include <cmath>
# include <vector>
# include <cstring>
using namespace std;

typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;


double rand_double(double min, double max)
{
	return min + (max - min) * (double)rand() / RAND_MAX ;
}



double dot(VECTOR &a, VECTOR &b)
{
	double sum = 0;
	for (int i = 0; i < a.size(); i++)
		sum += a[i] * b[i];
	return sum;
}


// D = diffusivity, ax0 = boundry for x = x0, bxL= boundry for x = xL, gt0 = boundry for t = t0
VECTOR Forward_Diffusion(string file_name,double(*ax0)(double),double(*bxL)(double),double(*gt0)(double),double D,double t0,double t_end,double x0,double xL,int Nt, int Nx)
{
	FILE *file = fopen(&file_name[0],"w"); // print the t,x,u(x,t) into this file
	fprintf(file,"# %s %s %s\n","x","t","u(x,t)");
	double dt = (t_end-t0)/Nt; // dt
	double dx = (xL-x0)/Nx; // dx
	double alpha = D*dt/dx/dx; // alpha
	if (alpha > 0.5) cout << "Warning: Numerical Instability" << endl;
    VECTOR V(Nx+1); 
	for (int i =0;i<Nx+1;i++){
		V[i] = gt0(x0+dx*i); // u(x,0) = gt0(x)
		fprintf(file,"%lg %lg %lg\n",x0+i*dx,t0,V[i]);
	}

	double vj_1; //to remember V[j]
	for (int i=0;i<Nt+1;i++){
		vj_1 = V[0];
		// In this iterations the values in vector V iteself changes to new t value
		for (int j = 0; j<Nx+1;j++){
			if (j == 0){ 
				V[j] = ax0(t0+i*dt); // 
				fprintf(file,"%lg %lg %lg\n",x0,t0+(i+1)*dt,V[j]);
			}
			else if(j == Nx) {
				V[j] = bxL(t0+i*dt);
				fprintf(file,"%lg %lg %lg\n",xL,t0+(i+1)*dt,V[j]);
			}
			else{
				double temp = V[j];
				V [j] = (1-2*alpha) * V[j] + alpha*(V[j+1]+vj_1);
				vj_1 = temp;
				fprintf(file,"%lg %lg %lg\n",x0+j*dx,t0+(i+1)*dt,V[j]);
			}
		}
	}
	fclose(file); // returns the u(x,t_end)
	return V;
}


// Returns the boundry condition u(0,t)
double ax0(double t){
	return 0;
}

// Returns the boundry condition u(L,t)
double bxL(double t){
	return 0;
}

// Returns the boundry condition u(L,t)
double gt0(double x){
	return 20*abs(sin(M_PI*x));
}


int main(){
	Forward_Diffusion("q3_forward.txt", &ax0,&bxL, &gt0, 1, 0, 4,  0,  2,  5000, 20);

}
/*
	The diffusion equation represents the heat at given space and time. There was intially a distribution of heat over the space given by 20|sin(pi*x)|. With time the heat. So clearly at x = 1 there is a sharp point in the function. However this sharpness vanishes over time (ref-q3_diffusion.png). Physically in case of diffusion the heat flows thowards low heat places. Hence over time the shapness is reduced. Also over time the system will be in equilibrium, so the local ups and downs starts to vanish.
*/