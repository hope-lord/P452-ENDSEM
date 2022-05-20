# include <iostream>
# include <vector>
# include <cmath>
using namespace std;

typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX; 

int _a = 572;
int _m = 16381;
int _seed = 98;


int my_randint(){
	int l = (_a*_seed)%_m;
	_seed = l;
	return l;
}

double my_randdouble(double min,double max){
	return min + (max-min)*(double)my_randint()/_m;
}


MATRIX RANDOM_WALK(int N){
	MATRIX xy(2,VECTOR(N+1));
	xy [0][0] = 0;
	xy [1][0] = 0;
	double r = 1;
	double dx,dy,theta;
	for (int i=1;i<=N;i++){
		theta = my_randdouble(0,2*M_PI);
		dy = r * sin(theta);
		dx = r* cos(theta);
		xy [0][i] = xy [0][i-1] + dx;
		xy [1][i] = xy [1][i-1] + dy;
	}
	return xy;
}



int main(){
	int nSims = 500;
	int N = 200;
	MATRIX XY;
	double disp2 = 0;
	for (int i=0;i<nSims;i++){
		XY = RANDOM_WALK(N);
		disp2 += pow(XY[0][N],2)+pow(XY[1][N],2); 
	}
	double rms = sqrt(disp2/nSims);
	cout << "RMS distance = "<< rms<<endl;
	cout << "SQRT(N) = "<< sqrt(N)<<endl;
	cout << "Absolute Difference between them =  "<< abs(sqrt(N)-rms)<<endl;
}
/*
RMS distance = 13.9414
SQRT(N) = 14.1421
Absolute Difference between them =  0.200727
*/