/*
    Assume that the wire is across the z axis, also assume cylindrical coordinate system. The center of the wire is at origin. So the end points of the wires are at -1 and 1. Let a point P at distance 1 radially from the wire. So potential due to length element dz at P is given by dz/(z^2+1). So to get the total potential we have to integrate from -1 to 1.
*/

#include<iostream>
#include<vector>
#include<cmath>
using namespace std;


typedef vector<double> VECTOR;


// Use the weight and abscissa to calculate the integration of the function
double Gauss_Qad(double (*f)(double),double x1,double x2,VECTOR &roots,VECTOR &weights){
	double sum = 0;
	double m = (x2-x1)/2; // m and c are used to transform the [-1,1] to any arbitary [a,b]
	double c = (x2+x1)/2;
	for(int i = 0; i < roots.size(); i++)
		sum += m*weights[i]*f(m*roots[i]+c);
	return sum;
}

double f(double z){
	return 1/sqrt(z*z+1);
}

int main(){
    double result4,result5,result6;
    VECTOR roots,weights;
	//Degree 4 quadrature
    roots = {0.861136311,0.339981043,-0.339981043,-0.861136311};
    weights = {0.347854845,0.652145154,0.652145154,0.347854845};
    result4 = Gauss_Qad(&f,-1,1,roots,weights);
    printf("Integration with 4th order quadrature, Result =  %0.9lf\n",result4);

    // Degree 5 Quadrature
    roots = {0.906179845,0.538469210,0,-0.538469210,-0.906179845};
    weights = {0.236926885,0.478628670,0.56888889,0.478628670,0.236926885};
    result5 = Gauss_Qad(&f,-1,1,roots,weights);
    printf("Integration with 5th order quadrature, Result =  %0.9lf\n",result5);
    
    // Degree 6 quadrature
    roots = {0.932469514,0.661209386,0.238619186,-0.238619186,-0.661209386,-0.932469514};
    weights = {0.171324492,0.360761573,0.467913934,0.467913934,0.360761573,0.171324492};
    result6 = Gauss_Qad(&f,-1,1,roots,weights);
    printf("Integration with 6th order quadrature, Result =  %0.9lf\n",result6);
    // double analytical = log ((sqrt(2)+1)/(sqrt(2)-1));
    cout << endl << "Relative error :"<<endl;
    cout << "Difference between 4th order and 5th order = "<< abs(result5 - result4)<<endl;
    cout << "Difference between 5th order and 6th order = "<< abs(result6 - result5)<<endl;
    cout << "Difference between 6th order and 4th order = "<< abs(result4 - result6)<<endl;
	return 0;
}

/*
Integration with 4th order quadrature, Result =  1.762054179
Integration with 5th order quadrature, Result =  1.762855332
Integration with 6th order quadrature, Result =  1.762730048

Relative error :
Difference between 4th order and 5th order = 0.000801153
Difference between 5th order and 6th order = 0.000125283
Difference between 6th order and 4th order = 0.00067587
*/