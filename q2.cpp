# include "Linear.cpp"
/* # include "cond.cpp" */

double phi(double x,int j){
	if (j==1) return x;
	if (j==2) return (3*x*x-1)/2;
    if (j==3) return (5*x*x*x-3*x)/2;
    if (j==4) return (35*x*x*x*x-30*x*x+3)/8;
    if (j==5) return (63*x*x*x*x*x-70*x*x*x+15*x)/8;
    if (j==6) return (231*x*x*x*x*x*x-315*x*x*x*x+105*x*x-5)/16;
	return 1;
}

void Legendre_poly_fit(const int n,VECTOR &xx,VECTOR &yy,VECTOR & ww){

    MATRIX A(n+1,VECTOR(n+1,0));
    const int N = xx.size();
    VECTOR x_vec(2*n+1);
    VECTOR yx_vec(n+1);
	double S_yy;
    // read data files
    for(int i=0;i<N;i++){
		if (ww[i]==0) ww[i] = 1;
        else ww[i] = 1/ww[i]/ww[i];
        /* sum_x += x;sum_y += y;sum_xy += x*y;sum_xx += x*x;sum_yy += y*y; */
		S_yy += yy[i]*yy[i]*ww[i];
        for(int j=0;j<n+1;j++){
            // x_vec[j] = pow(xx[i],j);
            yx_vec[j] += yy[i]*phi(xx[i],j)*ww[i];
            for(int k=j;k<n+1;k++){
                A[j][k] += ww[i]*phi(xx[i],k)*phi(xx[i],j);
                if(k!=j) A[k][j] = A[j][k];   
                // display_matrix(A); 
                // cout<<endl;
            }
            
        }
    }
    
	double S_xx = A[1][1];
	double S_xy = yx_vec[1];
    // solve linear equation
    A = inv(A);
    VECTOR coeff = matrix_mult(A,yx_vec);
    VECTOR coeff_err(n+1);
    for(int i=0;i<n+1;i++)
        coeff_err[i] = sqrt(A[i][i]);

    double chi2 = 0;
    for(int i=0;i<N;i++){
        double y_fit = 0;
		for(int l = 0;l<n+1;l++) y_fit += coeff[l]*phi(xx[i],l);
		/* polynomial(xx[i],coeff); */
        chi2 += (yy[i]-y_fit)*(yy[i]-y_fit)*ww[i];
    }
    double r = S_xy/S_xx/S_yy;

    // fprintf(stdout,"Pearson's r = %lg\n",r); 
    fprintf(stdout,"𝜒2/ν = %lg/%d\n",chi2,N-n-1);
    fprintf(stdout,"Coefficients:\n");
    for(int i=0;i<n+1;i++)
        fprintf(stdout,"a%d = %lg \n",i,coeff[i]);

    // // Covariance matrix
    // fprintf(stdout,"Covariance matrix:\n");
    // for(int i=0;i<n+1;i++){
    //     for(int j=0;j<n+1;j++){
    //         fprintf(stdout,"%lg ",A[i][j]);
    //     }
    //     fprintf(stdout,"\n");
// }

}

int main(){
	int N = 26;
    // read file
	FILE * file = fopen("esem4fit.txt","r");
	VECTOR xx(N),yy(N),ww(N);
	for(int i=0;i<N;i++){
		fscanf(file,"%lg %lg",&xx[i],&yy[i]);
		ww[i] = 1;
	}
    cout <<"The printed coefficients are the coefficents of Legendre polynomial for the fit.\n";
	cout << "Assuming 4th order polynomial :\n";
    Legendre_poly_fit(4,xx,yy,ww);
	return 0;
}

/*
Assuming 4th order polynomial :
𝜒2/ν = 0.00280429/21
Coefficients:
a0 = 0.0696578 
a1 = 0.00362402 
a2 = -0.0120826 
a3 = 0.0114262 
a4 = 0.110492 
*/