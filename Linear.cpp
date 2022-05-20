# include <iostream>
# include <vector>
# include <cmath>
using namespace std;
typedef vector<double> VECTOR;
typedef vector<VECTOR> MATRIX;

static VECTOR null_vector;
static MATRIX null_matrix;
/*---------------------------------------------------*/
bool partial_pivot(MATRIX &ar,VECTOR &br=null_vector,MATRIX &cr=null_matrix){
	int LIMIT=ar.size();
	for(int i=0;i<LIMIT;i++){
		if (ar[i][i]==0){
			bool isit=true;
			for(int t=0;t<LIMIT;t++){
				if(0!=ar[i][t])isit=false; 
			}
			if(isit) return (!isit);
			for(int j=i+1;i<LIMIT;i++){
				if(ar[j][i]!=0){
					for(int t=0;t<LIMIT;t++){
						swap(ar[j][t],ar[i][t]);
						if (cr!=null_matrix) swap(cr[j][t],cr[i][t]);
					}
					if (br!=null_vector)swap(br[j],br[i]);
					i=LIMIT;
				}
			}
		}
	}
	return true;
}

/*----------------------------------------------------------------*/


void gauss_jordan(MATRIX &ar, 
					VECTOR &br=null_vector,
					MATRIX &cr=null_matrix){
	int LIMIT=ar.size();

	for(int i=0;i<LIMIT;i++){
		if (ar[i][i]==0){
			bool c=partial_pivot(ar,br,cr);
			if(!c) {cout<<"WARNING: DETERMINANT IS 0.\n";return;}
		}
		double pivot=ar[i][i];
		if (br!=null_vector) br[i]/=pivot;
		for(int j=0;j<LIMIT;j++){
			ar[i][j]/=pivot;
			if (cr!=null_matrix) cr[i][j]/=pivot;
		}
		for(int j=0;j<LIMIT;j++){
			if(i!=j and ar[j][i]!=0){
				double fact=ar[j][i];
				for(int t=0;t<LIMIT;t++){
					ar[j][t]-=fact*ar[i][t];
					if (cr!=null_matrix) cr[j][t]-=fact*cr[i][t];
				}
				if(br!=null_vector)br[j]-=fact*br[i];
			}
		}
	}
}

/*-------------------------------------------*/


MATRIX inv(MATRIX &ar){
	int LIMIT=ar.size();
    MATRIX inv;
	for(int i=0;i<LIMIT;i++) {
		VECTOR m;
		for(int j=0;j<LIMIT;j++){
			if(i==j) m.push_back(1);
			else m.push_back(0);
		}
		inv.push_back(m);
	}
	gauss_jordan(ar,null_vector,inv);
	return inv;
}

double dot(VECTOR &ar,VECTOR &br){
	int len=ar.size();
	if(len!=br.size()){cout<<"Dimensions are not equal.\n";return 0;}
	double sum=0;
	for(int i=0;i<len;i++) sum+=ar[i]*br[i];
	return sum;
}

VECTOR matrix_mult(MATRIX &a,VECTOR &b){
	VECTOR axb;
	for(int i=0;i<a.size();i++) axb.push_back(dot(a[i],b));
	return axb;
}