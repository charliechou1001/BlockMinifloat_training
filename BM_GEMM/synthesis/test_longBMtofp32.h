#include <cmath>
#include "ap_int.h"
#include "typedef.h"

float BMToFP(ap_uint<W> a, ap_uint<2> modeOut, SEXP_T shared_exp){
	float two_exp;
	float frac;
	float sexp = std::pow(2,int(shared_exp));
	float res;

	if(modeOut == 0b00){//Res0
		ap_uint<ResM0> man = a(ResM0-1,0);
#if ResS0 !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if ResE0 !=0
		ap_uint<1> denorm = a(W-2,W-ResE0-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<ResE0+1> exp = a(W-2,W-ResE0-1) - (1<<(ResE0-1))+1 + denorm;
//		printf("exp: %d, a(W-2,W-ResE0-1): %d\n",exp.to_int(), a(W-2,W-ResE0-1).to_int());
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-ResM0) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-ResM0+1);
#endif
		res = sgn*frac*two_exp*sexp;
//		printf("frac: %f, two_exp: %f, sexp: %f\n",frac, two_exp, sexp);
	}else if(modeOut == 0b01){//Res1
		ap_uint<ResM1> man = a(ResM1-1,0);
#if ResS1 !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if ResE1 !=0
		ap_uint<1> denorm = a(W-2,W-ResE1-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<ResE1+1> exp = a(W-2,W-ResE1-1) - (1<<(ResE1-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-ResM1) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-ResM1+1);
#endif
		res = sgn*frac*two_exp*sexp;
	}else{//Res2
		ap_uint<ResM2> man = a(ResM2-1,0);
#if ResS2 !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if ResE2 !=0
		ap_uint<1> denorm = a(W-2,W-ResE2-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<ResE2+1> exp = a(W-2,W-ResE2-1) - (1<<(ResE2-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-ResM2) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-ResM2+1);
#endif
		res = sgn*frac*two_exp*sexp;
	}
	return res;
}


void CrossGEMMCompare(ap_uint<W> A[M][N], ap_uint<W> B[M][N], SEXP_T  betaA[BK_M][BK_N], SEXP_T betaB[BK_M][BK_N], ap_uint<2> modeOut){
	float sum =0;
	for(int i=0;i<BK_M;i++){
		for(int j=0;j<BK_N;j++){
			printf("betaA: %d, betaB: %d\n",betaA[i][j].to_int(), betaB[i][j].to_int() );
			for(int k1=0;k1<BLK_SIZE;k1++){
				for(int k2=0;k2<BLK_SIZE;k2++){
					float a = BMToFP(A[i*BLK_SIZE+k1][j*BLK_SIZE+k2],modeOut,betaA[i][j]);
					float b = BMToFP(B[i*BLK_SIZE+k1][j*BLK_SIZE+k2],modeOut,betaB[i][j]);
					printf("%s, %s\n", A[i*BLK_SIZE+k1][j*BLK_SIZE+k2].to_string(2).c_str(), B[i*BLK_SIZE+k1][j*BLK_SIZE+k2].to_string(2).c_str());
					if(a!=b){
						printf("x:%d, y:%d, A: %s, B: %s, a:%f, b:%f\n",i*BLK_SIZE+k1,j*BLK_SIZE+k2,
								A[i*BLK_SIZE+k1][j*BLK_SIZE+k2].to_string(2).c_str(), B[i*BLK_SIZE+k1][j*BLK_SIZE+k2].to_string(2).c_str(), a,b);
						sum += (a-b)*(a-b);
					}
//					if(A[i*TileSize+k1][j*TileSize+k2] !=B[i*TileSize+k1][j*TileSize+k2]){
//						printf("x:%d, y:%d, A: %s, B: %s\n",i*TileSize+k1,j*TileSize+k2,
//						A[i*TileSize+k1][j*TileSize+k2].to_string(2).c_str(), B[i*TileSize+k1][j*TileSize+k2].to_string(2).c_str());
//					}
				}
			}
		}
	}
	sum = sum/(M*N);
	printf("average MSE is: %f\n",sum);
}

void PrintResult(ap_uint<W> A[M][N], SEXP_T  betaA[BK_M][BK_N], ap_uint<2> modeOut){
	for(int i=0;i<BK_M;i++){
		for(int j=0;j<BK_N;j++){
			for(int k1=0;k1<BLK_SIZE;k1++){
				for(int k2=0;k2<BLK_SIZE;k2++){
					float a = BMToFP(A[i*BLK_SIZE+k1][j*BLK_SIZE+k2],modeOut,betaA[i][j]);
					printf("C[%d][%d]: %f,  ",i*BLK_SIZE+k1, j*BLK_SIZE+k2, a);
				}printf("\n");
			}
		}
	}
}
