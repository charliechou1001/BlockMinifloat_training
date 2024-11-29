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

void SingleMatrixCompare(ap_uint<W> A[SinM][SinN], ap_uint<W> B[SinM][SinN], SEXP_T  betaA, SEXP_T betaB, ap_uint<2> modeOut){
	float sum =0;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			float a = BMToFP(A[i][j],modeOut,betaA);
			float b = BMToFP(B[i][j],modeOut,betaB);
			if(a!=b){
				printf("i:%d, j:%d, A: %s, B: %s, a:%f, b:%f\n",i,j,A[i][j].to_string(2).c_str(), B[i][j].to_string(2).c_str(), a,b);
				sum += (a-b)*(a-b);
			}
		}
	}
	sum = sum/(SinM*SinN);
	printf("average MSE is: %f\n",sum);
}

void CrossGEMMCompare(ap_uint<W> A[M][N], ap_uint<W> B[M][N], SEXP_T  betaA[TILE_M][TILE_N], SEXP_T betaB[TILE_M][TILE_N], ap_uint<2> modeOut){
	float sum =0;
	for(int i=0;i<TILE_M;i++){
		for(int j=0;j<TILE_N;j++){
			for(int k1=0;k1<TileSize;k1++){
				for(int k2=0;k2<TileSize;k2++){
					float a = BMToFP(A[i*TileSize+k1][j*TileSize+k2],modeOut,betaA[i][j]);
					float b = BMToFP(B[i*TileSize+k1][j*TileSize+k2],modeOut,betaB[i][j]);
					if(a!=b){
						printf("x:%d, y:%d, A: %s, B: %s, a:%f, b:%f\n",i*TileSize+k1,j*TileSize+k2,
								A[i*TileSize+k1][j*TileSize+k2].to_string(2).c_str(), B[i*TileSize+k1][j*TileSize+k2].to_string(2).c_str(), a,b);
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

float InputBMToFP(ap_uint<W> a, ap_uint<2> mode, SEXP_T shared_exp){
	float two_exp;
	float frac;
	float sexp = std::pow(2,int(shared_exp));
	float res;

	if(mode == 0b0){//right mode 0
		ap_uint<M0R> man = a(M0R-1,0);
#if S0R !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if E0R !=0
		ap_uint<1> denorm = a(W-2,W-E0R-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<E0R> exp = a(W-2,W-E0R-1) - (1<<(E0R-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-M0R) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-M0R+1);
#endif
		res = sgn*frac*two_exp*sexp;
		printf("mode1 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
	}else if(mode == 0b01){//right mode 1
		ap_uint<M1R> man = a(M1R-1,0);
#if S1R !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if E1R !=0
		ap_uint<1> denorm = a(W-2,W-E1R-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<E1R> exp = a(W-2,W-E1R-1) - (1<<(E1R-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-M1R) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-M1R+1);
#endif
		res = sgn*frac*two_exp*sexp;
		printf("mode2 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
	}else if(mode == 0b10){//left mode 0
		ap_uint<M0L> man = a(M0L-1,0);
#if S0L !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if E0L !=0
		ap_uint<1> denorm = a(W-2,W-E0L-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<E0L> exp = a(W-2,W-E0L-1) - (1<<(E0L-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-M0L) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-M0L+1);
#endif
		res = sgn*frac*two_exp*sexp;
		printf("mode3 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
	}else{//left mode 1
		ap_uint<M1L> man = a(M1L-1,0);
#if S1L !=0
		ap_int<2> sgn = a[W-1]?-1:1;
#else
		ap_int<2> sgn = 1;
#endif
#if E1L !=0
		ap_uint<1> denorm = a(W-2,W-E1L-1)==0? ap_uint<1>(1):ap_uint<1>(0);
		ap_int<E1L> exp = a(W-2,W-E1L-1) - (1<<(E1L-1))+1 + denorm;
		two_exp = std::pow(2,exp.to_int());
		frac = man.to_float()*std::pow(2,-M1L) + (~denorm);
#else
		two_exp = 1;
		frac = man.to_float()*std::pow(2,-M1L+1);
#endif
		res = sgn*frac*two_exp*sexp;
		printf("mode4 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
	}
	return res;
}

void PrintFPInput(ap_uint<W> A[SinM][SinN], SEXP_T betaA,ap_uint<2> mode){
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			float a = InputBMToFP(A[i][j],mode,betaA);
			printf("%f, \n",a);
		}printf("\n");
	}
}

void PrintFPOutput(ap_uint<W> A[SinM][SinN], SEXP_T betaA,ap_uint<2> modeOut){
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			float a = BMToFP(A[i][j],modeOut,betaA);
			printf("%f, ",a);
		}printf("\n");
	}
}
