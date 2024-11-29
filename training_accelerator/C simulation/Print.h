#include "ap_int.h"


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
//		printf("mode1 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
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
//		printf("mode2 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
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
//		printf("mode3 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
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
//		printf("mode4 -- res: %f, frac:%f, two_exp:%f, sexp:%f, ",res,frac,two_exp,sexp);
	}
	return res;
}

//void PrintFPInput(ap_uint<W> A[SinM][SinN], SEXP_T betaA,ap_uint<2> mode){
//	for(int i=0;i<SinM;i++){
//		for(int j=0;j<SinN;j++){
//			float a = InputBMToFP(A[i][j],mode,betaA);
//			printf("%f, \n",a);
//		}printf("\n");
//	}
//}

void PrintMatrix(MemPack *A, SExpPack *betaA, int size, int betaA_addr, ap_uint<2> mode){
//	printf("--mode: %d\n",mode.to_int());
//	printf("size: %d\n",size);
	for(int i=0;i<size/(MAX_SIZE);i++){
		MemPack aa = A[i];
		for(int j=0;j<MAX_SIZE;j++){
//			printf("-----\n");
//			printf("betaA: %d, ",betaA[betaA_addr + int(i/BLK_SIZE)].range(SEXP*int(j/BLK_SIZE)+SEXP-1, SEXP*int(j/BLK_SIZE)).to_int());
			ap_uint<W> a = aa.range(j*W+W-1,j*W);
			SEXP_T beta = betaA[betaA_addr + int(i/BLK_SIZE)].range(SEXP*int(j/BLK_SIZE)+SEXP-1, SEXP*int(j/BLK_SIZE));
			float b = InputBMToFP(a,mode,beta);
			printf("%e, ", b );
//			printf("BM8: %s, float: %f\t",a.to_string(2).c_str(), b);
//			printf("betaA: %d, idx: %d, betaA_addr: %d\n",betaA[betaA_addr + int(i/MAX_SIZE)].to_int(), betaA_addr + int(i/MAX_SIZE), betaA_addr );
		}printf("\n");
	}
}


void PrintKadd(int size, XPack *in, SExpPack *beta_in){
	float mm,ll;
	for(int i=0;i<size;i++){
//		printf("beta_in: %d\n",beta_in[int(i/MAX_SIZE)].to_int());
		for(int j=0;j<MAX_SIZE;j++){
			ap_uint<Ws> aa = in[i].range(j*Ws+Ws-1,j*Ws);
			float sgn = aa[Ws-1]? -1 :1;
			ap_uint<Ws-1> frac = aa(Ws-2,0);
			SEXP_T beta = beta_in[int(i/BLK_SIZE)].range(SEXP*int(j/BLK_SIZE)+SEXP-1, SEXP*int(j/BLK_SIZE));

//			printf("beta: %f, fractional: %f, ", float(beta), frac.to_float()/std::pow(2,Ws-2) );
			mm = sgn*frac.to_float()/std::pow(2,Ws-2)*std::pow(2, float(beta));
			printf("%f, ",mm );
//			printf("%s, ",aa.to_string(2).c_str());
		}printf("\n");
	}
}

void PrintKaddValue(ap_int<Kadd> in, int expZeroPoint){
	float mm;
	mm = float(in)/std::pow(2,expZeroPoint);
	printf("%f, ",mm);
}
