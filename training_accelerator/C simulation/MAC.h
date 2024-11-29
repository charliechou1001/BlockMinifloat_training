#ifndef MACDEF_
#define MACDEF_

#include "ap_int.h"
#include "typedef.h"
#include "newCLZ.h"

//#define MAC_DEBUG

ap_int<Kadd> BMMAC(ap_uint<W> a, ap_uint<W> b, ap_int<Kadd> last, ap_uint<2> mode){
#pragma HLS INLINE
	const unsigned EA = (1<<EL)-1;
	const unsigned EB = (1<<ER)-1;

	ap_uint<1> left = mode[1];
	ap_uint<1> right = mode[0];

#if (S0L ==0 && S1L != 0)
	ap_uint<1> sgnA = left==0? ap_uint<1>(0) : a[W-1];
#elif(S1L == 0 && S0L >0)
	ap_uint<1> sgnA = left==0? a[W-1] : ap_uint<1>(0);
#elif (S1L != 0 && S0L != 0)
	ap_uint<1> sgnA = a[W-1];
#else
	ap_uint<1> sgnA = ap_uint<1>(0);
#endif

#if (S0R ==0 && S1R !=0)
	ap_uint<1> sgnB = right==0? ap_uint<1>(0) : b[W-1];
#elif(S1R == 0 && S0R !=0)
	ap_uint<1> sgnB = right==0? b[W-1] : ap_uint<1>(0);
#elif (S1R != 0 && S0R !=0)
	ap_uint<1> sgnB = b[W-1];
#else
	ap_uint<1> sgnB = ap_uint<1>(0);
#endif


#if (E0L == 0 && E1L != 0)
	ap_uint<EL> expA = left == 0? ap_uint<EL>(1) : a(W-1-S1L,W-1-S1L-E1L+1);
	ap_uint<1> denormA = (left == 1 && expA == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#elif (E0L != 0 && E1L == 0)
	ap_uint<EL> expA = left == 0? a(W-1-S0L,W-1-S0L-E0L+1) : ap_uint<EL>(1);
	ap_uint<1> denormA = (left == 0 && expA == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#elif (E0L !=0 && E1L != 0)
	ap_uint<EL> expA = left == 0? a(W-1-S0L,W-1-S0L-E0L+1) : a(W-1-S1L,W-1-S1L-E1L+1);
	ap_uint<1> denormA = (expA == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#endif

#if (E0R == 0 && E1R != 0)
	ap_uint<ER> expB = right == 0? ap_uint<ER>(1) : b(W-1-S1R,W-1-S1R-E1R+1);
	ap_uint<1> denormB = (right == 1 && expB == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#elif (E0R != 0 && E1R == 0)
	ap_uint<ER> expB = right == 0? b(W-1-S0R,W-1-S0R-E0R+1) : ap_uint<ER>(1);
	ap_uint<1> denormB = (right == 0 && expB == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#elif (E0R !=0 && E1R != 0)
	ap_uint<ER> expB = right == 0? b(W-1-S0R,W-1-S0R-E0R+1) : b(W-1-S1R,W-1-S1R-E1R+1);
	ap_uint<1> denormB = (expB == 0)? ap_uint<1>(1) : ap_uint<1>(0);
#endif

#if (E0L != 0 || E1L != 0)
	ap_uint<1> isBM_L0 = E0L != 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_uint<1> isBM_L1 = E1L != 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_uint<ML+1> fracA = left == 0? ( ~denormA & isBM_L0, a(M0L-1,0)) : ( ~denormA & isBM_L1, a(M1L-1,0));
	ap_int<1+EA+ML> A = ((1-2*sgnA)*ap_int<1+EA+ML>(fracA)) << (expA -1 + denormA);//(expA -1 + denormA) is the number for shift, -1 for e=0 & e=1 case
#else
	ap_int<ML+1> A = (left == 0)?  (1-2*sgnA)*ap_int<ML+1>(a(M0L-1,0)): (1-2*sgnA)*ap_int<ML+1>(a(M1L-1,0));
#endif

#if (E0R != 0 || E1R != 0)
	ap_uint<1> isBM_R0 = E0R != 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_uint<1> isBM_R1 = E1R != 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_uint<MR+1> fracB = right == 0? ( ~denormB & isBM_R0, b(M0R-1,0)) : ( ~denormB & isBM_R1, b(M1R-1,0));
	#if WT>0
	ap_int<1+EB+MR+WT> B = ( ((1-2*sgnB)*ap_int<1+EB+MR>(fracB)) << (expB -1 + denormB), ap_uint<WT>(0) );
	#else
	ap_int<1+EB+MR> B = ((1-2*sgnB)*ap_int<1+EB+MR>(fracB)) << (expB -1 + denormB);
	#endif
#else
	ap_int<MR+1> B = right == 0?  (1-2*sgnB)*ap_int<MR+1>(b(M0R-1,0)): (1-2*sgnB)*ap_int<MR+1>(b(M1R-1,0));
#endif

	ap_int<Kadd> Psum = (ap_int<Kadd-WT>(A*B), ap_uint<WT>(0)) + last;
#pragma HLS BIND_OP variable=Psum impl=dsp
	return Psum;
}

// template<int Ebit>
// void BMADD(ap_int<Kadd-WI-WT> mult, ap_int<Kadd> &Acc, ap_uint<Ebit> ebias, ap_uint<1> flag, ap_uint<1> flagtile){
// 	ap_int<Kadd> C = Acc;
// 	ap_int<Kadd-WI> mul = (mult,ap_uint<WT>(0));
// 	ap_uint<Ebit> Ka = flag? ebias : ap_uint<Ebit>(0);
// 	ap_uint<Ebit> Kb = (!flag & flagtile)? ebias : ap_uint<Ebit>(0);
// 	Acc = (mul >> Ka) + (C >> Kb);
// }

template<int Ebit>
void GetEbeta2(ap_uint<Ebit> &ebias, ap_uint<1> &flag, SEXP_T betaAB, SEXP_T &betaC, bool ksgn){
	ap_int<SEXP+1> betabias = ap_int<SEXP+1>(betaAB) - ap_int<SEXP+1>(betaC);
	ap_int<Ebit+1> diff;
	if(betabias > ap_int<Ebit+1>( (1<<Ebit)-1 )  ){
		diff = ap_int<Ebit+1>( (1<<Ebit)-1 );
	}else if(betabias < ap_int<Ebit+1>(  -(1<<Ebit)+1 ) ){
		diff = ap_int<Ebit+1>(  -(1<<Ebit)+1 );
	}else{
		diff = betabias;
	}

	ap_int<Ebit+1> bias = ksgn? ap_int<Ebit+1>(0) : diff;
	flag = bias[Ebit];
	ebias = bias>=0? bias(Ebit-1,0) : (~ap_uint<Ebit>(bias(Ebit-1,0))+1);
	betaC = (bias>=0 || ksgn)? betaAB : betaC;
}

template<int Ebit>
void BMInterAccum(bool rst, ap_int<Kadd> psum, ap_int<Kadd> &Acc, ap_uint<Ebit> ebias, ap_uint<1> flag){
	ap_int<Kadd> A = psum;
	ap_int<Kadd> B = Acc;
	ap_uint<Ebit> Ka = flag? ebias : ap_uint<Ebit>(0);
	ap_uint<Ebit> Kb = (!flag)? ebias : ap_uint<Ebit>(0);
	Acc = rst? A : ap_int<Kadd>( (A >> Ka) + (B >> Kb) );
}


template<int Zm=W, int Zsh=W>
SEXP_T BiasAdjustv3(SEXP_T betaIn, ap_uint<Zm> Z_min, ap_int<Zsh> &Zshift, ap_uint<2>modeIn, ap_uint<2> modeOut){

	//input parameter
		ap_uint<1> left = modeIn[1];
		ap_uint<1> right = modeIn[0];
	#if (E0L == 0 && E1L != 0)
		ap_uint<EL> biasA = left == 0? ap_uint<EL>(0) : ap_uint<EL>((1<<(EL-1))-1);
	#elif (E0L != 0 && E1L == 0)
		ap_uint<EL> biasA = left == 0? ap_uint<EL>((1<<(EL-1))-1) : ap_uint<EL>(0);
	#elif (E0L !=0 && E1L != 0)
		ap_uint<EL> biasA = left == 0? (1<<(E0L-1))-1 : (1<<(E1L-1))-1;
	#else
		ap_uint<EL> biasA = ap_uint<EL>(0);
	#endif

	#if (E0R == 0 && E1R != 0)
		ap_uint<ER> biasB = right == 0? ap_uint<ER>(0) : ap_uint<ER>((1<<(ER-1))-1);
	#elif (E0R != 0 && E1R == 0)
		ap_uint<ER> biasB = right == 0? ap_uint<ER>((1<<(ER-1))-1) : ap_uint<ER>(0);
	#elif (E0R !=0 && E1R != 0)
		ap_uint<ER> biasB = right == 0? (1<<(E0R-1))-1 : (1<<(E1R-1))-1;
	#else
		ap_uint<ER> biasB = ap_uint<ER>(0);
	#endif
		ap_uint<ML> manA = left == 0? M0L : M1L;
		ap_uint<MR> manB = right == 0? M0R : M1R;

		//output parameter
		ap_uint<OUTE> ResBias;
		ap_uint<OUTE> ResE;
		if(modeOut == 0b00){
		#if ResE0 !=0
			ResBias = ap_uint<OUTE>((1<<(ResE0-1))-1);
			ResE = ResE0;
		#else
			ResBias = 0;
			ResE = 0;
		#endif
		}else if(modeOut == 0b01){
		#if ResE1 !=0
			ResBias = ap_uint<OUTE>((1<<(ResE1-1))-1);
			ResE = ResE1;
		#else
			ResBias = 0;
			ResE = 0;
		#endif
		}else if(modeOut == 0b10){
		#if ResE2 != 0
			ResBias = ap_uint<OUTE>((1<<(ResE2-1))-1);
			ResE = ResE2;
		#else
			ResBias = 0;
			ResE = 0;
		#endif
		}else{
			ResBias = 0;
			ResE = 0;
		}
	ap_uint<8> CntOverflow = (Kadd-1) - (biasA-1+manA + biasB-1+manB) - ( (1<<ResE)-1 -ResBias + 1) -WT;
//---------------------------------------------------------------------
	ap_int<Zm> emax_shift = CntOverflow - Z_min;//1 for sign bit
//	printf("===CntOverflow:%d\t Z_min:%d\t emax_shift:%d\n",CntOverflow.to_int(),Z_min.to_int(), emax_shift.to_int());
	Zshift = emax_shift;
//	Zshift = emax_shift>0 ? ap_uint<Zsh>(emax_shift) : ap_uint<Zsh>(0);

	ap_int<SEXP+4> newbeta = ap_int<SEXP+4>(betaIn) + ap_int<SEXP+4>(Zshift);
	SEXP_T beta_out;
	if(newbeta > SEXP_T( (1<<(SEXP-1))-1 ) ){
		beta_out = SEXP_T( (1<<(SEXP-1))-1 );
	}else if(newbeta  < SEXP_T( -(1<<(SEXP-1))+1 ) ){
		beta_out = SEXP_T( -(1<<(SEXP-1))+1 );
	}else
		beta_out = newbeta;

	return beta_out;
}


template<int Zm=W>
ap_uint<Zm> CheckZmin(ap_uint<Kadd-1> sgn_rst){
#pragma HLS INLINE
	ap_uint<Zm> Z_min = CLZ64((sgn_rst,ap_uint<64-Kadd+1>(1)));
//	ap_uint<Zm> Z_min = CLZ32((sgn_rst,ap_uint<32-Kadd+1>(1)));
	return Z_min;
}

ap_uint<Kadd-1> GetUnsigned(ap_int<Kadd> mac_out){
	//sign bit and two's complement transform
	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);
	return sgn_rst;
}

//modeOut: 00: FW, 01:ER, 10:GD
template<int Zsh=W, int Zm=W, int mTR=4>
void Normalization(ap_int<Kadd> mac_out, ap_int<Zsh> Zshift, ap_uint<2> modeIn, ap_uint<2> modeOut, ap_uint<W> &R,
		ap_uint<1> block_out, ap_uint<Ws> &RB, ap_uint<1> relu_en, ap_uint<1> &relu){

	//input parameter
	ap_uint<1> left = modeIn[1];
	ap_uint<1> right = modeIn[0];
#if (E0L == 0 && E1L != 0)
	ap_uint<EL> biasA = left == 0? ap_uint<EL>(0) : ap_uint<EL>((1<<(EL-1))-1);
#elif (E0L != 0 && E1L == 0)
	ap_uint<EL> biasA = left == 0? ap_uint<EL>((1<<(EL-1))-1) : ap_uint<EL>(0);
#elif (E0L !=0 && E1L != 0)
	ap_uint<EL> biasA = left == 0? (1<<(E0L-1))-1 : (1<<(E1L-1))-1;
#else
	ap_uint<EL> biasA = ap_uint<EL>(0);
#endif

#if (E0R == 0 && E1R != 0)
	ap_uint<ER> biasB = right == 0? ap_uint<ER>(0) : ap_uint<ER>((1<<(ER-1))-1);
#elif (E0R != 0 && E1R == 0)
	ap_uint<ER> biasB = right == 0? ap_uint<ER>((1<<(ER-1))-1) : ap_uint<ER>(0);
#elif (E0R !=0 && E1R != 0)
	ap_uint<ER> biasB = right == 0? (1<<(E0R-1))-1 : (1<<(E1R-1))-1;
#else
	ap_uint<ER> biasB = ap_uint<ER>(0);
#endif
	ap_uint<ML> manA = left == 0? M0L : M1L;
	ap_uint<MR> manB = right == 0? M0R : M1R;

	//output parameter
	ap_uint<OUTE> ResBias;
	ap_uint<OUTE> ResE;
	if(modeOut == 0b00 ){
	#if ResE0 !=0
		ResBias = ap_uint<OUTE>((1<<(ResE0-1))-1);
		ResE = ResE0;
	#else
		ResBias = 0;
		ResE = 0;
	#endif
	}else if(modeOut == 0b01){
	#if ResE1 !=0
		ResBias = ap_uint<OUTE>((1<<(ResE1-1))-1);
		ResE = ResE1;
	#else
		ResBias = 0;
		ResE = 0;
	#endif
	}else if(modeOut == 0b10){
	#if ResE2 != 0
		ResBias = ap_uint<OUTE>((1<<(ResE2-1))-1);
		ResE = ResE2;
	#else
		ResBias = 0;
		ResE = 0;
	#endif
	}else{
		ResBias = 0;
		ResE = 0;
	}

	ap_uint<8> expZeroPoint = (Kadd -1) - (biasA -1 + manA + biasB -1 + manB) - WT;
	ap_uint<8> CntDenorm = (Kadd-1) - (biasA-1+manA + biasB-1+manB) +(ResBias-1) - WT;

	//sign bit and two's complement transform
	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);
	if(signR && relu_en){
		sgn_rst = 0;
		relu = 0;
	}else{
		sgn_rst = sgn_rst >> Zshift;
		relu = 1;
	}

#ifdef MAC_DEBUG
	printf("-------------------\n");
	printf("mac_out: %s,\n sgn_rst: %s,  Zshift:%d\n",mac_out.to_string(2).c_str(),sgn_rst.to_string(2).c_str(),Zshift.to_int());
#endif

	if(block_out == 0b0){

	//count leading zero and exponent
	ap_uint<Zm> count = CLZ64((sgn_rst,ap_uint<64-Kadd+1>(1)));

	ap_uint<OUTE> exp;
	if(count >= CntDenorm)//denorm and underflow
		exp = 0;
	else
		exp = expZeroPoint - (count + 1) + ResBias;//1 for implicit bit
#ifdef MAC_DEBUG
	printf("Kadd: %d, expZeroPoint: %d, count: %d, ResBias: %d\t",Kadd,expZeroPoint.to_int(),count.to_int(), ResBias.to_int());
	printf("exp: %s\n",exp.to_string(2).c_str());
#endif

	//mantissa and rounding
	ap_uint<Zm> count_adjust;
	ap_uint<OUTM> man;
	ap_uint<1> rand;
	if(count < CntDenorm && ResE != 0)
		count_adjust = count;
	else if(count >= CntDenorm && ResE != 0)//denorm and underflow
		count_adjust = ap_uint<8>(CntDenorm-1);
	else//BFP out
		count_adjust = expZeroPoint -2;//expZeroPoint -1, and another -1 to align to sgn_rst(Kadd-3,)
	sgn_rst = sgn_rst << count_adjust;
#ifdef MAC_DEBUG
	printf("count_adjust:%d\t", count_adjust.to_int());
	printf("sgn_rst after: %s\n",sgn_rst.to_string(2).c_str());
#endif

	if(modeOut == 0b00){
		man = sgn_rst(Kadd-3,Kadd-2-ResM0);
		ap_uint<1> guard = sgn_rst[Kadd-2-ResM0];
		ap_uint<1> round = sgn_rst[Kadd-3-ResM0];
		ap_uint<1> sticky = sgn_rst(Kadd-4-ResM0,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
		rand = round & sticky | guard&round&(~sticky);

		if(man < (1<<ResM0)-1)
			man = man + rand;
#if ResE0 != 0
		else if(rand==1 && exp < (1<<ResE0)-1){
			exp ++;
			man = 0;
		}
#endif
	}else if(modeOut == 0b01){
		man = sgn_rst(Kadd-3,Kadd-2-ResM1);
		ap_uint<1> guard = sgn_rst[Kadd-2-ResM1];
		ap_uint<1> round = sgn_rst[Kadd-3-ResM1];
		ap_uint<1> sticky = sgn_rst(Kadd-4-ResM1,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
		rand = round & sticky | guard&round&(~sticky);

		if(man < (1<<ResM1)-1)
			man = man + rand;
#if ResE1 != 0
		else if(rand==1 && exp < (1<<ResE1)-1){
			exp ++;
			man = 0;
		}
#endif
	}else{//(modeOut == 0b10)
		man = sgn_rst(Kadd-3,Kadd-2-ResM2);
		ap_uint<1> guard = sgn_rst[Kadd-2-ResM2];
		ap_uint<1> round = sgn_rst[Kadd-3-ResM2];
		ap_uint<1> sticky = sgn_rst(Kadd-4-ResM2,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
		rand = round & sticky | guard&round&(~sticky);

		if(man < (1<<ResM2)-1)
			man = man + rand;
#if ResE2 != 0
		else if(rand==1 && exp < (1<<ResE2)-1){
			exp ++;
			man = 0;
		}
#endif
	}
#ifdef MAC_DEBUG
	printf("sgn_rst(Kadd-3,Kadd-2-ResM2): %s\t",sgn_rst(Kadd-3,Kadd-2-ResM2).to_string(2).c_str() );
	printf("man:%s\n",man.to_string(2).c_str());
#endif

	//output
	if(modeOut == 0b00){
#if (ResS0 != 0 && ResE0 !=0)//BM
		R[W-1] = signR;
		R(W-2,W-1-ResE0) = exp;
		R(ResM0-1,0) = man;
#elif (ResS0 == 0 && ResE0 !=0)//unsigned BM
		R(W-1,W-ResE0) = exp;
		R(ResM0-1,0) = man;
#elif (ResS0 != 0 && ResE0 ==0)//BFP
		R[W-1] = signR;
		R(ResM0-1,0) = man;
#else //unsigned BFP
		R(ResM0-1,0) = man;
#endif
	}else if(modeOut == 0b01){
#if (ResS1 != 0 && ResE1 !=0)//BM
		R[W-1] = signR;
		R(W-2,W-1-ResE1) = exp;
		R(ResM1-1,0) = man;
#elif (ResS1 == 0 && ResE1 !=0)//unsigned BM
		R(W-1,W-ResE1) = exp;
		R(ResM1-1,0) = man;
#elif (ResS1 != 0 && ResE1 ==0)//BFP
		R[W-1] = signR;
		R(ResM1-1,0) = man;
#else //unsigned BFP
		R(ResM1-1,0) = man;
#endif
	}else{//(modeOut == 0b10)
#if (ResS2 != 0 && ResE2 !=0)//BM
		R[W-1] = signR;
		R(W-2,W-1-ResE2) = exp;
		R(ResM0-1,0) = man;
#elif (ResS2 == 0 && ResE2 !=0)//unsigned BM
		R(W-1,W-ResE2) = exp;
		R(ResM0-1,0) = man;
#elif (ResS2 != 0 && ResE2 ==0)//BFP
		R[W-1] = signR;
		R(ResM2-1,0) = man;
#else //unsigned BFP
		R(ResM2-1,0) = man;
#endif
	}

#ifdef MAC_DEBUG
	printf("R: %s\n",R.to_string(2).c_str());
#endif

	}else{//block_out == 0b1
		sgn_rst = sgn_rst << (expZeroPoint -2);
		ap_uint<Ws-1> man_out = sgn_rst(Kadd-3,Kadd-2-ResM3);
		RB(Ws-2,0) = man_out;
		RB[Ws-1] = signR;

#ifdef MAC_DEBUG
		printf("---\n");
	printf("expZeroPoint: %d, signR:%s\n",expZeroPoint.to_int(), signR.to_string(2).c_str() );
	printf("sgn_rst after shift:%s, man_out: %s \n", sgn_rst.to_string(2).c_str(), man_out.to_string(2).c_str() );
	printf("RB: %s\n",RB.to_string(2).c_str() );
#endif

	}

}


#endif
