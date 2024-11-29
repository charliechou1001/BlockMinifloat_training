#include "MAC.h"


void beta_SRL(int k, bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T &betaOut2){

	if((k & (BLK_SIZE-1))==0){
		ap_int<SEXP+4> betaplus = ap_int<SEXP+4>(betaA) + ap_int<SEXP+4>(betaB);
		SEXP_T betaAB;
		if(betaplus > SEXP_T( (1<<(SEXP-1))-1 ) ){
			betaAB = SEXP_T( (1<<(SEXP-1))-1 );
		}else if(betaplus < SEXP_T( -(1<<(SEXP-1))+1 ) ){
			betaAB = SEXP_T( -(1<<(SEXP-1))+1 );
		}else{
			betaAB = betaplus;
		}

		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );

	}

	if(k<0){
		betaC = 0;
	}

	if(flag==0b01){//load data
		betaOut1 = betaC;
	}else if(flag ==0b10){//shift data
		betaOut2 = betaOut1;
	}
}


void beta_SRL_2(int k,bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T betaOut2[MB], unsigned idx){

	if((k & (BLK_SIZE-1))==0){
		ap_int<SEXP+4> betaplus = ap_int<SEXP+4>(betaA) + ap_int<SEXP+4>(betaB);
		SEXP_T betaAB;
		if(betaplus > SEXP_T( (1<<(SEXP-1))-1 ) ){
			betaAB = SEXP_T( (1<<(SEXP-1))-1 );
		}else if(betaplus < SEXP_T( -(1<<(SEXP-1))+1 ) ){
			betaAB = SEXP_T( -(1<<(SEXP-1))+1 );
		}else{
			betaAB = betaplus;
		}

		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
	}

	if(k<0){
		betaC = 0;
	}

	if(flag==0b01){//load data
		betaOut1 = betaC;
	}else if(flag ==0b10){//shift data
		betaOut2[idx] = betaOut1;
	}
}

void PE2_SRL(int k, bool k_end, ap_uint<2> mode, int_t a, int_t b, ap_uint<EBIT> ebias,ap_uint<1> flg,
		ap_int<Kadd> &pSum, ap_int<Kadd> &accum, ap_int<Kadd> &out,
		ap_uint<2> flag, ap_int<Kadd> &sReg1, ap_int<Kadd> &sReg2   ){
#pragma HLS INLINE
	// Get previous sum
	ap_int<Kadd> last = (k & (BLK_SIZE-1))==0 ? ap_int<Kadd>(0) : pSum;

	// Update current sum
	// Handle boundary conditions
	int_t a_val = a;
	int_t b_val = b;
	pSum = BMMAC(a_val, b_val, last, mode);

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}


	if(k_end)
		out = accum;

	if(flag==0b01){//load data
		sReg1 = out;
	}else if(flag ==0b10){//shift data
		sReg2 = sReg1;
	}
}

template< int i, int j>
void PE3_SRL(int k, bool k_end, ap_uint<2> mode, int_t a, int_t b, ap_uint<EBIT> ebias,ap_uint<1> flg,
		ap_uint<2> flag, ap_int<Kadd> &sReg1, ap_int<Kadd> &sReg2   ){

	static ap_int<Kadd> pSum;
	static ap_int<Kadd> accum;
	static ap_int<Kadd> out;
	// Get previous sum
	ap_int<Kadd> last = (k & (BLK_SIZE-1))==0 ? ap_int<Kadd>(0) : pSum;

	// Update current sum
	// Handle boundary conditions
	int_t a_val = a;
	int_t b_val = b;
	pSum = BMMAC(a_val, b_val, last, mode);

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}


	if(k_end)
		out = accum;

	if(flag==0b01){//load data
		sReg1 = out;
	}else if(flag ==0b10){//shift data
		sReg2 = sReg1;
	}
}

template< int i, int j>
ap_int<Kadd> PE3(int k, bool k_end, ap_uint<2> mode, int_t a, int_t b, ap_uint<EBIT> ebias,ap_uint<1> flg){
#pragma HLS INLINE off
	static ap_int<Kadd> pSum;
	static ap_int<Kadd> accum;
	static ap_int <Kadd> out;
	// Get previous sum
	ap_int<Kadd> last = (k & (BLK_SIZE-1))==0 ? ap_int<Kadd>(0) : pSum;

	// Update current sum
	// Handle boundary conditions
	int_t a_val = a;
	int_t b_val = b;
	pSum = BMMAC(a_val, b_val, last, mode);

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}


	if(k_end)
		out = accum;

	return out;
}

