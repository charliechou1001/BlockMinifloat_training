#include "MAC.h"

void beta_SRL(int k, bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T &betaOut2){

	if((k & (BLK_SIZE-1))==0){
		SEXP_T betaAB = betaA + betaB;
//		printf("betaC before: %d, betaA: %d, betaB: %d,  ",betaC.to_int(),betaA.to_int(), betaB.to_int());
		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
//		printf("betaAB: %d, betaC: %d, ebias: %d, k_blk_start: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int(), k_blk_start);
	}

	if(flag==0b01){//load data
		betaOut1 = betaC;
	}else if(flag ==0b10){//shift data
		betaOut2 = betaOut1;
	}
}

void beta(int k, bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
		ap_uint<EBIT> &ebias, ap_uint<1> &flg){

	if((k & (BLK_SIZE-1))==0){
		SEXP_T betaAB = betaA + betaB;
//		printf("betaC before: %d, betaA: %d, betaB: %d,  ",betaC.to_int(),betaA.to_int(), betaB.to_int());
		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
//		printf("betaAB: %d, betaC: %d, ebias: %d, k_blk_start: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int(), k_blk_start);
	}

}

void beta_SRL_2(int k,bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T betaOut2[MB], unsigned idx){

	if((k & (BLK_SIZE-1))==0){
		SEXP_T betaAB = betaA + betaB;
//		printf("betaC before: %d, betaA: %d, betaB: %d,  ",betaC.to_int(),betaA.to_int(), betaB.to_int());
		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
//		printf("betaAB: %d, betaC: %d, ebias: %d, k_blk_start: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int(), k_blk_start);
	}

	if(flag==0b01){//load data
		betaOut1 = betaC;
	}else if(flag ==0b10){//shift data
		betaOut2[idx] = betaOut1;
//		printf("betaOut2[%d]: %d\n",idx, betaOut2[idx].to_int());
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
//	printf("a: %s, b:%s, last: %s, pSum: %s\n",a.to_string(2).c_str(), b.to_string(2).c_str(),
//			last.to_string(2).c_str(), pSum.to_string(2).c_str());

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}


//	if(int(k/BLK_SIZE)==0 && (((k+1) & (BLK_SIZE-1))==0 )){
//		accum = pSum;
//	}else if( int(k/BLK_SIZE)!=0 && (((k+1) & (BLK_SIZE-1))==0 ) ){
//		BMInterAccum<EBIT>(pSum, accum, ebias, flg);
////		printf("--accum: %s\n",accum.to_string(2).c_str());
//	}

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
//#pragma HLS INLINE
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
//	printf("a: %s, b:%s, last: %s, pSum: %s\n",a.to_string(2).c_str(), b.to_string(2).c_str(),
//			last.to_string(2).c_str(), pSum.to_string(2).c_str());

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}


//	if(int(k/BLK_SIZE)==0 && (((k+1) & (BLK_SIZE-1))==0 )){
//		accum = pSum;
//	}else if( int(k/BLK_SIZE)!=0 && (((k+1) & (BLK_SIZE-1))==0 ) ){
//		BMInterAccum<EBIT>(pSum, accum, ebias, flg);
////		printf("--accum: %s\n",accum.to_string(2).c_str());
//	}

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
//#pragma HLS INLINE
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
//	printf("a: %s, b:%s, last: %s, pSum: %s\n",a.to_string(2).c_str(), b.to_string(2).c_str(),
//			last.to_string(2).c_str(), pSum.to_string(2).c_str());

	if(((k+1) & (BLK_SIZE-1))==0 ){
		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
	}

//	if(i==0 && j==1)
//		printf("pSum: %s, accum:%s\t", pSum.to_string(2).c_str(), accum.to_string(2).c_str());

//	if(int(k/BLK_SIZE)==0 && (((k+1) & (BLK_SIZE-1))==0 )){
//		accum = pSum;
//	}else if( int(k/BLK_SIZE)!=0 && (((k+1) & (BLK_SIZE-1))==0 ) ){
//		BMInterAccum<EBIT>(pSum, accum, ebias, flg);
////		printf("--accum: %s\n",accum.to_string(2).c_str());
//	}

	if(k_end)
		out = accum;
//	else
//		out = 0;
//	printf("out:%d\n",out.to_int());
	return out;
}

void PEarray(int k, int size_k, ap_uint<2> flag, ap_uint<2> beta_flag, int ir, int ic,
			  int_t feederA[MAX_SIZE], int_t feederB[MAX_SIZE], ap_uint<2> mode, ap_int<Kadd> sReg_1[MAX_SIZE],
			  SExpPack *betaA, SExpPack *betaB, SEXP_T betasReg_1[MB][MB]){

#pragma HLS ARRAY_PARTITION variable=feederA dim=1 type=complete
#pragma HLS ARRAY_PARTITION variable=feederB dim=1 type=complete
#pragma HLS ARRAY_PARTITION variable=sReg_1 dim=1 type=complete
#pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 type=complete

	static ap_uint<EBIT> ebias[MB][MB];
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete
//	SEXP_T mul_beta[MB][MB];
//#pragma HLS ARRAY_PARTITION variable=mul_beta dim=0 complete
	static ap_uint<1> flg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete
	static SEXP_T betaRes[MB][MB], betasReg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete
    //------------------------------
	static ap_int<Kadd> pSum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = pSum dim = 0 complete
	static ap_int<Kadd> accum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = accum dim = 0 complete
	static ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete
	static ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete


    //------------------------------
    int TILE_K = size_k/MAX_SIZE;
    int kk =  k & (size_k-1);//k%(size_k-1)
    bool k_blk_start =  ( int(k/BLK_SIZE) & (size_k/BLK_SIZE-1) )==0;
    bool k_end = ( (k+1) & (size_k-1) ) ==0;

	//---------------shared exp----------------
	betax:for(int ip=0;ip<MB;ip++){
		#pragma HLS UNROLL
		betay:for(int iq=MB-1;iq>=0;iq--){
			#pragma HLS UNROLL
//					printf("ip: %d, iq: %d, b idx: %d\n",ip, iq,ic*TILE_K*MB);
			if(iq==MB-1)
				beta_SRL_2(k, k_blk_start, SEXP_T(betaA[ir*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip)),
						SEXP_T(betaB[ic*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*iq+SEXP-1,SEXP*iq)),
						betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg_1[ip], beta_flag[1]*( ( (k+1) & (size_k-1) ) -1) );
			else
				beta_SRL(k, k_blk_start, SEXP_T(betaA[ir*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip)),
						SEXP_T(betaB[ic*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*iq+SEXP-1,SEXP*iq)),
						betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg[ip][iq+1] );
		}
	}

	systolic2:
	for (int i = 0; i < MAX_SIZE; i++) {
		#pragma HLS UNROLL
	systolic3:
		for (int j = MAX_SIZE-1; j >=0; j--) {
			#pragma HLS UNROLL
			if(j == MAX_SIZE-1)
				PE2_SRL(kk, k_end, mode, feederA[i], feederB[j], ebias[int(i/BLK_SIZE)][int(j/BLK_SIZE)],
						 flg[int(i/BLK_SIZE)][int(j/BLK_SIZE)], pSum[i][j], accum[i][j], out[i][j], flag, sReg[i][j],sReg_1[i]);
			else
				PE2_SRL(kk, k_end, mode, feederA[i], feederB[j], ebias[int(i/BLK_SIZE)][int(j/BLK_SIZE)],
						 flg[int(i/BLK_SIZE)][int(j/BLK_SIZE)], pSum[i][j], accum[i][j], out[i][j], flag, sReg[i][j],sReg[i][j+1]);
		}
	}

}
