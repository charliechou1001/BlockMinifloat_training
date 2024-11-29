#include "MAC.h"
#include "PEarray.h"

void feeder2( MemPack *a, int_t feeder[MAX_SIZE], int idx, bool PE_en){
#pragma HLS INLINE
	MemPack rev = a[idx*PE_en];
//	printf("idx*PE_en: %d, idx:%d\n",idx*PE_en, idx);
	for(int ii=0;ii<MAX_SIZE;ii++){
		#pragma HLS UNROLL
		feeder[ii] = int_t(rev.range(ii*W+W-1,ii*W));
	}

}

//void PE2_SRL(int k, bool k_end, ap_uint<2> mode, int_t a, int_t b, ap_uint<EBIT> ebias,ap_uint<1> flg,
//		ap_int<Kadd> &pSum, ap_int<Kadd> &accum, ap_int<Kadd> &out,
//		ap_uint<2> flag, ap_int<Kadd> &sReg1, ap_int<Kadd> &sReg2   ){
//#pragma HLS INLINE
//	// Get previous sum
//	ap_int<Kadd> last = (k & (BLK_SIZE-1))==0 ? ap_int<Kadd>(0) : pSum;
//
//	// Update current sum
//	// Handle boundary conditions
//	int_t a_val = a;
//	int_t b_val = b;
//	pSum = BMMAC(a_val, b_val, last, mode);
////	printf("a: %s, b:%s, last: %s, pSum: %s\n",a.to_string(2).c_str(), b.to_string(2).c_str(),
////			last.to_string(2).c_str(), pSum.to_string(2).c_str());
//
//	if(((k+1) & (BLK_SIZE-1))==0 ){
//		BMInterAccum<EBIT>((int(k/BLK_SIZE)==0 ),pSum, accum, ebias, flg);
//	}
//
//
////	if(int(k/BLK_SIZE)==0 && (((k+1) & (BLK_SIZE-1))==0 )){
////		accum = pSum;
////	}else if( int(k/BLK_SIZE)!=0 && (((k+1) & (BLK_SIZE-1))==0 ) ){
////		BMInterAccum<EBIT>(pSum, accum, ebias, flg);
//////		printf("--accum: %s\n",accum.to_string(2).c_str());
////	}
//
//	if(k_end)
//		out = accum;
//
//	if(flag==0b01){//load data
//		sReg1 = out;
//	}else if(flag ==0b10){//shift data
//		sReg2 = sReg1;
//	}
//}



//void beta_SRL(int k, bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
//		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
//		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T &betaOut2){
//
//	if((k & (BLK_SIZE-1))==0){
//		SEXP_T betaAB = betaA + betaB;
////		printf("betaC before: %d, betaA: %d, betaB: %d,  ",betaC.to_int(),betaA.to_int(), betaB.to_int());
//		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
////		printf("betaAB: %d, betaC: %d, ebias: %d, k_blk_start: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int(), k_blk_start);
//	}
//
//	if(flag==0b01){//load data
//		betaOut1 = betaC;
//	}else if(flag ==0b10){//shift data
//		betaOut2 = betaOut1;
//	}
//}
//
//void beta_SRL_2(int k,bool k_blk_start, SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC,
//		ap_uint<EBIT> &ebias, ap_uint<1> &flg,
//		ap_uint<2> flag, SEXP_T &betaOut1, SEXP_T betaOut2[MB], unsigned idx){
//
//	if((k & (BLK_SIZE-1))==0){
//		SEXP_T betaAB = betaA + betaB;
////		printf("betaC before: %d, betaA: %d, betaB: %d,  ",betaC.to_int(),betaA.to_int(), betaB.to_int());
//		GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, k_blk_start );
////		printf("betaAB: %d, betaC: %d, ebias: %d, k_blk_start: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int(), k_blk_start);
//	}
//
//	if(flag==0b01){//load data
//		betaOut1 = betaC;
//	}else if(flag ==0b10){//shift data
//		betaOut2[idx] = betaOut1;
////		printf("betaOut2[%d]: %d\n",idx, betaOut2[idx].to_int());
//	}
//}


void CheckMax(int tt, ap_int<Kadd> sRegOut, ap_uint<Kadd-1> &ResAND,
//		ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_uint<Kadd-1> max_shift1[2], ap_uint<Kadd-1> max_shift2[2]){
#pragma HLS INLINE
	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	ap_uint<1> pp_idx = int(tt/BLK_SIZE) & 0b1;

	if(tt >=0 && tt < MAX_SIZE ){//data shift out

		//check max number of each column in PE array
		if( ((tt+1)&(BLK_SIZE-1)) == 0  ){//tt == BLK_SIZE*m-1
			max_shift1[~pp_idx] = ResAND | GetUnsigned(sRegOut);
//			max_shift1 = ResAND | GetUnsigned(sRegOut);
//			printf(".___%d___.",pp_idx.to_int());
			ResAND = 0b0;//clear register
		}else
			ResAND = ResAND | GetUnsigned(sRegOut);

//		printf("ResAND: %s, unsigned sRegOut:%s\n",ResAND.to_string(2).c_str(), GetUnsigned(sRegOut).to_string(2).c_str());
	}

	if(tt >= BLK_SIZE && tt<MAX_SIZE+BLK_SIZE){//max number shift vertically
//		printf("\\ ");
		max_shift2[pp_idx] = max_shift1[pp_idx];
//		max_shift2 = max_shift1;
	}

}

//void WriteBuffer(int tt, int tile_idx, ap_int<Kadd> sRegOut, ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2]){
//#pragma HLS INLINE
//	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
//	if(tt >=0 && tt < MAX_SIZE ){//data shift out
//		//write the PE array output to accbf buffer
////		int idx = ( int(tt/BLK_SIZE+tile_idx) + int(tt/BLK_SIZE+tile_idx)/3 )&3;//tile_idx =ir*TILE_C*MB+ic*MB;
////		printf("accbf idx: %d, %d\n",idx, BLK_SIZE*idx + BLK_SIZE -1 - (tt_mod));
////		accbf[BLK_SIZE*idx + BLK_SIZE -1 - (tt_mod)] = sRegOut;
//
//		ap_uint<2> idx = int(tt/BLK_SIZE+tile_idx) & 3;
//		if(idx[1] == 0b0){
//			accbf1[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
////			printf("test_accbf[1] : %d, idx:%d, tt:%d, tile_idx:%d\n", BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod), idx.to_int() , tt, tile_idx);
//		}else{
//			accbf2[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
////			printf("test_accbf[2] : %d, idx:%d, tt:%d, tile_idx:%d\n", BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod), idx.to_int() , tt, tile_idx);
//		}
//	}
//
//}

int_t WriteBuffer2(int tt, int tile_idx, int wr_tt, int wr_tile_idx,	ap_uint<2> mode, ap_uint<2> modeout,
		ap_int<Kadd> sRegOut, ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_int<W> Zshift[MB]){
#pragma HLS INLINE

	int_t R;
	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	if(tt >=0 && tt < MAX_SIZE ){//data shift out
		//write the PE array output to accbf buffer
//		int idx = ( int(tt/BLK_SIZE+tile_idx) + int(tt/BLK_SIZE+tile_idx)/3 )&3;//tile_idx =ir*TILE_C*MB+ic*MB;
//		printf("accbf idx: %d, %d\n",idx, BLK_SIZE*idx + BLK_SIZE -1 - (tt_mod));
//		accbf[BLK_SIZE*idx + BLK_SIZE -1 - (tt_mod)] = sRegOut;

		ap_uint<2> idx = int(tt/BLK_SIZE+tile_idx) & 3;
//		printf("idx: %d, tt:%d, tile idx:%d\n",idx.to_int(), tt, tile_idx);
		if(idx[1] == 0b0){
			accbf1[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
//			printf("test_accbf[1] : %d, idx:%d, tt:%d, tile_idx:%d\n", BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod), idx.to_int() , tt, tile_idx);
		}else{
			accbf2[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
//			printf("test_accbf[2] : %d, idx:%d, tt:%d, tile_idx:%d\n", BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod), idx.to_int() , tt, tile_idx);
		}
	}

	if(wr_tt >= 0 && wr_tt < MAX_SIZE){
		ap_uint<2> wr_idx = int(wr_tt/BLK_SIZE + wr_tile_idx*MB) & 3;//////
//				ap_uint<2> idx = int((wr_tt + tile_idx*MAX_SIZE)/(BLK_SIZE*2));
		ap_uint<Kadd> acc;
		if(wr_idx[1] == 0b0){
			acc = accbf1[BLK_SIZE*wr_idx[0] + ( wr_tt & (BLK_SIZE-1) ) ];
//					printf(" norm accbf1, acc:%s \n",acc.to_string(2).c_str());
		}else{
			acc = accbf2[BLK_SIZE*wr_idx[0] + ( wr_tt & (BLK_SIZE-1) ) ];
//					printf(" norm accbf2, acc:%s \n",acc.to_string(2).c_str());
		}

//				int idx =( int(wr_tt/BLK_SIZE+wr_tile_idx*MB) + int(wr_tt/BLK_SIZE+wr_tile_idx*MB)/3 )&3;
		Normalization<W,W,4>(acc, Zshift[MB-1-int(wr_tt/BLK_SIZE)], mode, modeout, R);

	}
	return R;

}

void AccumNorm(int tt, int tile_idx, int wr_tt, int wr_tile_idx, ap_uint<2> mode, ap_uint<2> modeout,
		ap_int<Kadd> sReg_1[MAX_SIZE], ap_int<W> Zshift[MB][MB], MemPack *c){

#pragma HLS INLINE

#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = Zshift dim = 1 complete

	static ap_int<Kadd> accbf1[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf1 dim = 1 complete
	static ap_int<Kadd> accbf2[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf2 dim = 1 complete
    int_t R[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=R dim=1 complete

    MemPack RR;

	for(int i = 0; i < MAX_SIZE; i++){
		#pragma HLS UNROLL
			R[i] = WriteBuffer2(tt, tile_idx*MB, wr_tt, wr_tile_idx, mode, modeout,
			sReg_1[i], accbf1[i], accbf2[i], Zshift[int(i/BLK_SIZE)]);
			RR.range(W*i+W-1,W*i) = R[i];
	}

	if(wr_tt >= 0 && wr_tt < MAX_SIZE){
		c[wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) )] = RR;
	}

}


void CheckMax_Beta(int tt, int tile_idx, ap_int<Kadd> sRegOut, ap_uint<Kadd-1> &ResAND,
//		ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_uint<Kadd-1> max_shift1[2], ap_uint<Kadd-1> &max_shift2,
		SEXP_T betaOut[MB], SEXP_T betaC[MB], ap_int<W> Zshift[MB], ap_uint<2> mode,  ap_uint<2> modeout){
#pragma HLS INLINE

	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	ap_uint<1> pp_idx = int(tt/BLK_SIZE) & 0b1;

	if(tt >=0 && tt < MAX_SIZE ){//data shift out

		//check max number of each column in PE array
		if( ((tt+1)&(BLK_SIZE-1)) == 0  ){//tt == BLK_SIZE*m-1
			max_shift1[~pp_idx] = ResAND | GetUnsigned(sRegOut);
//			max_shift1 = ResAND | GetUnsigned(sRegOut);
//			printf(".___%d___.",pp_idx.to_int());
			ResAND = 0b0;//clear register
		}else
			ResAND = ResAND | GetUnsigned(sRegOut);

//		printf("ResAND: %s, unsigned sRegOut:%s\n",ResAND.to_string(2).c_str(), GetUnsigned(sRegOut).to_string(2).c_str());

	}

	if(tt >= BLK_SIZE && tt<MAX_SIZE+BLK_SIZE){//max number shift vertically
//		printf("max_shift1: %s, ", max_shift1[pp_idx].to_string(2).c_str());
		if( tt_mod == 0 )
			max_shift2 = max_shift1[pp_idx];
		else if( (tt_mod>=0) && (tt_mod<=BLK_SIZE-1) )
			max_shift2 = max_shift2 | max_shift1[pp_idx];

		if( tt_mod == BLK_SIZE-1){
			ap_uint<W> Zmin = CheckZmin<W>(max_shift2);
			int idx = int(tt/BLK_SIZE)-1;
	//		printf("idx:%d\n",idx);
//			printf("Zshift idx: %d, tt:%d,  %d\n",MB-1-int((tt-BLK_SIZE)/BLK_SIZE), tt, int((tt-BLK_SIZE)/BLK_SIZE));
			betaC[idx] = BiasAdjustv3<W,W>( betaOut[idx], Zmin, Zshift[MB-1-idx], mode, modeout );
//			printf("int(tt/BLK_SIZE):%d,ResSum: %s,  Zshift: %d, betaOut: %d, betaC:%d\n", idx, max_shift2.to_string(2).c_str(),
//					Zshift[MB-1-idx].to_int(), betaOut[idx].to_int(),betaC[idx].to_int() );

		}
	}

}

void bmGemmv23(MemPack *a, // Read-Only Matrix A
        MemPack *b, // Read-Only Matrix B
		   MemPack *c,
		   SExpPack *betaA,
		   SExpPack *betaB,
		   SExpPack *betaC,
		   int size_m,    // The other Dimension of Matrix A
		   int size_n,    // Shared Dimension of A*B, or A.T*B, A*B.T, A.T*B.T
		   int size_k,     // The other Dimension of Matrix B
		   ap_uint<8> log_sizek,
		   ap_uint<8> log_sizek_tilec,
		   ap_uint<2> mode,
		   ap_uint<2> modeout)
{//assume that size_m, size_k, size_n can be divided by MAX_SIZE
#pragma HLS INTERFACE m_axi port=a offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=c offset=slave bundle=gmem2

    int TILE_R = size_m/MAX_SIZE;
    int TILE_C = size_n/MAX_SIZE;
    int TILE_K = size_k/MAX_SIZE;
    int BLK_K = size_k/BLK_SIZE;

    int_t feederA[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 complete
    int_t feederB[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederB dim=1 complete
    //------------------------------
    ap_int<Kadd> pSum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = pSum dim = 0 complete
    ap_int<Kadd> accum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = accum dim = 0 complete
    ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete
    ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete
    ap_int<Kadd> sReg_1[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 0 complete

    //------------------------------
	ap_uint<EBIT> ebias[MB][MB];
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete
//	SEXP_T mul_beta[MB][MB];
//#pragma HLS ARRAY_PARTITION variable=mul_beta dim=0 complete
	ap_uint<1> flg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete
	ap_int<W> Zshift[MB][MB];
#pragma HLS ARRAY_PARTITION variable=Zshift dim=1 complete
	SEXP_T betaRes[MB][MB], betasReg[MB][MB], betasReg_1[MB][MB], betaOutReg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 complete
#pragma HLS ARRAY_PARTITION variable=betaOutReg dim=1 complete
//------------------------------
	ap_uint<Kadd-1> ResAND[MAX_SIZE], ResAND_1[MAX_SIZE][2];
	ap_uint<Kadd-1> ResAND_2[MB];
#pragma HLS ARRAY_PARTITION variable = ResAND dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_2 dim = 1 complete
	ap_int<Kadd> accbf1[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf1 dim = 1 complete
	ap_int<Kadd> accbf2[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf2 dim = 1 complete
	//------------------------------
    int_t R[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=R dim=1 complete



	tile_loop:for(int k=0;k<TILE_R*TILE_C*size_k + MAX_SIZE + 2*BLK_SIZE+1;k++){
//		for(int k=0;k<size_k + MAX_SIZE + 2*BLK_SIZE +1;k++){
#pragma HLS PIPELINE

		//------------control and idx---------------
//		bool k_0 = ( k & (size_k-1) )==0;
    	bool k_end = ( (k+1) & (size_k-1) ) ==0;
//    	bool k_blk = (k+1) & ( BLK_SIZE-1 );
    	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;
//    	printf("k_0: %d, k_end: %d, ",k_0,k_end);

    	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;
    	bool SRL_en =  ( k & (size_k-1) ) >=0 && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k;
//    	int tile_idx = ( k & (size_k-1) );
    	ap_uint<2> flag;
    	if( ( k & (size_k-1) )==0 && SRL_en )
    		flag = 0b01;//load data
    	else if( SRL_en && ( (k & (size_k-1) )>0) && ((k & (size_k-1) ) <= MAX_SIZE) )
    		flag = 0b10;//shift
    	else
    		flag = 0b00;
//    	printf("k: %d, PE_en: %d, SRL_en:%d, flag: %d\n",k,feeder_en, SRL_en, flag.to_int());

//    	int ir = int(k/size_k) / TILE_C;
//    	int ic = int(k/size_k) & (TILE_C-1);
    	int ir = ( k >> log_sizek_tilec );
    	int ic = ( k >> log_sizek ) & (TILE_C-1);

    	int kk = k & (size_k-1);//k%(size_k-1)
//    	int tt1 = k % (size_k + MAX_SIZE + 2*BLK_SIZE) - size_k -1;
    	int tt = k > size_k-1 ? ( k &(size_k-1) )-1 : -1;
    	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;
//    	printf("k: %d, kk: %d, tt:%d, wr_tt:%d\t",k,kk,tt, wr_tt);

		int tile_idx = ir*TILE_C+ic-1;
//		int wr_tile_idx = ( int( (k  -1)/size_k) / TILE_C ) *TILE_C + ( int( (k -1)/size_k) & (TILE_C-1)) -1;
		int wr_tile_idx = ( (k  -1) >> log_sizek_tilec ) *TILE_C + ( ( (k -1)>> log_sizek) & (TILE_C-1)) -1;
//		printf("tile_idx: %d, wr_tile_idx: %d, %d\n",tile_idx,wr_tile_idx, int((k - 2*BLK_SIZE)/size_k));
    	//---------------beta idx------------------
    	ap_uint<2> beta_flag;
//    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k < TILE_R*TILE_C*size_k;
    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k >= size_k -1;
    	if( ( (k+1) & (size_k-1) )==0 && beta_SRL_en )
    		beta_flag = 0b01;
    	else if( beta_SRL_en && ( ((k+1) & (size_k-1) )>0) && (((k+1) & (size_k-1) ) <= MB) )
    		beta_flag = 0b10;
    	else
    		beta_flag = 0b00;

//    	printf("----ir: %d, ic: %d, beta_flag: %d\n",ir,ic, beta_flag.to_int());
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
//			if(ip==0 && iq==0)
//				printf("betaA: %d, betaB:%d, betaC:%d\n",SEXP_T(betaA[ir*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip)).to_int(),
//						SEXP_T(betaB[ic*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*iq+SEXP-1,SEXP*iq)).to_int(),  betaRes[ip][iq].to_int());
			}
		}

    	//---------------PE array-------------------
    	feeder2(a, feederA, (k & (size_k-1))+ir*size_k, feeder_en);
    	feeder2(b, feederB, (k & (size_k-1))+ic*size_k, feeder_en);
//    	printf("feederA idx: %d, feederB idx: %d\n",(k & (size_k-1))+ir*size_k, (k & (size_k-1))+ic*size_k);

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

//		for(int i = 0; i < MAX_SIZE; i++){
//		#pragma HLS UNROLL
//				R[i] = WriteBuffer2(tt, tile_idx*MB, wr_tt, wr_tile_idx, mode, modeout,
//				sReg_1[i], accbf1[i], accbf2[i], Zshift[int(i/BLK_SIZE)]);
//		}


		for(int i = 0; i < MAX_SIZE; i++){
#pragma HLS UNROLL
			//----------------Check Max block----------
			if( (i&(BLK_SIZE-1)) !=0 ){
//				printf("ResAND_1[%d] before:%s, ",i-1,ResAND_1[i-1].to_string(2).c_str());
				CheckMax(tt, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_1[i-1]  );
//				printf("ResAND_1[%d]:%s, \n",i-1,ResAND_1[i-1].to_string(2).c_str());
			}else{
//				printf("ResAND_2[%d] before:%s, ",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
				CheckMax_Beta(tt, tile_idx*MB, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_2[int(i/BLK_SIZE)],
						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], mode, modeout );
//				printf("ResAND_2[%d]:%s, \n",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
			}
		}

		AccumNorm(tt, tile_idx, wr_tt, wr_tile_idx, mode, modeout, sReg_1, Zshift, c);

		//-----beta write out-----
		for(int ip=0; ip<MB;ip++){
#pragma HLS UNROLL
			if(tt>=BLK_SIZE && tt <MAX_SIZE+BLK_SIZE && ( tt & (BLK_SIZE-1) )==BLK_SIZE-1 ){
				betaC[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
//				printf("k:%d, tt: %d, beta: %d, betaC idx:%d\n",k, tt, betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)].to_int(), tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE));
			}
		}

//		for(int i = 0; i < MAX_SIZE; i++){
//			#pragma HLS UNROLL
//			//---------------norm --------------------
//			if(wr_tt >= 0 && wr_tt < MAX_SIZE){
//				c[wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) )].range(W*i+W-1,W*i) = R[i];
////				printf("norm idx: %d, idx:%d, tt:%d, wr_tile_idx:%d\n", BLK_SIZE*idx[0] + ( wr_tt & (BLK_SIZE-1) ),idx.to_int(),tt,wr_tile_idx);
////				printf("wr_tile_idx: %d, R[%d][%d]: %s, accbf idx: %d, tt:%d, k:%d, wr_tt: %d\n ",wr_tile_idx, wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) ),
////						i, R[i].to_string(2).c_str(),BLK_SIZE*idx + ( wr_tt & (BLK_SIZE-1) ) , tt , k, wr_tt);
//
//			}
//		}

	}//k end


}

void bmGemmv24(MemPack *a, // Read-Only Matrix A
        MemPack *b, // Read-Only Matrix B
		   MemPack *c,
		   SExpPack *betaA,
		   SExpPack *betaB,
		   SExpPack *betaC,
		   int size_m,    // The other Dimension of Matrix A
		   int size_n,    // Shared Dimension of A*B, or A.T*B, A*B.T, A.T*B.T
		   int size_k,     // The other Dimension of Matrix B
		   ap_uint<8> log_sizek,
		   ap_uint<8> log_sizek_tilec,
		   ap_uint<2> mode,
		   ap_uint<2> modeout)
{//assume that size_m, size_k, size_n can be divided by MAX_SIZE
#pragma HLS INTERFACE m_axi port=a offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=c offset=slave bundle=gmem2

    int TILE_R = size_m/MAX_SIZE;
    int TILE_C = size_n/MAX_SIZE;
    int TILE_K = size_k/MAX_SIZE;
    int BLK_K = size_k/BLK_SIZE;

    int_t feederA[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 complete
    int_t feederB[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederB dim=1 complete
    //------------------------------
    ap_int<Kadd> pSum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = pSum dim = 0 complete
    ap_int<Kadd> accum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = accum dim = 0 complete
    ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete
    ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete
    ap_int<Kadd> sReg_1[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 0 complete

    //------------------------------
	ap_uint<EBIT> ebias[MB][MB];
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete
//	SEXP_T mul_beta[MB][MB];
//#pragma HLS ARRAY_PARTITION variable=mul_beta dim=0 complete
	ap_uint<1> flg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete
	ap_int<W> Zshift[MB][MB];
#pragma HLS ARRAY_PARTITION variable=Zshift dim=1 complete
	SEXP_T betaRes[MB][MB], betasReg[MB][MB], betasReg_1[MB][MB], betaOutReg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 complete
#pragma HLS ARRAY_PARTITION variable=betaOutReg dim=1 complete
//------------------------------
	ap_uint<Kadd-1> ResAND[MAX_SIZE], ResAND_1[MAX_SIZE][2];
	ap_uint<Kadd-1> ResAND_2[MB];
#pragma HLS ARRAY_PARTITION variable = ResAND dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_2 dim = 1 complete
	ap_int<Kadd> accbf1[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf1 dim = 1 complete
	ap_int<Kadd> accbf2[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf2 dim = 1 complete
	//------------------------------
    int_t R[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=R dim=1 complete



	tile_loop:for(int k=0;k<TILE_R*TILE_C*size_k + MAX_SIZE + 2*BLK_SIZE+1;k++){
//		for(int k=0;k<size_k + MAX_SIZE + 2*BLK_SIZE +1;k++){
#pragma HLS PIPELINE

		//------------control and idx---------------
//		bool k_0 = ( k & (size_k-1) )==0;
    	bool k_end = ( (k+1) & (size_k-1) ) ==0;
//    	bool k_blk = (k+1) & ( BLK_SIZE-1 );
    	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;
//    	printf("k_0: %d, k_end: %d, ",k_0,k_end);

    	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;
    	bool SRL_en =  ( k & (size_k-1) ) >=0 && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k;
//    	int tile_idx = ( k & (size_k-1) );
    	ap_uint<2> flag;
    	if( ( k & (size_k-1) )==0 && SRL_en )
    		flag = 0b01;//load data
    	else if( SRL_en && ( (k & (size_k-1) )>0) && ((k & (size_k-1) ) <= MAX_SIZE) )
    		flag = 0b10;//shift
    	else
    		flag = 0b00;
//    	printf("k: %d, PE_en: %d, SRL_en:%d, flag: %d\n",k,feeder_en, SRL_en, flag.to_int());

//    	int ir = int(k/size_k) / TILE_C;
//    	int ic = int(k/size_k) & (TILE_C-1);
    	int ir = ( k >> log_sizek_tilec );
    	int ic = ( k >> log_sizek ) & (TILE_C-1);

    	int kk = k & (size_k-1);//k%(size_k-1)
//    	int tt1 = k % (size_k + MAX_SIZE + 2*BLK_SIZE) - size_k -1;
    	int tt = k > size_k-1 ? ( k &(size_k-1) )-1 : -1;
    	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;
//    	printf("k: %d, kk: %d, tt:%d, wr_tt:%d\t",k,kk,tt, wr_tt);

		int tile_idx = ir*TILE_C+ic-1;
//		int wr_tile_idx = ( int( (k  -1)/size_k) / TILE_C ) *TILE_C + ( int( (k -1)/size_k) & (TILE_C-1)) -1;
		int wr_tile_idx = ( (k  -1) >> log_sizek_tilec ) *TILE_C + ( ( (k -1)>> log_sizek) & (TILE_C-1)) -1;
//		printf("tile_idx: %d, wr_tile_idx: %d, %d\n",tile_idx,wr_tile_idx, int((k - 2*BLK_SIZE)/size_k));
    	//---------------beta idx------------------
    	ap_uint<2> beta_flag;
//    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k < TILE_R*TILE_C*size_k;
    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k >= size_k -1;
    	if( ( (k+1) & (size_k-1) )==0 && beta_SRL_en )
    		beta_flag = 0b01;
    	else if( beta_SRL_en && ( ((k+1) & (size_k-1) )>0) && (((k+1) & (size_k-1) ) <= MB) )
    		beta_flag = 0b10;
    	else
    		beta_flag = 0b00;

    	//---------------PE array-------------------
    	feeder2(a, feederA, (k & (size_k-1))+ir*size_k, feeder_en);
    	feeder2(b, feederB, (k & (size_k-1))+ic*size_k, feeder_en);
//    	printf("feederA idx: %d, feederB idx: %d\n",(k & (size_k-1))+ir*size_k, (k & (size_k-1))+ic*size_k);

    	PEarray(k, size_k, flag, beta_flag, ir, ic, feederA, feederB, mode, sReg_1, betaA, betaB, betasReg_1);

		for(int i = 0; i < MAX_SIZE; i++){
#pragma HLS UNROLL
			//----------------Check Max block----------
			if( (i&(BLK_SIZE-1)) !=0 ){
//				printf("ResAND_1[%d] before:%s, ",i-1,ResAND_1[i-1].to_string(2).c_str());
				CheckMax(tt, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_1[i-1]  );
//				printf("ResAND_1[%d]:%s, \n",i-1,ResAND_1[i-1].to_string(2).c_str());
			}else{
//				printf("ResAND_2[%d] before:%s, ",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
				CheckMax_Beta(tt, tile_idx*MB, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_2[int(i/BLK_SIZE)],
						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], mode, modeout );
//				printf("ResAND_2[%d]:%s, \n",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
			}
		}

		AccumNorm(tt, tile_idx, wr_tt, wr_tile_idx, mode, modeout, sReg_1, Zshift, c);

		//-----beta write out-----
		for(int ip=0; ip<MB;ip++){
#pragma HLS UNROLL
			if(tt>=BLK_SIZE && tt <MAX_SIZE+BLK_SIZE && ( tt & (BLK_SIZE-1) )==BLK_SIZE-1 ){
				betaC[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
//				printf("k:%d, tt: %d, beta: %d, betaC idx:%d\n",k, tt, betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)].to_int(), tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE));
			}
		}


	}//k end


}


void bmGemmv25(MemPack *a, // Read-Only Matrix A
        MemPack *b, // Read-Only Matrix B
		   MemPack *c,
		   SExpPack *betaA,
		   SExpPack *betaB,
		   SExpPack *betaC,
		   int size_m,    // The other Dimension of Matrix A
		   int size_n,    // Shared Dimension of A*B, or A.T*B, A*B.T, A.T*B.T
		   int size_k,     // The other Dimension of Matrix B
		   ap_uint<8> log_sizek,
		   ap_uint<8> log_sizek_tilec,
		   ap_uint<2> mode,
		   ap_uint<2> modeout)
{//assume that size_m, size_k, size_n can be divided by MAX_SIZE
#pragma HLS INTERFACE m_axi port=a offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=c offset=slave bundle=gmem2

    int TILE_R = size_m/MAX_SIZE;
    int TILE_C = size_n/MAX_SIZE;
    int TILE_K = size_k/MAX_SIZE;
    int BLK_K = size_k/BLK_SIZE;

    int_t feederA[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 complete
    int_t feederB[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=feederB dim=1 complete
    //------------------------------
    ap_int<Kadd> pSum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = pSum dim = 0 complete
    ap_int<Kadd> accum[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = accum dim = 0 complete
    ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete
    ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete
    ap_int<Kadd> sReg_1[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 0 complete

    //------------------------------
	ap_uint<EBIT> ebias[MB][MB];
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete
//	SEXP_T mul_beta[MB][MB];
//#pragma HLS ARRAY_PARTITION variable=mul_beta dim=0 complete
	ap_uint<1> flg[MB][MB];
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete
	ap_int<W> Zshift[MB][MB];
#pragma HLS ARRAY_PARTITION variable=Zshift dim=1 complete
	SEXP_T betaRes[MB][MB], betasReg[MB][MB], betasReg_1[MB][MB], betaOutReg[MB][MB],betasReg_11[MB];
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete
#pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 complete
#pragma HLS ARRAY_PARTITION variable=betaOutReg dim=1 complete
#pragma HLS ARRAY_PARTITION variable=betasReg_11 dim=1 complete
//------------------------------
	ap_uint<Kadd-1> ResAND[MAX_SIZE], ResAND_1[MAX_SIZE][2];
	ap_uint<Kadd-1> ResAND_2[MB];
#pragma HLS ARRAY_PARTITION variable = ResAND dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_2 dim = 1 complete
	ap_int<Kadd> accbf1[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf1 dim = 1 complete
	ap_int<Kadd> accbf2[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf2 dim = 1 complete
	//------------------------------
    int_t R[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=R dim=1 complete



	tile_loop:for(int k=0;k<TILE_R*TILE_C*size_k + MAX_SIZE + 2*BLK_SIZE+1;k++){
//		for(int k=0;k<size_k + MAX_SIZE + 2*BLK_SIZE +1;k++){
#pragma HLS PIPELINE

		//------------control and idx---------------
//		bool k_0 = ( k & (size_k-1) )==0;
    	bool k_end = ( (k+1) & (size_k-1) ) ==0;
//    	bool k_blk = (k+1) & ( BLK_SIZE-1 );
    	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;
//    	printf("k_0: %d, k_end: %d, ",k_0,k_end);

    	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;
    	bool SRL_en =  ( k & (size_k-1) ) >=0 && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k;
//    	int tile_idx = ( k & (size_k-1) );
    	ap_uint<2> flag;
    	if( ( k & (size_k-1) )==0 && SRL_en )
    		flag = 0b01;//load data
    	else if( SRL_en && ( (k & (size_k-1) )>0) && ((k & (size_k-1) ) <= MAX_SIZE) )
    		flag = 0b10;//shift
    	else
    		flag = 0b00;
//    	printf("k: %d, PE_en: %d, SRL_en:%d, flag: %d\n",k,feeder_en, SRL_en, flag.to_int());

//    	int ir = int(k/size_k) / TILE_C;
//    	int ic = int(k/size_k) & (TILE_C-1);
    	int ir = ( k >> log_sizek_tilec );
    	int ic = ( k >> log_sizek ) & (TILE_C-1);

    	int kk = k & (size_k-1);//k%(size_k-1)
//    	int tt1 = k % (size_k + MAX_SIZE + 2*BLK_SIZE) - size_k -1;
    	int tt = k > size_k-1 ? ( k &(size_k-1) )-1 : -1;
    	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;
//    	printf("k: %d, kk: %d, tt:%d, wr_tt:%d\t",k,kk,tt, wr_tt);

		int tile_idx = ir*TILE_C+ic-1;
//		int wr_tile_idx = ( int( (k  -1)/size_k) / TILE_C ) *TILE_C + ( int( (k -1)/size_k) & (TILE_C-1)) -1;
		int wr_tile_idx = ( (k  -1) >> log_sizek_tilec ) *TILE_C + ( ( (k -1)>> log_sizek) & (TILE_C-1)) -1;
//		printf("tile_idx: %d, wr_tile_idx: %d, %d\n",tile_idx,wr_tile_idx, int((k - 2*BLK_SIZE)/size_k));
    	//---------------beta idx------------------
    	ap_uint<2> beta_flag;
//    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k < TILE_R*TILE_C*size_k;
    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k >= size_k -1;
    	if( ( (k+1) & (size_k-1) )==0 && beta_SRL_en )
    		beta_flag = 0b01;
    	else if( beta_SRL_en && ( ((k+1) & (size_k-1) )>0) && (((k+1) & (size_k-1) ) <= MB) )
    		beta_flag = 0b10;
    	else
    		beta_flag = 0b00;

    	//---------------PE array-------------------
    	feeder2(a, feederA, (k & (size_k-1))+ir*size_k, feeder_en);
    	feeder2(b, feederB, (k & (size_k-1))+ic*size_k, feeder_en);
//    	printf("feederA idx: %d, feederB idx: %d\n",(k & (size_k-1))+ir*size_k, (k & (size_k-1))+ic*size_k);

    	PEarray3(k, size_k, flag, beta_flag, ir, ic, feederA, feederB, mode, sReg_1, betaA, betaB, betasReg_1);

		for(int i = 0; i < MAX_SIZE; i++){
#pragma HLS UNROLL
			//----------------Check Max block----------
			if( (i&(BLK_SIZE-1)) !=0 ){
//				printf("ResAND_1[%d] before:%s, ",i-1,ResAND_1[i-1].to_string(2).c_str());
				CheckMax(tt, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_1[i-1]  );
//				printf("ResAND_1[%d]:%s, \n",i-1,ResAND_1[i-1].to_string(2).c_str());
			}else{
//				printf("ResAND_2[%d] before:%s, ",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
				CheckMax_Beta(tt, tile_idx*MB, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_2[int(i/BLK_SIZE)],
						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], mode, modeout );
//				printf("ResAND_2[%d]:%s, \n",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
			}
		}

		AccumNorm(tt, tile_idx, wr_tt, wr_tile_idx, mode, modeout, sReg_1, Zshift, c);

		//-----beta write out-----
		for(int ip=0; ip<MB;ip++){
#pragma HLS UNROLL
			if(tt>=BLK_SIZE && tt <MAX_SIZE+BLK_SIZE && ( tt & (BLK_SIZE-1) )==BLK_SIZE-1 ){
				betaC[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
//				printf("k:%d, tt: %d, beta: %d, betaC idx:%d\n",k, tt, betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)].to_int(), tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE));
			}
		}


	}//k end


}

