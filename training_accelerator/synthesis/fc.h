#include "typedef.h"
#include "modeldef.h"
#include "MAC.h"
#include "paraconfig.h"
#include "feeder.h"
#include "PEarray.h"

void CheckMax(int tt, ap_int<Kadd> sRegOut, ap_uint<Kadd-1> &ResAND,
//		ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_uint<Kadd-1> max_shift1[2], ap_uint<Kadd-1> max_shift2[2]){
//#pragma HLS INLINE
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

void WriteBuffer2(int tt, int tile_idx, int wr_tt, int wr_tile_idx, ap_uint<2> mode, ap_uint<2> modeout,
		ap_int<Kadd> sRegOut, ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_int<W> Zshift[MB], int_t &R, ap_uint<Ws> &RB, ap_uint<1> out_cond, ap_uint<2> TrainMode, ap_uint<1> relu){
#pragma HLS INLINE

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
		ap_uint<1> relu_en = TrainMode ==0b00 && !out_cond;
		Normalization<W,W,4>(acc, Zshift[MB-1-int(wr_tt/BLK_SIZE)], mode, modeout, R, out_cond, RB, relu_en, relu);

	}
//	return R;

}

void AccumNorm(int tt, int tile_idx, int wr_tt, int wr_tile_idx, ap_uint<2> mode, ap_uint<2> modeout,
		ap_int<Kadd> sReg_1[MAX_SIZE], ap_int<W> Zshift[MB][MB], ap_uint<1> out_cond, MemPack *c, XPack *out,
		ap_uint<2> TrainMode, ap_uint<1> DRELU[MAX_SIZE]){

#pragma HLS INLINE

#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 1 complete
#pragma HLS ARRAY_PARTITION variable = Zshift dim = 1 complete

	static ap_int<Kadd> accbf1[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf1 dim = 1 complete
	static ap_int<Kadd> accbf2[MAX_SIZE][BLK_SIZE*2];
#pragma HLS ARRAY_PARTITION variable = accbf2 dim = 1 complete
    int_t R[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=R dim=1 complete
    ap_uint<Ws> RB[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=RB dim=1 complete
    ap_uint<1> relu[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=relu dim=1 complete

    MemPack RR;
    XPack RBR;

	for(int i = 0; i < MAX_SIZE; i++){
		#pragma HLS UNROLL
			 WriteBuffer2(tt, tile_idx*MB, wr_tt, wr_tile_idx, mode, modeout,
			sReg_1[i], accbf1[i], accbf2[i], Zshift[int(i/BLK_SIZE)], R[i], RB[i], out_cond, TrainMode, relu);
			RR.range(W*i+W-1,W*i) = R[i];
			RBR.range(Ws*i+Ws-1,Ws*i) = RB[i];
	}

	if(wr_tt >= 0 && wr_tt < MAX_SIZE){
		if(out_cond)
			out[wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) )] = RBR;
		else
			c[wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) )] = RR;

		if(TrainMode==0b00 && !out_cond){
			for(int i = 0; i < MAX_SIZE; i++ ){
				#pragma HLS UNROLL
				DRELU[i] = relu[i];
			}
		}
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


void gemm( MemPack *a, // on-chip buffer
           MemPack *b, // from off-chip
		   MemPack *c,       // Output Result
		   BlkPack *ot,
		   SExpPack *betaA,//on-chip
		   SExpPack *betaB,//off-chip
		   SExpPack *betaC,
		   SExpPack *betaOT,
		   unsigned block,
		   unsigned layer,
		   ap_uint<2> TrainMode,
//		   unsigned W_addr,
		   unsigned Drelu_addr,
		   unsigned a_offset,
//		   unsigned betaA_addr,
		   unsigned betaB_addr
) {
//#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem

	int size_m, size_n, size_k;
	unsigned log_sizek, log_sizek_tilec;
	ap_uint<1> InMode;
	GEMMConfig(layer, TrainMode,size_m,size_n,size_k, log_sizek, log_sizek_tilec);
	ap_uint<2> MACmode;
	ap_uint<2> MACmodeOut;
	MACmodeOut = TrainMode;
	if(TrainMode==0b00)
		MACmode = 0b00;
	else if(TrainMode ==0b01)
		MACmode = 0b10;
	else
		MACmode = 0b11;

//    int TILE_R = size_m/MAX_SIZE;
//    int TILE_C = size_n/MAX_SIZE;
	int TILE_R = size_m >= MAX_SIZE? size_m/MAX_SIZE : 1;
	int TILE_C = size_n >= MAX_SIZE? size_n/MAX_SIZE : 1;
    int TILE_K = size_k >= MAX_SIZE? size_k/MAX_SIZE : 1;
    int BLK_K = size_k/BLK_SIZE;

    int_t localA[MAX_SIZE];
   #pragma HLS ARRAY_PARTITION variable=localA dim=0 complete
    int_t localB[MAX_SIZE];
   #pragma HLS ARRAY_PARTITION variable=localB dim=0 complete

    static ap_uint<1> DRELU[NBEATS_BLK][DRELU_LAYER_NUM/MAX_SIZE][MAX_SIZE];
	#pragma HLS ARRAY_PARTITION variable=DRELU dim=3 complete

    ap_int<Kadd> sReg_1[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 0 complete
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
#pragma HLS ARRAY_PARTITION variable = ResAND_1 dim = 0 complete
#pragma HLS ARRAY_PARTITION variable = ResAND_2 dim = 1 complete
//---------------------------

	SEXP_T betaAbf[MB][BETASIZE], betaBbf[MB][BETASIZE];

	ap_uint<1> TransA, TransB;
	if(TrainMode == 0b00 ){
		TransA = 0; TransB = 0b1;
	}else if(TrainMode == 0b10){
		TransA = 0b1; TransB = 0b1;
	}else{
		TransA = 0; TransB = 0;
	}
	BetaFeeder(betaA, TransA, betaAbf, TILE_K, 0 , (size_m/BLK_SIZE)*(size_k/BLK_SIZE) );
	BetaFeeder(betaB, TransB, betaBbf, TILE_K, betaB_addr , (size_n/BLK_SIZE)*(size_k/BLK_SIZE) );

	int aidx=0, bidx=0;
	ap_uint<1> feederA_mode, feederB_mode; ap_uint<1> shiftarray_outen;
    unsigned b_idx;
    unsigned a_idx=0;
    MemPack bIn, aIn;

    int begin = MAX_SIZE;
    int end = TILE_R*TILE_C*size_k + MAX_SIZE + 2*BLK_SIZE+1;
	tile_loop:for(int k=-begin;k<end;k++){
#pragma HLS PIPELINE

		//------------control and idx
    	bool k_end = ( (k+1) & (size_k-1) ) ==0;
//    	bool k_blk = (k+1) & ( BLK_SIZE-1 );
    	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;
//    	printf("k_0: %d, k_end: %d, ",k_0,k_end);

    	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;

    	bool SRL_en =  ( k & (size_k-1) ) >=0 && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k;
//    	int tile_idx = ( k & (size_k-1) );
    	ap_uint<2> flag;
    	if( ( k & (size_k-1) )==0 && k >= size_k )
    		flag = 0b01;//load data
    	else if( ( (k & (size_k-1) )>0) && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k )
    		flag = 0b10;//shift
    	else
    		flag = 0b00;

//    	int ir = int(k/size_k) / TILE_C;
//    	int ic = int(k/size_k) & (TILE_C-1);
//    	int ir = ( k >> log_sizek_tilec );
//    	int ic = ( k >> log_sizek ) & (TILE_C-1);
//    	int kk = k & (size_k-1);//k%(size_k-1)

		int ir = k>=0? ( k >> log_sizek_tilec ) : 0;
		int ic = k>=0? (( k >> log_sizek ) & (TILE_C-1)) :0;
		int kk = k>=0? k & (size_k-1) : k+MAX_SIZE;//k%(size_k-1)
    	int ik = int(kk/MAX_SIZE);

//    	int tt1 = k % (size_k + MAX_SIZE + 2*BLK_SIZE) - size_k -1;
    	int tt = k > size_k-1 ? ( k &(size_k-1) )-1 : -1;
    	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;
//    	printf("k: %d, kk: %d, tt:%d, wr_tt:%d\t",k,kk,tt, wr_tt);

		int tile_idx = ir*TILE_C+ic-1;
		if(tile_idx <0)
			tile_idx = 0;
		int wr_tile_idx = ( (k  -1) >> log_sizek_tilec ) *TILE_C + ( ( (k -1)>> log_sizek) & (TILE_C-1)) -1;
		if(wr_tile_idx < 0)
			wr_tile_idx = 0;
    	//---------------beta idx
    	ap_uint<2> beta_flag;
    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k >= size_k -1;
    	if( ( (k+1) & (size_k-1) )==0 && beta_SRL_en )
    		beta_flag = 0b01;
    	else if( beta_SRL_en && ( ((k+1) & (size_k-1) )>0) && (((k+1) & (size_k-1) ) <= MB) )
    		beta_flag = 0b10;
    	else
    		beta_flag = 0b00;

		if((k &(MAX_SIZE-1)) == 0){
			TileIdx(TrainMode, ir, ic, ik, size_m, size_n, size_k, TILE_R, TILE_C, TILE_K, aidx, bidx, k);
		}

		//------------------data feeder--------------------
		a_idx = TrainMode == 0b10? aidx + (kk & (MAX_SIZE-1)) : size_k*ir + kk + a_offset;
		if(TrainMode == 0b00)
			b_idx = (kk & (MAX_SIZE-1)) + bidx;
		else if(TrainMode == 0b10)
			b_idx = (kk & (MAX_SIZE-1)) + bidx + a_offset;
		else
			b_idx = size_k*ic + kk;
	//        printf("a_idx: %d, b_idx:%d\n", a_idx, b_idx);

		feederA_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && TrainMode==0b10 ? 1 : 0;
		feederB_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && (TrainMode == 0b00 || TrainMode == 0b10)? 1 : 0;
	//        printf("feederB_mode: %d\n",feederB_mode.to_int());
		shiftarray_outen= k>=0? 1:0;
		bIn = b[b_idx];
		aIn = a[a_idx];
		feederA(aIn,bIn,localA,TrainMode, feederA_mode,shiftarray_outen, a_idx, DRELU[block][Drelu_addr + size_k*ir + kk], block, layer);
		feederB(aIn,bIn,localB, TrainMode,feederB_mode,shiftarray_outen, b_idx, DRELU[block][Drelu_addr + bidx + (kk & (MAX_SIZE-1))], block, layer);


		PEarray3(k, size_k, flag, beta_flag, ir, ic, localA, localB, MACmode, sReg_1, betaAbf, betaBbf, betasReg_1);

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
						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], MACmode,MACmodeOut );
//				printf("ResAND_2[%d]:%s, \n",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
			}
		}

		ap_uint<1> out_cond = (layer == 7 && TrainMode ==0b00) || (layer == 5 && TrainMode ==0b00) ||
				(layer == 0 && TrainMode == 0b01);
		int index = (wr_tt >= 0 && wr_tt < MAX_SIZE)? wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) ) : 0;
		AccumNorm(tt, tile_idx, wr_tt, wr_tile_idx, MACmode,MACmodeOut, sReg_1, Zshift,out_cond, c, ot, TrainMode, DRELU[block][Drelu_addr+index]);

		//-----beta write out-----
		for(int ip=0; ip<MB;ip++){
#pragma HLS UNROLL
			if(tt>=BLK_SIZE && tt <MAX_SIZE+BLK_SIZE && ( tt & (BLK_SIZE-1) )==BLK_SIZE-1 ){
				if(out_cond)
					betaOT[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
				else
					betaC[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
//				printf("k:%d, tt: %d, beta: %d, betaC idx:%d\n",k, tt, betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)].to_int(), tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE));
			}
		}


	}//tile loop end


}

// void gemm( MemPack *a, // on-chip buffer
//            MemPack *b, // from off-chip
// 		   MemPack *c,       // Output Result
// 		   XPack *ot,
// 		   SExpPack *betaA,//on-chip
// 		   SExpPack *betaB,//off-chip
// 		   SExpPack *betaC,
// 		   SExpPack *betaOT,
// 		   unsigned block,
// 		   unsigned layer,
// 		   ap_uint<2> TrainMode,
// 		   unsigned W_addr,
// 		   unsigned Drelu_addr,
// 		   unsigned a_offset,
// //		   unsigned betaA_addr,
// 		   unsigned betaB_addr
// ) {

// 	int size_m, size_n, size_k;
// 	unsigned log_sizek, log_sizek_tilec;
// 	ap_uint<1> InMode;
// 	GEMMConfig(layer, TrainMode,size_m,size_n,size_k, log_sizek, log_sizek_tilec);
// 	ap_uint<2> MACmode;
// 	ap_uint<2> MACmodeOut;
// 	MACmodeOut = TrainMode;
// 	if(TrainMode==0b00)
// 		MACmode = 0b00;
// 	else if(TrainMode ==0b01)
// 		MACmode = 0b10;
// 	else
// 		MACmode = 0b11;

//     int TILE_R = size_m/MAX_SIZE;
//     int TILE_C = size_n/MAX_SIZE;
//     int TILE_K = size_k/MAX_SIZE;
//     int BLK_K = size_k/BLK_SIZE;

// //    printf("TILE_R: %d, TILE_C: %d, TILE_K: %d\n", TILE_R, TILE_C, TILE_K);
// //    printf("size_m: %d, size_n: %d, size_k: %d, log_sizek: %d, log_sizek_tilec:%d \n",size_m,size_n,size_k, log_sizek, log_sizek_tilec);

//     int_t localA[MAX_SIZE];
//    #pragma HLS ARRAY_PARTITION variable=localA dim=0 complete
//     int_t localB[MAX_SIZE];
//    #pragma HLS ARRAY_PARTITION variable=localB dim=0 complete

//     static ap_uint<1> DRELU[NBEATS_BLK][DRELU_LAYER_NUM/MAX_SIZE][MAX_SIZE];
// 	#pragma HLS ARRAY_PARTITION variable=DRELU dim=3 complete

//     ap_int<Kadd> sReg_1[MAX_SIZE];
// #pragma HLS ARRAY_PARTITION variable = sReg_1 dim = 0 complete
// 	ap_int<W> Zshift[MB][MB];
// #pragma HLS ARRAY_PARTITION variable=Zshift dim=1 complete
// 	SEXP_T betaRes[MB][MB], betasReg[MB][MB], betasReg_1[MB][MB], betaOutReg[MB][MB];
// #pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete
// #pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete
// #pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 complete
// #pragma HLS ARRAY_PARTITION variable=betaOutReg dim=1 complete
//     //------------------------------
// 	ap_uint<Kadd-1> ResAND[MAX_SIZE], ResAND_1[MAX_SIZE][2];
// 	ap_uint<Kadd-1> ResAND_2[MB];
// #pragma HLS ARRAY_PARTITION variable = ResAND dim = 1 complete
// #pragma HLS ARRAY_PARTITION variable = ResAND_1 dim = 1 complete
// #pragma HLS ARRAY_PARTITION variable = ResAND_2 dim = 1 complete
// //---------------------------


// 	if(TrainMode == 0b10){//GD, transposeA
// 		transA: for(int k = 0; k < MAX_SIZE; k++){
// 			ShiftArray1(localA, b, 0,0, k + W_addr);
// 		}
// 	}
// 	if(TrainMode == 0b00 || TrainMode == 0b10){//FW, GD, transposeB
// 		transB: for(int k = 0; k < MAX_SIZE; k++){
// 			unsigned sft_idx = TrainMode == 0b00? k + W_addr : k + a_offset;//W : a
// 			ShiftArray2(localB, a,b,DRELU[block][Drelu_addr+ k], sft_idx, 0, TrainMode, layer, 0);
// //			if(k==0)
// //				printf("%s\n",b[k+W_addr].to_string(2).c_str());
// 		}
// 	}
// //	printf("layer: %d TrainMode: %d\n",layer,TrainMode.to_int());
// //	printf("size_m: %d, size_k: %d, size_n: %d\n",size_m, size_k,size_n);

// 	SEXP_T betaAbf[MB][BETASIZE], betaBbf[MB][BETASIZE];
// 	ap_uint<1> TransA, TransB;
// 	if(TrainMode == 0b00 ){
// 		TransA = 0; TransB = 0b1;
// 	}else if(TrainMode == 0b10){
// 		TransA = 0b1; TransB = 0b1;
// 	}else{
// 		TransA = 0; TransB = 0;
// 	}
// 	BetaFeeder(betaA, TransA, betaAbf, TILE_K, 0 , (size_m/BLK_SIZE)*(size_k/BLK_SIZE) );
// 	BetaFeeder(betaB, TransB, betaBbf, TILE_K, betaB_addr , (size_n/BLK_SIZE)*(size_k/BLK_SIZE) );

// 	int aidx=0, bidx=0;
// 	tile_loop:for(int k=0;k<TILE_R*TILE_C*size_k + MAX_SIZE + 2*BLK_SIZE+1;k++){
// #pragma HLS PIPELINE

// 		//------------control and idx
//     	bool k_end = ( (k+1) & (size_k-1) ) ==0;
// //    	bool k_blk = (k+1) & ( BLK_SIZE-1 );
//     	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;
// //    	printf("k_0: %d, k_end: %d, ",k_0,k_end);

//     	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;
//     	bool SRL_en =  ( k & (size_k-1) ) >=0 && ((k & (size_k-1) ) <= MAX_SIZE) && k >= size_k;
// //    	int tile_idx = ( k & (size_k-1) );
//     	ap_uint<2> flag;
//     	if( ( k & (size_k-1) )==0 && SRL_en )
//     		flag = 0b01;//load data
//     	else if( SRL_en && ( (k & (size_k-1) )>0) && ((k & (size_k-1) ) <= MAX_SIZE) )
//     		flag = 0b10;//shift
//     	else
//     		flag = 0b00;
// //    	printf("k: %d, PE_en: %d, SRL_en:%d, flag: %d\n",k,feeder_en, SRL_en, flag.to_int());

// //    	int ir = int(k/size_k) / TILE_C;
// //    	int ic = int(k/size_k) & (TILE_C-1);
//     	int ir = ( k >> log_sizek_tilec );
//     	int ic = ( k >> log_sizek ) & (TILE_C-1);

//     	int kk = k & (size_k-1);//k%(size_k-1)
//     	int ik = int(kk/MAX_SIZE);
// //    	int tt1 = k % (size_k + MAX_SIZE + 2*BLK_SIZE) - size_k -1;
//     	int tt = k > size_k-1 ? ( k &(size_k-1) )-1 : -1;
//     	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;
// //    	printf("k: %d, kk: %d, tt:%d, wr_tt:%d\t",k,kk,tt, wr_tt);

// 		int tile_idx = ir*TILE_C+ic-1;
// 		if(tile_idx <0)
// 			tile_idx = 0;
// //		int wr_tile_idx = ( int( (k  -1)/size_k) / TILE_C ) *TILE_C + ( int( (k -1)/size_k) & (TILE_C-1)) -1;
// 		int wr_tile_idx = ( (k  -1) >> log_sizek_tilec ) *TILE_C + ( ( (k -1)>> log_sizek) & (TILE_C-1)) -1;
// 		if(wr_tile_idx < 0)
// 			wr_tile_idx = 0;
// //		printf("tile_idx: %d, wr_tile_idx: %d, %d\n",tile_idx,wr_tile_idx, int((k - 2*BLK_SIZE)/size_k));
//     	//---------------beta idx
//     	ap_uint<2> beta_flag;
// //    	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k < TILE_R*TILE_C*size_k;
//     	bool beta_SRL_en =  ( (k+1) & (size_k-1) ) >=0 && (( (k+1) & (size_k-1) ) <= MB) && k >= size_k -1;
//     	if( ( (k+1) & (size_k-1) )==0 && beta_SRL_en )
//     		beta_flag = 0b01;
//     	else if( beta_SRL_en && ( ((k+1) & (size_k-1) )>0) && (((k+1) & (size_k-1) ) <= MB) )
//     		beta_flag = 0b10;
//     	else
//     		beta_flag = 0b00;

// 		if((k &(MAX_SIZE-1)) == 0){
// 			if(TrainMode == 0b10){//GD, transposeA, transposeB
// 				ShiftArray1(localA, b, 1,0, W_addr);//A->B,outen=0,mode=1
// 			}
// 			if(TrainMode == 0b00 || TrainMode == 0b10 ){
// 				ShiftArray2(localB, a,b,DRELU[block][Drelu_addr], 0, 1, TrainMode, layer, 0);//A->B,outen=0,mode=1
// 			}
// 			TileIdx(TrainMode, ir, ic, ik, size_m, size_n, size_k, TILE_R, TILE_C, TILE_K, aidx, bidx);
// //			printf("ir: %d, ic: %d, ik:%d\n",ir, ic, ik);
// 		}

// 		//------------------data feeder--------------------
// 		if(TrainMode == 0b10){//GD
// 			ShiftArray1(localA, b, 0, 1, aidx + (kk & (MAX_SIZE-1)) + W_addr);
// //			printf("GD a: %d, aidx:%d", aidx + (kk & (MAX_SIZE-1)), aidx);
// 		}else{
// 			load_a:for(int ii=0;ii < MAX_SIZE; ii++){
// 				#pragma HLS UNROLL
// 				localA[ii] = int_t(a[ size_k*ir + kk + a_offset].range(ii*W+W-1,ii*W));
// 				if(TrainMode == 0b01){//ER
// 					localA[ii] = DReLu(localA[ii], DRELU[block][Drelu_addr + size_k*ir + kk][ii], layer);
// 				}
// //				else//FW
// //					localA[ii] = localA[ii];
// 			}
// //			printf("ER/FW a: %d, ", size_k*ir + kk + a_offset );
// 		}

// 		if(TrainMode == 0b00 || TrainMode == 0b10){
// //							ShiftArray2(localB, b, 0, 1, bidx + k);
// 			unsigned sft_idx = TrainMode == 0b00? (kk & (MAX_SIZE-1)) + bidx + W_addr : (kk & (MAX_SIZE-1)) + bidx + a_offset;//W : a
// 			ShiftArray2(localB, a,b,DRELU[block][Drelu_addr + bidx + (kk & (MAX_SIZE-1))], sft_idx, 0, TrainMode, layer, 1);
// //			printf("FW/GD b: %d\n", sft_idx);
// 		}else{
// //		    	        	printf("localB:\n");
// 			load_b:for(int jj=0;jj<MAX_SIZE;jj++){
// 				#pragma HLS UNROLL
// 				localB[jj] = int_t(b[ size_k*ic + kk + W_addr].range(jj*W+W-1,jj*W));
// 			}
// //			printf("ER b: %d\n", size_k*ic + kk + W_addr);
// 		}

// 		PEarray4(k, size_k, flag, beta_flag, ir, ic, localA, localB, MACmode, sReg_1, betaAbf, betaBbf, betasReg_1);

// 		for(int i = 0; i < MAX_SIZE; i++){
// #pragma HLS UNROLL
// 			//----------------Check Max block----------
// 			if( (i&(BLK_SIZE-1)) !=0 ){
// //				printf("ResAND_1[%d] before:%s, ",i-1,ResAND_1[i-1].to_string(2).c_str());
// 				CheckMax(tt, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_1[i-1]  );
// //				printf("ResAND_1[%d]:%s, \n",i-1,ResAND_1[i-1].to_string(2).c_str());
// 			}else{
// //				printf("ResAND_2[%d] before:%s, ",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
// 				CheckMax_Beta(tt, tile_idx*MB, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_2[int(i/BLK_SIZE)],
// 						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], MACmode,MACmodeOut );
// //				printf("ResAND_2[%d]:%s, \n",int(i/BLK_SIZE),ResAND_2[int(i/BLK_SIZE)].to_string(2).c_str());
// 			}
// 		}

// 		ap_uint<1> out_cond = (layer == 7 && TrainMode ==0b00) || (layer == 5 && TrainMode ==0b00) ||
// 				(layer == 0 && TrainMode == 0b01);
// 		AccumNorm(tt, tile_idx, wr_tt, wr_tile_idx, MACmode,MACmodeOut, sReg_1, Zshift,out_cond, c, ot);
// //		printf("k:%d, ir: %d, ic:%d, tt: %d, tile_idx: %d,  wr_tt: %d, wr_tile_idx: %d\n",k,ir,ic, tt, tile_idx, wr_tt, wr_tile_idx);
// //		if(wr_tt >= 0 && wr_tt < MAX_SIZE)
// //			printf("wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) ): %d\n",wr_tile_idx*MAX_SIZE + (MB-1-int(wr_tt/BLK_SIZE))*BLK_SIZE + ( wr_tt & (BLK_SIZE-1) ));

// 		//-----beta write out-----
// 		for(int ip=0; ip<MB;ip++){
// #pragma HLS UNROLL
// 			if(tt>=BLK_SIZE && tt <MAX_SIZE+BLK_SIZE && ( tt & (BLK_SIZE-1) )==BLK_SIZE-1 ){
// 				if(out_cond)
// 					betaOT[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
// 				else
// 					betaC[tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip) = betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)];
// //				printf("k:%d, tt: %d, beta: %d, betaC idx:%d\n",k, tt, betaOutReg[ip][int((tt-BLK_SIZE)/BLK_SIZE)].to_int(), tile_idx*MB+MB-1-int((tt-BLK_SIZE)/BLK_SIZE));
// 			}
// 		}


// 	}//tile loop end


// }

