#include "typedef.h"
#include "modeldef.h"
#include "MAC.h"
#include "paraconfig.h"
#include "feeder.h"
#include "PEarray.h"

//#define CM_DEBUG

void CheckMax(int tt, ap_int<Kadd> sRegOut, ap_uint<Kadd-1> &ResAND,
		ap_uint<Kadd-1> max_shift1[2], ap_uint<Kadd-1> max_shift2[2]){
#pragma HLS INLINE
	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	ap_uint<1> pp_idx = int(tt/BLK_SIZE) & 0b1;

	if(tt >=0 && tt < MAX_SIZE ){//data shift out

#ifdef CM_DEBUG
		printf("ResAND before: %s, ", ResAND.to_string(2).c_str() );
#endif

		//check max number of each column in PE array
		if( ((tt+1)&(BLK_SIZE-1)) == 0  ){//tt == BLK_SIZE*m-1
			max_shift1[~pp_idx] = ResAND | GetUnsigned(sRegOut);
			ResAND = 0b0;//clear register
		}else
			ResAND = ResAND | GetUnsigned(sRegOut);

#ifdef CM_DEBUG
		printf("ResAND: %s, unsigned sRegOut:%s, max_shift1[%d]: %s\n",ResAND.to_string(2).c_str(), GetUnsigned(sRegOut).to_string(2).c_str(),
				(~pp_idx).to_int(), max_shift1[~pp_idx].to_string(2).c_str() );
#endif
	}

	if(tt >= BLK_SIZE && tt<MAX_SIZE+BLK_SIZE){//max number shift vertically
		max_shift2[pp_idx] = max_shift1[pp_idx];
	}

}

void WriteBuffer2(int tt, int tile_idx, int wr_tt, int wr_tile_idx, ap_uint<2> mode, ap_uint<2> modeout,
		ap_int<Kadd> sRegOut, ap_int<Kadd> accbf1[BLK_SIZE*2],ap_int<Kadd> accbf2[BLK_SIZE*2],
		ap_int<W> Zshift[MB], int_t &R, ap_uint<Ws> &RB, ap_uint<1> out_cond, ap_uint<2> TrainMode, ap_uint<1> &relu){
#pragma HLS INLINE

	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	if(tt >=0 && tt < MAX_SIZE ){//data shift out

		ap_uint<2> idx = int(tt/BLK_SIZE+tile_idx) & 3;
		if(idx[1] == 0b0){
			accbf1[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
		}else{
			accbf2[BLK_SIZE*idx[0] + BLK_SIZE -1 - (tt_mod)] = sRegOut;
		}
	}

	if(wr_tt >= 0 && wr_tt < MAX_SIZE){
		ap_uint<2> wr_idx = int(wr_tt/BLK_SIZE + wr_tile_idx*MB) & 3;//////
		ap_uint<Kadd> acc;
		if(wr_idx[1] == 0b0){
			acc = accbf1[BLK_SIZE*wr_idx[0] + ( wr_tt & (BLK_SIZE-1) ) ];
		}else{
			acc = accbf2[BLK_SIZE*wr_idx[0] + ( wr_tt & (BLK_SIZE-1) ) ];
		}

		ap_uint<1> relu_en = TrainMode ==0b00 && !out_cond;

		Normalization<W,W,4>(acc, Zshift[MB-1-int(wr_tt/BLK_SIZE)], mode, modeout, R, out_cond, RB, relu_en, relu);
	}

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
			sReg_1[i], accbf1[i], accbf2[i], Zshift[int(i/BLK_SIZE)], R[i], RB[i], out_cond, TrainMode, relu[i]);
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
		ap_uint<Kadd-1> max_shift1[2], ap_uint<Kadd-1> &max_shift2,
		SEXP_T betaOut[MB], SEXP_T betaC[MB], ap_int<W> Zshift[MB], ap_uint<2> mode,  ap_uint<2> modeout){
#pragma HLS INLINE

	int tt_mod = tt & (BLK_SIZE-1);//tt%BLK_SIZE
	ap_uint<1> pp_idx = int(tt/BLK_SIZE) & 0b1;


	if(tt >=0 && tt < MAX_SIZE ){//data shift out
#ifdef CM_DEBUG
		printf("=ResAND before: %s, ",ResAND.to_string(2).c_str() );
#endif
		//check max number of each column in PE array
		if( ((tt+1)&(BLK_SIZE-1)) == 0  ){//tt == BLK_SIZE*m-1
			max_shift1[~pp_idx] = ResAND | GetUnsigned(sRegOut);
			ResAND = 0b0;//clear register
		}else
			ResAND = ResAND | GetUnsigned(sRegOut);
#ifdef CM_DEBUG
		printf("ResAND: %s, unsigned sRegOut:%s, max_shift1[%d]: %s\n",ResAND.to_string(2).c_str(), GetUnsigned(sRegOut).to_string(2).c_str(),
				(~pp_idx).to_int(), max_shift1[~pp_idx].to_string(2).c_str() );
#endif
	}

	if(tt >= BLK_SIZE && tt<MAX_SIZE+BLK_SIZE){//max number shift vertically
		if( tt_mod == 0 )
			max_shift2 = max_shift1[pp_idx];
		else if( (tt_mod>=0) && (tt_mod<=BLK_SIZE-1) )
			max_shift2 = max_shift2 | max_shift1[pp_idx];

		if( tt_mod == BLK_SIZE-1){
			ap_uint<W> Zmin = CheckZmin<W>(max_shift2);
			int idx = int(tt/BLK_SIZE)-1;
			betaC[idx] = BiasAdjustv3<W,W>( betaOut[idx], Zmin, Zshift[MB-1-idx], mode, modeout );

		}
	}

}


void gemm( MemPack *a, // on-chip buffer
           MemPack *b, // from off-chip, weight
		   MemPack *b2,// activation for BP
		   MemPack *c,       // Output Result
		   BlkPack *ot,
		   SExpPack *betaA,//on-chip
		   SExpPack *betaB,//off-chip
		   SExpPack *betaC,
		   SExpPack *betaOT,
		   unsigned block,
		   unsigned layer,
		   ap_uint<2> TrainMode,
		   unsigned Drelu_addr,
		   unsigned a_offset,
		   unsigned betaB_addr
) {

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
	BetaFeederA(betaA, betaB, TransA, betaAbf, size_m, size_k, 0, betaB_addr );
	BetaFeederB(betaA, betaB, TransB, betaBbf,  size_k, size_n, 0, betaB_addr , TrainMode );

	for(int i=0;i<MAX_SIZE;i++){
		ResAND[i] = 0;
	}


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
    	bool k_blk_start =  ( int(k/BLK_SIZE) & (BLK_K-1) )==0;

    	bool feeder_en = k>=0 && k < size_k*TILE_R*TILE_C;

    	bool SRL_en =  k >= size_k-1;
    	ap_uint<2> flag;

    		if( ( (k+1) & (size_k-1) )== 0 && SRL_en )//k = size_k -1
				flag = 0b01;//load data
			else if( ( (k & (size_k-1) )>=0) && ((k & (size_k-1) ) < MAX_SIZE) && SRL_en )
				flag = 0b10;//shift
			else
				flag = 0b00;

		int ir = k>=0? ( k >> log_sizek_tilec ) : 0;
		int ic = k>=0? (( k >> log_sizek ) & (TILE_C-1)) :0;
		int kk = k>=0? k & (size_k-1) : k+MAX_SIZE;//k%(size_k-1)
    	int ik = int(kk/MAX_SIZE);

    	int tt = k > size_k-1 ? ( k &(size_k-1) ) : -1;


    	int wr_tt = k > size_k ? ( (k - size_k - 2*BLK_SIZE -1) &(size_k-1) ) : -1 ;

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
		a_idx = TrainMode == 0b10? bidx + (kk & (MAX_SIZE-1)) + a_offset : size_k*ir + kk + a_offset;//new2

		if(TrainMode == 0b00)
			b_idx = (kk & (MAX_SIZE-1)) + bidx;
		else if(TrainMode == 0b10)
			b_idx = (kk & (MAX_SIZE-1)) + aidx;//new2
		else
			b_idx = size_k*ic + kk;

		feederA_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && TrainMode==0b10 ? 1 : 0;
		feederB_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && (TrainMode == 0b00 || TrainMode == 0b10)? 1 : 0;
		shiftarray_outen= k>=0? 1:0;
		bIn = TrainMode == 0b10? b2[b_idx] : b[b_idx];
		aIn = a[a_idx];
		feederA(aIn,bIn,localA, TrainMode, feederA_mode,shiftarray_outen, a_idx, DRELU[block][Drelu_addr + size_k*ir + kk], block, layer);
		feederB(aIn,bIn,localB, TrainMode, feederB_mode,shiftarray_outen, b_idx, DRELU[block][Drelu_addr + bidx + (kk & (MAX_SIZE-1))], block, layer);

		PEarray3(k, size_k, flag, beta_flag, ir, ic, localA, localB, MACmode, sReg_1, betaAbf, betaBbf, betasReg_1);

		for(int i = 0; i < MAX_SIZE; i++){
#pragma HLS UNROLL
			//----------------Check Max block----------
#ifdef CM_DEBUG
			printf("[%d] ", i);
#endif
			if( (i&(BLK_SIZE-1)) !=0 ){
				CheckMax(tt, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_1[i-1]  );
			}else{
				CheckMax_Beta(tt, tile_idx*MB, sReg_1[i], ResAND[i], ResAND_1[i], ResAND_2[int(i/BLK_SIZE)],
						betasReg_1[int(i/BLK_SIZE)], betaOutReg[int(i/BLK_SIZE)],  Zshift[int(i/BLK_SIZE)], MACmode,MACmodeOut );
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
			}
		}


	}//tile loop end


}


