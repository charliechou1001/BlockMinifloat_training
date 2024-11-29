#include "ap_int.h" 
#include "PE.h" 

void PEarray3(int k, int size_k, ap_uint<2> flag, ap_uint<2> beta_flag, int ir, int ic, 
                   int_t feederA[MAX_SIZE], int_t feederB[MAX_SIZE], ap_uint<2> mode, ap_int<Kadd> sReg_1[MAX_SIZE], 
                   SExpPack *betaA, SExpPack *betaB, SEXP_T betasReg_1[MB][MB]){
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=feederB dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=sReg_1 dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 type=complete

static ap_uint<EBIT> ebias[MB][MB];
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete 
static ap_uint<1> flg[MB][MB]; 
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete 
static SEXP_T betaRes[MB][MB], betasReg[MB][MB]; 
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete 
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete

static ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete
    static ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete

    int TILE_K = size_k/MAX_SIZE;
    int kk =  k & (size_k-1);//k%(size_k-1)
    bool k_blk_start =  ( int(k/BLK_SIZE) & (size_k/BLK_SIZE-1) )==0;
    bool k_end = ( (k+1) & (size_k-1) ) ==0;

    //---------------shared exp----------------
    betax:for(int ip=0;ip<MB;ip++){
	#pragma HLS UNROLL
	betay:for(int iq=MB-1;iq>=0;iq--){
		#pragma HLS UNROLL
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

    //------------------PE array-----------------
    out[0][15] = PE3<0,15>(kk, k_end, mode, feederA[0], feederB[15], ebias[0][1], flg[0][1]);
    out[0][14] = PE3<0,14>(kk, k_end, mode, feederA[0], feederB[14], ebias[0][1], flg[0][1]);
    out[0][13] = PE3<0,13>(kk, k_end, mode, feederA[0], feederB[13], ebias[0][1], flg[0][1]);
    out[0][12] = PE3<0,12>(kk, k_end, mode, feederA[0], feederB[12], ebias[0][1], flg[0][1]);
    out[0][11] = PE3<0,11>(kk, k_end, mode, feederA[0], feederB[11], ebias[0][1], flg[0][1]);
    out[0][10] = PE3<0,10>(kk, k_end, mode, feederA[0], feederB[10], ebias[0][1], flg[0][1]);
    out[0][9] = PE3<0,9>(kk, k_end, mode, feederA[0], feederB[9], ebias[0][1], flg[0][1]);
    out[0][8] = PE3<0,8>(kk, k_end, mode, feederA[0], feederB[8], ebias[0][1], flg[0][1]);
    out[0][7] = PE3<0,7>(kk, k_end, mode, feederA[0], feederB[7], ebias[0][0], flg[0][0]);
    out[0][6] = PE3<0,6>(kk, k_end, mode, feederA[0], feederB[6], ebias[0][0], flg[0][0]);
    out[0][5] = PE3<0,5>(kk, k_end, mode, feederA[0], feederB[5], ebias[0][0], flg[0][0]);
    out[0][4] = PE3<0,4>(kk, k_end, mode, feederA[0], feederB[4], ebias[0][0], flg[0][0]);
    out[0][3] = PE3<0,3>(kk, k_end, mode, feederA[0], feederB[3], ebias[0][0], flg[0][0]);
    out[0][2] = PE3<0,2>(kk, k_end, mode, feederA[0], feederB[2], ebias[0][0], flg[0][0]);
    out[0][1] = PE3<0,1>(kk, k_end, mode, feederA[0], feederB[1], ebias[0][0], flg[0][0]);
    out[0][0] = PE3<0,0>(kk, k_end, mode, feederA[0], feederB[0], ebias[0][0], flg[0][0]);
    out[1][15] = PE3<1,15>(kk, k_end, mode, feederA[1], feederB[15], ebias[0][1], flg[0][1]);
    out[1][14] = PE3<1,14>(kk, k_end, mode, feederA[1], feederB[14], ebias[0][1], flg[0][1]);
    out[1][13] = PE3<1,13>(kk, k_end, mode, feederA[1], feederB[13], ebias[0][1], flg[0][1]);
    out[1][12] = PE3<1,12>(kk, k_end, mode, feederA[1], feederB[12], ebias[0][1], flg[0][1]);
    out[1][11] = PE3<1,11>(kk, k_end, mode, feederA[1], feederB[11], ebias[0][1], flg[0][1]);
    out[1][10] = PE3<1,10>(kk, k_end, mode, feederA[1], feederB[10], ebias[0][1], flg[0][1]);
    out[1][9] = PE3<1,9>(kk, k_end, mode, feederA[1], feederB[9], ebias[0][1], flg[0][1]);
    out[1][8] = PE3<1,8>(kk, k_end, mode, feederA[1], feederB[8], ebias[0][1], flg[0][1]);
    out[1][7] = PE3<1,7>(kk, k_end, mode, feederA[1], feederB[7], ebias[0][0], flg[0][0]);
    out[1][6] = PE3<1,6>(kk, k_end, mode, feederA[1], feederB[6], ebias[0][0], flg[0][0]);
    out[1][5] = PE3<1,5>(kk, k_end, mode, feederA[1], feederB[5], ebias[0][0], flg[0][0]);
    out[1][4] = PE3<1,4>(kk, k_end, mode, feederA[1], feederB[4], ebias[0][0], flg[0][0]);
    out[1][3] = PE3<1,3>(kk, k_end, mode, feederA[1], feederB[3], ebias[0][0], flg[0][0]);
    out[1][2] = PE3<1,2>(kk, k_end, mode, feederA[1], feederB[2], ebias[0][0], flg[0][0]);
    out[1][1] = PE3<1,1>(kk, k_end, mode, feederA[1], feederB[1], ebias[0][0], flg[0][0]);
    out[1][0] = PE3<1,0>(kk, k_end, mode, feederA[1], feederB[0], ebias[0][0], flg[0][0]);
    out[2][15] = PE3<2,15>(kk, k_end, mode, feederA[2], feederB[15], ebias[0][1], flg[0][1]);
    out[2][14] = PE3<2,14>(kk, k_end, mode, feederA[2], feederB[14], ebias[0][1], flg[0][1]);
    out[2][13] = PE3<2,13>(kk, k_end, mode, feederA[2], feederB[13], ebias[0][1], flg[0][1]);
    out[2][12] = PE3<2,12>(kk, k_end, mode, feederA[2], feederB[12], ebias[0][1], flg[0][1]);
    out[2][11] = PE3<2,11>(kk, k_end, mode, feederA[2], feederB[11], ebias[0][1], flg[0][1]);
    out[2][10] = PE3<2,10>(kk, k_end, mode, feederA[2], feederB[10], ebias[0][1], flg[0][1]);
    out[2][9] = PE3<2,9>(kk, k_end, mode, feederA[2], feederB[9], ebias[0][1], flg[0][1]);
    out[2][8] = PE3<2,8>(kk, k_end, mode, feederA[2], feederB[8], ebias[0][1], flg[0][1]);
    out[2][7] = PE3<2,7>(kk, k_end, mode, feederA[2], feederB[7], ebias[0][0], flg[0][0]);
    out[2][6] = PE3<2,6>(kk, k_end, mode, feederA[2], feederB[6], ebias[0][0], flg[0][0]);
    out[2][5] = PE3<2,5>(kk, k_end, mode, feederA[2], feederB[5], ebias[0][0], flg[0][0]);
    out[2][4] = PE3<2,4>(kk, k_end, mode, feederA[2], feederB[4], ebias[0][0], flg[0][0]);
    out[2][3] = PE3<2,3>(kk, k_end, mode, feederA[2], feederB[3], ebias[0][0], flg[0][0]);
    out[2][2] = PE3<2,2>(kk, k_end, mode, feederA[2], feederB[2], ebias[0][0], flg[0][0]);
    out[2][1] = PE3<2,1>(kk, k_end, mode, feederA[2], feederB[1], ebias[0][0], flg[0][0]);
    out[2][0] = PE3<2,0>(kk, k_end, mode, feederA[2], feederB[0], ebias[0][0], flg[0][0]);
    out[3][15] = PE3<3,15>(kk, k_end, mode, feederA[3], feederB[15], ebias[0][1], flg[0][1]);
    out[3][14] = PE3<3,14>(kk, k_end, mode, feederA[3], feederB[14], ebias[0][1], flg[0][1]);
    out[3][13] = PE3<3,13>(kk, k_end, mode, feederA[3], feederB[13], ebias[0][1], flg[0][1]);
    out[3][12] = PE3<3,12>(kk, k_end, mode, feederA[3], feederB[12], ebias[0][1], flg[0][1]);
    out[3][11] = PE3<3,11>(kk, k_end, mode, feederA[3], feederB[11], ebias[0][1], flg[0][1]);
    out[3][10] = PE3<3,10>(kk, k_end, mode, feederA[3], feederB[10], ebias[0][1], flg[0][1]);
    out[3][9] = PE3<3,9>(kk, k_end, mode, feederA[3], feederB[9], ebias[0][1], flg[0][1]);
    out[3][8] = PE3<3,8>(kk, k_end, mode, feederA[3], feederB[8], ebias[0][1], flg[0][1]);
    out[3][7] = PE3<3,7>(kk, k_end, mode, feederA[3], feederB[7], ebias[0][0], flg[0][0]);
    out[3][6] = PE3<3,6>(kk, k_end, mode, feederA[3], feederB[6], ebias[0][0], flg[0][0]);
    out[3][5] = PE3<3,5>(kk, k_end, mode, feederA[3], feederB[5], ebias[0][0], flg[0][0]);
    out[3][4] = PE3<3,4>(kk, k_end, mode, feederA[3], feederB[4], ebias[0][0], flg[0][0]);
    out[3][3] = PE3<3,3>(kk, k_end, mode, feederA[3], feederB[3], ebias[0][0], flg[0][0]);
    out[3][2] = PE3<3,2>(kk, k_end, mode, feederA[3], feederB[2], ebias[0][0], flg[0][0]);
    out[3][1] = PE3<3,1>(kk, k_end, mode, feederA[3], feederB[1], ebias[0][0], flg[0][0]);
    out[3][0] = PE3<3,0>(kk, k_end, mode, feederA[3], feederB[0], ebias[0][0], flg[0][0]);
    out[4][15] = PE3<4,15>(kk, k_end, mode, feederA[4], feederB[15], ebias[0][1], flg[0][1]);
    out[4][14] = PE3<4,14>(kk, k_end, mode, feederA[4], feederB[14], ebias[0][1], flg[0][1]);
    out[4][13] = PE3<4,13>(kk, k_end, mode, feederA[4], feederB[13], ebias[0][1], flg[0][1]);
    out[4][12] = PE3<4,12>(kk, k_end, mode, feederA[4], feederB[12], ebias[0][1], flg[0][1]);
    out[4][11] = PE3<4,11>(kk, k_end, mode, feederA[4], feederB[11], ebias[0][1], flg[0][1]);
    out[4][10] = PE3<4,10>(kk, k_end, mode, feederA[4], feederB[10], ebias[0][1], flg[0][1]);
    out[4][9] = PE3<4,9>(kk, k_end, mode, feederA[4], feederB[9], ebias[0][1], flg[0][1]);
    out[4][8] = PE3<4,8>(kk, k_end, mode, feederA[4], feederB[8], ebias[0][1], flg[0][1]);
    out[4][7] = PE3<4,7>(kk, k_end, mode, feederA[4], feederB[7], ebias[0][0], flg[0][0]);
    out[4][6] = PE3<4,6>(kk, k_end, mode, feederA[4], feederB[6], ebias[0][0], flg[0][0]);
    out[4][5] = PE3<4,5>(kk, k_end, mode, feederA[4], feederB[5], ebias[0][0], flg[0][0]);
    out[4][4] = PE3<4,4>(kk, k_end, mode, feederA[4], feederB[4], ebias[0][0], flg[0][0]);
    out[4][3] = PE3<4,3>(kk, k_end, mode, feederA[4], feederB[3], ebias[0][0], flg[0][0]);
    out[4][2] = PE3<4,2>(kk, k_end, mode, feederA[4], feederB[2], ebias[0][0], flg[0][0]);
    out[4][1] = PE3<4,1>(kk, k_end, mode, feederA[4], feederB[1], ebias[0][0], flg[0][0]);
    out[4][0] = PE3<4,0>(kk, k_end, mode, feederA[4], feederB[0], ebias[0][0], flg[0][0]);
    out[5][15] = PE3<5,15>(kk, k_end, mode, feederA[5], feederB[15], ebias[0][1], flg[0][1]);
    out[5][14] = PE3<5,14>(kk, k_end, mode, feederA[5], feederB[14], ebias[0][1], flg[0][1]);
    out[5][13] = PE3<5,13>(kk, k_end, mode, feederA[5], feederB[13], ebias[0][1], flg[0][1]);
    out[5][12] = PE3<5,12>(kk, k_end, mode, feederA[5], feederB[12], ebias[0][1], flg[0][1]);
    out[5][11] = PE3<5,11>(kk, k_end, mode, feederA[5], feederB[11], ebias[0][1], flg[0][1]);
    out[5][10] = PE3<5,10>(kk, k_end, mode, feederA[5], feederB[10], ebias[0][1], flg[0][1]);
    out[5][9] = PE3<5,9>(kk, k_end, mode, feederA[5], feederB[9], ebias[0][1], flg[0][1]);
    out[5][8] = PE3<5,8>(kk, k_end, mode, feederA[5], feederB[8], ebias[0][1], flg[0][1]);
    out[5][7] = PE3<5,7>(kk, k_end, mode, feederA[5], feederB[7], ebias[0][0], flg[0][0]);
    out[5][6] = PE3<5,6>(kk, k_end, mode, feederA[5], feederB[6], ebias[0][0], flg[0][0]);
    out[5][5] = PE3<5,5>(kk, k_end, mode, feederA[5], feederB[5], ebias[0][0], flg[0][0]);
    out[5][4] = PE3<5,4>(kk, k_end, mode, feederA[5], feederB[4], ebias[0][0], flg[0][0]);
    out[5][3] = PE3<5,3>(kk, k_end, mode, feederA[5], feederB[3], ebias[0][0], flg[0][0]);
    out[5][2] = PE3<5,2>(kk, k_end, mode, feederA[5], feederB[2], ebias[0][0], flg[0][0]);
    out[5][1] = PE3<5,1>(kk, k_end, mode, feederA[5], feederB[1], ebias[0][0], flg[0][0]);
    out[5][0] = PE3<5,0>(kk, k_end, mode, feederA[5], feederB[0], ebias[0][0], flg[0][0]);
    out[6][15] = PE3<6,15>(kk, k_end, mode, feederA[6], feederB[15], ebias[0][1], flg[0][1]);
    out[6][14] = PE3<6,14>(kk, k_end, mode, feederA[6], feederB[14], ebias[0][1], flg[0][1]);
    out[6][13] = PE3<6,13>(kk, k_end, mode, feederA[6], feederB[13], ebias[0][1], flg[0][1]);
    out[6][12] = PE3<6,12>(kk, k_end, mode, feederA[6], feederB[12], ebias[0][1], flg[0][1]);
    out[6][11] = PE3<6,11>(kk, k_end, mode, feederA[6], feederB[11], ebias[0][1], flg[0][1]);
    out[6][10] = PE3<6,10>(kk, k_end, mode, feederA[6], feederB[10], ebias[0][1], flg[0][1]);
    out[6][9] = PE3<6,9>(kk, k_end, mode, feederA[6], feederB[9], ebias[0][1], flg[0][1]);
    out[6][8] = PE3<6,8>(kk, k_end, mode, feederA[6], feederB[8], ebias[0][1], flg[0][1]);
    out[6][7] = PE3<6,7>(kk, k_end, mode, feederA[6], feederB[7], ebias[0][0], flg[0][0]);
    out[6][6] = PE3<6,6>(kk, k_end, mode, feederA[6], feederB[6], ebias[0][0], flg[0][0]);
    out[6][5] = PE3<6,5>(kk, k_end, mode, feederA[6], feederB[5], ebias[0][0], flg[0][0]);
    out[6][4] = PE3<6,4>(kk, k_end, mode, feederA[6], feederB[4], ebias[0][0], flg[0][0]);
    out[6][3] = PE3<6,3>(kk, k_end, mode, feederA[6], feederB[3], ebias[0][0], flg[0][0]);
    out[6][2] = PE3<6,2>(kk, k_end, mode, feederA[6], feederB[2], ebias[0][0], flg[0][0]);
    out[6][1] = PE3<6,1>(kk, k_end, mode, feederA[6], feederB[1], ebias[0][0], flg[0][0]);
    out[6][0] = PE3<6,0>(kk, k_end, mode, feederA[6], feederB[0], ebias[0][0], flg[0][0]);
    out[7][15] = PE3<7,15>(kk, k_end, mode, feederA[7], feederB[15], ebias[0][1], flg[0][1]);
    out[7][14] = PE3<7,14>(kk, k_end, mode, feederA[7], feederB[14], ebias[0][1], flg[0][1]);
    out[7][13] = PE3<7,13>(kk, k_end, mode, feederA[7], feederB[13], ebias[0][1], flg[0][1]);
    out[7][12] = PE3<7,12>(kk, k_end, mode, feederA[7], feederB[12], ebias[0][1], flg[0][1]);
    out[7][11] = PE3<7,11>(kk, k_end, mode, feederA[7], feederB[11], ebias[0][1], flg[0][1]);
    out[7][10] = PE3<7,10>(kk, k_end, mode, feederA[7], feederB[10], ebias[0][1], flg[0][1]);
    out[7][9] = PE3<7,9>(kk, k_end, mode, feederA[7], feederB[9], ebias[0][1], flg[0][1]);
    out[7][8] = PE3<7,8>(kk, k_end, mode, feederA[7], feederB[8], ebias[0][1], flg[0][1]);
    out[7][7] = PE3<7,7>(kk, k_end, mode, feederA[7], feederB[7], ebias[0][0], flg[0][0]);
    out[7][6] = PE3<7,6>(kk, k_end, mode, feederA[7], feederB[6], ebias[0][0], flg[0][0]);
    out[7][5] = PE3<7,5>(kk, k_end, mode, feederA[7], feederB[5], ebias[0][0], flg[0][0]);
    out[7][4] = PE3<7,4>(kk, k_end, mode, feederA[7], feederB[4], ebias[0][0], flg[0][0]);
    out[7][3] = PE3<7,3>(kk, k_end, mode, feederA[7], feederB[3], ebias[0][0], flg[0][0]);
    out[7][2] = PE3<7,2>(kk, k_end, mode, feederA[7], feederB[2], ebias[0][0], flg[0][0]);
    out[7][1] = PE3<7,1>(kk, k_end, mode, feederA[7], feederB[1], ebias[0][0], flg[0][0]);
    out[7][0] = PE3<7,0>(kk, k_end, mode, feederA[7], feederB[0], ebias[0][0], flg[0][0]);
    out[8][15] = PE3<8,15>(kk, k_end, mode, feederA[8], feederB[15], ebias[1][1], flg[1][1]);
    out[8][14] = PE3<8,14>(kk, k_end, mode, feederA[8], feederB[14], ebias[1][1], flg[1][1]);
    out[8][13] = PE3<8,13>(kk, k_end, mode, feederA[8], feederB[13], ebias[1][1], flg[1][1]);
    out[8][12] = PE3<8,12>(kk, k_end, mode, feederA[8], feederB[12], ebias[1][1], flg[1][1]);
    out[8][11] = PE3<8,11>(kk, k_end, mode, feederA[8], feederB[11], ebias[1][1], flg[1][1]);
    out[8][10] = PE3<8,10>(kk, k_end, mode, feederA[8], feederB[10], ebias[1][1], flg[1][1]);
    out[8][9] = PE3<8,9>(kk, k_end, mode, feederA[8], feederB[9], ebias[1][1], flg[1][1]);
    out[8][8] = PE3<8,8>(kk, k_end, mode, feederA[8], feederB[8], ebias[1][1], flg[1][1]);
    out[8][7] = PE3<8,7>(kk, k_end, mode, feederA[8], feederB[7], ebias[1][0], flg[1][0]);
    out[8][6] = PE3<8,6>(kk, k_end, mode, feederA[8], feederB[6], ebias[1][0], flg[1][0]);
    out[8][5] = PE3<8,5>(kk, k_end, mode, feederA[8], feederB[5], ebias[1][0], flg[1][0]);
    out[8][4] = PE3<8,4>(kk, k_end, mode, feederA[8], feederB[4], ebias[1][0], flg[1][0]);
    out[8][3] = PE3<8,3>(kk, k_end, mode, feederA[8], feederB[3], ebias[1][0], flg[1][0]);
    out[8][2] = PE3<8,2>(kk, k_end, mode, feederA[8], feederB[2], ebias[1][0], flg[1][0]);
    out[8][1] = PE3<8,1>(kk, k_end, mode, feederA[8], feederB[1], ebias[1][0], flg[1][0]);
    out[8][0] = PE3<8,0>(kk, k_end, mode, feederA[8], feederB[0], ebias[1][0], flg[1][0]);
    out[9][15] = PE3<9,15>(kk, k_end, mode, feederA[9], feederB[15], ebias[1][1], flg[1][1]);
    out[9][14] = PE3<9,14>(kk, k_end, mode, feederA[9], feederB[14], ebias[1][1], flg[1][1]);
    out[9][13] = PE3<9,13>(kk, k_end, mode, feederA[9], feederB[13], ebias[1][1], flg[1][1]);
    out[9][12] = PE3<9,12>(kk, k_end, mode, feederA[9], feederB[12], ebias[1][1], flg[1][1]);
    out[9][11] = PE3<9,11>(kk, k_end, mode, feederA[9], feederB[11], ebias[1][1], flg[1][1]);
    out[9][10] = PE3<9,10>(kk, k_end, mode, feederA[9], feederB[10], ebias[1][1], flg[1][1]);
    out[9][9] = PE3<9,9>(kk, k_end, mode, feederA[9], feederB[9], ebias[1][1], flg[1][1]);
    out[9][8] = PE3<9,8>(kk, k_end, mode, feederA[9], feederB[8], ebias[1][1], flg[1][1]);
    out[9][7] = PE3<9,7>(kk, k_end, mode, feederA[9], feederB[7], ebias[1][0], flg[1][0]);
    out[9][6] = PE3<9,6>(kk, k_end, mode, feederA[9], feederB[6], ebias[1][0], flg[1][0]);
    out[9][5] = PE3<9,5>(kk, k_end, mode, feederA[9], feederB[5], ebias[1][0], flg[1][0]);
    out[9][4] = PE3<9,4>(kk, k_end, mode, feederA[9], feederB[4], ebias[1][0], flg[1][0]);
    out[9][3] = PE3<9,3>(kk, k_end, mode, feederA[9], feederB[3], ebias[1][0], flg[1][0]);
    out[9][2] = PE3<9,2>(kk, k_end, mode, feederA[9], feederB[2], ebias[1][0], flg[1][0]);
    out[9][1] = PE3<9,1>(kk, k_end, mode, feederA[9], feederB[1], ebias[1][0], flg[1][0]);
    out[9][0] = PE3<9,0>(kk, k_end, mode, feederA[9], feederB[0], ebias[1][0], flg[1][0]);
    out[10][15] = PE3<10,15>(kk, k_end, mode, feederA[10], feederB[15], ebias[1][1], flg[1][1]);
    out[10][14] = PE3<10,14>(kk, k_end, mode, feederA[10], feederB[14], ebias[1][1], flg[1][1]);
    out[10][13] = PE3<10,13>(kk, k_end, mode, feederA[10], feederB[13], ebias[1][1], flg[1][1]);
    out[10][12] = PE3<10,12>(kk, k_end, mode, feederA[10], feederB[12], ebias[1][1], flg[1][1]);
    out[10][11] = PE3<10,11>(kk, k_end, mode, feederA[10], feederB[11], ebias[1][1], flg[1][1]);
    out[10][10] = PE3<10,10>(kk, k_end, mode, feederA[10], feederB[10], ebias[1][1], flg[1][1]);
    out[10][9] = PE3<10,9>(kk, k_end, mode, feederA[10], feederB[9], ebias[1][1], flg[1][1]);
    out[10][8] = PE3<10,8>(kk, k_end, mode, feederA[10], feederB[8], ebias[1][1], flg[1][1]);
    out[10][7] = PE3<10,7>(kk, k_end, mode, feederA[10], feederB[7], ebias[1][0], flg[1][0]);
    out[10][6] = PE3<10,6>(kk, k_end, mode, feederA[10], feederB[6], ebias[1][0], flg[1][0]);
    out[10][5] = PE3<10,5>(kk, k_end, mode, feederA[10], feederB[5], ebias[1][0], flg[1][0]);
    out[10][4] = PE3<10,4>(kk, k_end, mode, feederA[10], feederB[4], ebias[1][0], flg[1][0]);
    out[10][3] = PE3<10,3>(kk, k_end, mode, feederA[10], feederB[3], ebias[1][0], flg[1][0]);
    out[10][2] = PE3<10,2>(kk, k_end, mode, feederA[10], feederB[2], ebias[1][0], flg[1][0]);
    out[10][1] = PE3<10,1>(kk, k_end, mode, feederA[10], feederB[1], ebias[1][0], flg[1][0]);
    out[10][0] = PE3<10,0>(kk, k_end, mode, feederA[10], feederB[0], ebias[1][0], flg[1][0]);
    out[11][15] = PE3<11,15>(kk, k_end, mode, feederA[11], feederB[15], ebias[1][1], flg[1][1]);
    out[11][14] = PE3<11,14>(kk, k_end, mode, feederA[11], feederB[14], ebias[1][1], flg[1][1]);
    out[11][13] = PE3<11,13>(kk, k_end, mode, feederA[11], feederB[13], ebias[1][1], flg[1][1]);
    out[11][12] = PE3<11,12>(kk, k_end, mode, feederA[11], feederB[12], ebias[1][1], flg[1][1]);
    out[11][11] = PE3<11,11>(kk, k_end, mode, feederA[11], feederB[11], ebias[1][1], flg[1][1]);
    out[11][10] = PE3<11,10>(kk, k_end, mode, feederA[11], feederB[10], ebias[1][1], flg[1][1]);
    out[11][9] = PE3<11,9>(kk, k_end, mode, feederA[11], feederB[9], ebias[1][1], flg[1][1]);
    out[11][8] = PE3<11,8>(kk, k_end, mode, feederA[11], feederB[8], ebias[1][1], flg[1][1]);
    out[11][7] = PE3<11,7>(kk, k_end, mode, feederA[11], feederB[7], ebias[1][0], flg[1][0]);
    out[11][6] = PE3<11,6>(kk, k_end, mode, feederA[11], feederB[6], ebias[1][0], flg[1][0]);
    out[11][5] = PE3<11,5>(kk, k_end, mode, feederA[11], feederB[5], ebias[1][0], flg[1][0]);
    out[11][4] = PE3<11,4>(kk, k_end, mode, feederA[11], feederB[4], ebias[1][0], flg[1][0]);
    out[11][3] = PE3<11,3>(kk, k_end, mode, feederA[11], feederB[3], ebias[1][0], flg[1][0]);
    out[11][2] = PE3<11,2>(kk, k_end, mode, feederA[11], feederB[2], ebias[1][0], flg[1][0]);
    out[11][1] = PE3<11,1>(kk, k_end, mode, feederA[11], feederB[1], ebias[1][0], flg[1][0]);
    out[11][0] = PE3<11,0>(kk, k_end, mode, feederA[11], feederB[0], ebias[1][0], flg[1][0]);
    out[12][15] = PE3<12,15>(kk, k_end, mode, feederA[12], feederB[15], ebias[1][1], flg[1][1]);
    out[12][14] = PE3<12,14>(kk, k_end, mode, feederA[12], feederB[14], ebias[1][1], flg[1][1]);
    out[12][13] = PE3<12,13>(kk, k_end, mode, feederA[12], feederB[13], ebias[1][1], flg[1][1]);
    out[12][12] = PE3<12,12>(kk, k_end, mode, feederA[12], feederB[12], ebias[1][1], flg[1][1]);
    out[12][11] = PE3<12,11>(kk, k_end, mode, feederA[12], feederB[11], ebias[1][1], flg[1][1]);
    out[12][10] = PE3<12,10>(kk, k_end, mode, feederA[12], feederB[10], ebias[1][1], flg[1][1]);
    out[12][9] = PE3<12,9>(kk, k_end, mode, feederA[12], feederB[9], ebias[1][1], flg[1][1]);
    out[12][8] = PE3<12,8>(kk, k_end, mode, feederA[12], feederB[8], ebias[1][1], flg[1][1]);
    out[12][7] = PE3<12,7>(kk, k_end, mode, feederA[12], feederB[7], ebias[1][0], flg[1][0]);
    out[12][6] = PE3<12,6>(kk, k_end, mode, feederA[12], feederB[6], ebias[1][0], flg[1][0]);
    out[12][5] = PE3<12,5>(kk, k_end, mode, feederA[12], feederB[5], ebias[1][0], flg[1][0]);
    out[12][4] = PE3<12,4>(kk, k_end, mode, feederA[12], feederB[4], ebias[1][0], flg[1][0]);
    out[12][3] = PE3<12,3>(kk, k_end, mode, feederA[12], feederB[3], ebias[1][0], flg[1][0]);
    out[12][2] = PE3<12,2>(kk, k_end, mode, feederA[12], feederB[2], ebias[1][0], flg[1][0]);
    out[12][1] = PE3<12,1>(kk, k_end, mode, feederA[12], feederB[1], ebias[1][0], flg[1][0]);
    out[12][0] = PE3<12,0>(kk, k_end, mode, feederA[12], feederB[0], ebias[1][0], flg[1][0]);
    out[13][15] = PE3<13,15>(kk, k_end, mode, feederA[13], feederB[15], ebias[1][1], flg[1][1]);
    out[13][14] = PE3<13,14>(kk, k_end, mode, feederA[13], feederB[14], ebias[1][1], flg[1][1]);
    out[13][13] = PE3<13,13>(kk, k_end, mode, feederA[13], feederB[13], ebias[1][1], flg[1][1]);
    out[13][12] = PE3<13,12>(kk, k_end, mode, feederA[13], feederB[12], ebias[1][1], flg[1][1]);
    out[13][11] = PE3<13,11>(kk, k_end, mode, feederA[13], feederB[11], ebias[1][1], flg[1][1]);
    out[13][10] = PE3<13,10>(kk, k_end, mode, feederA[13], feederB[10], ebias[1][1], flg[1][1]);
    out[13][9] = PE3<13,9>(kk, k_end, mode, feederA[13], feederB[9], ebias[1][1], flg[1][1]);
    out[13][8] = PE3<13,8>(kk, k_end, mode, feederA[13], feederB[8], ebias[1][1], flg[1][1]);
    out[13][7] = PE3<13,7>(kk, k_end, mode, feederA[13], feederB[7], ebias[1][0], flg[1][0]);
    out[13][6] = PE3<13,6>(kk, k_end, mode, feederA[13], feederB[6], ebias[1][0], flg[1][0]);
    out[13][5] = PE3<13,5>(kk, k_end, mode, feederA[13], feederB[5], ebias[1][0], flg[1][0]);
    out[13][4] = PE3<13,4>(kk, k_end, mode, feederA[13], feederB[4], ebias[1][0], flg[1][0]);
    out[13][3] = PE3<13,3>(kk, k_end, mode, feederA[13], feederB[3], ebias[1][0], flg[1][0]);
    out[13][2] = PE3<13,2>(kk, k_end, mode, feederA[13], feederB[2], ebias[1][0], flg[1][0]);
    out[13][1] = PE3<13,1>(kk, k_end, mode, feederA[13], feederB[1], ebias[1][0], flg[1][0]);
    out[13][0] = PE3<13,0>(kk, k_end, mode, feederA[13], feederB[0], ebias[1][0], flg[1][0]);
    out[14][15] = PE3<14,15>(kk, k_end, mode, feederA[14], feederB[15], ebias[1][1], flg[1][1]);
    out[14][14] = PE3<14,14>(kk, k_end, mode, feederA[14], feederB[14], ebias[1][1], flg[1][1]);
    out[14][13] = PE3<14,13>(kk, k_end, mode, feederA[14], feederB[13], ebias[1][1], flg[1][1]);
    out[14][12] = PE3<14,12>(kk, k_end, mode, feederA[14], feederB[12], ebias[1][1], flg[1][1]);
    out[14][11] = PE3<14,11>(kk, k_end, mode, feederA[14], feederB[11], ebias[1][1], flg[1][1]);
    out[14][10] = PE3<14,10>(kk, k_end, mode, feederA[14], feederB[10], ebias[1][1], flg[1][1]);
    out[14][9] = PE3<14,9>(kk, k_end, mode, feederA[14], feederB[9], ebias[1][1], flg[1][1]);
    out[14][8] = PE3<14,8>(kk, k_end, mode, feederA[14], feederB[8], ebias[1][1], flg[1][1]);
    out[14][7] = PE3<14,7>(kk, k_end, mode, feederA[14], feederB[7], ebias[1][0], flg[1][0]);
    out[14][6] = PE3<14,6>(kk, k_end, mode, feederA[14], feederB[6], ebias[1][0], flg[1][0]);
    out[14][5] = PE3<14,5>(kk, k_end, mode, feederA[14], feederB[5], ebias[1][0], flg[1][0]);
    out[14][4] = PE3<14,4>(kk, k_end, mode, feederA[14], feederB[4], ebias[1][0], flg[1][0]);
    out[14][3] = PE3<14,3>(kk, k_end, mode, feederA[14], feederB[3], ebias[1][0], flg[1][0]);
    out[14][2] = PE3<14,2>(kk, k_end, mode, feederA[14], feederB[2], ebias[1][0], flg[1][0]);
    out[14][1] = PE3<14,1>(kk, k_end, mode, feederA[14], feederB[1], ebias[1][0], flg[1][0]);
    out[14][0] = PE3<14,0>(kk, k_end, mode, feederA[14], feederB[0], ebias[1][0], flg[1][0]);
    out[15][15] = PE3<15,15>(kk, k_end, mode, feederA[15], feederB[15], ebias[1][1], flg[1][1]);
    out[15][14] = PE3<15,14>(kk, k_end, mode, feederA[15], feederB[14], ebias[1][1], flg[1][1]);
    out[15][13] = PE3<15,13>(kk, k_end, mode, feederA[15], feederB[13], ebias[1][1], flg[1][1]);
    out[15][12] = PE3<15,12>(kk, k_end, mode, feederA[15], feederB[12], ebias[1][1], flg[1][1]);
    out[15][11] = PE3<15,11>(kk, k_end, mode, feederA[15], feederB[11], ebias[1][1], flg[1][1]);
    out[15][10] = PE3<15,10>(kk, k_end, mode, feederA[15], feederB[10], ebias[1][1], flg[1][1]);
    out[15][9] = PE3<15,9>(kk, k_end, mode, feederA[15], feederB[9], ebias[1][1], flg[1][1]);
    out[15][8] = PE3<15,8>(kk, k_end, mode, feederA[15], feederB[8], ebias[1][1], flg[1][1]);
    out[15][7] = PE3<15,7>(kk, k_end, mode, feederA[15], feederB[7], ebias[1][0], flg[1][0]);
    out[15][6] = PE3<15,6>(kk, k_end, mode, feederA[15], feederB[6], ebias[1][0], flg[1][0]);
    out[15][5] = PE3<15,5>(kk, k_end, mode, feederA[15], feederB[5], ebias[1][0], flg[1][0]);
    out[15][4] = PE3<15,4>(kk, k_end, mode, feederA[15], feederB[4], ebias[1][0], flg[1][0]);
    out[15][3] = PE3<15,3>(kk, k_end, mode, feederA[15], feederB[3], ebias[1][0], flg[1][0]);
    out[15][2] = PE3<15,2>(kk, k_end, mode, feederA[15], feederB[2], ebias[1][0], flg[1][0]);
    out[15][1] = PE3<15,1>(kk, k_end, mode, feederA[15], feederB[1], ebias[1][0], flg[1][0]);
    out[15][0] = PE3<15,0>(kk, k_end, mode, feederA[15], feederB[0], ebias[1][0], flg[1][0]);

    for(int i=0;i<MAX_SIZE;i++){
		#pragma HLS UNROLL
    	for(int j=MAX_SIZE-1;j>=0;j--){
			#pragma HLS UNROLL
    		if(j!=MAX_SIZE-1){
        		if(flag==0b01)
        			sReg[i][j] = out[i][j];
        		else if(flag==0b10)
        			sReg[i][j+1] = sReg[i][j];
    		}else{
        		if(flag==0b01)
        			sReg[i][j] = out[i][j];
        		else if(flag==0b10)
        			sReg_1[i] = sReg[i][MAX_SIZE-1];
    		}
    	}
    }
}
