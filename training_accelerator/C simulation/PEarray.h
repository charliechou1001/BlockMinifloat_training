#include "ap_int.h" 
#include "PE.h" 
#include "Print.h"

void PEarray3(int k, int size_k, ap_uint<2> flag, ap_uint<2> beta_flag, int ir, int ic, 
                   int_t feederA[MAX_SIZE], int_t feederB[MAX_SIZE], ap_uint<2> mode, ap_int<Kadd> sReg_1[MAX_SIZE], 
                   SEXP_T betaAbf[MB][BETASIZE], SEXP_T betaBbf[MB][BETASIZE], SEXP_T betasReg_1[MB][MB]){
#pragma HLS INLINE
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=feederB dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=betaAbf dim=1 type=complete 
 #pragma HLS ARRAY_PARTITION variable=betaBbf dim=1 type=complete 
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

    bool k_tile_start = ( k & (size_k-1) ) ==0;

    //---------------shared exp----------------
    betax:for(int ip=0;ip<MB;ip++){
	#pragma HLS UNROLL
	betay:for(int iq=MB-1;iq>=0;iq--){
		#pragma HLS UNROLL
		if(iq==MB-1)
			beta_SRL_2(k, k_tile_start, betaAbf[ip][ir*TILE_K*MB+int(kk/BLK_SIZE)], betaBbf[iq][ic*TILE_K*MB+int(kk/BLK_SIZE)],
					betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg_1[ip], beta_flag[1]*( ( (k+1) & (size_k-1) ) -1) );
		else
			beta_SRL(k, k_tile_start, betaAbf[ip][ir*TILE_K*MB+int(kk/BLK_SIZE)], betaBbf[iq][ic*TILE_K*MB+int(kk/BLK_SIZE)],
					betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg[ip][iq+1] );

	}
   }

    //------------------PE array-----------------
    out[0][3] = PE3<0,3>(kk, k_end, mode, feederA[0], feederB[3], ebias[0][1], flg[0][1]);
    out[0][2] = PE3<0,2>(kk, k_end, mode, feederA[0], feederB[2], ebias[0][1], flg[0][1]);
    out[0][1] = PE3<0,1>(kk, k_end, mode, feederA[0], feederB[1], ebias[0][0], flg[0][0]);
    out[0][0] = PE3<0,0>(kk, k_end, mode, feederA[0], feederB[0], ebias[0][0], flg[0][0]);
    out[1][3] = PE3<1,3>(kk, k_end, mode, feederA[1], feederB[3], ebias[0][1], flg[0][1]);
    out[1][2] = PE3<1,2>(kk, k_end, mode, feederA[1], feederB[2], ebias[0][1], flg[0][1]);
    out[1][1] = PE3<1,1>(kk, k_end, mode, feederA[1], feederB[1], ebias[0][0], flg[0][0]);
    out[1][0] = PE3<1,0>(kk, k_end, mode, feederA[1], feederB[0], ebias[0][0], flg[0][0]);
    out[2][3] = PE3<2,3>(kk, k_end, mode, feederA[2], feederB[3], ebias[1][1], flg[1][1]);
    out[2][2] = PE3<2,2>(kk, k_end, mode, feederA[2], feederB[2], ebias[1][1], flg[1][1]);
    out[2][1] = PE3<2,1>(kk, k_end, mode, feederA[2], feederB[1], ebias[1][0], flg[1][0]);
    out[2][0] = PE3<2,0>(kk, k_end, mode, feederA[2], feederB[0], ebias[1][0], flg[1][0]);
    out[3][3] = PE3<3,3>(kk, k_end, mode, feederA[3], feederB[3], ebias[1][1], flg[1][1]);
    out[3][2] = PE3<3,2>(kk, k_end, mode, feederA[3], feederB[2], ebias[1][1], flg[1][1]);
    out[3][1] = PE3<3,1>(kk, k_end, mode, feederA[3], feederB[1], ebias[1][0], flg[1][0]);
    out[3][0] = PE3<3,0>(kk, k_end, mode, feederA[3], feederB[0], ebias[1][0], flg[1][0]);

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
    			if(flag==0b01 || flag==0b10 )
    			sReg_1[i] = sReg[i][MAX_SIZE-1];
    			if(flag==0b01)
    				sReg[i][j] = out[i][j];

    		}
    	}
    }


}
