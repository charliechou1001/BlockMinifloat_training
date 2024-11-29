#include "typedef.h"


int_t DReLu(int_t a, ap_uint<1> drelu, int layer){
	int_t res;
	if(layer != 5 && layer != 7){//ER
		res = drelu != ap_uint<1>(0)? a : int_t(0);
	}else
		res = a;
	return res;
}

int_t RELU(int_t c, ap_uint<1> &drelu, int layer){
	int_t im;
	if(layer != 5 && layer != 7){
		im = c > 0? c: int_t(0);
		drelu = c > 0? ap_uint<1>(1) : ap_uint<1>(0);
	}else{
		im = c;
	}
	return im;
}

void ShiftArray1(int_t feeder[MAX_SIZE], MemPack MM, ap_uint<1> mode, ap_uint<1> outen, unsigned idx){
#pragma HLS INLINE
	static int_t A1[MAX_SIZE][MAX_SIZE];//in
#pragma HLS ARRAY_PARTITION variable=A1 dim=0 type=complete
	static int_t B1[MAX_SIZE][MAX_SIZE];//out
#pragma HLS ARRAY_PARTITION variable=B1 dim=0 type=complete

//	MemPack MMIn;
//	MMIn = MM[idx];

	for(int i=0;i<MAX_SIZE;i++){
		#pragma HLS UNROLL
		for(int j=0;j<MAX_SIZE;j++){
			#pragma HLS UNROLL

			if(mode==0b1){//A->B
				B1[i][j] = A1[i][j];
			}

			//load
			if(i == MAX_SIZE -1)
				A1[MAX_SIZE-1][j] = MM.range(j*W+W-1,j*W);
//				A1[MAX_SIZE-1][j] = MM[idx].range(j*W+W-1,j*W);
			else
				A1[i][j] = A1[i+1][j];
//				printf("[shift] A[MAX_SIZE-1][%d]: %d ",j,A[MAX_SIZE-1][j].to_int());
			if(outen){//transpose out
				if(j==0)
					feeder[i] = B1[i][0];
				else
					B1[i][j-1] = B1[i][j];
			}

		}
//		printf("\n");
	}

}

void betaTrans(SExpPack *betaIn,
			   SEXP_T betaout[MB][BETASIZE],
			   unsigned in_idx,
			   unsigned out_idx,
			   ap_uint<1> mode,  ap_uint<1> outen){
#pragma HLS INLINE

	static SEXP_T A[MB][MB];//in
#pragma HLS ARRAY_PARTITION variable=A dim=0 type=complete
	static SEXP_T B[MB][MB];//out
#pragma HLS ARRAY_PARTITION variable=B dim=0 type=complete

	for(int i=0;i<MB;i++){
		#pragma HLS UNROLL
		for(int j=0;j<MB;j++){
			#pragma HLS UNROLL

			if(mode==0b1){//A->B
				B[i][j] = A[i][j];
			}else{
				//load
				if(i == MB -1)
					A[MB-1][j] = betaIn[in_idx].range(j*SEXP+SEXP-1,j*SEXP);
				else
					A[i][j] = A[i+1][j];
				if(outen){//transpose out
					if(j==0)
						betaout[i][out_idx] = B[i][0];
					else
						B[i][j-1] = B[i][j];
				}
			}
		}
	}

}

void BetaFeeder(SExpPack *ba,
				ap_uint<1> Trans,
				SEXP_T betabf[MB][BETASIZE],
				unsigned TILE_K,
				unsigned beta_addr,
				unsigned len
){
#pragma HLS ARRAY_PARTITION variable=betabf dim=1 type=complete

	ap_uint<1> mode;
	ap_uint<1> outen;

	if(Trans){
		for(int k=0;k<len/MB+1;k++){
			mode = 0b1;
			betaTrans(ba, betabf, 0, 0, mode, outen);
			for(int i=0;i<MB;i++){
				mode = 0b0;
				outen = k==0? 0: 0b1;
				betaTrans(ba, betabf, k*MB+i, (k-1)*MB+i, mode, outen);
			}
		}
	}else{
		for(int k=0;k<len/MB;k++){
			for(int i=0;i<MB;i++){
#pragma HLS UNROLL
				betabf[i][k] = ba[k+beta_addr].range(i*SEXP+SEXP-1, i*SEXP);
			}
		}
	}

}


void ShiftArray2(int_t feeder[MAX_SIZE], MemPack MEMInA, MemPack MEMInB, ap_uint<1> DRELU[MAX_SIZE],
		unsigned idx, ap_uint<1> mode, ap_uint<2> TrainMode, int layer, ap_uint<1> outen){
#pragma HLS INLINE
//	MemPack MEMInBB;

	static int_t A2[MAX_SIZE][MAX_SIZE];//in
#pragma HLS ARRAY_PARTITION variable=A2 dim=0 type=complete
	static int_t B2[MAX_SIZE][MAX_SIZE];//out
#pragma HLS ARRAY_PARTITION variable=B2 dim=0 type=complete
	int_t P[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=P dim=0 complete

//	MEMInBB = MEMInB[idx];

	for(int i=0;i<MAX_SIZE;i++){
		#pragma HLS UNROLL
		for(int j=0;j<MAX_SIZE;j++){
			#pragma HLS UNROLL

			if(mode==0b1){//A->B
				B2[i][j] = A2[i][j];
			}

			//load
			if(i == MAX_SIZE -1){
				if(TrainMode == 0b10){//GD, a.T
//					P[j] = MEMInA[idx].range(j*W+W-1,j*W);
					P[j] = MEMInA.range(j*W+W-1,j*W);
					A2[MAX_SIZE-1][j] = DReLu(P[j], DRELU[j], layer);
				}else{//FW, W.T
//					P[j] = MEMInB[idx].range(j*W+W-1,j*W);
					P[j] = MEMInB.range(j*W+W-1,j*W);
					A2[MAX_SIZE-1][j] = P[j];
				}
			}else
				A2[i][j] = A2[i+1][j];

			if(outen){//transpose out
				if(j==0)
					feeder[i] = B2[i][0];
				else
					B2[i][j-1] = B2[i][j];
			}

		}
	}

//	for(int i=0;i<MAX_SIZE;i++){
//		for(int j=0;j<MAX_SIZE;j++){
//			printf("A2[%d][%d]: %d, ",i,j,A2[i][j].to_int());
//		}printf("\n");
//	}printf("------\n");

//
//	for(int i=0;i<MAX_SIZE;i++){
//		printf("P[%d]: %d, ",i,P[i].to_int());
//	}printf("\n------\n");
}

void feederA(MemPack a,
			 MemPack b,
			 int_t localA[MAX_SIZE],
			 ap_uint<2> TrainMode,
			 ap_uint<1> shiftarray_mode, ap_uint<1> shiftarray_outen, unsigned idx,
			 ap_uint<1> drelu[MAX_SIZE], int block, int layer
){
//#pragma HLS INLINE
//	if(TrainMode==0b10){//GD
//		ShiftArray1( localA, b[idx], shiftarray_mode, shiftarray_outen, idx );
//	}else{
//		load_a:for(int ii=0;ii < MAX_SIZE; ii++){
//			#pragma HLS UNROLL
//			localA[ii] = int_t(a[ idx].range(ii*W+W-1,ii*W));
//			if(TrainMode == 0b01){//ER
//				localA[ii] = DReLu(localA[ii], drelu[ii], layer);
//			}
//		}
//	}

	if(TrainMode==0b10){//GD
//		bIn = b[idx];
		ShiftArray1( localA, b, shiftarray_mode, shiftarray_outen, idx );
	}else{
		load_a:for(int ii=0;ii < MAX_SIZE; ii++){
			#pragma HLS UNROLL
			localA[ii] = int_t(a.range(ii*W+W-1,ii*W));
			if(TrainMode == 0b01){//ER
				localA[ii] = DReLu(localA[ii], drelu[ii], layer);
			}
		}
	}
}

void feederB(MemPack a,
			MemPack b,
			int_t localB[MAX_SIZE],
			ap_uint<2> TrainMode,
			ap_uint<1> shiftarray_mode, ap_uint<1> shiftarray_outen, unsigned idx,
			ap_uint<1> drelu[MAX_SIZE], int block, int layer
){
//#pragma HLS INLINE
//	if(TrainMode == 0b00 || TrainMode == 0b10){
//		ShiftArray2( localB, a, b[idx], drelu, idx, shiftarray_mode, TrainMode, layer, shiftarray_outen );
//	}else{
//		load_b:for(int jj=0;jj<MAX_SIZE;jj++){
//			#pragma HLS UNROLL
//			localB[jj] = int_t(b[idx].range(jj*W+W-1,jj*W));
////		    	        		printf("%d\t",localB[jj].to_int());
//		}
//	}

	if(TrainMode == 0b00 || TrainMode == 0b10){
		ShiftArray2( localB, a, b, drelu, idx, shiftarray_mode, TrainMode, layer, shiftarray_outen );
	}else{
		load_b:for(int jj=0;jj<MAX_SIZE;jj++){
			#pragma HLS UNROLL
			localB[jj] = int_t(b.range(jj*W+W-1,jj*W));
//		    	        		printf("%d\t",localB[jj].to_int());
		}
	}

}


void TileIdx(ap_uint<2> TrainMode, int ir, int ic, int ik, int size_m, int size_n, int size_k,
		int TILE_R, int TILE_C, int TILE_K, int &aidx, int &bidx, int k
){
	if(TrainMode == 0b10 && k>=0){//GD, transposeA, transposeB
		if(ic == TILE_C-1 && ik == TILE_K-1){
			aidx = ((ir+1)%TILE_R)*MAX_SIZE;
		}else
			aidx = size_m*((ik+1)%TILE_K) + ir*MAX_SIZE;

		if(ik == TILE_K-1)
			bidx = ((ic+1)%TILE_C)*MAX_SIZE;
		else
			bidx = size_n*((ik+1)%TILE_K) + ic*MAX_SIZE;
//    				printf("aidx: %d, bidx: %d, ir: %d, ik:%d, ic: %d\n",aidx, bidx, ir, ik, ic);
	}else if(TrainMode == 0b00  && k>=0){//FW, transposeB
		if(ik == TILE_K-1){
			bidx = ((ic+1)%TILE_C)*MAX_SIZE;
		}else
			bidx = size_n*((ik+1)%TILE_K) + ic*MAX_SIZE;

		aidx=0;
	}else{
		aidx=0;bidx=0;
	}
}


//void FeederTest1(MemPack *a, // this function hasn't been tested ,but only used for synthesis for exploration
//        MemPack *b,
//		unsigned W_addr,
//		unsigned a_offset){
//
//#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem
//
//    int size_k=16, size_m=16, size_n=16;
//	int TILE_R = size_m/MAX_SIZE; int TILE_C = size_n/MAX_SIZE; int TILE_K = size_k/MAX_SIZE;
//	unsigned log_sizek=4, log_sizek_tilec=5;
//
//    int_t localA[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localA dim=0 complete
//    int_t localB[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localB dim=0 complete
//
//    static ap_uint<1> DRELU[MAX_SIZE];
//	#pragma HLS ARRAY_PARTITION variable=DRELU dim=0 complete
//
//    int layer=0; int block =0;int Drelu_idx=0;
//    ap_uint<2> TrainMode=0b00;
//    ap_uint<1> shiftarray_mode; ap_uint<1> shiftarray_outen; unsigned idx;
//    int aidx, bidx;
//    unsigned b_idx;
//    unsigned a_idx=0;
//
//    int begin = MAX_SIZE;
//    int end = 2*2*size_k + MAX_SIZE + 2*BLK_SIZE+1;
//	for(int k=-begin;k<end;k++){
//		int ir = k>=0? ( k >> log_sizek_tilec ) : 0;
//		int ic = k>=0? (( k >> log_sizek ) & (TILE_C-1)) :0;
//
//		int kk = k>=0? k & (size_k-1) : k+MAX_SIZE;//k%(size_k-1)
//		int ik = int(kk/MAX_SIZE);
////		printf("ir:%d, ic:%d, kk:%d, ik:%d, (k &(MAX_SIZE-1)):%d\n",ir,ic, kk,ik,(k &(MAX_SIZE-1)));
//
//		if((k &(MAX_SIZE-1)) == 0){
//			TileIdx(TrainMode, ir, ic, ik, size_m, size_n, size_k, TILE_R, TILE_C, TILE_K, aidx, bidx, k);
////			printf("aidx: %d, bidx:%d, \n", aidx,bidx);
//			shiftarray_mode=1;shiftarray_outen=0;
//			if(k>=0 && TrainMode==0b10)
//				feederA(a,b,localA,TrainMode, shiftarray_mode,shiftarray_outen, a_idx, DRELU, block, layer);
//			if(k>=0 && (TrainMode == 0b00 || TrainMode == 0b10))
//				feederB(a,b,localB, TrainMode,shiftarray_mode,shiftarray_outen, b_idx, DRELU, block, Drelu_idx, layer);
//		}
//
//        a_idx = TrainMode == 0b10? aidx + (kk & (MAX_SIZE-1)) + W_addr : size_k*ir + kk + a_offset;
//        if(TrainMode == 0b00)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + W_addr;
//        else if(TrainMode == 0b10)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + a_offset;
//        else
//        	b_idx = size_k*ic + kk + W_addr;
////        printf("a_idx: %d, b_idx:%d\n", a_idx, b_idx);
//
//        shiftarray_mode=0;shiftarray_outen= k>=0? 1:0;
//		feederA(a,b,localA,TrainMode, shiftarray_mode,shiftarray_outen, a_idx, DRELU, block, layer);
//		feederB(a,b,localB, TrainMode,shiftarray_mode,shiftarray_outen, b_idx, DRELU, block, Drelu_idx, layer);
//
//
////		printf("localA[0]: %d, localB[0]:%d\n",localA[0].to_int(), localB[0].to_int());
//	}
//
//}

//void FeederTest2(MemPack *a,
//        MemPack *b,
//		unsigned W_addr,
//		unsigned a_offset){
//
//#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem
//
//    int size_k=16, size_m=16, size_n=16;
//	int TILE_R = size_m/MAX_SIZE; int TILE_C = size_n/MAX_SIZE; int TILE_K = size_k/MAX_SIZE;
//	unsigned log_sizek=4, log_sizek_tilec=5;
//
//    int_t localA[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localA dim=0 complete
//    int_t localB[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localB dim=0 complete
//
//    static ap_uint<1> DRELU[MAX_SIZE];
//	#pragma HLS ARRAY_PARTITION variable=DRELU dim=0 complete
//
//    int layer=0; int block =0;int Drelu_idx=0;
//    ap_uint<2> TrainMode=0b10;
//    ap_uint<1> feederA_mode, feederB_mode; ap_uint<1> shiftarray_outen; unsigned idx;
//    int aidx, bidx;
//    unsigned b_idx;
//    unsigned a_idx=0;
//
//    int begin = MAX_SIZE;
//    int end = 2*2*size_k + MAX_SIZE + 2*BLK_SIZE+1;
//	for(int k=-begin;k<end;k++){
//		int ir = k>=0? ( k >> log_sizek_tilec ) : 0;
//		int ic = k>=0? (( k >> log_sizek ) & (TILE_C-1)) :0;
//
//		int kk = k>=0? k & (size_k-1) : k+MAX_SIZE;//k%(size_k-1)
//		int ik = int(kk/MAX_SIZE);
////		printf("ir:%d, ic:%d, kk:%d, ik:%d, (k &(MAX_SIZE-1)):%d\n",ir,ic, kk,ik,(k &(MAX_SIZE-1)));
//
//		if((k &(MAX_SIZE-1)) == 0){
//			TileIdx(TrainMode, ir, ic, ik, size_m, size_n, size_k, TILE_R, TILE_C, TILE_K, aidx, bidx, k);
////			printf("aidx: %d, bidx:%d, \n", aidx,bidx);
////			shiftarray_mode=1;shiftarray_outen=0;
////			if(k>=0 && TrainMode==0b10)
////				feederA(a,b,localA,TrainMode, shiftarray_mode,shiftarray_outen, a_idx, DRELU, block, layer);
////			if(k>=0 && (TrainMode == 0b00 || TrainMode == 0b10))
////				feederB(a,b,localB, TrainMode,shiftarray_mode,shiftarray_outen, b_idx, DRELU, block, Drelu_idx, layer);
//		}
//
//        a_idx = TrainMode == 0b10? aidx + (kk & (MAX_SIZE-1)) + W_addr : size_k*ir + kk + a_offset;
//        if(TrainMode == 0b00)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + W_addr;
//        else if(TrainMode == 0b10)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + a_offset;
//        else
//        	b_idx = size_k*ic + kk + W_addr;
////        printf("a_idx: %d, b_idx:%d\n", a_idx, b_idx);
//
//        feederA_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && TrainMode==0b10 ? 1 : 0;
//        feederB_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && (TrainMode == 0b00 || TrainMode == 0b10)? 1 : 0;
////        printf("feederB_mode: %d\n",feederB_mode.to_int());
//        shiftarray_outen= k>=0? 1:0;
//		feederA(a,b,localA,TrainMode, feederA_mode,shiftarray_outen, a_idx, DRELU, block, layer);
//		feederB(a,b,localB, TrainMode,feederB_mode,shiftarray_outen, b_idx, DRELU, block, Drelu_idx, layer);
//
////		for(int i=0;i<MAX_SIZE;i++){
////			printf("localB[%d]:%d, ", i, localB[i].to_int());
////		}printf("\n");
//
////		printf("localA[0]: %d, localB[0]:%d\n",localA[0].to_int(), localB[0].to_int());
//	}
//
//}


//void FeederTest21( MemPack *a, // on-chip buffer
//        MemPack *b, // from off-chip
////		   MemPack *c,       // Output Result
////		   BlkPack *ot,
////		   SExpPack *betaA,//on-chip
////		   SExpPack *betaB,//off-chip
////		   SExpPack *betaC,
////		   SExpPack *betaOT,
//		   unsigned block,
//		   unsigned layer,
//		   ap_uint<2> TrainMode,
//		   unsigned W_addr,
//		   unsigned Drelu_addr,
//		   unsigned a_offset,
////		   unsigned betaA_addr,
//		   unsigned betaB_addr
//){
//
//#pragma HLS INTERFACE m_axi port=b offset=slave bundle=gmem
//
//#pragma HLS INTERFACE mode=s_axilite port=block bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=layer bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=TrainMode bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=W_addr bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=Drelu_addr bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=a_offset bundle=control
//#pragma HLS INTERFACE mode=s_axilite port=betaB_addr bundle=control
//
////#pragma HLS INTERFACE m_axi port=betaB offset=slave bundle=gmem1
//
//    int size_k=16, size_m=16, size_n=16;
//	int TILE_R = size_m/MAX_SIZE; int TILE_C = size_n/MAX_SIZE; int TILE_K = size_k/MAX_SIZE;
//	unsigned log_sizek=4, log_sizek_tilec=5;
//
//    int_t localA[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localA dim=0 complete
//    int_t localB[MAX_SIZE];
//   #pragma HLS ARRAY_PARTITION variable=localB dim=0 complete
//
//    static ap_uint<1> DRELU[5][100][MAX_SIZE];
//	#pragma HLS ARRAY_PARTITION variable=DRELU dim=3 complete
//
////    int layer=0; int block =0;
//    int Drelu_idx=0;
////    ap_uint<2> TrainMode=0b10;
//    ap_uint<1> feederA_mode, feederB_mode; ap_uint<1> shiftarray_outen; unsigned idx;
//    int aidx, bidx;
//    unsigned b_idx;
//    unsigned a_idx=0;
////    unsigned Drelu_addr=0;
//
//    int begin = MAX_SIZE;
//    int end = 2*2*size_k + MAX_SIZE + 2*BLK_SIZE+1;
//	for(int k=-begin;k<end;k++){
//#pragma HLS PIPELINE
//		int ir = k>=0? ( k >> log_sizek_tilec ) : 0;
//		int ic = k>=0? (( k >> log_sizek ) & (TILE_C-1)) :0;
//
//		int kk = k>=0? k & (size_k-1) : k+MAX_SIZE;//k%(size_k-1)
//		int ik = int(kk/MAX_SIZE);
////		printf("ir:%d, ic:%d, kk:%d, ik:%d, (k &(MAX_SIZE-1)):%d\n",ir,ic, kk,ik,(k &(MAX_SIZE-1)));
//
//		if((k &(MAX_SIZE-1)) == 0){
//			TileIdx(TrainMode, ir, ic, ik, size_m, size_n, size_k, TILE_R, TILE_C, TILE_K, aidx, bidx, k);
////			printf("aidx: %d, bidx:%d, \n", aidx,bidx);
////			shiftarray_mode=1;shiftarray_outen=0;
////			if(k>=0 && TrainMode==0b10)
////				feederA(a,b,localA,TrainMode, shiftarray_mode,shiftarray_outen, a_idx, DRELU, block, layer);
////			if(k>=0 && (TrainMode == 0b00 || TrainMode == 0b10))
////				feederB(a,b,localB, TrainMode,shiftarray_mode,shiftarray_outen, b_idx, DRELU, block, Drelu_idx, layer);
//		}
//
//        a_idx = TrainMode == 0b10? aidx + (kk & (MAX_SIZE-1)) + W_addr : size_k*ir + kk + a_offset;
//        if(TrainMode == 0b00)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + W_addr;
//        else if(TrainMode == 0b10)
//        	b_idx = (kk & (MAX_SIZE-1)) + bidx + a_offset;
//        else
//        	b_idx = size_k*ic + kk + W_addr;
////        printf("a_idx: %d, b_idx:%d\n", a_idx, b_idx);
//
//        feederA_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && TrainMode==0b10 ? 1 : 0;
//        feederB_mode = ((k &(MAX_SIZE-1)) == 0) && k>=0 && (TrainMode == 0b00 || TrainMode == 0b10)? 1 : 0;
////        printf("feederB_mode: %d\n",feederB_mode.to_int());
//        shiftarray_outen= k>=0? 1:0;
////		feederA(a,b,localA,TrainMode, feederA_mode,shiftarray_outen, a_idx, DRELU[0][Drelu_addr + size_k*ir + kk], block, layer);
////		feederB(a,b,localB, TrainMode,feederB_mode,shiftarray_outen, b_idx, DRELU[0][Drelu_addr + bidx + (kk & (MAX_SIZE-1))], block, Drelu_idx, layer);
//
//		if(TrainMode==0b10){//GD
//	//		bIn = b[idx];
//			ShiftArray1( localA, b[idx], feederA_mode, shiftarray_outen, idx );
//		}else{
//			load_a:for(int ii=0;ii < MAX_SIZE; ii++){
//				#pragma HLS UNROLL
//				localA[ii] = int_t(a[ idx].range(ii*W+W-1,ii*W));
//				if(TrainMode == 0b01){//ER
//					localA[ii] = DReLu(localA[ii], DRELU[0][Drelu_addr + size_k*ir + kk][ii], layer);
//				}
//			}
//		}
//
//		if(TrainMode == 0b00 || TrainMode == 0b10){
//	//		bIn = b[idx];
//			ShiftArray2( localB, a, b[idx], DRELU[0][Drelu_addr + bidx + (kk & (MAX_SIZE-1))], idx, feederB_mode, TrainMode, layer, shiftarray_outen );
//		}else{
//			load_b:for(int jj=0;jj<MAX_SIZE;jj++){
//				#pragma HLS UNROLL
//				localB[jj] = int_t(b[idx].range(jj*W+W-1,jj*W));
//	//		    	        		printf("%d\t",localB[jj].to_int());
//			}
//		}
//
////		for(int i=0;i<MAX_SIZE;i++){
////			printf("localB[%d]:%d, ", i, localB[i].to_int());
////		}printf("\n");
//
////		printf("localA[0]: %d, localB[0]:%d\n",localA[0].to_int(), localB[0].to_int());
//	}
//
//}
