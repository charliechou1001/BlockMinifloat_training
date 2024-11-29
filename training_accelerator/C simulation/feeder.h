#include "typedef.h"


int_t DReLu(int_t a, ap_uint<1> drelu, int layer){
	int_t res;
	if(layer != 5 && layer != 7){//ER
		res = drelu != ap_uint<1>(0)? a : int_t(0);
	}else
		res = a;
	return res;
}


void ShiftArray1(int_t feeder[MAX_SIZE], MemPack MM, ap_uint<1> mode, ap_uint<1> outen, unsigned idx){
#pragma HLS INLINE
	static int_t A1[MAX_SIZE][MAX_SIZE];//in
#pragma HLS ARRAY_PARTITION variable=A1 dim=0 type=complete
	static int_t B1[MAX_SIZE][MAX_SIZE];//out
#pragma HLS ARRAY_PARTITION variable=B1 dim=0 type=complete

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
			else
				A1[i][j] = A1[i+1][j];
			if(outen){//transpose out
				if(j==0)
					feeder[i] = B1[i][0];
				else
					B1[i][j-1] = B1[i][j];
			}

		}
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
			}

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

void betaTransB(SExpPack *betaIn1,
				SExpPack *betaIn2,
			   SEXP_T betaout[MB][BETASIZE],
			   unsigned in_idx_a,
			   unsigned in_idx_w,
			   unsigned out_idx,
			   ap_uint<1> mode,  ap_uint<1> outen, ap_uint<2> TrainMode){
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
			}
			//load
			if(i == MB -1){
				if(TrainMode == 0b10){//GD, a.T
					A[MB-1][j] = betaIn1[in_idx_a].range(j*SEXP+SEXP-1,j*SEXP);//from betaA
				}else//FW, W.T
					A[MB-1][j] = betaIn2[in_idx_w].range(j*SEXP+SEXP-1,j*SEXP);//from betaB
			}else
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

void BetaFeederA(SExpPack *ba,
				SExpPack *bb,
				ap_uint<1> Trans,
				SEXP_T betabf[MB][BETASIZE],
				unsigned size_m,
				unsigned size_k,
				unsigned betaA_addr,
				unsigned betaB_addr
){
#pragma HLS ARRAY_PARTITION variable=betabf dim=1 type=complete

	ap_uint<1> mode;
	ap_uint<1> outen;
	unsigned len = (size_m*size_k)/(MAX_SIZE*BLK_SIZE);
	if(Trans){
		unsigned blk_n = size_m /MAX_SIZE;
		unsigned blk_k = size_k / MAX_SIZE;
		for(int k=0;k<len/MB+1;k++){
			mode = 0b1;
			for(int i=0;i<MB;i++){

				mode = i==0? 0b1: 0b0;
				outen = k==0? 0: 0b1;
				betaTrans(bb, betabf, betaB_addr + MB*blk_k*(k%blk_n ) + int(k/blk_n)*MB +i, (k-1)*MB+i, mode, outen);
			}
		}
	}else{

		for(int k=0;k<len;k++){
			for(int i=0;i<MB;i++){
#pragma HLS UNROLL
				betabf[i][k] = ba[k+betaA_addr].range(i*SEXP+SEXP-1, i*SEXP);
			}
		}
	}

}

void BetaFeederB(SExpPack *ba,
				SExpPack *bb,
				ap_uint<1> Trans,
				SEXP_T betabf[MB][BETASIZE],
				unsigned size_k,
				unsigned size_n,
				unsigned betaA_addr,
				unsigned betaB_addr,
				ap_uint<2> TrainMode
){
#pragma HLS ARRAY_PARTITION variable=betabf dim=1 type=complete

	ap_uint<1> mode;
	ap_uint<1> outen;
	unsigned blk_n = size_n /MAX_SIZE;
	unsigned blk_k = size_k / MAX_SIZE;
	unsigned len = (size_n*size_k)/(MAX_SIZE*BLK_SIZE);

	if(Trans){
		for(int k=0;k<len/MB+1;k++){
			mode = 0b1;
			for(int i=0;i<MB;i++){

				mode = i==0? 0b1: 0b0;
				outen = k==0? 0: 0b1;
				betaTransB(ba, bb, betabf, betaA_addr + MB*blk_k*(k%blk_n ) + int(k/blk_n)*MB +i,
						betaB_addr + MB*blk_k*(k%blk_n ) + int(k/blk_n)*MB +i, (k-1)*MB+i, mode, outen, TrainMode);

			}
		}
	}else{
		for(int k=0;k<len;k++){
			for(int i=0;i<MB;i++){
#pragma HLS UNROLL
				betabf[i][k] = bb[k+betaB_addr].range(i*SEXP+SEXP-1, i*SEXP);
			}
		}
	}

}

void ShiftArray2(int_t feeder[MAX_SIZE], MemPack MEMInA, MemPack MEMInB, ap_uint<1> DRELU[MAX_SIZE],
		unsigned idx, ap_uint<1> mode, ap_uint<2> TrainMode, int layer, ap_uint<1> outen){
#pragma HLS INLINE

	static int_t A2[MAX_SIZE][MAX_SIZE];//in
#pragma HLS ARRAY_PARTITION variable=A2 dim=0 type=complete
	static int_t B2[MAX_SIZE][MAX_SIZE];//out
#pragma HLS ARRAY_PARTITION variable=B2 dim=0 type=complete
	int_t P[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=P dim=0 complete


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
					P[j] = MEMInA.range(j*W+W-1,j*W);
					A2[MAX_SIZE-1][j] = DReLu(P[j], DRELU[j], layer);
				}else{//FW, W.T
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

}

void feederA(MemPack a,
			 MemPack b,
			 int_t localA[MAX_SIZE],
			 ap_uint<2> TrainMode,
			 ap_uint<1> shiftarray_mode, ap_uint<1> shiftarray_outen, unsigned idx,
			 ap_uint<1> drelu[MAX_SIZE], int block, int layer
){


	if(TrainMode==0b10){//GD
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


	if(TrainMode == 0b00 || TrainMode == 0b10){
		ShiftArray2( localB, a, b, drelu, idx, shiftarray_mode, TrainMode, layer, shiftarray_outen );
	}else{
		load_b:for(int jj=0;jj<MAX_SIZE;jj++){
			#pragma HLS UNROLL
			localB[jj] = int_t(b.range(jj*W+W-1,jj*W));
		}
	}

}


void TileIdx(ap_uint<2> TrainMode, int ir, int ic, int ik, int size_m, int size_n, int size_k,
		int TILE_R, int TILE_C, int TILE_K, int &aidx, int &bidx, int k
){
	if(TrainMode == 0b10 && k>=0){//GD, transposeA, transposeB
		if(ic == TILE_C-1 && ik == TILE_K-1){
			aidx = ((ir+1)&(TILE_R-1))*MAX_SIZE;
		}else
			aidx = size_m*((ik+1)&(TILE_K-1)) + ir*MAX_SIZE;

		if(ik == TILE_K-1)
			bidx = ((ic+1)&(TILE_C-1))*MAX_SIZE;
		else
			bidx = size_n*((ik+1)&(TILE_K-1)) + ic*MAX_SIZE;
	}else if(TrainMode == 0b00  && k>=0){//FW, transposeB
		if(ik == TILE_K-1){
			bidx = ((ic+1)&(TILE_C-1))*MAX_SIZE;
		}else
			bidx = size_n*((ik+1)&(TILE_K-1)) + ic*MAX_SIZE;

		aidx=0;
	}else{
		aidx=0;bidx=0;
	}
}
