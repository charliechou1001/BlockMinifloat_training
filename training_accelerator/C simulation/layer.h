#ifndef LAYERDEF_
#define LAYERDEF_

#include "ap_int.h"
#include "typedef.h"
//#include "Cvt.h"
#include "MAC.h"
#include "paraconfig.h"

void Load(XPack in[], XPack IM[], SExpPack *betaIn, SExpPack *betaIM, int in_addr, int im_addr, int betaIn_addr, int betaIM_addr){
#pragma HLS INLINE
	load:for(int t=0;t < B*LOOKBACK/MAX_SIZE; t++){
		IM[t+ im_addr] = in[t + in_addr];
	}

	for(int k=0; k < B*LOOKBACK/(BLK_SIZE*MAX_SIZE);k++){//beta
		betaIM[k+betaIM_addr] = betaIn[k+betaIn_addr];
	}
}

void LoadWeight(MemPack *WA, MemPack Wbf[], int layer, int W_addr, ap_uint<2> TrainMode){
	int len;
	WAReadConfig(layer, TrainMode, len);
	load_weight: for(int t=0;t<len;t++){
		Wbf[t] = WA[t+W_addr];
	}
}

void WriteWeight(MemPack *WA, MemPack Wbf[], int layer, int W_addr, ap_uint<2> TrainMode){
	int Wlen;
	WAReadConfig(layer, TrainMode, Wlen);
	write_weight: for(int t=0;t<Wlen;t++){
		WA[t+W_addr] = Wbf[t];
	}
}

void ResidualToBM(XPack *res, MemPack *In, MemPack *Act, SExpPack *betaIM, SExpPack *betaA, SExpPack *betaWA,
		unsigned size, unsigned offset, unsigned A_offset, unsigned Act_addr, unsigned betaAct_addr,ap_uint<2> TrainMode){

	ap_uint<Ws> resdata[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=resdata dim=1 type=complete
	ap_uint<Ws> resdata_2[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=resdata_2 dim=1 type=complete
	XPack res_in, res_in_2;
	MemPack Indata, Indata_2;
	ap_uint<W> BMdata[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=BMdata dim=1 type=complete
	ap_uint<W> BMdata_2[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=BMdata_2 dim=1 type=complete

#pragma HLS DATAFLOW
	res_to_bm:for(int t=0;t < size; t++){
#pragma HLS PIPELINE
		res_in = res[t+offset];
		for(int j=0;j<MAX_SIZE;j++){
#pragma HLS UNROLL
			resdata[j] = res_in.range(j*Ws+Ws-1,j*Ws);
#if W == 4
			if(TrainMode==0b00)
				BMdata[j] =  resdata[j].range(Ws-2,Ws-W-1);//activation, unsigned 4 bit integer
			else
				BMdata[j] =  resdata[j].range(Ws-1,Ws-W);//error, BM<0,3>
#else
			ap_uint<1> guard = resdata[j][Ws-W];
			ap_uint<1> round = resdata[j][Ws-W-1];
			ap_uint<1> sticky = resdata[j](Ws-W-2,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
			ap_uint<1> randn = round & sticky | guard&round&(~sticky);

			BMdata[j] =  resdata[j].range(Ws-1,Ws-W);//activation and error, BM<0,7>
			BMdata[j] = BMdata[j] + randn;
#endif
		}


		for(int j=0;j<MAX_SIZE;j++){
#pragma HLS UNROLL
			Indata.range(j*W+W-1,W*j) = BMdata[j];
		};
		In[t+A_offset] = Indata;
	}

	if(TrainMode == 0b00){
		res_to_act:for(int t=0;t < size; t++){
	#pragma HLS PIPELINE
			res_in_2 = res[t+offset];
			for(int j=0;j<MAX_SIZE;j++){
	#pragma HLS UNROLL
				resdata_2[j] = res_in_2.range(j*Ws+Ws-1,j*Ws);
	#if W == 4
				if(TrainMode==0b00)
					BMdata_2[j] =  resdata_2[j].range(Ws-2,Ws-W-1);//activation, unsigned 4 bit integer
				else
					BMdata_2[j] =  resdata_2[j].range(Ws-1,Ws-W);//error, BM<0,3>
	#else
				BMdata_2[j] =  resdata_2[j].range(Ws-1,Ws-W);//activation and error, BM<0,7>
	#endif
			}
			for(int j=0;j<MAX_SIZE;j++){
	#pragma HLS UNROLL
				Indata_2.range(j*W+W-1,W*j) = BMdata_2[j];
			};
			Act[t+Act_addr] = Indata_2;
		}

	}

	res_to_bm_beta:for(int i=0;i<size/BLK_SIZE;i++){
		betaA[i] = betaIM[i+offset/BLK_SIZE];

		SEXP_T aa = betaIM[i+offset/BLK_SIZE].range(W-1,0);
		SEXP_T bb = betaA[i].range(W-1,0);
 	}


	if(TrainMode==0b00){
		res_to_act_beta:for(int i=0;i<size/BLK_SIZE;i++){
			betaWA[i+betaAct_addr] = betaIM[i+offset/(MAX_SIZE*BLK_SIZE)];
	 	}
	}
}

void WriteACT(MemPack *A, MemPack *WA, SExpPack *betaA, SExpPack *betaWA, int layer, int A_offset,int Act_addr, int betaAct_addr){
#pragma HLS DATAFLOW
	int size;
	switch(layer+1){
		case 1:
			size = LOOKBACK*B /MAX_SIZE;
			break;
		case 2:
		case 3:
		case 4:
		case 5://5 and 7 share the same BP activation
		case 7:
			size = LK*B /MAX_SIZE;
			break;
		case 6:
		case 8:
			size = THETA*B/MAX_SIZE;
			break;
	}

	wr_act:for(int t=0; t < size; t++){//data
		WA[t+ Act_addr] = A[t + A_offset];
	}

	wr_act_beta:for(int k=0; k < size/BLK_SIZE;k++){//beta
		betaWA[k + betaAct_addr] = betaA[k];
	}
}

void WriteBack(MemPack *A, MemPack *C, MemPack *Act, SExpPack *betaA, SExpPack *betaC, SExpPack *betaWA, int layer, int A_offset,int Act_addr, int betaAct_addr){
#pragma HLS DATAFLOW
	int size;
	switch(layer+1){
		case 1:
		case 2:
		case 3:
		case 4:
			size = LK*B;
			break;
		case 5:
		case 7:
			size = THETA*B;
			break;
	}

	wr_back:for(int t=0; t < size/MAX_SIZE; t++){//data
		A[t+ A_offset] = C[t];
		Act[t+ Act_addr] = C[t];////
	}

	wr_back_beta:for(int k=0;k < size/(MAX_SIZE*BLK_SIZE); k++){//beta
		betaA[k] = betaC[k];
		betaWA[k + betaAct_addr] = betaC[k];////
	}
}

void WriteBackBP(MemPack *A, MemPack *C, SExpPack *betaA, SExpPack *betaC, int layer, int A_offset){
#pragma HLS DATAFLOW
	int size;
	switch(layer+1){
		case 1:
		case 2:
		case 3:
		case 4:
			size = LK*B;
			break;
		case 5:
		case 7:
			size = THETA*B;
			break;
	}

	wr_back:for(int t=0; t < size/MAX_SIZE; t++){//data
		A[t+ A_offset] = C[t];
	}

	wr_back_beta:for(int k=0;k < size/(MAX_SIZE*BLK_SIZE); k++){//beta
		betaA[k] = betaC[k];
	}
}

void WriteGDbuffer(MemPack *in, MemPack *out, SExpPack *betain, SExpPack *betaout, int layer){
#pragma HLS DATAFLOW
	int size;
	switch(layer+1){
		case 1:
			size = LOOKBACK*LK;
			break;
		case 2:
		case 3:
		case 4:
			size = LK*LK;
			break;
		case 5:
		case 7:
			size = THETA*LK;
			break;
		case 6:
			size = THETA*LOOKBACK;
			break;
		case 8:
			size = THETA*FORECAST;
			break;
	}

	wr_gd:for(int t=0; t < size/MAX_SIZE; t++){//data
		out[t] = in[t];

	}

	wr_gd_beta:for(int k=0;k < size/(MAX_SIZE*BLK_SIZE); k++){//beta
		betaout[k] = betain[k];
	}


}

void WriteBKRes(XPack *IM, XPack *HP, SExpPack *betaIM, SExpPack *betaHP){

	BKRes:for(int t=0; t < B*LOOKBACK/MAX_SIZE; t++){//data
		IM[t] = HP[t];
	}

	BKRes_beta:for(int k=0;k < B*LOOKBACK/(MAX_SIZE*BLK_SIZE); k++){//beta
		betaIM[k] = betaHP[k];
	}

}


#endif
