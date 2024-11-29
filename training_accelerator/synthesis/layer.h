#ifndef LAYERDEF_
#define LAYERDEF_

#include "ap_int.h"
#include "typedef.h"
//#include "Cvt.h"
#include "MAC.h"
//#include "test_longBMtofp32.h"


//void Load(X_T in[], XPack IM[], SExpPack *betaIn, SExpPack *betaIM, int in_addr, int im_addr, int betaIn_addr, int betaIM_addr){
//	load:for(int t=0;t < B*LOOKBACK/MAX_SIZE; t++){
//		XPack imdata;
//		X_T d_in;
//		for(int j=0;j < MAX_SIZE/BANDWIDTH;j++){
//
//		d_in = in[t%LOOKBACK+ int(t/LOOKBACK)*LOOKBACK*(MAX_SIZE/BANDWIDTH) + j*LOOKBACK + in_addr];
//		imdata.range(j*BANDWIDTH*Ws + BANDWIDTH*Ws -1,j*BANDWIDTH*Ws) = d_in;
//		}
//		IM[t+ im_addr] = imdata;
//	}
//
//	for(int k=0; k < B*LOOKBACK/(BLK_SIZE*MAX_SIZE);k++){//beta
//		betaIM[k+betaIM_addr] = betaIn[k+betaIn_addr];
//	}
//}

void Load(XPack in[], XPack IM[], SExpPack *betaIn, SExpPack *betaIM, int in_addr, int im_addr, int betaIn_addr, int betaIM_addr){
#pragma HLS INLINE
	load:for(int t=0;t < B*LOOKBACK/MAX_SIZE; t++){
		IM[t+ im_addr] = in[t + in_addr];
	}

	for(int k=0; k < B*LOOKBACK/(BLK_SIZE*MAX_SIZE);k++){//beta
		betaIM[k+betaIM_addr] = betaIn[k+betaIn_addr];
	}
}

//void Load(X_T in[], XPack IM[], SExpPack *betaIn, SExpPack *betaIM, int in_addr, int im_addr, int betaIn_addr, int betaIM_addr){
//#pragma HLS DATAFLOW
//	X_T d_in[B*LOOKBACK/BANDWIDTH];
//
//	load:for(int t=0;t < B*LOOKBACK/BANDWIDTH; t++){
//		d_in[t] = in[t + in_addr];
//		}
//	load_cvt:for(int t=0;t < B*LOOKBACK/MAX_SIZE; t++){
//		XPack imdata;
//		for(int j=0;j < MAX_SIZE/BANDWIDTH;j++){
//
//		imdata.range(j*BANDWIDTH*Ws + BANDWIDTH*Ws -1,j*BANDWIDTH*Ws) =
//				d_in[t%LOOKBACK+ int(t/LOOKBACK)*LOOKBACK*(MAX_SIZE/BANDWIDTH) + j*LOOKBACK];
//		}
//		IM[t+ im_addr] = imdata;
//	}
//
//
////	for(int k=0; k < B*(LOOKBACK+FORECAST)/(MAX_SIZE*MAX_SIZE);k++){//beta
//	load_betaIM:for(int k=0; k < B*LOOKBACK/(BLK_SIZE*MAX_SIZE);k++){//beta
//		betaIM[k+betaIM_addr] = betaIn[k+betaIn_addr];
////		printf("betaIn: %d, betaIM:%d\n",betaIn[k+betaIn_addr].to_int(), betaIM[k+betaIM_addr].to_int());
//	}
//}

void LoadWeight(MemPack *WA, MemPack Wbf[], int W_addr, int len){

	load_weight: for(int t=0;t<len/MAX_SIZE;t++){
		Wbf[t] = WA[t+W_addr];
	}
}

void WriteWeight(MemPack *WA, MemPack Wbf[], int W_addr, int len){
	write_weight: for(int t=0;t<len/MAX_SIZE;t++){
		WA[t+W_addr] = Wbf[t];
	}
}

void ResidualToBM(XPack *res, MemPack *In, SExpPack *betaIM, SExpPack *betaA, unsigned size, unsigned offset, unsigned A_offset){
#pragma HLS DATAFLOW
	res_to_bm:for(int t=0;t < size; t++){
		MemPack Indata;
		for(int j=0;j<MAX_SIZE;j++){
#pragma HLS UNROLL
			ap_uint<Ws> resdata = res[t+offset].range(j*Ws+Ws-1,j*Ws);
//				printf("res_idx: %d,  %s\n",t+offset, resdata.to_string(2).c_str());
//			int_t Rdata;
//			Normalization<W,W,4>(resdata, ap_uint<W>(0), ap_uint<2>(0), ap_uint<2>(0), Rdata);
//				printf("resdata: %s, Rdata: %s\n",resdata.to_string(2).c_str(), Rdata.to_string(2).c_str());
			Indata.range( j*W+W-1, j*W ) =  resdata.range(Ws-1,Ws-W);
//				printf("i*BANDWIDTH+j: %d\n",i*BANDWIDTH+j);
		}
		In[t+A_offset] = Indata;
//		printf("Indata: %s, In: %s\n",Indata.to_string(2).c_str(), In[t].to_string(2).c_str());
	}

	res_to_bm_beta:for(int i=0;i<size/BLK_SIZE;i++){
		betaA[i+A_offset/(MAX_SIZE*BLK_SIZE)] = betaIM[i+offset/(MAX_SIZE*BLK_SIZE)];
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
//	printf("size: %d\n",size);
//	printf("WRITE_ACT: \n");
	wr_act:for(int t=0; t < size; t++){//data
		WA[t+ Act_addr] = A[t + A_offset];
	}

	wr_act_beta:for(int k=0; k < size/BLK_SIZE;k++){//beta
		betaWA[k + betaAct_addr] = betaA[k];
	}
}

void WriteBack(MemPack *A, MemPack *C, SExpPack *betaA, SExpPack *betaC, int layer, int A_offset){
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
//		case 6:
//		case 8:
//			size = THETA*B;
//			break;
	}
//	if(layer==7)
//		printf("size: %d, THETA: %d, B:%d,  A_offset: %d\n",size,THETA, B, A_offset);

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
//	if(layer==7)
//		printf("size: %d, THETA: %d, B:%d,  A_offset: %d\n",size,THETA, B, A_offset);

	wr_gd:for(int t=0; t < size/MAX_SIZE; t++){//data
		out[t] = in[t];

	}

	wr_gd_beta:for(int k=0;k < size/(MAX_SIZE*BLK_SIZE); k++){//beta
		betaout[k] = betain[k];
	}


}

#endif
