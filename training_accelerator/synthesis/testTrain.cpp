
#include "train.h"
#include "modeldef.h"
#include "ap_int.h"

//int main(){
//	X_T In[B*(LOOKBACK + FORECAST)/BANDWIDTH];
//	ap_uint<Ws> XX[B][LOOKBACK+FORECAST];
//	SExpPack betaX[B/MAX_SIZE*(LOOKBACK+FORECAST)/MAX_SIZE];
//	ap_uint<W> wt[ACT_NUM + W_NUM];
//
//	for(int i=0;i<B*(LOOKBACK)/BANDWIDTH;i++){
//		for(int j=0;j<BANDWIDTH;j++){
////			printf("idx: %d, %d\t",j+ int(i/LOOKBACK)*BANDWIDTH,i%LOOKBACK);
//			In[i].range(j*Ws+Ws-1,j*Ws) = XX[j+ int(i/LOOKBACK)*BANDWIDTH][i%LOOKBACK];
//		}
//
//	}
//	for(int i=0;i<B*FORECAST/BANDWIDTH;i++){
//		for(int j=0;j<BANDWIDTH;j++){
//			In[i+B*LOOKBACK/BANDWIDTH].range(j*Ws+Ws-1,j*Ws) = XX[j+ int(i/FORECAST)*BANDWIDTH][i%FORECAST+LOOKBACK];
////			printf("%d,  %d, XX: %s, \n",j+ int(i/FORECAST)*BANDWIDTH,i%FORECAST+LOOKBACK, XX[j+ int(i/FORECAST)*BANDWIDTH][i%FORECAST+LOOKBACK].to_string(2).c_str());
//		}
//	}
////	for(int i=0;i<B*(LOOKBACK+FORECAST)/BANDWIDTH;i++){
////		printf("In[%d]: %s\n",i, In[i].to_string(2).c_str());
////	}
//
//	MemPack Weight[(ACT_NUM + W_NUM)/MAX_SIZE];
////	printf("Weight\n");
////	for(int i=0;i<(ACT_NUM + W_NUM)/MAX_SIZE;i++){
////		Weight[i]=MemPack(0);
////	}
//
////	int size;
////	for(int i=0;i<16;i++){//two nbeats block have 16 layers
////		if ((i+1)%8==0)//layer 7
////			size = FORECAST;
////		else
////			size = LK;
////
//////		printf("idx: %d\n",i);
////		for(int j=0;j<LK/MAX_SIZE;j++){
////			for(int k=0;k<size;k++){
////				for(int ii=0;ii<MAX_SIZE;ii++){
////					Weight[(i - int(i/8))*LK*LK/MAX_SIZE + int(i/8)*FORECAST*LK/MAX_SIZE + k + size*j ].range(ii*W+W-1,ii*W)
////							= wt[i*LK*LK + LK*k + MAX_SIZE*j + ii];
////				}
//////				printf("Weight idx: %d, wt idx: %d\n",(i - int(i/8))*LK*LK/MAX_SIZE + int(i/8)*FORECAST*LK/MAX_SIZE + k + size*j,i*LK*LK + LK*k + MAX_SIZE*j);
////			}
////		}
////	}
//
////
////	printf("Weight\n");
////	for(int i=0;i<(ACT_NUM + W_NUM)/MAX_SIZE;i++){
////		printf("%s\n",Weight[i].to_string(2).c_str());
////	}
//
////	SEXP_T betaIn[EPOCH*BATCH_NUM][INSIZE/MAX_SIZE];
//	int X_addr = 0;
//	int betaIn_addr=0;
//
//
//	Train(Weight, In, X_addr, betaX, 0);
//}

int main(){
//	printf("in:%n, w:%d\n",EPOCH * B*(LOOKBACK + FORECAST)/MAX_SIZE, (ACT_NUM + W_NUM)/MAX_SIZE);
	XPack In[EPOCH * B*(LOOKBACK + FORECAST)/MAX_SIZE];
	MemPack Weight[(ACT_NUM + W_NUM)/MAX_SIZE];
	SExpPack betaX[EPOCH * B/MAX_SIZE*(LOOKBACK+FORECAST)/MAX_SIZE];


	Train(Weight, In, betaX);

	return 0;
}
