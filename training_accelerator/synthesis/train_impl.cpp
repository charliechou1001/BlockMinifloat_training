#include "train.h"
#include "modeldef.h"

//void train_impl(MemPack *WA, X_T *X, SExpPack *betaX){
//#pragma HLS INTERFACE m_axi port=WA offset=slave bundle=gmem
//#pragma HLS INTERFACE m_axi port=X offset=slave bundle=gmem1
//#pragma HLS INTERFACE m_axi port=betaX offset=slave bundle=gmem2
//
//	int X_addr;
//	int betaX_addr;
//
//	epoch:for(int i=0;i<EPOCH;i++){
//		Train(WA, X, X_addr, betaX, betaX_addr);
//		X_addr += B*(LOOKBACK + FORECAST)/BANDWIDTH;
//		betaX_addr += B/MAX_SIZE*(LOOKBACK+FORECAST)/MAX_SIZE;
//	}
//
//}

int impl(){
	XPack In[EPOCH * B*(LOOKBACK + FORECAST)/BANDWIDTH];
	MemPack Weight[(ACT_NUM + W_NUM)/MAX_SIZE];
	SExpPack betaX[EPOCH * B/MAX_SIZE*(LOOKBACK+FORECAST)/MAX_SIZE];

	Train(Weight, In, betaX);

	return 0;
}
