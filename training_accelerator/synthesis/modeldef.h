#ifndef MODELDEF_
#define MODELDEF_
#include "typedef.h"

const unsigned LAYER_FC = 4;//FC stack
const unsigned LAYER_THETA = 2;//backcast, forecast stack
const unsigned LAYER_NUM = LAYER_FC + 2*LAYER_THETA;

//model size parameter
const unsigned NBEATS_BLK = 2;
//layer size parameter
const unsigned B = 128;//batch size
const unsigned LK = 64;//LK length

const unsigned INPUT = 12;//input length, M4 yearly
const unsigned OUTPUT = 6;//output length
const unsigned THA = INPUT + OUTPUT;
//#if ((INPUT % MAX_SIZE) !=0)//padding
	const unsigned LOOKBACK = ( int(INPUT/MAX_SIZE) + 1 )*MAX_SIZE;
//#else
//	const unsigned LOOKBACK = INPUT;
//#endif
//#if (OUTPUT % MAX_SIZE) !=0//padding
	const unsigned FORECAST = ( int(OUTPUT/MAX_SIZE) + 1 )*MAX_SIZE;
//#else
//	const unsigned FORECAST = OUTPUT;
//#endif
//#if (THA % MAX_SIZE) != 0
//	const unsigned THETA = ( int(THA/MAX_SIZE) + 1 )*MAX_SIZE;
//#else
	const unsigned THETA = THA;
//#endif

const unsigned LOG_LOOKBACK = 5;//log2(LOOKBACK)
const unsigned LOG_LK = 9;
const unsigned LOG_THETA = 5;
const unsigned LOG_B = 10;

const unsigned LOG_LB_MS = 1;//log2(LOOKBACK/MAX_SIZE)
const unsigned LOG_LK_MS = 4;//log2(LK/MAX_SIZE)
const unsigned LOG_THETA_MS = 1;//log2(THETA/MAX_SIZE)
const unsigned LOG_B_MS = 5;//log2(B/MAX_SIZE)

const unsigned DRELU_LAYER_NUM = (4*B*LK + 2*B*THETA);//total de_relu number
const unsigned W_block = LK*LOOKBACK + 3*LK*LK + 2*THETA*LK + THETA*LOOKBACK + THETA*FORECAST;
const unsigned Act_block = B*LOOKBACK + 4*B*LK + 2*B*THETA;
const unsigned RELU_block = 4*B*LK + 2*B*THETA;
const unsigned W_NUM = NBEATS_BLK * W_block;//total weight number
const unsigned ACT_NUM = NBEATS_BLK * (B*LOOKBACK + 4*B*LK + 2*B*THETA);

const unsigned BFSIZE = B*LK / MAX_SIZE;
const unsigned WSIZE = LK*LK / MAX_SIZE;
const unsigned INSIZE = B*LOOKBACK / MAX_SIZE;
const unsigned OTSIZE = B*FORECAST / MAX_SIZE;

#endif
