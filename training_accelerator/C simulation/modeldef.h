#ifndef MODELDEF_
#define MODELDEF_
#include "typedef.h"

const unsigned LAYER_FC = 4;//FC stack
const unsigned LAYER_THETA = 2;//backcast, forecast stack
const unsigned LAYER_NUM = LAYER_FC + 2*LAYER_THETA;

//model size parameter
const unsigned NBEATS_BLK = 2;
//layer size parameter
const unsigned B = 8;//batch size, could be divided by MAX_SIZE
const unsigned LK = 8;//LK length

const unsigned INPUT = 5;//input length, M4 yearly
const unsigned OUTPUT = 3;//output length
const unsigned BSIZE = 6;//btach size without padding
const unsigned THA = INPUT + OUTPUT;

const unsigned LOOKBACK = ( int(INPUT/MAX_SIZE) + 1 )*MAX_SIZE;
const unsigned FORECAST = 2*MAX_SIZE;
const unsigned THETA = THA;


const unsigned LOG_LOOKBACK = 3;//log2(LOOKBACK)
const unsigned LOG_FORECAST = 3;//log2(FORECAST)/
const unsigned LOG_LK = 3;
const unsigned LOG_THETA = 3;
const unsigned LOG_B = 3;

const unsigned LOG_LB_MS = 1;//log2(LOOKBACK/MAX_SIZE)
const unsigned LOG_LK_MS = 1;//log2(LK/MAX_SIZE)
const unsigned LOG_THETA_MS = 1;//log2(THETA/MAX_SIZE)
const unsigned LOG_B_MS = 1;//log2(B/MAX_SIZE)
const unsigned LOG_FORE_MS = 1;//log2(FORECAST/MAX_SIZE)

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

#define DATA_NUM 16

#endif
