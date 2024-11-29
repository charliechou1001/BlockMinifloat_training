#ifndef PARACONF_
#define PARACONF_

#include "modeldef.h"

void WAReadConfig(int layer, ap_uint<2> TrainMode, int &len){
#pragma HLS INLINE
	if(TrainMode==0b00 || TrainMode==0b01){//W, FW&ER
		switch((layer+1)){
		case 1:
			len = LOOKBACK * LK /MAX_SIZE;
			break;
		case 2:
		case 3:
		case 4:
			len = LK * LK /MAX_SIZE;
			break;
		case 5:
		case 7:
			len = LK * THETA /MAX_SIZE;
			break;
		case 6:
			len = THETA * LOOKBACK /MAX_SIZE;
			break;
		case 8:
			len = THETA * FORECAST /MAX_SIZE;
			break;
		}
	}else{//GD
		switch((layer+1)){
		case 1:
			len = B * LOOKBACK /MAX_SIZE;//block input
			break;
		case 2:
		case 3:
		case 4:
			len = B * LK /MAX_SIZE;//fc1-3
			break;
		case 5:
		case 7:
			len = B * LK /MAX_SIZE;//fc4
			break;
		case 6:
//			len = B * LOOKBACK /MAX_SIZE;
			len = B * THETA /MAX_SIZE; //b1
			break;
		case 8:
//			len = B * FORECAST /MAX_SIZE;
			len = B * THETA /MAX_SIZE;//f1
			break;
		}
	}
}

void GEMMConfig(int layer, ap_uint<2> TrainMode, int &size_m, int &size_n, int &size_k, unsigned &log_sizek, unsigned &log_sizek_tilec){

	if(TrainMode==0b00){//FW
		switch((layer+1)){
		case 1:
			size_m = B; size_k = LOOKBACK; size_n = LK;
			log_sizek = LOG_LOOKBACK; log_sizek_tilec = LOG_LOOKBACK + LOG_LK_MS;
			break;
		case 2:
		case 3:
		case 4:
			size_m = B; size_k = LK; size_n = LK;
			log_sizek = LOG_LK; log_sizek_tilec = LOG_LK + LOG_LK_MS;
			break;
		case 5:
		case 7:
			size_m = B; size_k = LK; size_n = THETA;
			log_sizek = LOG_LK; log_sizek_tilec = LOG_LK + LOG_THETA_MS;
			break;
		case 6:
			size_m = B; size_k = THETA; size_n = LOOKBACK;
			log_sizek = LOG_THETA; log_sizek_tilec = LOG_THETA + LOG_LB_MS;
			break;
		case 8:
			size_m = B; size_k = THETA; size_n = FORECAST;
			log_sizek = LOG_THETA; log_sizek_tilec = LOG_THETA + LOG_FORE_MS;
			break;
		}
	}else if(TrainMode==0b01){//ER

		switch((layer+1)){
		case 1:
			size_m = B; size_k = LK; size_n = LOOKBACK;
			log_sizek = LOG_LK; log_sizek_tilec = LOG_LK + LOG_LB_MS;
			break;
		case 2:
		case 3:
		case 4:
			size_m = B; size_k = LK; size_n = LK;
			log_sizek = LOG_LK; log_sizek_tilec = LOG_LK + LOG_LK_MS;
			break;
		case 5:
		case 7:
			size_m = B; size_k = THETA; size_n = LK;
			log_sizek = LOG_THETA; log_sizek_tilec = LOG_THETA + LOG_LK_MS;
			break;
		case 6:
			size_m = B; size_k = LOOKBACK; size_n = THETA;
			log_sizek = LOG_LOOKBACK; log_sizek_tilec = LOG_LOOKBACK + LOG_THETA_MS;
			break;
		case 8:
			size_m = B; size_k = FORECAST; size_n = THETA;
			log_sizek = LOG_FORECAST; log_sizek_tilec = LOG_FORECAST + LOG_THETA_MS;
			break;
		}
	}else{//GD
		switch((layer+1)){
		case 1:
			size_n = LK; size_k = B; size_m = LOOKBACK;
			log_sizek = LOG_B; log_sizek_tilec = LOG_B + LOG_LK_MS;
			break;
		case 2:
		case 3:
		case 4:
			size_n = LK; size_k = B; size_m = LK;
			log_sizek = LOG_B; log_sizek_tilec = LOG_B + LOG_LK_MS;
			break;
		case 5:
		case 7:
			size_n = THETA; size_k = B; size_m = LK;
			log_sizek = LOG_B; log_sizek_tilec = LOG_B + LOG_THETA_MS;
			break;
		case 6:
			size_n = LOOKBACK; size_k = B; size_m = THETA;
			log_sizek = LOG_B; log_sizek_tilec = LOG_B + LOG_LB_MS;
			break;
		case 8:
			size_n = FORECAST; size_k = B; size_m = THETA;
			log_sizek = LOG_B; log_sizek_tilec = LOG_B + LOG_FORE_MS;

			break;
		}
	}
}

void AddrConfig(int layer, int block, int &W_addr, int &Act_addr, int &Drelu_addr,
		int &betaAct_addr, int &betaW_addr){

	//weight
	unsigned Wlayer;
	switch((layer+1)){
		case 1:
			Wlayer = 0;
			break;
		case 2:
			Wlayer = LK*LOOKBACK;
			break;
		case 3:
			Wlayer = LK*LOOKBACK + LK*LK;
			break;
		case 4:
			Wlayer = LK*LOOKBACK + 2*LK*LK;
			break;
		case 5:
			Wlayer = LK*LOOKBACK + 3*LK*LK;
			break;
		case 6:
			Wlayer = LK*LOOKBACK + 3*LK*LK + THETA*LK;
			break;
		case 7:
			Wlayer = LK*LOOKBACK + 3*LK*LK + THETA*LK + THETA*LOOKBACK;
			break;
		case 8:
			Wlayer = LK*LOOKBACK + 3*LK*LK + 2*THETA*LK + THETA*LOOKBACK;
			break;
	}
	W_addr = (W_block * block + Wlayer) / MAX_SIZE;
	betaW_addr = Wlayer / (MAX_SIZE*BLK_SIZE);

	//activation
	unsigned Actlayer;
	switch((layer+1)){
		case 1:
			Actlayer = 0;
			break;
		case 2:
			Actlayer = B*LOOKBACK;
			break;
		case 3:
			Actlayer = B*LOOKBACK + B*LK;
			break;
		case 4:
			Actlayer = B*LOOKBACK + 2*B*LK;
			break;
		case 5:
		case 7:
			Actlayer = B*LOOKBACK + 3*B*LK;
			break;
		case 6:
			Actlayer = B*LOOKBACK + 4*B*LK ;
			break;
		case 8:
			Actlayer = B*LOOKBACK + 4*B*LK + B*THETA;
			break;
	}
	Act_addr = (Act_block * block + Actlayer) / MAX_SIZE;
	betaAct_addr =  W_block / (MAX_SIZE*BLK_SIZE) + Actlayer / (MAX_SIZE*BLK_SIZE);//W_block / (MAX_SIZE*MAX_SIZE): offset for weight beta

	unsigned dRELUlayer;
	switch((layer+1)){
		case 1:
			dRELUlayer = 0;
			break;
		case 2:
			dRELUlayer = B*LK;
			break;
		case 3:
			dRELUlayer = 2*B*LK;
			break;
		case 4:
			dRELUlayer = 3*B*LK;
			break;
		case 5:
			dRELUlayer = 4*B*LK;
			break;
		case 6:
		case 7:///
		case 8:
			dRELUlayer = 4*B*LK + B*THETA;
			break;
	}

	Drelu_addr =  dRELUlayer / MAX_SIZE;
}



#endif








