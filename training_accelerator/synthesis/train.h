#ifndef TRAINDEF_
#define TRAINDEF_

#include "fc.h"
#include "layer.h"
#include "Vadd.h"
#include "MAPELoss.h"
#include "SGD.h"
#include "er_op.h"
//#include "test_longBMtofp32.h"
//#include "data.h"

//SEXP_T betawt[NBEATS_BLK][W_block/(MAX_SIZE*MAX_SIZE)] =
//{-2., -2.,-3., -2.,-2., -2., -2., -2.,-2., -2., -2., -2.,-2., -2., -2., -2., -2., -2., -2., -2.,-2., -2., -2., -2.,-2., -2., -2., -2.,-2., -2.,
//-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2};


void Train(MemPack *WA, XPack *X, SExpPack *betaX){
#pragma HLS INTERFACE m_axi port=WA offset=slave bundle=WA_gmem
#pragma HLS INTERFACE m_axi port=X offset=slave bundle=X1
#pragma HLS INTERFACE m_axi port=betaX offset=slave bundle=gmem2


#pragma HLS allocation function instances=WriteWeight limit=1
#pragma HLS allocation function instances=WriteGDbuffer limit=1
#pragma HLS allocation function instances=WriteACT limit=1

	int X_addr=0;
	int betaX_addr=0;

	//	printf("train start\n");
	int W_addr, Act_addr, Drelu_addr, A_offset;
	int betaAct_addr, betaW_addr, betaIM_addr;
	int Act_offset = W_NUM / MAX_SIZE;

	MemPack A[BFSIZE*2], LP[BFSIZE], GD_BF[LK*LK/MAX_SIZE], Wbf[BFSIZE];
	XPack HP[BFSIZE];
	static XPack IM[INSIZE + OTSIZE];//addr 0: input and backcast, addr INSIZE: label, addr INSIZE+OTSIZE: forecast

	SExpPack betaA[2][BFSIZE/BLK_SIZE], betaLP[INSIZE/BLK_SIZE], betaHP[INSIZE/BLK_SIZE];
	SExpPack betaIM[(INSIZE + OTSIZE)/BLK_SIZE];
	SExpPack betaWA[NBEATS_BLK][(Act_block + W_block)/(MAX_SIZE*BLK_SIZE)];
	SExpPack betaGD[LK*LK/(MAX_SIZE*BLK_SIZE)];

	ap_uint<MB> shift[(INSIZE+OTSIZE)/BLK_SIZE];
//	SEXP_T betaW[8*NBEATS_BLK][W_block/(MAX_SIZE*MAX_SIZE)];
//	SEXP_T betaAct[15*NBEATS_BLK][(B*LK)/(MAX_SIZE*MAX_SIZE)];
	int betaAct_offset = 8*NBEATS_BLK*MB;
	int Wlen;

	epoch:for(int i=0;i<EPOCH;i++){

		ini_sft:for(int i=0;i<(INSIZE+OTSIZE)/BLK_SIZE;i++){
			shift[i] = 0;
	//		betaIM[i] = SEXP_T(-20);
		}

		ini_im:for(int i=0;i<(INSIZE+OTSIZE);i++){
			IM[i] = XPack(0);
		}

		//	printf("W_NUM: %d, W_NUM/MAX_SIZE:%d\n",W_NUM, W_NUM/MAX_SIZE);
		betaIM_addr = 0;
		Load(X, IM, betaX, betaIM, X_addr, 0, betaX_addr, betaIM_addr);

		ap_uint<2> TrainMode;int layer;
	//	Forward
		TrainMode = 0b00;
		FW:for(int blk=0; blk < NBEATS_BLK; blk++){//nbeats blocks
	//		FW:for(int blk=0; blk < 1; blk++){//nbeats blocks
	//		printf("IM input: \n" );
	//		PrintKadd(INSIZE, 15, IM, betaIM);
			A_offset = 0;
			ResidualToBM(IM,A, betaIM, betaA[0], B*LOOKBACK/MAX_SIZE,0, A_offset);
	//		printf("input: \n");
	//		PrintMatrix(A, betaA[0], INSIZE*MAX_SIZE, 0,0, ap_uint<2>(0));

			layer = 0;

			FC4:for(int jj=0;jj < LAYER_FC; jj++){//4 FC layer
	//			FC4:for(int jj=0;jj < 1; jj++){//4 FC layer
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
	//			printf("[ADDR]: layer:%d, blk: %d, W_addr:%d, Act_addr: %d, Drelu_addr: %d,\n betaAct_addr: %d, betaW_addr:%d\n",
	//											layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WriteACT(A,WA, betaA[A_offset/BFSIZE], betaWA[blk], layer, A_offset, Act_addr,betaAct_addr);

				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);
				layer ++;
			}

			if(blk != NBEATS_BLK -1){

				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				A_offset = BFSIZE;
				WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);
				layer ++;

				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
	//			printf("[ADDR]: layer:%d, blk: %d, W_addr:%d, Act_addr: %d, Drelu_addr: %d,\n betaAct_addr: %d, betaW_addr:%d\n",
	//								layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WriteACT(A,WA, betaA[A_offset/BFSIZE], betaWA[blk], layer, A_offset, Act_addr,betaAct_addr);
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				layer ++;
				//backcast residual
				VecADD(IM, HP, betaIM, betaHP, shift, B, LOOKBACK, 0, ap_uint<1>(1),TrainMode);
			}

			layer = 6;
			A_offset = 0;
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			WriteACT(A,WA, betaA[A_offset/BFSIZE], betaWA[blk], layer, A_offset, Act_addr,betaAct_addr);
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, W_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
			WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);
			layer ++;

			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			WriteACT(A,WA, betaA[A_offset/BFSIZE], betaWA[blk], layer, A_offset, Act_addr,betaAct_addr);
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, W_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

			//forecast sum
			VecADD(IM, HP, betaIM, betaHP, shift, B, FORECAST, INSIZE, ap_uint<1>(0),TrainMode);
		}

		MAPELoss(IM,X,betaIM,betaX,B*LOOKBACK/MAX_SIZE,X_addr, betaX_addr);



		BP:for(int blk= NBEATS_BLK-1;blk >=0;blk--){
			A_offset = 0;
			ResidualToBM(IM,A, betaIM, betaA[0], B*FORECAST/MAX_SIZE,B*LOOKBACK/MAX_SIZE,A_offset);

	//		printf("loss input: \n");
	//		PrintMatrix(A, betaA[0], B*LK, A_offset,A_offset/BFSIZE, ap_uint<2>(0));

			//forecast 2
			layer = 7;
			TrainMode = 0b10;//GD
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, Act_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
	//		printf("gradient:\n");
	//		PrintKadd(THETA*FORECAST/MAX_SIZE, 15, HP, betaHP);
			WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

			TrainMode =0b01;//ER
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, W_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
			SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer);//weight update
			WriteWeight(WA, Wbf, W_addr, Wlen);

			layer--;
			WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);
	//		printf("A: \n");
	//		PrintMatrix(A, betaA[0], B*LK, A_offset,A_offset/BFSIZE, ap_uint<2>(0));

			//forecast 1
			TrainMode = 0b10;//GD
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, Act_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
			WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

			TrainMode =0b01;//ER
			WAReadConfig(layer, TrainMode, Wlen);
			LoadWeight(WA, Wbf, W_addr, Wlen);
			gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
			SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer);//weight update
//			SGD(Wbf, LP, betaWA[blk], betaLP, layer);
			WriteWeight(WA, Wbf, W_addr, Wlen);

			layer=3;
			WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

			if(blk != NBEATS_BLK -1){
				A_offset = BFSIZE;
				ResidualToBM(IM,A, betaIM, betaA[0], B*LOOKBACK/MAX_SIZE,0, A_offset);
				ER_OP( A, betaA[0], betaA[1], B, LOOKBACK, 0,0,ap_uint<1>(0) );//negative

				//backcast 2
				layer = 5;
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, Act_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer);//weight update
				WriteWeight(WA, Wbf, W_addr, Wlen);

				layer--;
				WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

				//backcast 1
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, Act_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer);//weight update
				WriteWeight(WA, Wbf, W_addr, Wlen);

				layer--;
				WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

				//branch error addition
				ER_OP( A, betaA[0], betaA[1], B, LK, 0,BFSIZE,ap_uint<1>(1) );//addition
			}

			layer = 3;
			A_offset = 0;
			BPFC4:for(int jj = LAYER_FC-1; jj >= 0; jj--){
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, Act_addr, Wlen);
				gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				WAReadConfig(layer, TrainMode, Wlen);
				LoadWeight(WA, Wbf, W_addr, Wlen);
				if(layer !=0 || blk != 0){
					gemm(A,Wbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
					WriteBack(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);
				}
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer);//weight update
				WriteWeight(WA, Wbf, W_addr, Wlen);
				layer--;
			}

			if(blk!=0){//backcast residual
				TrainMode =0b01;//ER
				ap_uint<1> op = blk == NBEATS_BLK -1? ap_uint<1>(0) : ap_uint<1>(1);
				VecADD(IM, HP, betaIM, betaHP, shift, B, LOOKBACK, 0, op,TrainMode);
			}

		}

	X_addr += B*(LOOKBACK + FORECAST)/MAX_SIZE;
	betaX_addr += B/MAX_SIZE*(LOOKBACK+FORECAST)/MAX_SIZE;
	}
}

#endif
