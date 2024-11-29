#ifndef TRAINDEF_
#define TRAINDEF_

#include "fc.h"
#include "layer.h"
#include "Vadd.h"
#include "MAPELoss.h"
#include "SGD.h"
#include "er_op.h"
//#include "Print.h"
//#include "data.h"

void Train(MemPack *WA, MemPack *Act, XPack *X, SExpPack *betaX){
#pragma HLS INTERFACE m_axi port=WA offset=slave bundle=WA_gmem
#pragma HLS INTERFACE m_axi port=Act offset=slave bundle=Act_gmem
#pragma HLS INTERFACE m_axi port=X offset=slave bundle=X1
#pragma HLS INTERFACE m_axi port=betaX offset=slave bundle=gmem2


#pragma HLS allocation function instances=WriteWeight limit=1
#pragma HLS allocation function instances=WriteGDbuffer limit=1
#pragma HLS allocation function instances=WriteACT limit=1

	float sum_loss;

	int X_addr=0;
	int betaX_addr=0;

	//	printf("train start\n");
	int W_addr, Act_addr, Drelu_addr, A_offset;
	int betaAct_addr, betaW_addr, betaIM_addr;
	int Act_offset = W_NUM / MAX_SIZE;

	MemPack A[BFSIZE*2], LP[BFSIZE];
	MemPack GD_BF[LK*LK/MAX_SIZE];
	MemPack  Wbf[LK*LK/MAX_SIZE], Actbf[BFSIZE];
	XPack HP[BFSIZE];
	static XPack IM[INSIZE + OTSIZE];//addr 0: input and backcast, addr INSIZE: label, addr INSIZE+OTSIZE: forecast

	SExpPack betaA[2][BFSIZE/BLK_SIZE], betaLP[INSIZE/BLK_SIZE], betaHP[INSIZE/BLK_SIZE];
	SExpPack betaIM[(INSIZE + OTSIZE)/BLK_SIZE];
	SExpPack betaWA[NBEATS_BLK][(Act_block + W_block)/(MAX_SIZE*BLK_SIZE)];
	SExpPack betaGD[LK*LK/(MAX_SIZE*BLK_SIZE)];

	for(int ii=0; ii<NBEATS_BLK;ii++){
		for(int jj=0; jj<W_block/(MAX_SIZE*BLK_SIZE); jj++ ){
			betaWA[ii][jj] = betaX[jj+ W_block/(MAX_SIZE*BLK_SIZE)*ii + DATA_NUM * B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)];
		}
	}

	ap_uint<MB> shift[(INSIZE+OTSIZE)/BLK_SIZE];
	int betaAct_offset = 8*NBEATS_BLK*MB;
	int Wlen;

	epoch:for(int ep=0;ep<EPOCH;ep++){

		printf("epoch: %d, iter: %d, label: %d\t ", int(ep/16), ep, ep%DATA_NUM );


		ini_sft:for(int i=0;i<(INSIZE+OTSIZE)/BLK_SIZE;i++){
			shift[i] = 0;
		}

		for(int i=0;i<(INSIZE+OTSIZE)/BLK_SIZE;i++){
			for(int j=0;j<MB;j++){
				betaIM[i].range(SEXP*j+SEXP-1,SEXP*j) = SEXP_T(-90);
			}
		}

		ini_im:for(int i=0;i<(INSIZE+OTSIZE);i++){
			IM[i] = XPack(0);
		}

		betaIM_addr = 0;
		X_addr = (ep % DATA_NUM)*B*(LOOKBACK + FORECAST)/MAX_SIZE;

		betaX_addr = X_addr/BLK_SIZE;

		Load(X, IM, betaX, betaIM, X_addr, 0, betaX_addr, betaIM_addr);

		ap_uint<2> TrainMode;int layer;
	//	Forward
		TrainMode = 0b00;
		FW:for(int blk=0; blk < NBEATS_BLK; blk++){//nbeats blocks

			A_offset = 0;
			layer = 0;
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			ResidualToBM(IM,A, Act, betaIM, betaA[0], betaWA[blk],B*LOOKBACK/MAX_SIZE,0, A_offset, Act_addr,betaAct_addr,TrainMode);

			FC4:for(int jj=0;jj < LAYER_FC; jj++){

				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

				layer ++;
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WriteBack(A,LP,Act,betaA[A_offset/BFSIZE], betaLP, betaWA[blk], layer-1, A_offset, Act_addr,betaAct_addr);

			}

			if(blk != NBEATS_BLK -1){

				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);

				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				A_offset = BFSIZE;


				layer ++;
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				WriteBack(A,LP,Act,betaA[A_offset/BFSIZE], betaLP, betaWA[blk], layer-1, A_offset, Act_addr,betaAct_addr);
				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

				//backcast residual
				VecADD(IM, HP, betaIM, betaHP, shift, B, LOOKBACK, 0, ap_uint<1>(1),TrainMode);

			}

			layer = 6;
			A_offset = 0;
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

			layer ++;
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			WriteBack(A,LP,Act,betaA[A_offset/BFSIZE], betaLP, betaWA[blk], layer-1, A_offset, Act_addr,betaAct_addr);
			LoadWeight(WA, Wbf, layer, W_addr, TrainMode);

			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

			//forecast sum
			VecADD(IM, HP, betaIM, betaHP, shift, B, FORECAST, INSIZE, ap_uint<1>(0),TrainMode);
		}


		X_addr = B*LOOKBACK/MAX_SIZE + (ep%DATA_NUM )*B*(LOOKBACK + FORECAST)/MAX_SIZE;
		betaX_addr = X_addr/BLK_SIZE;

		float res = MAPEPrint(IM, X, betaIM, betaX, B*LOOKBACK/MAX_SIZE,X_addr, betaX_addr );
		MAPELoss(IM,X,betaIM,betaX,B*LOOKBACK/MAX_SIZE,X_addr, betaX_addr);


		if( (ep+1) % 16 ==0 ){
			sum_loss += res;
			printf("===average loss: %f\n", sum_loss/16 );
			sum_loss = 0;
		}else
			sum_loss += res;


		BP:for(int blk= NBEATS_BLK-1;blk >=0;blk--){
			A_offset = 0;
			ResidualToBM(IM,A, Act, betaIM, betaA[0], betaWA[blk], B*FORECAST/MAX_SIZE,B*LOOKBACK/MAX_SIZE,A_offset,Act_addr,betaAct_addr,0b10 );


			//forecast 2
			layer = 7;
			TrainMode = 0b10;//GD
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			LoadWeight(Act, Actbf, layer, Act_addr, TrainMode);


			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
			WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

			TrainMode =0b01;//ER
			LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);


			SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer, betaW_addr);//weight update
			WriteWeight(WA, Wbf, layer, W_addr, TrainMode);

			layer--;
			WriteBackBP(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

			//forecast 1
			TrainMode = 0b10;//GD
			AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
			LoadWeight(Act, Actbf, layer, Act_addr, TrainMode);
			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
			WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

			TrainMode =0b01;//ER
			LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
			gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);

			SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer, betaW_addr);//weight update
			WriteWeight(WA, Wbf, layer, W_addr, TrainMode);

			layer=3;
			WriteBackBP(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

			if(blk != NBEATS_BLK -1){
				A_offset = BFSIZE;
				ResidualToBM(IM,A, Act, betaIM, betaA[1], betaWA[blk], B*LOOKBACK/MAX_SIZE,0,A_offset,Act_addr,betaAct_addr,TrainMode );
				NEG( A, B, LOOKBACK, A_offset );//negative

				//backcast 2
				layer = 5;
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				LoadWeight(Act, Actbf, layer, Act_addr, TrainMode);

				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer, betaW_addr);//weight update
				WriteWeight(WA, Wbf, layer, W_addr, TrainMode);

				layer--;
				WriteBackBP(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

				//backcast 1
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				LoadWeight(Act, Actbf, layer, Act_addr, TrainMode);
				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer, betaW_addr);//weight update
				WriteWeight(WA, Wbf, layer, W_addr, TrainMode);

				layer--;
				WriteBackBP(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

				//branch error addition
				ER_OP( A, betaA[0], betaA[1], B, LK, 0,BFSIZE );//addition
			}

			layer = 3;
			A_offset = 0;
			BPFC4:for(int jj = LAYER_FC-1; jj >= 0; jj--){
				TrainMode = 0b10;//GD
				AddrConfig(layer, blk, W_addr, Act_addr, Drelu_addr, betaAct_addr, betaW_addr);
				LoadWeight(Act, Actbf, layer, Act_addr, TrainMode);

				gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaAct_addr);
				WriteGDbuffer(LP, GD_BF, betaLP, betaGD, layer);

				TrainMode =0b01;//ER
				LoadWeight(WA, Wbf, layer, W_addr, TrainMode);
				if(layer !=0 || blk != 0){
					gemm(A,Wbf,Actbf,LP,HP, betaA[A_offset/BFSIZE], betaWA[blk],betaLP,betaHP,blk, layer, TrainMode, Drelu_addr, A_offset, betaW_addr);
					WriteBackBP(A,LP,betaA[A_offset/BFSIZE], betaLP, layer, A_offset);

				}
				SGD(Wbf, GD_BF, betaWA[blk], betaGD, layer, betaW_addr);//weight update
				WriteWeight(WA, Wbf, layer, W_addr, TrainMode);
				layer--;
			}

			if(blk <NBEATS_BLK -1 && blk!=0){//backcast residual
				TrainMode =0b01;//ER
				ap_uint<1> op = blk == NBEATS_BLK -1? ap_uint<1>(0) : ap_uint<1>(1);

				VecADD(IM, HP, betaIM, betaHP, shift, B, LOOKBACK, 0, op,TrainMode);
			}else if(blk == NBEATS_BLK -1){
				WriteBKRes(IM, HP, betaIM, betaHP);

			}

		}

	}
}

#endif

