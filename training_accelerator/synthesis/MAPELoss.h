//MAPE loss derivation reference:  https://stats.stackexchange.com/questions/316076/gradient-and-hessian-of-the-mape
#ifndef MAPEDEF_
#define MAPEDEF_

#include "typedef.h"
#include "newCLZ.h"
//#include "Cvt.h"
#include "modeldef.h"

void MAPECore(ap_int<Ws> &y, ap_int<Ws> label, SEXP_T beta_y, SEXP_T beta_label){
#pragma HLS INLINE
	const ap_int<Ws> coff = 1/(B*OUTPUT);
	ap_int<Ws> a,b,c;
	ap_int<W> ebias = beta_y - beta_label;
	ap_int<W> Ka = ebias >= 0? ap_int<W>(0) : ap_int<W>(-ebias);
	ap_int<W> Kb = ebias < 0? ap_int<W>(0) : ap_int<W>(ebias);

//	printf("Ka: %d, Kb:%d\n",Ka.to_int(), Kb.to_int());

	a(Ws-1,0) = y(Ws-1,0);
	b(Ws-1,0) = label!=ap_int<Ws>(0)? label(Ws-1,0) : ap_int<Ws>(1);
//	printf("a: %f, b: %f, label:%s\t",float(a), float(b),label.to_string(2).c_str());

	ap_int<Ws> div_res = coff*(1/b);
	c = ( (a >> Ka) - (b >> Kb) )>=0? div_res : ap_int<Ws>(-div_res);

//	c = ( (a >> Ka) - (b >> Kb) )>=0? ap_int<Ws>(1)/(B*OUTPUT*b) : ap_int<Ws>(-1)/(B*OUTPUT*b);

//	printf("c before: %f\n",float(c));
//	c = c >>(fw_frac - er_frac);
//	printf("fw_frac - er_frac: %d, c after: %f\n",fw_frac - er_frac,float(c));
	y(Ws-1,0) = c(Ws-1,0);
}

//void MAPELoss(XPack *Y, X_T *LABEL, SExpPack *betaY, SExpPack *betaL, int y_addr, int X_addr, int betaX_addr){
//	ap_int<Ws> Abuf[HALF_SIZE],Bbuf[HALF_SIZE];
//#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
//#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
//
//	X_T label_bf[B*FORECAST/MAX_SIZE];
//	for(int i=0;i<B*FORECAST/MAX_SIZE;i++ ){
//		label_bf[i] = LABEL[i];
//	}
//
//	unsigned TILE_M = B/MAX_SIZE;
//	unsigned TILE_N = FORECAST/BLK_SIZE;
//
//	mape:for(int m=0;m < TILE_M;m++){
//		for(int n=0;n< TILE_N;n++){
//
//			for(int i=0;i<BLK_SIZE;i++){
//				for(int j=0;j<MAX_SIZE/BANDWIDTH;j++){
////#pragma HLS PIPELINE
//					for(int k=0;k<BANDWIDTH;k++){
//						#pragma HLS UNROLL
//						if(n!=TILE_N-1 && i<BLK_SIZE - (FORECAST-OUTPUT) ){
//	//						printf("===i:%d, j:%d,k:%d\n",i,j,k);
//							Abuf[k] = Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range((k+j*BANDWIDTH)*Ws+Ws-1 ,(k+j*BANDWIDTH)*Ws);
//							Bbuf[k] = label_bf[i*(MAX_SIZE/BANDWIDTH)+j+(m*TILE_N+n)*BLK_SIZE*(MAX_SIZE/BANDWIDTH) + X_addr].range(k*Ws+Ws-1 ,k*Ws);
//	//						printf("A idx: %d, k+j*BANDWIDTH: %d, B idx: %d, k: %d\n",i+(m*TILE_N+n)*MAX_SIZE,k+j*BANDWIDTH,
//	//								i+j*MAX_SIZE+(m*TILE_N+n)*MAX_SIZE*(MAX_SIZE/BANDWIDTH),  k);
////							Bbuf[k] = InToResidual(indata);
//							MAPECore(Abuf[k], Bbuf[k], betaY[m*TILE_N+n + y_addr/BLK_SIZE].range(int((k+j*BANDWIDTH)/BLK_SIZE)*SEXP + SEXP-1,int((k+j*BANDWIDTH)/BLK_SIZE)*SEXP),
//									betaL[m*TILE_N+n+betaX_addr].range(int(k/BLK_SIZE)*SEXP + SEXP-1,int(k/BLK_SIZE)*SEXP));
//							Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range((k+j*BANDWIDTH)*Ws+Ws-1 ,(k+j*BANDWIDTH)*Ws) = Abuf[k];
//						}
//					}
//				}
//			}
//
//			for(int i=0;i<MAX_SIZE/BANDWIDTH;i++){
//				for(int j=0; j<MM;j++){
//#pragma HLS UNROLL
//					betaY[m*TILE_N+n + y_addr/BLK_SIZE].range(SEXP*(i*MM+j)+SEXP-1, SEXP*(i*MM+j))
//							= -betaL[m*TILE_N+n+betaX_addr].range(SEXP*j+SEXP-1, SEXP*j);
//				}
//			}
//
////			printf("betaY: %d, idx: %d\n", betaY[m*TILE_N+n + y_addr/MAX_SIZE].to_int(),m*TILE_N+n + y_addr/MAX_SIZE);
//
//		}
//	}
//
//}

void MAPELoss(XPack *Y, XPack *LABEL, SExpPack *betaY, SExpPack *betaL, int y_addr, int X_addr, int betaX_addr){
	ap_int<Ws> Abuf[MAX_SIZE],Bbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1

	XPack label_bf[B*FORECAST/MAX_SIZE];
	for(int i=0;i<B*FORECAST/MAX_SIZE;i++ ){
		label_bf[i] = LABEL[i + X_addr];
	}

	unsigned TILE_M = B/MAX_SIZE;
	unsigned TILE_N = FORECAST/BLK_SIZE;

	mape:for(int m=0;m < TILE_M;m++){
		for(int n=0;n< TILE_N;n++){

			for(int i=0;i<BLK_SIZE;i++){
				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
					if(n!=TILE_N-1 && i<BLK_SIZE - (FORECAST-OUTPUT) ){
//						printf("===i:%d, j:%d,k:%d\n",i,j,k);
						Abuf[j] = Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range(j*Ws+Ws-1,j*Ws);
						Bbuf[j] = label_bf[i+(m*TILE_N+n)*BLK_SIZE].range(j*Ws+Ws-1 ,j*Ws);
						MAPECore(Abuf[j], Bbuf[j], betaY[m*TILE_N+n + y_addr/BLK_SIZE].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP),
								betaL[m*TILE_N+n+betaX_addr].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP));
						Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range(j*Ws+Ws-1,j*Ws) = Abuf[j];
					}
				}
			}

			betaY[m*TILE_N+n + y_addr/BLK_SIZE] = -betaL[m*TILE_N+n+betaX_addr];

		}
	}

}

#endif

