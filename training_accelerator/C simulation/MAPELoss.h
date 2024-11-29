//MAPE loss derivation reference:  https://stats.stackexchange.com/questions/316076/gradient-and-hessian-of-the-mape
#ifndef MAPEDEF_
#define MAPEDEF_

#include "typedef.h"
#include "newCLZ.h"
//#include "Cvt.h"
#include "modeldef.h"
#include <cmath>

using FT = ap_fixed<Ws, 2, AP_RND, AP_SAT>;

void MAPECore(ap_uint<Ws> &y, ap_uint<Ws> label, SEXP_T beta_y, SEXP_T beta_label){
#pragma HLS INLINE
	FT a=0,b=0,c=0;
	ap_int<W> ebias = beta_y - beta_label;
	ap_int<W> Ka = ebias >= 0? ap_int<W>(0) : ap_int<W>(-ebias);
	ap_int<W> Kb = ebias < 0? ap_int<W>(0) : ap_int<W>(ebias);

	a(Ws-2,0) = y(Ws-2,0);
	a = y[Ws-1] == 0? a : FT(-a);

	b(Ws-2,0) = label(Ws-2,0);

	b = label[Ws-1] == 0? b : FT(-b);
	b = b!=0? b : FT(1);

	c = ( (a >> Ka) - (b >> Kb) )>=0? FT(1)/(BSIZE*OUTPUT*b) : FT(-1)/(BSIZE*OUTPUT*b);

	y[Ws-1] = c[Ws-1];
	y(Ws-2,0) = c[Ws-1]==0? c(Ws-2,0) : (~c(Ws-2,0) +1 );

}


void MAPELoss(XPack *Y, XPack *LABEL, SExpPack *betaY, SExpPack *betaL, int y_addr, int X_addr, int betaX_addr){
	ap_uint<Ws> Abuf[MAX_SIZE],Bbuf[MAX_SIZE], Dbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Dbuf type=complete dim=1

	ap_uint<Ws> Cbuf[BLK_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Cbuf type=complete dim=2

	ap_uint<Ws> check[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=check type=complete dim=1
	ap_uint<Ws> cc_max[MB];
#pragma HLS ARRAY_PARTITION variable=cc_max type=complete dim=1
	ap_uint<W> Z_min[MB];
#pragma HLS ARRAY_PARTITION variable=Z_min type=complete dim=1

	XPack label_bf[B*FORECAST/MAX_SIZE];
	for(int i=0;i<B*FORECAST/MAX_SIZE;i++ ){
		label_bf[i] = LABEL[i + X_addr];
	}

	unsigned TILE_M = B/MAX_SIZE;
	unsigned TILE_N = FORECAST/BLK_SIZE;

	mape:for(int m=0;m < TILE_M;m++){
		for(int n=0;n< TILE_N;n++){

			for(int i=0;i<MAX_SIZE;i++){
#pragma HLS UNROLL
				check[i] = 0;
			}

			for(int i=0;i<BLK_SIZE;i++){
				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
					bool conda = n < TILE_N/2-1 || ( n == TILE_N/2-1 && i<BLK_SIZE - (FORECAST/2-OUTPUT) );
					bool condb = m < TILE_M-1 || ( m == TILE_M-1 && j<MAX_SIZE - (B - BSIZE) );

						Abuf[j] = Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range(j*Ws+Ws-1,j*Ws);
						Bbuf[j] = label_bf[i+(m*TILE_N+n)*BLK_SIZE].range(j*Ws+Ws-1 ,j*Ws);
						MAPECore(Abuf[j], Bbuf[j], betaY[m*TILE_N+n + y_addr/BLK_SIZE].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP),
								betaL[m*TILE_N+n+betaX_addr].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP));
						Cbuf[i][j] = conda != 0 && condb != 0? Abuf[j] : ap_uint<Ws>(0);
						check[j] = check[j] | Cbuf[i][j];

				}
			}

			for(int j=0;j<MB;j++){
		#pragma HLS UNROLL
				cc_max[j] = 0;
				for(int i=0;i<BLK_SIZE;i++){
		#pragma HLS UNROLL
					cc_max[j] = check[i+j*BLK_SIZE] | cc_max[j];
				}

				Z_min[j] = CLZ16( (cc_max[j](Ws-2,0),ap_uint<1>(1)) );
			}

			SEXP_T betaC;
			SExpPack betaA;
			SExpPack betaB = betaL[m*TILE_N+n+betaX_addr];
			for(int j=0;j<MB;j++){
				#pragma HLS UNROLL
				bool condc = n <= TILE_N/2-1 &&( m < TILE_M-1 || (m == TILE_M-1 && j < MB - int( (B - BSIZE)/BLK_SIZE ) ) );
				SEXP_T bb = betaB.range(SEXP*j+SEXP-1,SEXP*j);
				betaC = condc? SEXP_T(-bb- Z_min[j]) : SEXP_T(-90);
				betaA.range(SEXP*j+SEXP-1,SEXP*j) = betaC;
			}
			betaY[m*TILE_N+n + y_addr/BLK_SIZE] = betaA;

			for(int i=0;i<BLK_SIZE;i++){
				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
					Dbuf[j](Ws-2,0) = Cbuf[i][j](Ws-2,0) << Z_min[int(j/BLK_SIZE)];
					Dbuf[j][Ws-1] =  Cbuf[i][j][Ws-1];
					Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range(j*Ws+Ws-1,j*Ws) = Dbuf[j];
				}
			}


		}
	}

}

float MAPEPrint_Core(ap_uint<Ws> Abuf, ap_uint<Ws> Bbuf, SEXP_T betaA, SEXP_T betaB){
	FT a=0,b=0,c=0;

	a(Ws-2,0) = Abuf(Ws-2,0);
	a = Abuf[Ws-1] == 0? a : FT(-a);

	b(Ws-2,0) = Bbuf(Ws-2,0);
	b = Bbuf[Ws-1] == 0? b : FT(-b);

	float aa = a.to_float() * std::pow(2,betaA.to_int() );
	float bb = b.to_float() * std::pow(2,betaB.to_int() );

	float div;
	if(b!=0){
		div = (aa-bb)/bb;
		if(div < 0)
			div = -div;
	}else
		div =0;

	return div;

}

float MAPEPrint(XPack *Y, XPack *LABEL, SExpPack *betaY, SExpPack *betaL, int y_addr, int X_addr, int betaX_addr){

	ap_uint<Ws> Abuf[MAX_SIZE],Bbuf[MAX_SIZE];

	XPack label_bf[B*FORECAST/MAX_SIZE];
	for(int i=0;i<B*FORECAST/MAX_SIZE;i++ ){
		label_bf[i] = LABEL[i + X_addr];
	}

	unsigned TILE_M = B/MAX_SIZE;
	unsigned TILE_N = FORECAST/BLK_SIZE;

    float sum=0;

//    printf("LossPrint\n");
	mape:for(int m=0;m < TILE_M;m++){
		for(int n=0;n< TILE_N;n++){

			for(int i=0;i<BLK_SIZE;i++){
				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
						Abuf[j] = Y[i+(m*TILE_N+n)*BLK_SIZE + y_addr].range(j*Ws+Ws-1,j*Ws);
						Bbuf[j] = label_bf[i+(m*TILE_N+n)*BLK_SIZE].range(j*Ws+Ws-1 ,j*Ws);

						SEXP_T betaA = betaY[m*TILE_N+n + y_addr/BLK_SIZE].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP);
						SEXP_T betaB = betaL[m*TILE_N+n+betaX_addr].range(int(j/BLK_SIZE)*SEXP + SEXP-1,int(j/BLK_SIZE)*SEXP);


						bool conda = n < TILE_N/2-1 || ( n == TILE_N/2-1 && i<BLK_SIZE - (FORECAST/2-OUTPUT) );
						bool condb = m < TILE_M-1 || ( m == TILE_M-1 && j<MAX_SIZE - (B - BSIZE) );

						float res = conda != 0 && condb != 0? MAPEPrint_Core(Abuf[j], Bbuf[j], betaA, betaB) : 0;
						sum += res;

				}
			}

		}
	}

	printf("loss: %f\n", sum/(BSIZE*OUTPUT) );

	return sum/(BSIZE*OUTPUT);
}

#endif

