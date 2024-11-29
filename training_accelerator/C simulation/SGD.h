#include "newCLZ.h"
#include "typedef.h"
#include "modeldef.h"


using ST = ap_fixed<Ws, 2, AP_RND, AP_SAT>;

void SGDCore(ap_uint<W> w, ap_uint<W> g, ap_uint<Ws> &fracC, SEXP_T betaW, SEXP_T betaG){

	ap_uint<1> sgnW = w[W-1];
	ap_uint<1> sgnG = g[W-1];

	ap_int<EBIT> ebias = betaW - betaG;
	ap_int<EBIT> Kw = ebias >= 0? ap_int<EBIT>(0) : ap_int<EBIT>(-ebias);
	ap_int<EBIT> Kg = ebias < 0? ap_int<EBIT>(0) : ap_int<EBIT>(ebias);

	ST a = 0, b =0, c=0;
	a(Ws-2,Ws-W) = w(W-2,0);
	a = sgnW == 0? a : ST(-a);
	b(Ws-2,Ws-W) = g(W-2,0);
	b = sgnG == 0? b : ST(-b);

	ST learning_rate = 0.05;

	c = (a >> Kw) - learning_rate*(b >> Kg);
	fracC[Ws-1] = c[Ws-1];
	fracC(Ws-2,0) = c[Ws-1]==0? c(Ws-2,0) : (~c(Ws-2,0) +1 );


}


void SGD(MemPack *WG, MemPack *G, SExpPack *betaW, SExpPack *betaG, unsigned layer, unsigned betaW_addr){

	ap_uint<W> Abuf[MAX_SIZE], Bbuf[MAX_SIZE], Wbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Wbuf type=complete dim=1
	ap_uint<Ws> fracC[BLK_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=fracC type=complete dim=2
	ap_int<5> Zshift[MB];
	ap_uint<W> Z_min[MB];
#pragma HLS ARRAY_PARTITION variable=Zshift type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Z_min type=complete dim=1
	ap_uint<Ws> check[MAX_SIZE], cc_max[MB], cc1[MAX_SIZE], cc2[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=check type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc1 type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc_max type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc2 type=complete dim=1
	SEXP_T betaAA[MB], betaBB[MB];
#pragma HLS ARRAY_PARTITION variable=betaAA type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=betaBB type=complete dim=1
	SExpPack betaW_WR;
	static ap_uint<9> state[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=state type=complete dim=1

	for(int i=0;i<MAX_SIZE;i++){
#pragma HLS UNROLL
		state[i] = i+1;
	}

	unsigned Bsize; unsigned Len; unsigned size;
	switch((layer+1)){
		case 1:
			size = LK*LOOKBACK; Bsize = LOOKBACK; Len = LK;
			break;
		case 2:
		case 3:
		case 4:
			size = LK * LK; Bsize = LK; Len = LK;
			break;
		case 5:
		case 7:
			size = THETA * LK; Bsize = THETA; Len = LK;
			break;
		case 6:
			size = LOOKBACK * THETA; Bsize = LOOKBACK; Len = THETA;
			break;
		case 8:
			size = FORECAST * THETA; Bsize = FORECAST; Len = THETA;
			break;
	}

	unsigned TILE_M = Bsize/MAX_SIZE;
	unsigned TILE_N = Len/BLK_SIZE;

	outer:for(int m=0;m<TILE_M;m++){
		for(int n=0;n<TILE_N;n++){

			for(int j=0;j<MB;j++){
#pragma HLS UNROLL
				betaAA[j] = betaW[m*TILE_N+n + betaW_addr].range(j*SEXP+SEXP-1, j*SEXP);
				betaBB[j] = betaG[m*TILE_N+n].range(j*SEXP+SEXP-1, j*SEXP);
			}

			for(int j=0;j<MAX_SIZE;j++){
				#pragma HLS UNROLL
				check[j] = 0;
			}

			for(int t=0;t<BLK_SIZE;t++){
//#pragma HLS PIPELINE
				for(int i=0; i<MAX_SIZE;i++){
#pragma HLS UNROLL
					Abuf[i] = WG[t + (m*TILE_N+n)*BLK_SIZE].range(W*i+W-1,W*i);
					Bbuf[i] = G[t + (m*TILE_N+n)*BLK_SIZE].range(W*i+W-1,W*i);
					SGDCore(Abuf[i], Bbuf[i], fracC[t][i],betaAA[int(i/BLK_SIZE)],betaBB[int(i/BLK_SIZE)] );
					cc1[i] = cc1[i] | fracC[t][i];
				}
			}
		//---------------------------------------------------------------------

			for(int j=0;j<MAX_SIZE;j++){
				#pragma HLS UNROLL
				cc2[j] = cc1[j];
			}

			for(int j=0;j<MB;j++){
		#pragma HLS UNROLL
				cc_max[j] = 0;
				for(int i=0;i<BLK_SIZE;i++){
		#pragma HLS UNROLL
					cc_max[j] = cc2[i+BLK_SIZE*j] | cc_max[j];
				}

				ap_uint<16> tt;
				tt(Ws-1,1) = cc_max[j](Ws-2,0);
				tt[0] = 1;
				Z_min[j] = CLZ16( tt );
				betaAA[j] = betaAA[j] >= betaBB[j]? (betaAA[j] - Z_min[j]) : (betaBB[j] - Z_min[j]);
				betaW_WR.range(j*SEXP+SEXP-1, j*SEXP) = betaAA[j];

			}
			betaW[m*TILE_N+n + betaW_addr] = betaW_WR;

		//---------------------------------------------------------------------
			for(int i=0;i<BLK_SIZE;i++){
				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
					ap_uint<Ws> aa = fracC[i][j] << Z_min[int(j/BLK_SIZE)];

					ap_uint<1> guard = aa[Ws-W];
					ap_uint<1> round = aa[Ws-W-1];
					ap_uint<1> sticky = aa(Ws-W-2,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
					ap_uint<1> randm = round & sticky | guard&round&(~sticky);

					Wbuf[j](W-2,0) = aa(Ws-2,Ws-W) + randm;
					Wbuf[j][W-1] = fracC[i][j][Ws-1];
					WG[i + (m*TILE_N+n)*BLK_SIZE].range(W*j+W-1,W*j) = Wbuf[j];
				}
			}

	}}




}

