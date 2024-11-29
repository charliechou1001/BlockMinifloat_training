#include "newCLZ.h"
#include "typedef.h"
#include "modeldef.h"

void SGDCore(ap_uint<W> w, ap_uint<W> g, ADOUT_T &fracC,
		SEXP_T betaW, SEXP_T betaG, ap_uint<ADINT+ADDEC+1> &check){

	//	printf("w: %s,  g: %s\n",w.to_string(2).c_str(), g.to_string(2).c_str());
	ap_uint<1> sgnW = w[W-1];
	ap_uint<1> sgnG = g[W-1];

	ap_int<EBIT> ebias = betaW - betaG + ALPHA;
	ap_int<EBIT> Kw = ebias >= 0? ap_int<EBIT>(0) : ap_int<EBIT>(-ebias);
	ap_int<EBIT> Kg = ebias < 0? ap_int<EBIT>(0) : ap_int<EBIT>(ebias);

#if (E0R != 0)
	ap_uint<E0R> etaW = (1<<(E0R-1)) -1;//minifloat bias
	ap_uint<1> denormW = w(W-2,W-1-E0R) == 0? ap_uint<1>(1) : ap_uint<1>(0);
	ADIN_T fracW=0;
	fracW(ADDEC,ADDEC-M0R) = (~denormW, w(M0R-1,0));
	fracW = (1-2*sgnW)*( fracW << ( w(W-2,W-1-E0R) - etaW - Kw + denormW) );
//	printf("w(W-2,W-1-E0R): %d, etaW: %d, Kw: %d\n",w(W-2,W-1-E0R).to_int(), etaW.to_int(), Kw.to_int());
//	printf("( w(W-2,W-1-E0R) - etaW - Kw): %d\t",int( w(W-2,W-1-E0R) - etaW - Kw));
//	printf("fracW: %s\n",fracW.to_string(2).c_str());
#else
	ADIN_T fracW=0;
	fracW(ADDEC,ADDEC-M0R+1) = w(M0R-1,0);
	fracW = (1-2*sgnW)* (fracW >> Kw);
#endif

#if (ResE2 != 0)
	ap_uint<ResE2> etaG = (1<<(ResE2-1)) -1;//minifloat bias
	ap_uint<1> denormG = g(W-2,W-1-ResE2) == 0? ap_uint<1>(1) : ap_uint<1>(0);
	ADIN_T fracG=0;
	fracG(ADDEC,ADDEC-ResM2) = (~denormG, g(ResM2-1,0));
	fracG = (1-2*sgnG)*( fracG << ( g(W-2,W-1-ResE2) - etaG -Kg + denormG) );
//	printf("g(W-2,W-1-ResE2): %d, etaG: %d, Kg: %d\n",g(W-2,W-1-ResE2).to_int(), etaG.to_int(), Kg.to_int());
//	printf("( g(W-2,W-1-ResE2) - etaG -Kg) ): %d\t", int( g(W-2,W-1-ResE2) - etaG -Kg) );
//	printf("fracG: %s\n",fracG.to_string(2).c_str());
#else
	ADIN_T fracG=0;
	fracG(ADDEC,ADDEC-ResM2+1) = g(ResM2-1,0);
	fracG = (1-2*sgnG)* (fracG >> Kg);
#endif

//	printf("fracW float: %f,   fracG float: %f\n",float(fracW), float(fracG));

	fracC = fracW - fracG;
//	printf("Kw: %d, Kg: %d\n",Kw.to_int(),Kg.to_int());
//	 printf("fracC: %s fracC float: %f\n",fracC.to_string(2).c_str(), float(fracC));

	check = fracC>=0? fracC(ADINT+ADDEC, 0) : (~fracC(ADINT+ADDEC, 0)+1);
}

void betaW_Adjust(ap_uint<4> Z_min, SEXP_T &betaW, SEXP_T betaG, ap_int<5> &Zshift){
	Zshift = LE - LEW - Z_min;
	betaW = betaW >= betaG? (betaW + Zshift) : (betaG + Zshift);
}

//template<int mTR=4>
//ap_uint<1> st_round(ap_uint<mTR> tail){
//	ap_uint<mTR> random_bits;
//	ap_uint<1> m;
//
//	//LFSR 9 & st rounding
//	ap_uint<9> state = INIT;
//	m = state[8] ^~ state[4];
//	state(8,1) = state(7,0);
//	state[0] = m;
//	random_bits = state >> (9-mTR);
//
//	return (tail >= random_bits)? ap_uint<1>(1) : ap_uint<1>(0);
//}

template<int mTR=4>
void SGDNorm(ADOUT_T fracC, ap_uint<W> &w, ap_int<5> &Zshift){
#pragma HLS INLINE
//	printf("fracC input: %s  FP: %f\n", fracC.to_string(2).c_str(), float(fracC));
	fracC = fracC >> Zshift;
	ap_uint<1> sgnW = fracC[ADINT+ADDEC+1];
	ap_uint<ADINT+ADDEC> sgn_res = sgnW == 0? fracC(ADINT+ADDEC-1,0) : ( ~fracC(ADINT+ADDEC-1,0)+1 );

//	printf("fracC: %s\n  sgn_res: %s\n",fracC.to_string(2).c_str(), sgn_res.to_string(2).c_str());

#if E0R!=0
	#if W==8
		ap_uint<W> count = CLZ64((sgn_res, ap_uint<64-ADINT-ADDEC>(0)));
	#else
		ap_uint<W> count = CLZ16((sgn_res, ap_uint<16-ADINT-ADDEC>(0)));
	#endif

	count = count - ADINT + LEW;//align to the exp max value
	ap_uint<W> count_adjust;

//#if E0R!=0
	ap_uint<E0R> expW;
	if(count >= ((1<<E0R)-1) ){
		expW = 0;
		count_adjust = ((1<<E0R)-1) -1 + ADINT - LEW;
	}else{
		expW = (1<<E0R)-1 - count;
		count_adjust = count + ADINT - LEW;
	}
//	printf("count: %d, count_adjust:%d, expW: %d\n", count.to_int(), count_adjust.to_int(), expW.to_int());
	sgn_res = sgn_res << count_adjust;
	ap_uint<M0R> manW = sgn_res(ADINT+ADDEC-2,ADINT+ADDEC-M0R-1);//more -1 for removing implicit bit
//	printf("sgn_res: %s, manW: %s\n", sgn_res.to_string(2).c_str(), manW.to_string(2).c_str());
	ap_uint<1> rand = st_round<mTR>(sgn_res(ADINT+ADDEC-M0R-1,ADINT+ADDEC-M0R-1 - mTR+1) );
	if(manW < (1<<M0R)-1)
		manW = manW + rand;

	w[W-1] = sgnW;
	w(W-2,W-2-E0R+1) = expW;
	w(M0R-1,0) = manW;
#else
	sgn_res = sgn_res << (ADINT - LEW);
	ap_uint<M0R> manW = sgn_res(ADINT+ADDEC-1,ADINT+ADDEC-M0R);
//	printf("sgn_res: %s, manW: %s\n", sgn_res.to_string(2).c_str(), manW.to_string(2).c_str());
//	ap_uint<1> rand = st_round<mTR>(sgn_res(ADINT+ADDEC-M0R,ADINT+ADDEC-M0R- mTR+1) );
//	if(manW < (1<<M0R)-1)
//		manW = manW + rand;

	w[W-1] = sgnW;
	w(M0R-1,0) = manW;
#endif
}

void SGD(MemPack *WG, MemPack *G, SExpPack *betaW, SExpPack *betaG, unsigned layer){

	ap_uint<W> Abuf[MAX_SIZE], Bbuf[MAX_SIZE], Wbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Wbuf type=complete dim=1
	ADOUT_T fracC[BLK_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=fracC type=complete dim=2
	ap_int<5> Zshift[MB];
	ap_uint<W> Z_min[MB];
#pragma HLS ARRAY_PARTITION variable=Zshift type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Z_min type=complete dim=1
	ap_uint<ADINT+ADDEC+1> check[MAX_SIZE], cc_max[MB], cc1[MAX_SIZE], cc2[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=check type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc1 type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc_max type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=cc2 type=complete dim=1
	SEXP_T betaAA[MB], betaBB[MB];
#pragma HLS ARRAY_PARTITION variable=betaAA type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=betaBB type=complete dim=1
	SExpPack betaW_WR;

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
				cc_max[j] = 0;
				betaAA[j] = betaW[m*TILE_N+n].range(j*SEXP+SEXP-1, j*SEXP);
				betaBB[j] = betaG[m*TILE_N+n].range(j*SEXP+SEXP-1, j*SEXP);
			}

			for(int t=0;t<BLK_SIZE;t++){
#pragma HLS PIPELINE
				if(t==0){
					for(int j=0;j<MAX_SIZE;j++){
						#pragma HLS UNROLL
						check[j] = 0;
					}
				}

				for(int i=0; i<MAX_SIZE;i++){
#pragma HLS UNROLL
					Abuf[i] = WG[t].range(W*i+W-1,W*i);
					Bbuf[i] = G[t].range(W*i+W-1,W*i);
		//			printf("------------[%d][%d]----------------\n",t,i);
					SGDCore(Abuf[i], Bbuf[i], fracC[t][i],betaAA[int(i/BLK_SIZE)],
							betaBB[int(i/BLK_SIZE)], check[i]);
					cc1[i] = cc1[i] | check[i];
				}
			}
		//---------------------------------------------------------------------

			for(int j=0;j<MAX_SIZE;j++){
				#pragma HLS UNROLL
				cc2[j] = cc1[j];
			}

			for(int j=0;j<MB;j++){
		#pragma HLS UNROLL
				for(int i=0;i<BLK_SIZE;i++){
		#pragma HLS UNROLL
					cc_max[j] = cc2[i] | cc_max[j];
				}

				#if W==8
					Z_min[j] = CLZ64((cc_max[j], ap_uint<64-ADINT-ADDEC>(0)));
				#else
					Z_min[j] = CLZ16((cc_max[j], ap_uint<16-ADINT-ADDEC>(0)));
				#endif

				betaW_Adjust(Z_min[j], betaAA[j] , betaBB[j], Zshift[j]);
				betaW_WR.range(j*SEXP+SEXP-1, j*SEXP) = betaAA[j];
		//		printf("Zmin: %d, Zshift: %d, betaW: %d\n",Z_min[j].to_int(),Zshift[j].to_int(), betaW.to_int());
			}
			betaW[m*TILE_N+n] = betaW_WR;

		//---------------------------------------------------------------------

			for(int t=0;t<BLK_SIZE;t++){
#pragma HLS PIPELINE
				for(int i=0;i<MAX_SIZE;i++){
#pragma HLS UNROLL
//					printf("------------[%d][%d]----------------\n",t,i);
					SGDNorm<4>(fracC[t][i], Wbuf[i], Zshift[int(i/BLK_SIZE)]);
					WG[t].range(W*i+W-1,W*i) = Wbuf[i];
//					printf("W: %s\n",Wbuf[i].to_string(2).c_str());
				}
			}
	}}

//	for(int i=0;i<TILE_M*TILE_N;i++){
//#pragma HLS UNROLL
//		betaW[i] = betaW_bf[i];
//	}



}

