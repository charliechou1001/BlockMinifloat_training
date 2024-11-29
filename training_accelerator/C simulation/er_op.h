#ifndef ER_OP_DEF_
#define ER_OP_DEF_

#include "newCLZ.h"
#include "typedef.h"
//#include <cmath>

//#define ER_DEBUG
//#define MAC_DEBUG

using ET = ap_fixed<Ws, 2, AP_RND, AP_SAT>;

void ERCore(ap_uint<W> w, ap_uint<W> g, ap_int<Kadd> &fracC,
SEXP_T betaW, SEXP_T betaG){

		ap_uint<1> sgnW = w[W-1];
		ap_uint<1> sgnG = g[W-1];

		ap_int<EBIT> ebias = betaW - betaG;
		ap_int<EBIT> Kw = ebias >= 0? ap_int<EBIT>(0) : ap_int<EBIT>(-ebias);
		ap_int<EBIT> Kg = ebias < 0? ap_int<EBIT>(0) : ap_int<EBIT>(ebias);

		ET a = 0, b =0, c=0;
		a(Ws-2,Ws-W) = w(W-2,0);
		a = sgnW == 0? a : ST(-a);
		b(Ws-2,Ws-W) = g(W-2,0);
		b = sgnG == 0? b : ST(-b);

		c = (a >> Kw) - (b >> Kg);
		fracC[Ws-1] = c[Ws-1];
		fracC(Ws-2,0) = c[Ws-1]==0? c(Ws-2,0) : (~c(Ws-2,0) +1 );

}

template<int Zm=W, int Zsh=W>
void ERBias(SEXP_T &betaW, SEXP_T betaG, ap_uint<Zm> Z_min, ap_int<Zsh> &Zshift){

#if E1L!=0
	const int eL = (1<<(E1L-1))-1;
#else
	const int eL = 0;
#endif

	ap_uint<8> CntOverflow = (Kadd-1) - (eL-1+M1L)*2 - ( (1<<E1L)-1 -eL + 1) -WT;
//---------------------------------------------------------------------
	ap_int<Zm> emax_shift = CntOverflow - Z_min;//1 for sign bit
	Zshift = emax_shift;
	betaW = betaW >= betaG? (betaW+Zshift) : (betaG+Zshift);
}


template<int Zsh=W, int Zm=W, int mTR=4>
void ERNorm(ap_int<Kadd> mac_out, ap_int<Zsh> Zshift,ap_uint<W> &R){

#if E1L!=0
	const int eL = (1<<(E1L-1))-1;
#else
	const int eL = 0;
#endif

	ap_uint<8> expZeroPoint = (Kadd -1) - (eL-1+M1L)*2 - WT;//used for validation
	ap_uint<8> CntDenorm = (Kadd-1) - (eL-1+M1L)*2 +(eL-1) - WT;//used for validation


	//sign bit and two's complement transform
	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);
	sgn_rst = sgn_rst >> Zshift;
#ifdef ER_DEBUG
	printf("mac_out: %s,\n sgn_rst: %s\n",mac_out.to_string(2).c_str(),sgn_rst.to_string(2).c_str());
#endif

	ap_uint<Zm> count;

#if E1L != 0
	count = CLZ64((sgn_rst,ap_uint<64-Kadd+1>(1)));

	ap_uint<E1L> exp;
	if(count >= CntDenorm)//denorm and underflow
		exp = 0;
	else
		exp = expZeroPoint - (count + 1) + eL;//1 for implicit bit
#endif

#ifdef ER_DEBUG
//	printf("count_adjust:%d\t", count_adjust.to_int());
	printf("Kadd: %d, expZeroPoint: %d, count: %d, ResBias: %d\t",Kadd,expZeroPoint.to_int(),count.to_int(), eL);
//	printf("exp: %s\n",exp.to_string(2).c_str());
#endif

	//mantissa and rounding
	ap_uint<Zm> count_adjust;
	ap_uint<M1L> man;
	ap_uint<1> rand;
	if(count < CntDenorm && E1L != 0)
		count_adjust = count;
	else if(count >= CntDenorm && E1L != 0)//denorm and underflow
		count_adjust = ap_uint<8>(CntDenorm-1);
	else//BFP out
		count_adjust = expZeroPoint -2;//expZeroPoint -1, and another -1 to align to sgn_rst(Kadd-3,)
	sgn_rst = sgn_rst << count_adjust;
#ifdef ER_DEBUG
	printf("count_adjust:%d\t", count_adjust.to_int());
	printf("sgn_rst after: %s\n",sgn_rst.to_string(2).c_str());
#endif

	man = sgn_rst(Kadd-3,Kadd-2-M1L);
	ap_uint<1> guard = sgn_rst[Kadd-2-M1L];
	ap_uint<1> round = sgn_rst[Kadd-3-M1L];
	ap_uint<1> sticky = sgn_rst(Kadd-4-M1L,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
	rand = round & sticky | guard&round&(~sticky);
	if(man < (1<<M1L)-1)
		man = man + rand;
#if E1L != 0
	else if(rand==1 && exp < (1<<E1L)-1){
		exp ++;
		man = 0;
	}
#endif

#ifdef ER_DEBUG
	printf("man:%s\n",man.to_string(2).c_str());
#endif

	//output
#if (S1L != 0 && E1L !=0)//BM
		R[W-1] = signR;
		R(W-2,W-1-E1L) = exp;
		R(M1L-1,0) = man;
#elif (S1L == 0 && E1L !=0)//unsigned BM
		R(W-1,W-E1L) = exp;
		R(M1L-1,0) = man;
#elif (S1L != 0 && E1L ==0)//BFP
		R[W-1] = signR;
		R(M1L-1,0) = man;
#else //unsigned BFP
		R(M1L-1,0) = man;
#endif


#ifdef ER_DEBUG
	printf("R: %s\n",R.to_string(2).c_str());
#endif

}

void ER_OP(MemPack *A, SExpPack *betaA1, SExpPack *betaA2,
		unsigned Bsize, unsigned Len,int addr1, int addr2){

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

	unsigned TILE_M = Bsize/MAX_SIZE;
	unsigned TILE_N = Len/BLK_SIZE;

	outer:for(int m=0;m<TILE_M;m++){
			for(int n=0;n<TILE_N;n++){

				for(int j=0;j<MB;j++){
	#pragma HLS UNROLL
					betaAA[j] = betaA1[m*TILE_N+n].range(j*SEXP+SEXP-1, j*SEXP);
					betaBB[j] = betaA2[m*TILE_N+n].range(j*SEXP+SEXP-1, j*SEXP);
				}

				for(int j=0;j<MAX_SIZE;j++){
					#pragma HLS UNROLL
					check[j] = 0;
				}

				for(int t=0;t<BLK_SIZE;t++){
	//#pragma HLS PIPELINE
					for(int i=0; i<MAX_SIZE;i++){
	#pragma HLS UNROLL
						Abuf[i] = A[t+(m*TILE_N+n)*BLK_SIZE+addr1].range(i*W+W-1,i*W);
						Bbuf[i] = A[t+(m*TILE_N+n)*BLK_SIZE+addr2].range(i*W+W-1,i*W);
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
				betaA1[m*TILE_N+n] = betaW_WR;

			//---------------------------------------------------------------------
				for(int i=0;i<BLK_SIZE;i++){
					for(int j=0;j<MAX_SIZE;j++){
						#pragma HLS UNROLL
						ap_uint<Ws> aa = fracC[i][j] << Z_min[int(j/BLK_SIZE)];
						ap_uint<1> guard = aa[Ws-W];
						ap_uint<1> round = aa[Ws-W-1];
						ap_uint<1> sticky = aa(Ws-W-2,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
						ap_uint<1> randn = round & sticky | guard&round&(~sticky);

						Wbuf[j](W-2,0) = aa(Ws-2,Ws-W) + randn;
						Wbuf[j][W-1] = fracC[i][j][Ws-1];
						A[i+(m*TILE_N+n)*BLK_SIZE+addr1].range(W*j+W-1,W*j) = Wbuf[j];
					}
				}

		}}

}

void NEG(MemPack *A, unsigned Bsize, unsigned Len,int addr1){
	ap_uint<W> Abuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1

	smaller:for(int t=0;t<Bsize*Len/MAX_SIZE;t++){
#pragma HLS PIPELINE
		for(int i=0; i<MAX_SIZE;i++){
			#pragma HLS UNROLL
			Abuf[i] = A[t+addr1].range(W*i+W-1,W*i);

			A[t+addr1].range(W*i+W-1,W*i) = (~Abuf[i][W-1] , Abuf[i](W-2,0));
		}
	}
}


#endif
