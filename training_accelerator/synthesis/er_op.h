#ifndef ER_OP_DEF_
#define ER_OP_DEF_

#include "newCLZ.h"
#include "typedef.h"
//#include <cmath>

//#define ER_DEBUG
//#define MAC_DEBUG

void ERCore(ap_uint<W> w, ap_uint<W> g, ap_int<Kadd> &fracC,
		SEXP_T betaW, SEXP_T betaG, ap_uint<Kadd-1> &check){

//	printf("w: %s,  g: %s\n",w.to_string(2).c_str(), g.to_string(2).c_str());
	ap_uint<1> sgnW = w[W-1];
	ap_uint<1> sgnG = g[W-1];

	ap_int<EBIT> ebias = betaW - betaG;
	ap_int<EBIT> Kw = ebias >= 0? ap_int<EBIT>(0) : ap_int<EBIT>(-ebias);
	ap_int<EBIT> Kg = ebias < 0? ap_int<EBIT>(0) : ap_int<EBIT>(ebias);
//	printf("Kw: %d, Kg: %d\n",Kw.to_int(), Kg.to_int());

#if E1L!=0
	const int eL = (1<<(E1L-1))-1;
#else
	const int eL = 0;
#endif
	int shift =  2*(M1L -1+eL) + WT;
//	printf("expZeroPoint: %d\n",expZeroPoint.to_int());

#if (E1L != 0)
	ap_uint<E1L> etaW = (1<<(E1L-1)) -1;//minifloat bias
	ap_uint<1> denormW = w(W-2,W-1-E1L) == 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_int<Kadd> fracW = (1-2*sgnW)*( ap_int<Kadd>( (~denormW, w(M1L-1,0)) ) << ( (w(W-2,W-1-E1L) - 1 + denormW) - Kw + (shift-M1L)) );

	ap_uint<E1L> etaG = (1<<(E1L-1)) -1;//minifloat bias
	ap_uint<1> denormG = g(W-2,W-1-E1L) == 0? ap_uint<1>(1) : ap_uint<1>(0);
	ap_int<Kadd> fracG = (1-2*sgnG)*( ap_int<Kadd>( (~denormG, g(M1L-1,0)) ) << ( (g(W-2,W-1-E1L) - 1 + denormG) - Kg + (shift-M1L)) );

//	printf("fracW:%s\n", fracW.to_string(2).c_str());
	#ifdef ER_DEBUG
//	printf("w(W-2,W-1-E1L): %d, etaW: %d, Kw: %d\n",w(W-2,W-1-E1L).to_int(), etaW.to_int(), Kw.to_int());
	printf("( w(W-2,W-1-E1L) - 1 - Kw): %d\n",int( w(W-2,W-1-E1L) - 1 - Kw));
	printf("fracW: %s, float fracW: %f\n",fracW.to_string(2).c_str(), float(fracW)/std::pow(2,shift));
	printf("g: %s\n",g.to_string(2).c_str());
	#endif
#else

	ap_int<Kadd> fracW = (1-2*sgnW)*( ap_int<Kadd>( w(M1L-1,0) ) <<( (shift-M1L+1) - Kw)  );
	ap_int<Kadd> fracG = (1-2*sgnG)*( ap_int<Kadd>( g(M1L-1,0) ) <<( (shift-M1L+1) - Kg)  );
#endif

	#ifdef ER_DEBUG
//	int shift =  2*(M1L -1) + WT;
	printf("fracG after >>Kg: %s\n", g.to_string(2).c_str());
	printf("fracW float: %f,   fracG float: %f\n",float(fracW)/std::pow(2,shift), float(fracG)/std::pow(2,shift));
	#endif
	fracC = fracW + fracG;
	#ifdef ER_DEBUG
//	printf("Kw: %d, Kg: %d\n",Kw.to_int(),Kg.to_int());
	 printf("fracC: %s fracC float: %f\n",fracC.to_string(2).c_str(), float(fracC)/std::pow(2,shift));
	#endif
	check = fracC>=0? fracC(Kadd-2, 0) : (~fracC(Kadd-2, 0)+1);

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
//	printf("CntOverflow:%d\t Z_min:%d\t emax_shift:%d\n",CntOverflow.to_int(),Z_min.to_int(), emax_shift.to_int());
	Zshift = emax_shift;
//	Zshift = emax_shift>0 ? ap_uint<Zsh>(emax_shift) : ap_uint<Zsh>(0);
//	beta =  beta + Zshift;
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

//	printf("expZeroPoint: %d,  CntDenorm: %d\n",expZeroPoint.to_int(),CntDenorm.to_int());

	//sign bit and two's complement transform
	ap_uint<1> signR = mac_out[Kadd-1];
	ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(mac_out(Kadd-2,0))+1) : mac_out(Kadd-2,0);
	sgn_rst = sgn_rst >> Zshift;
#ifdef MAC_DEBUG
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

#ifdef MAC_DEBUG
//	printf("count_adjust:%d\t", count_adjust.to_int());
	printf("Kadd: %d, expZeroPoint: %d, count: %d, ResBias: %d\t",Kadd,expZeroPoint.to_int(),count.to_int(), eL);
	printf("exp: %s\n",exp.to_string(2).c_str());
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
#ifdef MAC_DEBUG
	printf("count_adjust:%d\t", count_adjust.to_int());
	printf("sgn_rst after: %s\n",sgn_rst.to_string(2).c_str());
#endif

	man = sgn_rst(Kadd-3,Kadd-2-M1L);
	ap_uint<1> guard = sgn_rst[Kadd-2-M1L];
	ap_uint<1> round = sgn_rst[Kadd-3-M1L];
	ap_uint<1> sticky = sgn_rst(Kadd-4-M1L,0)==0? ap_uint<1>(0) : ap_uint<1>(1);
//		printf("sgn_rst: %s, guard: %s, round: %s, sticky: %s, sgn_rst(Kadd-4-ResM0,0): %s\n",
//				sgn_rst.to_string(2).c_str(),guard.to_string(2).c_str(), round.to_string(2).c_str(), sticky.to_string(2).c_str(), sgn_rst(Kadd-4-ResM0,0).to_string(2).c_str());
	rand = round & sticky | guard&round&(~sticky);
//		rand = st_round<mTR>(sgn_rst(Kadd-2-ResM0,Kadd-2-ResM0-mTR+1));
	if(man < (1<<M1L)-1)
		man = man + rand;
#if E1L != 0
	else if(rand==1 && exp < (1<<E1L)-1){
		exp ++;
		man = 0;
	}
#endif

#ifdef MAC_DEBUG
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


#ifdef MAC_DEBUG
	printf("R: %s\n",R.to_string(2).c_str());
#endif

}

void ER_OP(MemPack *A, SExpPack *betaA1, SExpPack *betaA2,
		unsigned Bsize, unsigned Len,int addr1, int addr2, ap_uint<1> is_add){

	unsigned TILE_M = Bsize/MAX_SIZE;
	unsigned TILE_N = Len/BLK_SIZE;

	ap_uint<W> Abuf[MAX_SIZE];
	ap_uint<W> Bbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
	ap_uint<W> Wbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Wbuf type=complete dim=1

	ap_int<Kadd> fracC[BLK_SIZE][MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=fracC type=complete dim=2
	ap_uint<Kadd-1> AddAnd[MAX_SIZE], AddAnd2[MAX_SIZE], check[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=AddAnd type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=AddAnd2 type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=check type=complete dim=1
	ap_int<Kadd> CheckZ[MB];
	ap_uint<W> Z[MB];
	ap_int<W> Zshift[MB];
#pragma HLS ARRAY_PARTITION variable=CheckZ type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Z type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Zshift type=complete dim=1
	SEXP_T betaAA[MB], betaBB[MB];
#pragma HLS ARRAY_PARTITION variable=betaAA type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=betaBB type=complete dim=1
	SExpPack betaW_WR;
//	SExpPack betaW_bf[LK*LK/(MAX_SIZE*BLK_SIZE)];

	if(is_add == 0b1){
		outer:for(int m=0;m<TILE_M;m++){
			for(int n=0;n<TILE_N;n++){

			for(int i=0;i<MAX_SIZE;i++){
				#pragma HLS UNROLL
				check[i] = 0;
			}
			for(int j=0;j<MB;j++){
				#pragma HLS UNROLL
				betaAA[j] = betaA1[m*TILE_N+n].range(j*SEXP+SEXP-1,j*SEXP);
				betaBB[j] = betaA2[m*TILE_N+n].range(j*SEXP+SEXP-1,j*SEXP);
			}
			first:for(int t=0;t<BLK_SIZE;t++){
	#pragma HLS PIPELINE
				loadloop:for(int i=0; i<MAX_SIZE;i++){
					#pragma HLS UNROLL
					Abuf[i] = A[t+(m*TILE_N+n)*BLK_SIZE+addr1].range(i*W+W-1,i*W);
					Bbuf[i] = A[t+(m*TILE_N+n)*BLK_SIZE+addr2].range(i*W+W-1,i*W);
	//				printf("---------[%d][%d]--------\n",t+(m*TILE_N+n)*MAX_SIZE,i);
					ERCore(Abuf[i], Bbuf[i], fracC[t][i], betaAA[int(i/BLK_SIZE)], betaBB[int(i/BLK_SIZE)], check[i]);
					AddAnd[i] = AddAnd[i] | check[i];
				}
			}
		//---------------------------------------------------------------------
			for(int j=0;j<MAX_SIZE;j++){
#pragma HLS UNROLL
				AddAnd2[j] = AddAnd[j];
			}

			for(int j=0;j<MB;j++){
#pragma HLS UNROLL
				for(int i=0;i<BLK_SIZE;i++){
					#pragma HLS UNROLL
					CheckZ[j] = CheckZ[j] | AddAnd2[i];
				}
			}

			for(int i=0;i<MB;i++){
#pragma HLS UNROLL
				Z[i] = CLZ64((CheckZ[i],ap_uint<64-Kadd+1>(1)));
//				SEXP_T betaA = betaA1[m*TILE_N+n].range(i*SEXP+SEXP-1,i*SEXP);
//				SEXP_T betaB = betaA2[m*TILE_N+n].range(i*SEXP+SEXP-1,i*SEXP);
		//		printf("betaA1 before: %d, ",betaA1[m*TILE_N+n+W_addr/MAX_SIZE].to_int());
				ERBias<W,W>(betaAA[i], betaBB[i], Z[i], Zshift[i]);
				betaW_WR.range(i*SEXP+SEXP-1,i*SEXP) = betaAA[i];
			}
			betaA1[m*TILE_N+n] = betaW_WR;
	//		printf("betaA1: %d, betaA2: %d, Z:%d, Zshift:%d\n",betaA1[m*TILE_N+n+W_addr/MAX_SIZE].to_int(), betaA2[m*TILE_N+n].to_int(), Z.to_int(), Zshift.to_int());
		//---------------------------------------------------------------------

			secondone:for(int t=0;t<BLK_SIZE;t++){
	#pragma HLS PIPELINE
				norm:for(int i=0;i<MAX_SIZE;i++){
					#pragma HLS UNROLL
	//				printf("===========[%d][%d]=========\n",t+m*MAX_SIZE,i+n*MAX_SIZE);
					ERNorm<W,W,4>(fracC[t][i], Zshift[int(i/BLK_SIZE)], Wbuf[i]);
					A[t+(m*TILE_N+n)*BLK_SIZE+addr1].range(i*W+W-1,i*W) = Wbuf[i];
				}
			}

		}}

//		for(int i=0;i<TILE_M*TILE_N;i++){
//			betaA1[i] = betaW_bf[i];
//		}

	}else{
		smaller:for(int t=0;t<Bsize*Len/MAX_SIZE;t++){
			for(int i=0; i<MAX_SIZE;i++){
				#pragma HLS UNROLL
				Abuf[i] = A[t+addr1].range(W*i+W-1,W*i);
				ap_uint<W> res = Abuf[i];
				A[t+addr1].range(W*i+W-1,W*i) = (~res[W-1] , res(W-2,0));
			}
		}
	}
}


#endif
