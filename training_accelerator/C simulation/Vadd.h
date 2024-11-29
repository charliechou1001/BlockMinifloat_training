#ifndef VADDDEF_
#define VADDDEF_

#include "typedef.h"
#include "newCLZ.h"

//used in both activation / error, a is residual, b is NBEATS block output
//op: 0 add, 1 sub

template<int Ebit>
void VADDCore(ap_uint<Ws> &a, ap_uint<Ws> b, SEXP_T betaA, SEXP_T betaB, ap_uint<1> shift, ap_uint<1> shift_o, ap_uint<1> op){
#pragma HLS INLINE
//	printf("--------------------\n");
	ap_int<Ebit> ebias = betaA - betaB;
//	printf("ebias: %d, betaA: %d, betaB: %d\n",ebias.to_int(),betaA.to_int(),betaB.to_int());
	ap_int<Ebit> Ka = ebias >= 0? ap_int<Ebit>(0) : ap_int<Ebit>(-ebias);
	ap_int<Ebit> Kb = ebias < 0? ap_int<Ebit>(0) : ap_int<Ebit>(ebias);
//	printf("Ka: %d, Kb:%d, shift:%d\n",Ka.to_int(),Kb.to_int(),shift.to_int());

	ap_uint<1> sgnA = a[Ws-1];
	ap_uint<1> sgnB = ( ~b[Ws-1] & op ) | (b[Ws-1] & ~op);
	ap_uint<Ws> A = a(Ws-2,0) >> ( Ka + shift);
	ap_uint<Ws> B = b(Ws-2,0) >> Kb;

	ap_uint<Ws> fracA = sgnA == 0? A : ap_uint<Ws>(~A+1);
	ap_uint<Ws> fracB = sgnB == 0? B : ap_uint<Ws>(~B+1);
	ap_uint<Ws> fracC = fracA + fracB;

	ap_uint<1> signR = fracC[Ws-1];
	ap_uint<Ws-1> sgn_rst = signR? (~fracC(Ws-2,0)+1) : fracC(Ws-2,0);
	ap_uint<Ws+1> c = (signR, sgn_rst);

//	printf("a: %s, b:%s,\n A:%s, B:%s,\n fracA:%s, fracB: %s, fracC:%s,\n c:%s\n",
//			a.to_string(2).c_str(), b.to_string(2).c_str(), A.to_string(2).c_str(),B.to_string(2).c_str(),
//			fracA.to_string(2).c_str(),fracB.to_string(2).c_str(),fracC.to_string(2).c_str(),
//			c.to_string(2).c_str() );

//	const int ss = 14;
//	printf("a: %s, %f\n b: %s, %f\n c:%s, %f\n",a.to_string(2).c_str(), float(a)/std::pow(2,ss),
//			b.to_string(2).c_str(), float(b)/std::pow(2,ss), c.to_string(2).c_str(), float(c)/std::pow(2,ss));

	if(c > ((1<<Ws)-1) ){
		a = ((1<<Ws)-1);
		shift_o = 0b1;
	}else if( c < -(1<<Ws)+1 ){
		a = -(1<<Ws)+1;
		shift_o = 0b1;
	}else{
		a = c;
//		shift = 0b0;
	}
//	ap_uint<Ws-1> check = c>=0? c(Ws-2, 0) : (~c(Ws-2, 0)+1);
//	AddAND = AddAND | check;
////	printf("AddAND: %s, check: %s\n",AddAND.to_string(2).c_str(), check.to_string(2).c_str());
//
//	a = c;
}

void VecADD(XPack *res, XPack *out, SExpPack *betaR, SExpPack *betaT, ap_uint<MB> *shift,
		unsigned Bsize, unsigned Len, int res_addr, ap_uint<1> op, ap_uint<2> modeIn){

	unsigned TILE_M = Bsize/MAX_SIZE;
	unsigned TILE_N = Len/BLK_SIZE;

	ap_uint<Ws> Abuf[MAX_SIZE], Bbuf[MAX_SIZE];
#pragma HLS ARRAY_PARTITION variable=Abuf type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=Bbuf type=complete dim=1
	ap_uint<1> shift_o[MAX_SIZE];
	ap_uint<1> shift_o1[MB];
#pragma HLS ARRAY_PARTITION variable=shift_o type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=shift_o1 type=complete dim=1
	SEXP_T bt1[MB], bt2[MB];
#pragma HLS ARRAY_PARTITION variable=bt1 type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=bt2 type=complete dim=1

	outer:for(int m=0;m<TILE_M;m++){
		for(int n=0;n<TILE_N;n++){
//#pragma HLS DATAFLOW
			ap_uint<MB> cc = shift[m*TILE_N+n];
//			SExpPack bt1 = betaR[m*TILE_N+n+res_addr/BLK_SIZE];
//			SExpPack bt2 = betaT[m*TILE_N+n];

			for(int i=0;i<MB;i++){
				#pragma HLS UNROLL
				bt1[i] = betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP) + cc[i];
				bt2[i] = betaT[m*TILE_N+n].range(i*SEXP+SEXP-1, i*SEXP);
//				betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP) = shift[m*TILE_N+n][i] + betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP);
//				shift[m*TILE_N+n][i] = 0b0;
				shift_o1[i] = 0b0;
			}
			add:for(int t=0;t<BLK_SIZE;t++){
				#pragma HLS PIPELINE

				for(int i=0; i<MAX_SIZE;i++){
					#pragma HLS UNROLL
					Abuf[i] = res[t+(m*TILE_N+n)*BLK_SIZE+res_addr].range(i*Ws+Ws-1 ,i*Ws);
					Bbuf[i] = out[t+(m*TILE_N+n)*BLK_SIZE].range(i*Ws+Ws-1 ,i*Ws);
					ap_uint<1> sh = cc[int(i/BLK_SIZE)];
//					printf("----\n Abuf before: %s, Bbuf:%s, bt1:%d, bt2:%d, ", Abuf[i].to_string(2).c_str(), Bbuf[i].to_string(2).c_str(),
//							bt1[int(i/BLK_SIZE)].to_int(), bt2[int(i/BLK_SIZE)].to_int() );
					VADDCore<W>(Abuf[i], Bbuf[i], bt1[int(i/BLK_SIZE)], bt2[int(i/BLK_SIZE)], sh, shift_o[i], op);
//					printf("Abuf after: %s\n", Abuf[i].to_string(2).c_str() );
					res[t+(m*TILE_N+n)*BLK_SIZE+res_addr].range(i*Ws+Ws-1 ,i*Ws) = Abuf[i];
//					ss[int(i/BLK_SIZE)] = sh;
				}
//				shift[m*TILE_N+n] = ss;
			}

			for(int i=0;i<MB;i++){
				for(int j=0;j<BLK_SIZE;j++){
					#pragma HLS UNROLL
					shift_o1[i] = shift_o1[i] | shift_o[j+i*BLK_SIZE];
				}
			}
			ap_uint<MB> dd;
			for(int i=0;i<MB;i++){
#pragma HLS UNROLL
				dd[i] = shift_o1[i];
			}
			shift[m*TILE_N+n] = dd;

			for(int i=0;i<MB;i++){
				#pragma HLS UNROLL
//				SEXP_T betaA = betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP);
//				SEXP_T betaB = betaT[m*TILE_N+n].range(i*SEXP+SEXP-1, i*SEXP);
//				betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP) = betaA >= betaB ? betaA : betaB;
				betaR[m*TILE_N+n+res_addr/BLK_SIZE].range(i*SEXP+SEXP-1, i*SEXP) = bt1[i] >= bt2[i]? bt1[i] : bt2[i];
			}

		}
	}

}




#endif
