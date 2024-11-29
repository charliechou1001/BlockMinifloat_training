#ifndef TYPEDEF_
#define TYPEDEF_

#include "ap_int.h"

//W = S0L+E0L+M0L = S1L+E1L+M1L = S0R+E0R+M0R = S1R+E1R+M1R
//#define W 8
const unsigned W=8;
//bit-length of data, 8 or 4 bit
#define WI 5
#define WT 3
//extra bit for accumulator to avoid overflow

//left input mode 0(FW activation)
#define S0L 1
#define E0L 2
#define M0L 5

//left input mode 1(error)
#define S1L 1
#define E1L 0
#define M1L 7

//right input mode 0(weight)
#define S0R 1
#define E0R 4
#define M0R 3

//right input mode 1(BP activation for gradient)
#define S1R 1
#define E1R 0
#define M1R 7

//output mode 0(forward)
#define ResS0 1
#define ResE0 3
#define ResM0 4

//output mode 1(error)
#define ResS1 1
#define ResE1 0
#define ResM1 7

//output mode 2(gradient)
#define ResS2 0
#define ResE2 4
#define ResM2 4

//#define M 16
//#define N 16
//#define K 16

#define MAX_SIZE 16

const unsigned L = 64;
const unsigned M = L;
const unsigned N = L;
const unsigned K = L;

const unsigned TILE_M = M/MAX_SIZE;//used for test_longBMtofp32.h
const unsigned TILE_N = N/MAX_SIZE;
const unsigned TILE_K = K/MAX_SIZE;
//const unsigned TileSize = MAX_SIZE;

#define BLK_SIZE 8
//const unsigned BLK_SIZE = 4;
const unsigned MB = MAX_SIZE/BLK_SIZE;
const unsigned BK_M = M/BLK_SIZE;
const unsigned BK_N = N/BLK_SIZE;
const unsigned BK_K = K/BLK_SIZE;


#if E0L==0 && E1L==0
	#define EL 1
#else
	#if E0L > E1L
		#define EL E0L
	#else
		#define EL E1L
	#endif
#endif

#if E0R==0 && E1R==0
	#define ER 1
#else
	#if(E0R > E1R)
		#define ER E0R
	#else
		#define ER E1R
	#endif
#endif

#if M0L==0 & M1L==0
	#define ML 1
#else
	#if M0L > M1L
		#define ML M0L
	#else
		#define ML M1L
	#endif
#endif

#if M0R==0 && M1R==0
	#define MR 1
#else
	#if(M0R > M1R)
		#define MR M0R
	#else
		#define MR M1R
	#endif
#endif

#if ResE0==0 && ResE1==0 && ResE2==0
	#define OUTE 1
#else
	#if(ResE0 >= ResE1 && ResE0 >= ResE2)
		#define OUTE ResE0
	#elif(ResE1 > ResE0 && ResE1 >= ResE2)
		#define OUTE ResE1
	#elif(ResE2 > ResE0 && ResE2 > ResE1)
		#define OUTE ResE2
	#else
		#define OUTE 1
	#endif
#endif

//#if(ResE0 >= ResE1 && ResE0 >= ResE2)
//	#define OUTE ResE0
//#elif(ResE1 > ResE0 && ResE1 >= ResE2)
//	#define OUTE ResE1
//#elif(ResE2 > ResE0 && ResE2 > ResE1)
//	#define OUTE ResE2
//#else
//	#define OUTE 1
//#endif

#if ResM0 ==0 && ResM1==0 && ResM2==0
	#define OUTM 1
#else
	#if(ResM0 >= ResM1 && ResM0 >= ResM2)
		#define OUTM ResM0
	#elif(ResM1 > ResM0 && ResM1 >= ResM2)
		#define OUTM ResM1
	#elif(ResM2 > ResM0 && ResM2 > ResM1)
		#define OUTM ResM2
	#else
		#define OUTM 1
	#endif
#endif

//#define BIT W
using int_t = ap_uint<W>;
using MemPack = ap_uint<W*MAX_SIZE>;

//const unsigned Kadd = (1+(std::pow(2,EL)-1+ML)+(std::pow(2,ER)-1+MR)+WI);
const unsigned Kadd = (1+((1<<EL)-1+ML)+((1<<ER)-1+MR)+WI+WT);
#define EBIT 4

using KaddPack = ap_uint<Kadd*MAX_SIZE>;

const unsigned SEXP= 8;
typedef ap_int<SEXP> SEXP_T;
using SExpPack = ap_uint<SEXP*MB>;

#define INIT 101001100
//for st_round


#endif
