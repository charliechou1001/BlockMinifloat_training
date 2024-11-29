#ifndef TYPEDEF_
#define TYPEDEF_

#include "ap_int.h"

//W = S0L+E0L+M0L = S1L+E1L+M1L = S0R+E0R+M0R = S1R+E1R+M1R
#define W 8
//const unsigned W=8;
#define WO 16
#define Ws 16
//bit-length of data, 8 or 4 bit
//#define WI 5
//#define WT 3
#define WI 5
#define WT 3
//extra bit for accumulator to avoid overflow

//left input mode 0(FW activation)
#define S0L 1
#define E0L 0
#define M0L 7

//left input mode 1(error)
#define S1L 1
#define E1L 0
#define M1L 7

//right input mode 0(weight)
#define S0R 1
#define E0R 0
#define M0R 7

//right input mode 1(BP activation for gradient)
#define S1R 1
#define E1R 0
#define M1R 7

//output mode 0(forward)
#define ResS0 1
#define ResE0 0
#define ResM0 7

//output mode 1(error)
#define ResS1 1
#define ResE1 0
#define ResM1 7

//output mode 2(gradient)
#define ResS2 1
#define ResE2 0
#define ResM2 7

//output mode 3(block out)
#define ResS3 0
#define ResM3 15

//#define MAX_SIZE 4
#define MAX_SIZE 4
const unsigned LOG_MAX_SIZE = 2;

const unsigned BLK_SIZE = 2;
const unsigned MB = MAX_SIZE/BLK_SIZE;

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

//-----------sgd-------------
#if E0R ==0
#define ETAW 0
#else
#define ETAW ( (1<<(E0R-1)) -1 )
#endif

#if ResE2 ==0
#define ETAG 0
#else
#define ETAG ( (1<<(ResE2-1)) -1 )
#endif

#if E0R >= ResE2//extended addend input integer part
#define ADINT ( (1<<E0R)-1 - ETAW + 1 )
#else
#define ADINT ( (1<<ResE2)-1 - ETAG + 1 )
#endif

#define LEW ( (1<<E0R)-1 - ETAW + 1 )
#define LE ( ADINT + 1 )

//extended addend input decimal part
#if ( M0R + ETAW -1) > ( ResM2 + ETAG -1)
#define ADDEC ( M0R + ETAW -1)
#else
#define ADDEC ( ResM2 + ETAG -1)
#endif

using ADIN_T = ap_fixed<ADINT+ADDEC+1,ADINT+1,AP_RND,AP_SAT>;
using ADOUT_T = ap_fixed<ADINT+ADDEC+2,ADINT+2,AP_RND,AP_SAT>;

using int_t = ap_uint<W>;
using MemPack = ap_uint<W*MAX_SIZE>;
using BlkPack = ap_uint<WO*MAX_SIZE>;
using XPack = ap_uint<Ws*MAX_SIZE>;

const unsigned BANDWIDTH = MAX_SIZE/2;//BANDWIDTH < MAX_SIZE
const unsigned HALF_SIZE = MAX_SIZE/2;
using X_T = ap_uint<Ws*BANDWIDTH>;
const unsigned MM = BANDWIDTH/BLK_SIZE;

const unsigned Kadd = (1+((1<<EL)-1+ML)+((1<<ER)-1+MR)+WI+WT);
#define EBIT 8
#define ALPHA 3//learning rate, equivalent to 1/8 = 0.125

using KaddPack = ap_uint<Kadd*MAX_SIZE>;
using accum_t = ap_int<Kadd>;

using AccMemPack = ap_uint<Kadd*HALF_SIZE>;


const unsigned SEXP= 8;
typedef ap_int<SEXP> SEXP_T;
using SExpPack = ap_uint<SEXP*MB>;
using SExpIn = ap_uint<SEXP*BANDWIDTH/BLK_SIZE>;
const unsigned BETASIZE = 8;

#define INIT 101001100

#define EPOCH 16*80

#endif
