#include "bmgemm_hw.h"

void impl(){
	unsigned tileM = M/MAX_SIZE;
	unsigned tileN = N/MAX_SIZE;
	unsigned tileK = K/MAX_SIZE;
	ap_uint<2> mode, modeout;
	ap_uint<8> log_sizek = 4;
	ap_uint<8> log_sizek_tilec = 5;
	MemPack A[tileM*tileK*MAX_SIZE],B[tileK*tileN*MAX_SIZE],C[tileM*tileN*MAX_SIZE];
//	KaddPack D[tileM*tileN*MAX_SIZE];

	SExpPack betaA[BK_M*BK_K/MB], betaB[BK_K*BK_N/MB],betaC[BK_M*BK_K];

	bmGemmv25(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);

//	bmGemmv1(A,B,C, betaA, betaB, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);
}

