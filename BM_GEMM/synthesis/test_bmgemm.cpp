#include "bmgemm_hw.h"
#include "test_longBMtofp32.h"
//#include "pydata.h"//sued for 16*16*16 gemm
#include "pydata1.h"

unsigned tileM = M/MAX_SIZE;
unsigned tileN = N/MAX_SIZE;
unsigned tileK = K/MAX_SIZE;

void Test2(){
		ap_uint<2> mode, modeout;
		MemPack A[tileM*tileK*MAX_SIZE],B[tileK*tileN*MAX_SIZE],C[tileM*tileN*MAX_SIZE];

		mode = 0b11;
		modeout = 0b0;
		for(int i=0;i<tileM*tileK*MAX_SIZE;i++){
			for(int j=0;j<MAX_SIZE;j++){
				A[i].range(W*j+W-1,W*j) = CTA30[j+ int(i/(tileK*MAX_SIZE))*MAX_SIZE][i%(tileK*MAX_SIZE)];
//				if(i==0)
//					printf("%s, ",CTA30[j+ int(i/(tileK*MAX_SIZE))*MAX_SIZE][i%(tileK*MAX_SIZE)].to_string(2).c_str());
			}
		}
//		printf("\n");

		for(int k=0;k<tileK*tileN*MAX_SIZE;k++){
			for(int j=0;j<MAX_SIZE;j++){
				B[k].range(W*j+W-1,W*j) = CTB30[k%(tileK*MAX_SIZE)][j+int(k/(tileK*MAX_SIZE))*MAX_SIZE];
//				if(k==0)
//					printf("%s, ", CTB30[k%(tileK*MAX_SIZE)][j+int(k/(tileK*MAX_SIZE))*MAX_SIZE].to_string(2).c_str());
			}
		}
//		printf("\n");

		SExpPack betaA[BK_M*BK_K/MB], betaB[BK_K*BK_N/MB],betaC[BK_M*BK_K];
		for(int i=0;i<BK_M*BK_K/MB;i++){
			for(int j=0;j<MB;j++){
				betaA[i].range(SEXP*j+SEXP-1,j*SEXP) =  betaCTA30[j+ int(i/(BK_K))*MB][i%(BK_K)];
//				printf("%d,  %d,  %d, ||\t",j+ int(i/(BK_K))*BLK_SIZE,i%(BK_K),betaA[i].range(SEXP*j+SEXP-1,j*SEXP).to_int());
			}
//			printf("\n");
		}
//		printf("\n");
		for(int k=0;k<BK_K*BK_N/MB;k++){
			for(int j=0;j<MB;j++){
				betaB[k].range(SEXP*j+SEXP-1,j*SEXP) = betaCTB30[k%(BK_K)][j+int(k/(BK_K))*MB];
//				printf("betaB: %d %d   %d,||\t ",k%(BK_K), j+int(k/(BK_K))*MB, betaB[k].range(SEXP*j+SEXP-1,j*SEXP).to_int());
			}
//			printf("\n");
		}


		ap_uint<W> resultC[M][N];
		SEXP_T betaRC[BK_M][BK_N];
		printf("===A*B===\n");
//		printf("-----------------------------------------------------------\n");
//		bmGemmv2(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);

		ap_uint<8> log_sizek = 0;
		int kk =K;
		while(kk!=1){
			log_sizek++;
			kk = kk/2;
		}//log2(size_k)

		ap_uint<8> log_sizek_tilec = 0;
		kk = K*N/MAX_SIZE;
		while(kk!=1){
			log_sizek_tilec++;
			kk = kk/2;
		}
//		printf(" log_sizek:%d, log_sizek_tilec:%d\n ",log_sizek.to_int(), log_sizek_tilec.to_int());
//		bmGemmv23(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);
//		printf("==============\n");

		bmGemmv25(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);// DON"T USE LOG FUNCTION!!!!! THE RESULT WILL BE WRONG!

//		bmGemmv3(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);

//		printf("C:\n");
		for(int i=0;i<tileM*tileN*MAX_SIZE;i++){
			for(int j=0;j<MAX_SIZE;j++){
				resultC[j+ int(i/(tileN*MAX_SIZE))*MAX_SIZE][i%(tileN*MAX_SIZE)] = C[i].range(W*j+W-1,W*j);
//				printf("%s, ",resultC[j+ int(i/(tileN*MAX_SIZE))*MAX_SIZE][i%(tileN*MAX_SIZE)].to_string(2).c_str());
			}
//			printf("\n");
		}
//		printf("betaC: \n");
		for(int i=0;i<BK_M*BK_N/MB;i++){
			for(int j=0;j<MB;j++){
				betaRC[j+ int(i/(BK_N))*MB][i%(BK_N)] = betaC[i].range(SEXP*j+SEXP-1,j*SEXP);
//				printf("%d,  %d,  %d, ||\t",j+ int(i/(BK_N))*BLK_SIZE,i%(BK_N),SEXP_T(betaRC[j+ int(i/(BK_N))*MB][i%(BK_N)] ).to_int());
			}
//			printf("\n");
		}

		CrossGEMMCompare(resultC,CTC1, betaRC, betaCTC1, modeout);
}

void Test_blank(){//to find any seg fault problem
	ap_uint<2> mode, modeout;
	MemPack A[tileM*tileK*MAX_SIZE],B[tileK*tileN*MAX_SIZE],C[tileM*tileN*MAX_SIZE];
	SExpPack betaA[BK_M*BK_K/MB], betaB[BK_K*BK_N/MB],betaC[BK_M*BK_K];
	mode = 0b11;
	modeout = 0b0;

	ap_uint<8> log_sizek = 0;
	int kk =K;
	while(kk!=1){
		log_sizek++;
		kk = kk/2;
	}//log2(size_k)

	ap_uint<8> log_sizek_tilec = 0;
	kk = K*N/MAX_SIZE;
	while(kk!=1){
		log_sizek_tilec++;
		kk = kk/2;
	}

	bmGemmv25(A,B,C, betaA, betaB, betaC, M,N,K,log_sizek, log_sizek_tilec, mode, modeout);

}

void TestPE(){
	ap_int<W> a[16] ={0b10000000, 0b10000111, 0b00001111, 0b10001111, 0b10100101, 0b10000000, 0b10010101, 0b10000000,
					0b10101001, 0b10000101, 0b00000011, 0b00000001, 0b01011100, 0b10111011, 0b10000111, 0b01001101};
	SEXP_T betaA[4] = { -1, 0, -1, -1};
	ap_int<W> b[16] = {0b10000000, 0b10000000, 0b00010110, 0b10001010, 0b00101000, 0b10000010, 0b00000101, 0b10000000,
						0b11100000, 0b00000111, 0b10010110, 0b11000111, 0b10101000, 0b11000010, 0b10010000, 0b00100010};
	SEXP_T betaB[4] = {-1, -1, -1, -1};

	ap_uint<2> mode = 0b11;
	ap_uint<EBIT> ebias; ap_uint<1> flg;SEXP_T betaAB, betaC;
	ap_int<Kadd> pSum, accum, out;
	ap_uint<2> flag;
	ap_int<Kadd> sReg1, sReg2;
	SEXP_T beta_sReg1, beta_sReg2;
	for(unsigned i=0;i<16;i++){
//		printf("i: %d\n",i);
		bool k_end = i==15? 1: 0;
		bool k_blk_start =  ( int(i/BLK_SIZE) & (16/BLK_SIZE-1) )==0;
		beta_SRL(i, k_blk_start, betaA[int(i/BLK_SIZE)], betaB[int(i/BLK_SIZE)], betaC, ebias, flg, flag, beta_sReg1, beta_sReg2 );

//		if(i%BLK_SIZE==0){
//			betaAB = betaA[int(i/BLK_SIZE)] + betaB[int(i/BLK_SIZE)];
//			printf("betaC before: %d,  ",betaC.to_int());
//			GetEbeta2<EBIT>( ebias, flg, betaAB, betaC, i );
//			printf("betaAB: %d, betaC: %d, ebias: %d\n", betaAB.to_int(), betaC.to_int(), ebias.to_int());
//		}
		PE2_SRL(i, k_end, mode, a[i], b[i], ebias, flg, pSum, accum, out, flag, sReg1, sReg2);
	}

	printf("out: %s\n",out.to_string(2).c_str());

}



int main(){
//	Test_blank();
	Test2();
//	TestPE();
	return 0;
}
