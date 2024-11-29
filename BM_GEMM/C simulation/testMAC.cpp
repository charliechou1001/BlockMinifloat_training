#include "ap_int.h"
//#include "test_data.h"
#include "pydata.h"
#include "crossTileGEMM.h"
//#include "MAC.h"
#include "test_longBMtofp32.h"

//void mult_test(){
//	ap_uint<W> a;
//	ap_uint<W> b;
//	ap_uint<2> mode;
//	ap_int<Kadd-WI> mult1;
//
////	a = 0b10011011;//(2,5)denorm,(3,4)
////	b = 0b01010010;//(1,6), (4,3)
//
//	a = 0b10011011;//(2,5)denorm,signed BFP
//	b = 0b11110010;//signed BFP, (4,3)
//
//	mode = 0;
//	BMMul(a,b,mult1,mode);
//	printf("--------------------------\n");
//	mode = 1;
//	BMMul(a,b,mult1,mode);
//	printf("--------------------------\n");
//	mode = 0b10;
//	BMMul(a,b,mult1,mode);
//	printf("--------------------------\n");
//	mode = 0b11;
//	BMMul(a,b,mult1,mode);
//}

void crossGEMM2_test(){
	ap_uint<2> mode, modeOut;
	ap_uint<W> C[M][N];
	SEXP_T betaC[TILE_M][TILE_N];

	mode = 0b11;//left, right ==0
	modeOut = 0b0;
	CrossTileGEMM2(CrA30, CrB30, betaCrA30, betaCrB30, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_3_0, betaC, betaCrC_3_0, modeOut);
	printf("------------------------\n");
	modeOut = 0b1;
	CrossTileGEMM2(CrA31, CrB31, betaCrA31, betaCrB31, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_3_1, betaC, betaCrC_3_1, modeOut);
	printf("------------------------\n");
//	modeOut = 0b10;
//	CrossTileGEMM2(CrA32, CrB32, betaCrA32, betaCrB32, betaC, C, mode, modeOut);//
//	CrossGEMMCompare(C, CrC_3_2, betaC, betaCrC_3_2, modeOut);
//	printf("=========================\n");
//
//	mode = 0b10;//left, right ==0
//	modeOut = 0b0;
//	CrossTileGEMM2(CrA20, CrB20, betaCrA20, betaCrB20, betaC, C, mode, modeOut);//
//	CrossGEMMCompare(C, CrC_2_0, betaC, betaCrC_2_0, modeOut);
//	printf("------------------------\n");
//	modeOut = 0b1;
//	CrossTileGEMM2(CrA21, CrB21, betaCrA21, betaCrB21, betaC, C, mode, modeOut);//
//	CrossGEMMCompare(C, CrC_2_1, betaC, betaCrC_2_1, modeOut);
//	printf("------------------------\n");
//	modeOut = 0b10;
//	CrossTileGEMM2(CrA22, CrB22, betaCrA22, betaCrB22, betaC, C, mode, modeOut);//
//	CrossGEMMCompare(C, CrC_2_2, betaC, betaCrC_2_2, modeOut);
//	printf("=========================\n");

}

void crossGEMM_test(){
	ap_uint<2> mode, modeOut;
	ap_uint<W> C[M][N];
	SEXP_T betaC[TILE_M][TILE_N];

	mode = 0b11;//left, right ==0
	modeOut = 0b0;
	CrossTileGEMM(CrA30, CrB30, betaCrA30, betaCrB30, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_3_0, betaC, betaCrC_3_0, modeOut);
	printf("------------------------\n");
	modeOut = 0b1;
	CrossTileGEMM(CrA31, CrB31, betaCrA31, betaCrB31, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_3_1, betaC, betaCrC_3_1, modeOut);
	printf("------------------------\n");
	modeOut = 0b10;
	CrossTileGEMM(CrA32, CrB32, betaCrA32, betaCrB32, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_3_2, betaC, betaCrC_3_2, modeOut);
	printf("=========================\n");

	mode = 0b10;//left, right ==0
	modeOut = 0b0;
	CrossTileGEMM(CrA20, CrB20, betaCrA20, betaCrB20, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_2_0, betaC, betaCrC_2_0, modeOut);
	printf("------------------------\n");
	modeOut = 0b1;
	CrossTileGEMM(CrA21, CrB21, betaCrA21, betaCrB21, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_2_1, betaC, betaCrC_2_1, modeOut);
	printf("------------------------\n");
	modeOut = 0b10;
	CrossTileGEMM(CrA22, CrB22, betaCrA22, betaCrB22, betaC, C, mode, modeOut);//
	CrossGEMMCompare(C, CrC_2_2, betaC, betaCrC_2_2, modeOut);
	printf("=========================\n");
}

void singleGEMM_test(){
	SEXP_T betaC;
	ap_uint<2> mode;
	ap_uint<W> C[SinM][SinN];
	ap_uint<2> modeOut;
//	(1,2,5) (1,0,7) || (1,4,3) (1,0,7)
	mode=0;
	modeOut=0b00;//BM
	SingleTileGEMM(SA00,SB00,betaSA00,betaSB00,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_0_0, betaC, betaSC_0_0, modeOut);
	printf("------------------------\n");

	modeOut=0b01;//BFP
	SingleTileGEMM(SA01,SB01,betaSA01,betaSB01,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_0_1, betaC, betaSC_0_1, modeOut);
	printf("------------------------\n");

	modeOut=0b10;//unsigned
	SingleTileGEMM(SA02,SB02,betaSA02,betaSB02,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_0_2, betaC, betaSC_0_2, modeOut);
	printf("=========================\n");

	mode = 0b1;
	modeOut=0b00;
	SingleTileGEMM(SA10,SB10,betaSA10,betaSB10,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_1_0, betaC, betaSC_1_0, modeOut);
	printf("------------------------\n");

	modeOut=0b1;
	SingleTileGEMM(SA11,SB11,betaSA11,betaSB11,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_1_1, betaC, betaSC_1_1, modeOut);
	printf("------------------------\n");

	modeOut=0b10;
	SingleTileGEMM(SA12,SB12,betaSA12,betaSB12,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_1_2, betaC, betaSC_1_2, modeOut);
	printf("=======================\n");

	mode = 0b10;
	modeOut=0b0;
	SingleTileGEMM(SA20,SB20,betaSA20,betaSB20,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_2_0, betaC, betaSC_2_0, modeOut);
	printf("------------------------\n");

	modeOut=0b1;
	SingleTileGEMM(SA21,SB21,betaSA21,betaSB21,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_2_1, betaC, betaSC_2_1, modeOut);
	printf("------------------------\n");

	modeOut=0b10;
	SingleTileGEMM(SA22,SB22,betaSA22,betaSB22,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_2_2, betaC, betaSC_2_2, modeOut);
	printf("=======================\n");

	mode = 0b11;
	modeOut=0b0;
	SingleTileGEMM(SA30,SB30,betaSA30,betaSB30,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_3_0, betaC, betaSC_3_0, modeOut);
	printf("------------------------\n");

	modeOut=0b1;
	SingleTileGEMM(SA31,SB31,betaSA31,betaSB31,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_3_1, betaC, betaSC_3_1, modeOut);
	printf("------------------------\n");

	modeOut=0b10;
	SingleTileGEMM(SA32,SB32,betaSA32,betaSB32,betaC,C,mode,modeOut);
	SingleMatrixCompare(C, SC_3_2, betaC, betaSC_3_2, modeOut);
	printf("========================\n");
}

//void PrintFPdata(){
//	ap_uint<2> mode, modeOut;
//
//	printf("SA30: \n");
//	PrintFPInput(SA30, betaSA30, 3);
//
//	printf("SB30: \n");
//	PrintFPInput(SB30, betaSB30, 1);
//
//	printf("SC_3_0: \n");
//	PrintFPOutput(SC_3_0, betaSC_3_0, 0);
//
//}

int main(){
//	singleGEMM_test();
//	crossGEMM_test();
	crossGEMM2_test();

//	ap_uint<W> aa=0b11111111;
//	BMToFP(aa, 0, -4);

//	PrintFPdata();
}
