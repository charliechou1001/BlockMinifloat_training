#include "Train.h"
#include "netdata1.h"
//#include "Print.h"


void TrainTest(){

	XPack Input[DATA_NUM * (INSIZE + OTSIZE)];
	for(int t=0;t<DATA_NUM;t++){
		for(int ii=0;ii<INSIZE;ii++){
			for(int jj=0;jj<MAX_SIZE;jj++){
				Input[ii + t*(INSIZE + OTSIZE)].range(Ws*jj+Ws-1,Ws*jj) = In[t*B + jj+int(ii/LOOKBACK)*MAX_SIZE][ii%LOOKBACK];
	//			printf("%d, %d, %d\n",ii,jj+int(ii/LOOKBACK)*MAX_SIZE, ii%LOOKBACK);
			}
		}
		for(int ii=0;ii<OTSIZE;ii++){
			for(int jj=0;jj<MAX_SIZE;jj++){
				Input[ii + INSIZE + t*(INSIZE + OTSIZE) ].range(Ws*jj+Ws-1,Ws*jj) = label[t*B + jj+int(ii/FORECAST)*MAX_SIZE][ii%FORECAST];
	//			printf("%d, %d, %d | \t",ii,jj+int(ii/FORECAST)*MAX_SIZE, ii%FORECAST);
			}
	//		printf("\n");
		}
	}

	static MemPack Wgt[W_NUM];
	for(int ii=0;ii<W_NUM;ii++ ){
		for(int jj=0;jj<MAX_SIZE;jj++){
			Wgt[ii].range(W*jj+W-1,W*jj) = Weight[ii][jj];
		}
	}
//	printf("Wgt[0]:%s\n",Wgt[0].to_string(2).c_str());

	SExpPack betaX[DATA_NUM * B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE) + W_NUM/(MAX_SIZE*BLK_SIZE)];
	for(int ii=0;ii<DATA_NUM * B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE);ii++ ){
		for(int jj=0;jj<MB;jj++){
			betaX[ii].range(SEXP*jj+SEXP-1,jj*SEXP) = betaIn[ii][jj];
		}
	}
	for(int ii=0;ii<W_NUM/(MAX_SIZE*BLK_SIZE);ii++ ){
		for(int jj=0;jj<MB;jj++){
			betaX[ii+DATA_NUM * B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)].range(SEXP*jj+SEXP-1,jj*SEXP) = betaW[ii][jj];
		}
//		printf("betaX: %s \n", betaX[ii+B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)].to_string(2).c_str());
	}

	MemPack Act[ACT_NUM/MAX_SIZE];
	Train(Wgt, Act, Input, betaX);


}

void PrintTest(){

	ap_uint<W> W1[4][4] = {
	0b00101110, 0b01101010, 0b10011010, 0b01100000,
	0b10010101, 0b00001110, 0b01110100, 0b11110111,
	0b11010001, 0b10100000, 0b10110010, 0b01101111,
	0b11010011, 0b10111011, 0b11011001, 0b11111000,
	};
	ap_int<8> betaW1[2][2] = {
	-2, -2,
	-2, -2,
	};
	SExpPack betaX[2];
	for(int ii=0;ii<2;ii++ ){
		for(int jj=0;jj<MB;jj++){
			betaX[ii].range(SEXP*jj+SEXP-1,jj*SEXP) = betaW1[ii][jj];
		}
	}
	MemPack Wgt[4];
	for(int ii=0;ii<4;ii++ ){
		for(int jj=0;jj<MAX_SIZE;jj++){
			Wgt[ii].range(W*jj+W-1,W*jj) = W1[ii][jj];
		}
	}

	PrintMatrix(Wgt, betaX, 4*MAX_SIZE, 0, 0);

}

void PrintKaddTest(){

	ap_uint<Ws> W1[4][4] = {
	0b0010111000000000, 0b0110101000000000, 0b1001101000000000, 0b0110000000000000,
	0b1001010100000000, 0b0000111000000000, 0b0111010000000000, 0b1111011100000000,
	0b1101000100000000, 0b1010000000000000, 0b1011001000000000, 0b0110111100000000,
	0b1101001100000000, 0b1011101100000000, 0b1101100100000000, 0b1111100000000000,
	};
	ap_int<8> betaW1[2][2] = {
	-2, -2,
	-2, -2,
	};
	SExpPack betaX[2];
	for(int ii=0;ii<2;ii++ ){
		for(int jj=0;jj<MB;jj++){
			betaX[ii].range(SEXP*jj+SEXP-1,jj*SEXP) = betaW1[ii][jj];
		}
	}
	XPack Wgt[4];
	for(int ii=0;ii<4;ii++ ){
		for(int jj=0;jj<MAX_SIZE;jj++){
			Wgt[ii].range(Ws*jj+Ws-1,Ws*jj) = W1[ii][jj];
		}
	}

	PrintKadd( 4, Wgt, betaX );

}

void TestUnsigned(){
	ap_int<Kadd> mac_out = 0b1111111111111111111111111;
	ap_uint<Kadd-1> sgn_rst = GetUnsigned(mac_out);
	printf("sgn_rst: %s\n",sgn_rst.to_string(2).c_str() );
	mac_out = 0b1111;
	sgn_rst = GetUnsigned(mac_out);
	printf("sgn_rst: %s\n",sgn_rst.to_string(2).c_str() );
}

void Test_ap_int(){
	ap_uint<8> a = 0b11000011;//-67
	ap_int<8> b,c;//ap_int accept 2's complement code
	b = a;
	c = ap_int<8>(a)>>2;
	printf("b: %d, %s, c:%d\n",b.to_int(), b.to_string(2).c_str(), c.to_int() );

	using XXT = ap_fixed<8, 2, AP_RND, AP_SAT>;
	XXT aa=0, bb=0, cc=0, dd=0, ee=0;
	aa(7,0) = a(7,0);

	bb(6,0) = a(6,0);
	bb = -bb;
	ee = bb >> 2;

	ap_uint<8> sdf = (ap_uint<1>(1), (~a(6,0)+1) );
	cc = sdf;

	dd(6,0) = a(6,0);
	dd[7] = 1;
	printf("aa: %f, %s, bb: %f, %s\n, cc: %f, %s,sdf:%s,  ee: %f, %s\n dd:%f, %s\n",aa.to_float(), aa.to_string(2).c_str(),bb.to_float(), bb.to_string(2).c_str(),
			cc.to_float(), cc.to_string(2).c_str(), sdf.to_string(2).c_str(), ee.to_float(), ee.to_string(2).c_str(),
			dd.to_float(), dd.to_string(2).c_str() );
}

void Test_beta(){
	SEXP_T a = SEXP_T( (1<<(SEXP-1))-1 );
	SEXP_T b = SEXP_T( -(1<<(SEXP-1))+1 );
	ap_int<EBIT+1> c = ap_int<EBIT+1>( (1<<(EBIT))-1 );
	ap_int<EBIT+1> d = ap_int<EBIT+1>(  -(1<<(EBIT+1))+1 );
	printf("a: %d, b:%d, c:%d, d:%d \n",a.to_int(), b.to_int(), c.to_int(), d.to_int() );

}

void Test_idx(){
//	(k+1) & (size_k-1) )== 0

	int size_k = 4;
	for(int k=0;k<12;k++){
		printf("k: %d, %d,  %d \n", k, ((k+1) & (size_k-1)), (k & (size_k-1)) );

	}

}

//void MAPEPrint(XPack *Y, XPack *LABEL, SExpPack *betaY, SExpPack *betaL, int y_addr, int X_addr, int betaX_addr)

void Test_MAPEPrint(){

	XPack in1[OTSIZE], in2[OTSIZE], in3[OTSIZE];
	for(int ii=0;ii<OTSIZE;ii++){
		for(int jj=0;jj<MAX_SIZE;jj++){
			in1[ii].range(Ws*jj+Ws-1,Ws*jj) = label[10*B + jj+int(ii/FORECAST)*MAX_SIZE][ii%FORECAST];//label 11
			in2[ii].range(Ws*jj+Ws-1,Ws*jj) = label[12*B + jj+int(ii/FORECAST)*MAX_SIZE][ii%FORECAST];//label 13
			in3[ii].range(Ws*jj+Ws-1,Ws*jj) = label[14*B + jj+int(ii/FORECAST)*MAX_SIZE][ii%FORECAST];//label 15
		}
	}

	SExpPack betaX1[B*(FORECAST)/(MAX_SIZE*BLK_SIZE)];
	SExpPack betaX2[B*(FORECAST)/(MAX_SIZE*BLK_SIZE)];
	SExpPack betaX3[B*(FORECAST)/(MAX_SIZE*BLK_SIZE)];

	for(int ii=0;ii< B*(FORECAST)/(MAX_SIZE*BLK_SIZE);ii++ ){
		for(int jj=0;jj<MB;jj++){
			betaX1[ii].range(SEXP*jj+SEXP-1,jj*SEXP) =
					betaIn[B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)*10+B*(LOOKBACK)/(MAX_SIZE*BLK_SIZE)+ii ][jj];
			betaX2[ii].range(SEXP*jj+SEXP-1,jj*SEXP) =
					betaIn[B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)*12+B*(LOOKBACK)/(MAX_SIZE*BLK_SIZE)+ii ][jj];
			betaX3[ii].range(SEXP*jj+SEXP-1,jj*SEXP) =
					betaIn[B*(LOOKBACK + FORECAST)/(MAX_SIZE*BLK_SIZE)*14+B*(LOOKBACK)/(MAX_SIZE*BLK_SIZE)+ii ][jj];
		}
	}

	MAPEPrint( in1, in2, betaX1, betaX2, 0,0,0 );
	MAPEPrint( in1, in3, betaX1, betaX3, 0,0,0 );
	MAPEPrint( in2, in3, betaX2, betaX3, 0,0,0 );

}


void Test_st(){
	static ap_uint<9> state=1;
	ap_uint<1> a;
	for(int i=0;i<50;i++){
		a = st_round<4>(state);
		printf("state: %d, a: %d\n", state.to_int(), a.to_int() );
	}

}

int main(){
//	PrintKaddTest();
//	Test_MAPEPrint();
	TrainTest();
//	Test_st();

//	Test_beta();
//	Test_idx();
//	Test_ap_int();
//	TestUnsigned();
//	PrintTest();
}
