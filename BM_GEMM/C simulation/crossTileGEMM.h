// this is the GEMM for arithmetic test
#include "MAC.h"


#define TILE_K K/TileSize
#define TILE_N N/TileSize
#define TILE_M M/TileSize
#define Ebit 4


void SingleTileGEMM(ap_uint<W> A[SinM][SinK], ap_uint<W> B[SinK][SinN],SEXP_T betaA, SEXP_T betaB, SEXP_T &betaC, ap_uint<W> C[SinM][SinN],
				ap_uint<2> mode, ap_uint<2> modeOut){
	ap_int<Kadd> accum[SinM][SinN];
	ap_int<Kadd-WI-WT> mult;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
//			printf("-------[%d][%d]----------\n",i,j);
			for(int k=0;k<SinK;k++){
				if(k==0)
					accum[i][j] = 0;
				BMMul(A[i][k], B[k][j], mult, mode);
				ap_int<Kadd-WI> mul = (mult,ap_uint<WT>(0));
				accum[i][j] += mul;
//				printf("A: %s, B: %s, mult: %s,  accum: %s\n",A[i][k].to_string(2).c_str(),B[k][j].to_string(2).c_str(), mult.to_string(2).c_str(), accum[i][j].to_string(2).c_str());
//				printf("---------------\n");
			}
//			printf("---------------\n");
		}
	}

	//print accumulation result for verfication
//	printf("accum out: \n");
//	float acc;
//	int decimal = 7;//(biasA -1 + manA + biasB -1 + manB)
//	for(int i=0;i<SinM;i++){
//		for(int j=0;j<SinN;j++){
//			acc = float(accum[i][j])/std::pow(2,decimal)*std::pow(2,int(betaA+betaB));
//			printf("accum[%d][%d]: %f\t",i,j,acc);
//
////			printf("accum[%d][%d]: %s\n",i,j,Kulisch[i][j].to_string(2).c_str());
//		}
//		printf("\n");
//	}

	//check max & convert
	ap_uint<W> Z_min=(1<<W)-1;
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
			CheckZmin<W>(accum[i][j],Z_min);
//			printf("---------------\n");
		}
	}
//	printf("Zmin: %d, ",Z_min.to_int());

	//adjust bias
//	printf("=============================\n");
//	SEXP_T beta;
	betaC = betaA + betaB;
	ap_int<W> Zshift=0;
	BiasAdjust<W,W>(betaC, Z_min, Zshift,mode,modeOut);
//	printf("betaC: %d, Zshift: %d\n",betaC.to_int(), Zshift.to_int());
	for(int i=0;i<SinM;i++){
		for(int j=0;j<SinN;j++){
//			printf("------------[%d][%d]---------------\n",i,j);
			Normalization<W,W,4>(accum[i][j],Zshift,mode,modeOut,C[i][j]);
		}
	}
}

void CrossTileGEMM2(ap_uint<W> A[M][K], ap_uint<W> B[K][N], SEXP_T betaA[TILE_M][TILE_K],
		SEXP_T betaB[TILE_K][TILE_N], SEXP_T betaC[TILE_M][TILE_N], ap_uint<W> C[M][N],
		ap_uint<2> mode, ap_uint<2> modeOut){

	ap_int<Kadd> psum[TileSize][TileSize];
	ap_int<Kadd> accumulator[M][N];
	ap_int<Kadd-WI-WT> mult;

	ap_uint<4> ebias;
	ap_uint<1> flag;
	SEXP_T mul_beta;

//	for(int i=0;i<1;i++){
//		for(int j=0;j<1;j++){
	for(int i=0;i<TILE_M;i++){
		for(int j=0;j<TILE_N;j++){
			for(int ii=0;ii<TileSize;ii++){
				for(int jj=0;jj<TileSize;jj++){

//					printf("-------[%d][%d]----------\n",i*TileSize+ii,j*TileSize+jj);
					for(int k=0;k<TILE_K;k++){
						//shared exponent bias
						mul_beta = betaA[i][k] + betaB[k][j];
//						printf("betaA: %d\t, betaB: %d\t, i: %d, j: %d, k:%d,  ",betaA[i][k].to_int(),betaB[k][j].to_int(),i,j,k);
//						printf("betaC before: %d\t",betaC[i][j].to_int());
						GetEbeta2<Ebit>(ebias,flag,mul_beta,betaC[i][j],k);
//						printf("betaC: %d\n",betaC[i][j].to_int());

						//intra-block partial sum accumulation
						for(int kk=0;kk<TileSize;kk++){
							if(kk==0)
								psum[ii][jj] = 0;
							BMMul(A[i*TileSize+ii][k*TileSize+kk], B[k*TileSize+kk][j*TileSize+jj], mult, mode);
							#if WT>0
							ap_int<Kadd-WI> mul = (mult,ap_uint<WT>(0));//in-block accumulation
							#else
							ap_int<Kadd-WI> mul = mult;//in-block accumulation
							#endif
							psum[ii][jj] += mul;

							if(ii+i*TileSize==0&&jj+j*TileSize==0)
								printf("k:%d, mul: %s, psum: %s\n", k*TileSize+kk, mul.to_string(2).c_str(), psum[ii][jj].to_string(2).c_str());
						}

						//inter-block accumulation
						if(k==0)
							accumulator[ii+i*TileSize][jj+j*TileSize] = psum[ii][jj];
						else
							BMInterAccum<Ebit>(psum[ii][jj], accumulator[ii+i*TileSize][jj+j*TileSize],ebias,flag);

						if(ii+i*TileSize==0&&jj+j*TileSize==0)
							printf("---accum: %s\n", accumulator[ii+i*TileSize][jj+j*TileSize].to_string(2).c_str());

					}
				}
			}
		}
	}

	printf("accum out: \n");
	float accum;
	int expZeroPoint = 11;
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
//			printf("---------[%d][%d]--------------\n",i,j);
//			printf("%s,\t", accumulator[i][j].to_string(2).c_str());
			accum = float(accumulator[i][j])/std::pow(2,expZeroPoint)*std::pow(2,betaC[int(i/TileSize)][int(j/TileSize)].to_int());
			printf("accum[%d][%d]: %f\t",i,j,accum);

		}
		printf("\n");
	}

	//====================================check max========================================================
		ap_uint<8> Z_min[TILE_M][TILE_N];
		//initialization
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
				Z_min[i][j]=(1<<W)-1;
			}
		}
//		printf("Zmin: \n");
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
				for(int k1=0;k1<TileSize;k1++){
					for(int k2=0;k2<TileSize;k2++){
						CheckZmin<8>(accumulator[i*TileSize+k1][j*TileSize+k2],Z_min[i][j]);
						if(i==0 && j==0){
//							printf("Z_min: %d\n",Z_min[i][j].to_int());
//							printf("accum[%d][%d]: %s\n", i*TileSize+k1, j*TileSize+k2, accumulator[i*TileSize+k1][j*TileSize+k2].to_string(2).c_str());
						}
					}

				}
//				printf("%d\t",Z_min[i][j].to_int());
			}
//			printf("\n");
		}
		//================================Bias Adjust========================================================
		ap_int<W> Zshift[TILE_M][TILE_N];
//		printf("betaC: \n");
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
//				printf("i:%d,j:%d,Z_min: %d, Z_shift: %d, betaC before: %d\t",i,j,Z_min[i][j].to_int(), Zshift[i][j].to_int(), betaC[i][j].to_int());
				BiasAdjust<W,W>(betaC[i][j], Z_min[i][j], Zshift[i][j],mode,modeOut);
//				printf("after %d\n",betaC[i][j].to_int());
			}
//			printf("\n");
		}
		//=================================Normalization======================================================
			for(int i=0;i<TILE_M;i++){
				for(int j=0;j<TILE_N;j++){
					for(int k1=0;k1<TileSize;k1++){
						for(int k2=0;k2<TileSize;k2++){
//							printf("------------[%d][%d]---------------\n",i,j);
							Normalization<W,W,4>(accumulator[i*TileSize+k1][j*TileSize+k2],Zshift[i][j],mode,modeOut,C[i*TileSize+k1][j*TileSize+k2]);
						}
					}
				}
			}
}

void CrossTileGEMM(ap_uint<W> A[M][K], ap_uint<W> B[K][N], SEXP_T betaA[TILE_M][TILE_K],
		SEXP_T betaB[TILE_K][TILE_N], SEXP_T betaC[TILE_M][TILE_N], ap_uint<W> C[M][N],
		ap_uint<2> mode, ap_uint<2> modeOut){

//	printf("start\n");

	ap_int<Kadd> accumulator[M][N];
	ap_int<Kadd-WI-WT> mult;
	ap_uint<4> ebias;
	ap_uint<1> flag, k_flg;
	SEXP_T mul_beta;
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
//	for(int i=0;i<TileSize;i++){
//		for(int j=0;j<TileSize;j++){
//			printf("-------[%d][%d]----------\n",i,j);
			for(int k=0;k<K;k++){
				if(k==0)
					accumulator[i][j] = 0;
				if(k%TileSize==0){
					mul_beta = betaA[int(i/TileSize)][int(k/TileSize)] + betaB[int(k/TileSize) ][ int(j/TileSize)];
//					printf("betaA: %d\t betaB:%d\t, i:%d   j:%d   k:%d ",betaA[int(i/TileSize)][int(k/TileSize)].to_int(),
//							betaB[int(k/TileSize) ][ int(j/TileSize)].to_int(), int(i/TileSize), int(j/TileSize), int(k/TileSize));
					k_flg = k==0? 1:0;
//					printf("betaC before: %d\t",betaC[int(i/TileSize)][int(j/TileSize)].to_int());
					GetEbeta<Ebit>(ebias,flag,mul_beta,betaC[int(i/TileSize)][int(j/TileSize)],k_flg);
//					printf("betaC: %d\n",betaC[int(i/TileSize)][int(j/TileSize)].to_int());
				}
				BMMul(A[i][k], B[k][j], mult, mode);
				ap_uint<1> flagtile = k%TileSize == 0? ap_uint<1>(1) : ap_uint<1>(0);
				BMADD<Ebit>(mult, accumulator[i][j], ebias, flag, flagtile);
//				if(i<TileSize && j<TileSize)
//					printf("k: %d, accumulator[%d][%d]: %s\n",k,i,j,accumulator[i][j].to_string(2).c_str());
//				printf("-------------------------------------\n");
			}
		}
	}

	printf("accum out: \n");
	float accum;
	int expZeroPoint = 12;
	for(int i=0;i<M;i++){
		for(int j=0;j<N;j++){
//			printf("---------[%d][%d]--------------\n",i,j);
//			printf("%s,\t", accumulator[i][j].to_string(2).c_str());
			accum = float(accumulator[i][j])/std::pow(2,expZeroPoint)*std::pow(2,betaC[int(i/TileSize)][int(j/TileSize)].to_int());
			printf("accum[%d][%d]: %f\t",i,j,accum);

//			ap_uint<Kadd> xx = accumulator[i][j];
//			ap_uint<1> signR = xx[Kadd-1];
//			ap_uint<Kadd-1> sgn_rst = signR? (~ap_uint<Kadd-1>(xx(Kadd-2,0))+1) : xx(Kadd-2,0);
//			ap_uint<W> ZZ = (1<<W)-1;
//			CheckZmin<W>(sgn_rst,ZZ);
//			printf("accum[%d][%d]: %s\n",i,j,accumulator[i][j].to_string(2).c_str());
//			printf("sgn_rst Z: %d, ", ZZ.to_int());

		}
		printf("\n");
	}

	//====================================check max========================================================
		ap_uint<W> Z_min[TILE_M][TILE_N];
		//initialization
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
				Z_min[i][j]=(1<<W)-1;
			}
		}
//		printf("Zmin: \n");
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
				for(int k1=0;k1<TileSize;k1++){
					for(int k2=0;k2<TileSize;k2++){
						CheckZmin<W>(accumulator[i*TileSize+k1][j*TileSize+k2],Z_min[i][j]);
//						if(i==0 && j==0)
//							printf("Z_min: %d\n",Z_min[i][j].to_int());
					}

				}
//				printf("%d\t",Z_min[i][j].to_int());
			}
//			printf("\n");
		}
		//================================Bias Adjust========================================================
		ap_int<W> Zshift[TILE_M][TILE_N];
//		printf("betaC: \n");
		for(int i=0;i<TILE_M;i++){
			for(int j=0;j<TILE_N;j++){
//				printf("i:%d,j:%d,Z_min: %d, Z_shift: %d, betaC before: %d\t",i,j,Z_min[i][j].to_int(), Zshift[i][j].to_int(), betaC[i][j].to_int());
				BiasAdjust<W,W>(betaC[i][j], Z_min[i][j], Zshift[i][j],mode,modeOut);
//				printf("after %d\n",betaC[i][j].to_int());
			}
//			printf("\n");
		}
		//=================================Normalization======================================================
			for(int i=0;i<TILE_M;i++){
				for(int j=0;j<TILE_N;j++){
					for(int k1=0;k1<TileSize;k1++){
						for(int k2=0;k2<TileSize;k2++){
//							printf("------------[%d][%d]---------------\n",i,j);
							Normalization<W,W,4>(accumulator[i*TileSize+k1][j*TileSize+k2],Zshift[i][j],mode,modeOut,C[i*TileSize+k1][j*TileSize+k2]);
						}
					}
				}
			}

//			printf("C: \n");
//			for(int i=0;i<M;i++){
//				for(int j=0;j<N;j++){
//					printf("%s,\t",C[i][j].to_string(2).c_str());
//				}printf("\n");
//			}

}
