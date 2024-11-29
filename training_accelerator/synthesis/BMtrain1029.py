import CppHeaderParser

cH = CppHeaderParser.CppHeader("typedef.h")
# print(cH.defines[24])
# print(cH.defines[25])
# MAX_SIZE = int(cH.defines[24].split(" ")[-1])
# BLK_SIZE = int(cH.defines[25].split(" ")[-1])
MAX_SIZE = 32
BLK_SIZE = 16

filename = 'PEarray.h'
with open(filename, 'w') as f:
    f.write("#include \"ap_int.h\" \n")
    f.write("#include \"PE.h\" \n\n")
    f.write("void PEarray3(int k, int size_k, ap_uint<2> flag, ap_uint<2> beta_flag, int ir, int ic, \n \
                  int_t feederA[MAX_SIZE], int_t feederB[MAX_SIZE], ap_uint<2> mode, ap_int<Kadd> sReg_1[MAX_SIZE], \n \
                  SEXP_T betaAbf[MB][BETASIZE], SEXP_T betaBbf[MB][BETASIZE], SEXP_T betasReg_1[MB][MB]){\n")
    f.write("#pragma HLS INLINE\n\
#pragma HLS ARRAY_PARTITION variable=feederA dim=1 type=complete \n \
#pragma HLS ARRAY_PARTITION variable=feederB dim=1 type=complete \n \
#pragma HLS ARRAY_PARTITION variable=betaAbf dim=1 type=complete \n \
#pragma HLS ARRAY_PARTITION variable=betaBbf dim=1 type=complete \n \
#pragma HLS ARRAY_PARTITION variable=sReg_1 dim=1 type=complete \n \
#pragma HLS ARRAY_PARTITION variable=betasReg_1 dim=1 type=complete\n\n")
            
    f.write("static ap_uint<EBIT> ebias[MB][MB];\n\
#pragma HLS ARRAY_PARTITION variable=ebias dim=0 complete \n\
static ap_uint<1> flg[MB][MB]; \n\
#pragma HLS ARRAY_PARTITION variable=flg dim=0 complete \n\
static SEXP_T betaRes[MB][MB], betasReg[MB][MB]; \n\
#pragma HLS ARRAY_PARTITION variable=betaRes dim=0 complete \n\
#pragma HLS ARRAY_PARTITION variable=betasReg dim=0 complete\n\n\
static ap_int<Kadd> sReg[MAX_SIZE][MAX_SIZE];\n\
#pragma HLS ARRAY_PARTITION variable = sReg dim = 0 complete\n\
    static ap_int<Kadd> out[MAX_SIZE][MAX_SIZE];\n\
#pragma HLS ARRAY_PARTITION variable = out dim = 0 complete\n\
\n")
    
    f.write("    int TILE_K = size_k/MAX_SIZE;\n\
    int kk =  k & (size_k-1);//k%(size_k-1)\n\
    bool k_blk_start =  ( int(k/BLK_SIZE) & (size_k/BLK_SIZE-1) )==0;\n\
    bool k_end = ( (k+1) & (size_k-1) ) ==0;\n\n")

##old version
    f.write("    //---------------shared exp----------------\n\
    betax:for(int ip=0;ip<MB;ip++){\n\
	#pragma HLS UNROLL\n\
	betay:for(int iq=MB-1;iq>=0;iq--){\n\
		#pragma HLS UNROLL\n\
		if(iq==MB-1)\n\
			beta_SRL_2(k, k_blk_start, betaAbf[ip][ir*TILE_K*MB+int(kk/BLK_SIZE)], betaBbf[ip][ic*TILE_K*MB+int(kk/BLK_SIZE)],\n\
					betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg_1[ip], beta_flag[1]*( ( (k+1) & (size_k-1) ) -1) );\n\
		else\n\
			beta_SRL(k, k_blk_start, betaAbf[ip][ir*TILE_K*MB+int(kk/BLK_SIZE)], betaBbf[ip][ic*TILE_K*MB+int(kk/BLK_SIZE)],\n\
					betaRes[ip][iq], ebias[ip][iq], flg[ip][iq], beta_flag, betasReg[ip][iq], betasReg[ip][iq+1] );\n\
	}\n\
   }\n\n")

#new version
 #    f.write("    //---------------shared exp----------------\n\
 #    betax:for(int ip=0;ip<MB;ip++){\n\
	# #pragma HLS UNROLL\n\
 #    	betay:for(int iq=MB-1;iq>=0;iq--){\n\
 #        #pragma HLS UNROLL\n\
 #            beta(k, k_blk_start, SEXP_T(betaA[ir*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*ip+SEXP-1,SEXP*ip)),\n\
 #            SEXP_T(betaB[ic*TILE_K*MB+int(kk/BLK_SIZE)].range(SEXP*iq+SEXP-1,SEXP*iq)),\n\
 #            betaRes[ip][iq], ebias[ip][iq], flg[ip][iq]);\n\
 #        }\n\
 #    }\n\n")

 #    f.write("\n\
 #    for(int i=0;i<MB;i++){\n\
	# #pragma HLS UNROLL\n\
 #    	for(int j=MB-1;j>=0;j--){\n\
	#     #pragma HLS UNROLL\n\
 #    	    if(j!=MB-1){\n\
 #                if(beta_flag==0b01){//load data\n\
 #        	    betasShift[i][j] = betaRes[i][j];\n\
 #    		}else if(beta_flag==0b10)//shift data\n\
 #        	    betasShift[i][j+1] = betasShift[i][j];\n\
 #    	    }else{\n\
 #        	if(beta_flag==0b01)\n\
 #        	     betasShift[i][j] = betaRes[i][j];\n\
 #        	else if(beta_flag==0b10)\n\
 #        	    betasReg_1[i][beta_flag[1]*( ( (k+1) & (size_k-1) ) -1)] = betasShift[i][MB-1];\n\
 #    	    }\n\
 #    	}\n\
 #    }\n\n")


    f.write("    //------------------PE array-----------------\n")
    for i in range(MAX_SIZE):
        for j in range(MAX_SIZE-1,-1,-1):
            f.write("    out[%d][%d] = PE3<%d,%d>(kk, k_end, mode, feederA[%d], feederB[%d], ebias[%d][%d], flg[%d][%d]);\n"
                    %(i,j, i,j, i,j, int(i/BLK_SIZE),int(j/BLK_SIZE), int(i/BLK_SIZE),int(j/BLK_SIZE)))
    
    f.write("\n    for(int i=0;i<MAX_SIZE;i++){\n\
		#pragma HLS UNROLL\n\
    	for(int j=MAX_SIZE-1;j>=0;j--){\n\
			#pragma HLS UNROLL\n\
    		if(j!=MAX_SIZE-1){\n\
        		if(flag==0b01)\n\
        			sReg[i][j] = out[i][j];\n\
        		else if(flag==0b10)\n\
        			sReg[i][j+1] = sReg[i][j];\n\
    		}else{\n\
        		if(flag==0b01)\n\
        			sReg[i][j] = out[i][j];\n\
        		else if(flag==0b10)\n\
        			sReg_1[i] = sReg[i][MAX_SIZE-1];\n\
    		}\n\
    	}\n\
    }\n")


    f.write("}\n")#end
    
