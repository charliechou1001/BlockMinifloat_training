# BlockMinifloat_training

__Abstract:__ Time series forecasting is the problem of predicting future data samples from historical information, and recent deep neural network (DNN) based techniques have achieved excellent results compared with conventional statistical approaches. Many applications at the  edge can utilise this technology, and most implementations have focused on inference; an ability to train at the edge would enable the DNN to adapt to changing conditions. Unfortunately, training requires approximately three times more memory and computation than inference. Moreover, edge applications are often constrained by energy efficiency. In this work, we implement a block minifloat (BM) training accelerator for a time series prediction network, N-BEATS. Our architecture involves a mixed precision GEMM accelerator that utilizes BM arithmetic. We use a 4-bit DSP packing scheme to optimize the implementation further, achieving a throughput of 779 Gops. The resulting power efficiency is 42.4 Gops/W, 3.1x better than a graphics processing unit in a similar technology. 

## 1.Directory
Each directory includes the "C simulation" and  "synthesis" sub-directory. "C simulation" is the code to verify computation correctness, and the "synthesis" is the code to generate the FPGA .bin file.

### BM_GEMM
It includes the BM GEMM kernel implementation that supports runtime configuration of precision. The top function is the function _bmGemmv25_ in bmgemm_hw.h.

### training_accelerator 
It includes the N-BEATS training implementation using the BM GEMM kernel. The top function is the function _Train_ in train.h.

## 2. Hardware specification
- Hardware Device: Xilinx Alveo U50
- Development Platform: xilinx_u50_gen3x16_xdma_201920_3
- Compiler: Vitis HLS  2021.2

## 3. I/O connectivity
sp=Train.WA: HBM[0]  
sp=Train.Act: HBM[1]  
sp=Train.X:HBM[2]  
sp=Train.betaX:HBM[3]  

