#ifndef CLZDEF_
#define CLZDEF_

#include "ap_int.h"

ap_uint<2> enc(ap_uint<2> a){
#pragma HLS INLINE
	ap_uint<2> out;
	out[1] = (~a[1])&(~a[0]);
	out[0] = (~a[1])&a[0];
	return out;
}

//len is original bit length
template<int len>
ap_uint<len> clzEncode(ap_uint<len> v){
#pragma HLS INLINE
	ap_uint<len> a;
	EncodeLoop:for(int i=0;i<len/2;i++){
#pragma HLS UNROLL
		a(2*i+1,2*i) = enc(v(2*i+1,2*i));
	}
	return a;
}

//n=i+2
//template<int n>
//ap_uint<n+1> clzi(ap_uint<2*n> v){
//	if(v[2*n-1]==0)
//		return (ap_uint<1>(v[2*n-1]&v[n-1]),ap_uint<1>(0),ap_uint<n-1>(v(2*n-2,n)));
//	else
//		return (ap_uint<1>(v[2*n-1]&v[n-1]),ap_uint<1>(~v[n-1]),ap_uint<n-1>(v(n-2,0)));
//}

//n=i+2
template<int n>
ap_uint<n+1> clzi(ap_uint<2*n> v){
#pragma HLS INLINE
	ap_uint<n+1> out;
	out[n] = ap_uint<1>(v[2*n-1]&v[n-1]);
	out(n-1,0) = (v[2*n-1]==0)? (ap_uint<1>(0),ap_uint<n-1>(v(2*n-2,n))) : (ap_uint<1>(~v[n-1]),ap_uint<n-1>(v(n-2,0)));
	return out;
}

//P = len/(2^(i+1)), P is pair number, len is original bit length
template<int i, int P>
ap_uint<(i+3)*P/2> clzOneloop(ap_uint<(i+2)*P> v){
#pragma HLS INLINE
	ap_uint<(i+3)*P/2> a;
	Oneloop:for(int k=0;k<P/2;k++){
#pragma HLS UNROLL
//		printf("idx: %d\t%d\t",2*(i+2)*k+2*i+3,2*(i+2)*k);
//		printf("%s\t",v(2*(i+2)*k+2*i+3,2*(i+2)*k).to_string(2).c_str());
//		printf("idx a:%d\t%d\n",(i+3)*k+i+2,(i+3)*k);
		a((i+3)*k+i+2,(i+3)*k) = clzi<i+2>(v(2*(i+2)*k+2*i+3,2*(i+2)*k));
//		printf("a: %s\n",a((i+3)*k+i+2,(i+3)*k).to_string(2).c_str());
	}
	return a;
}

//template<int len=8>
ap_uint<4> CLZ8(ap_uint<8> v){
#pragma HLS INLINE
	const int len=8;
	//encode
	ap_uint<len> e;
	e= clzEncode<len>(v);

	const int i1=0;
	const int i2=1;
	const int i3=2;
	const int P1=len>>1;//P1/2
	const int P2=len>>2;//P1/4
	const int P3=len>>3;
	ap_uint<(i2+2)*P2> clz1 = clzOneloop<i1,P1>(e);
//	printf("clz1:%s\n",clz1.to_string(2).c_str());
	ap_uint<(i3+2)*P3> clz2 = clzOneloop<i2,P2>(clz1);
//	printf("clz2:%s\n",clz2.to_string(2).c_str());
	return clz2;
}

//16bit CLZ,len=16
//template<int len=16>
ap_uint<5> CLZ16(ap_uint<16> v){
#pragma HLS INLINE
	const int len=16;
	//encode
	ap_uint<len> e;
	e= clzEncode<len>(v);

	const int i1=0;
	const int i2=1;
	const int i3=2;
	const int i4=3;
	const int P1=len>>1;//P1/2
	const int P2=len>>2;//P1/4
	const int P3=len>>3;//P1/8
	const int P4=len>>4;
	ap_uint<(i2+2)*P2> clz1 = clzOneloop<i1,P1>(e);
//	printf("clz1:%s\n",clz1.to_string(2).c_str());
	ap_uint<(i3+2)*P3> clz2 = clzOneloop<i2,P2>(clz1);
//	printf("clz2:%s\n",clz2.to_string(2).c_str());
	ap_uint<(i4+2)*P4> clz3 = clzOneloop<i3,P3>(clz2);
//	printf("clz3:%s\n",clz3.to_string(2).c_str());
	return clz3;
}

ap_uint<5> CLZ32(ap_uint<32> v){
#pragma HLS INLINE
	const int len=32;
	//encode
	ap_uint<len> e;
	e= clzEncode<len>(v);

	const int i1=0;
	const int i2=1;
	const int i3=2;
	const int i4=3;
	const int i5=4;
	const int P1=len>>1;//P1/2
	const int P2=len>>2;//P1/4
	const int P3=len>>3;//P1/8
	const int P4=len>>4;
	const int P5=len>>5;
	ap_uint<(i2+2)*P2> clz1 = clzOneloop<i1,P1>(e);
//	printf("clz1:%s\n",clz1.to_string(2).c_str());
	ap_uint<(i3+2)*P3> clz2 = clzOneloop<i2,P2>(clz1);
//	printf("clz2:%s\n",clz2.to_string(2).c_str());
	ap_uint<(i4+2)*P4> clz3 = clzOneloop<i3,P3>(clz2);
//	printf("clz3:%s\n",clz3.to_string(2).c_str());
	ap_uint<(i5+2)*P5> clz4 = clzOneloop<i4,P4>(clz3);
//	printf("clz4:%s\n",clz4.to_string(2).c_str());
	return clz4;
}

ap_uint<8> CLZ64(ap_uint<64> v){
#pragma HLS INLINE
	const int len=64;
	//encode
	ap_uint<len> e;
	e= clzEncode<len>(v);

	const int i1=0;
	const int i2=1;
	const int i3=2;
	const int i4=3;
	const int i5=4;
	const int i6=5;
	const int P1=len>>1;//P1/2
	const int P2=len>>2;//P1/4
	const int P3=len>>3;//P1/8
	const int P4=len>>4;
	const int P5=len>>5;
	const int P6=len>>6;
	ap_uint<(i2+2)*P2> clz1 = clzOneloop<i1,P1>(e);
//	printf("clz1:%s\n",clz1.to_string(2).c_str());
	ap_uint<(i3+2)*P3> clz2 = clzOneloop<i2,P2>(clz1);
//	printf("clz2:%s\n",clz2.to_string(2).c_str());
	ap_uint<(i4+2)*P4> clz3 = clzOneloop<i3,P3>(clz2);
//	printf("clz3:%s\n",clz3.to_string(2).c_str());
	ap_uint<(i5+2)*P5> clz4 = clzOneloop<i4,P4>(clz3);
//	printf("clz4:%s\n",clz4.to_string(2).c_str());
	ap_uint<(i6+2)*P6> clz5 = clzOneloop<i5,P5>(clz4);
//	printf("clz5:%s\n",clz5.to_string(2).c_str());
	return clz5;
}

ap_uint<8> CLZ128(ap_uint<128> v){
#pragma HLS INLINE
	const int len=128;
	//encode
	ap_uint<len> e;
	e= clzEncode<len>(v);

	const int i1=0;
	const int i2=1;
	const int i3=2;
	const int i4=3;
	const int i5=4;
	const int i6=5;
	const int i7=6;
	const int P1=len>>1;//P1/2
	const int P2=len>>2;//P1/4
	const int P3=len>>3;//P1/8
	const int P4=len>>4;
	const int P5=len>>5;
	const int P6=len>>6;
	const int P7=len>>7;
	ap_uint<(i2+2)*P2> clz1 = clzOneloop<i1,P1>(e);
//	printf("clz1:%s\n",clz1.to_string(2).c_str());
	ap_uint<(i3+2)*P3> clz2 = clzOneloop<i2,P2>(clz1);
//	printf("clz2:%s\n",clz2.to_string(2).c_str());
	ap_uint<(i4+2)*P4> clz3 = clzOneloop<i3,P3>(clz2);
//	printf("clz3:%s\n",clz3.to_string(2).c_str());
	ap_uint<(i5+2)*P5> clz4 = clzOneloop<i4,P4>(clz3);
//	printf("clz4:%s\n",clz4.to_string(2).c_str());
	ap_uint<(i6+2)*P6> clz5 = clzOneloop<i5,P5>(clz4);
//	printf("clz5:%s\n",clz5.to_string(2).c_str());
	ap_uint<(i7+2)*P7> clz6 = clzOneloop<i6,P6>(clz5);
	return clz6;
}

template<int mTR=4>
ap_uint<1> st_round(ap_uint<mTR> tail){
	ap_uint<mTR> random_bits;
	ap_uint<1> m;

	//LFSR 9 & st rounding
	static ap_uint<9> state = INIT;
	m = state[8] ^~ state[4];
	state(8,1) = state(7,0);
	state[0] = m;
	random_bits = state >> (9-mTR);

	return (tail >= random_bits)? ap_uint<1>(1) : ap_uint<1>(0);
}

#endif
