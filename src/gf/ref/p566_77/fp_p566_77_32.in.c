
//Automatically generated modular arithmetic C code
//Command line : python external/modarith/monty.py 32 18598086829520555783846380975110189094438587774212218431037944647610435872268509543080432501675502406444356152151604227666069441056269378093262747991060326103104220171337727
//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)

#include <stdio.h>
#include <stdint.h>

#define sspint int32_t
#define spint uint32_t
#define dpint uint64_t
#define sdpint int64_t
#define Wordlength 32
#define Nlimbs 20
#define Radix 29
#define Nbits 573
#define Nbytes 72

#define MONTGOMERY
//propagate carries
static inline spint prop(spint *n) {
	int i;
	spint mask=((spint)1<<29u)-(spint)1;
	sspint carry=(sspint)n[0];
	carry>>=29u;
	n[0]&=mask;
	for (i=1;i<19;i++) {
		carry+=(sspint)n[i];
		n[i] = (spint)carry & mask;
		carry>>=29u;
	}
	n[19]+=(spint)carry;
	return -((n[19]>>1)>>30u);
}

//propagate carries and add p if negative, propagate carries again
static spint flatten(spint *n) {
	spint carry=prop(n);
	n[0]-=(spint)1u&carry;
	n[19]+=((spint)0x268000u)&carry;
	(void)prop(n);
	return (carry&1);
}

//Montgomery final subtract
static spint modfsb(spint *n) {
	n[0]+=(spint)1u;
	n[19]-=(spint)0x268000u;
	return flatten(n);
}

//Modular addition - reduce less than 2p
static void modadd(const spint *a,const spint *b,spint *n) {
	spint carry;
	n[0]=a[0]+b[0];
	n[1]=a[1]+b[1];
	n[2]=a[2]+b[2];
	n[3]=a[3]+b[3];
	n[4]=a[4]+b[4];
	n[5]=a[5]+b[5];
	n[6]=a[6]+b[6];
	n[7]=a[7]+b[7];
	n[8]=a[8]+b[8];
	n[9]=a[9]+b[9];
	n[10]=a[10]+b[10];
	n[11]=a[11]+b[11];
	n[12]=a[12]+b[12];
	n[13]=a[13]+b[13];
	n[14]=a[14]+b[14];
	n[15]=a[15]+b[15];
	n[16]=a[16]+b[16];
	n[17]=a[17]+b[17];
	n[18]=a[18]+b[18];
	n[19]=a[19]+b[19];
	n[0]+=(spint)2u;
	n[19]-=(spint)0x4d0000u;
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[19]+=((spint)0x4d0000u)&carry;
	(void)prop(n);
}

//Modular subtraction - reduce less than 2p
static void modsub(const spint *a,const spint *b,spint *n) {
	spint carry;
	n[0]=a[0]-b[0];
	n[1]=a[1]-b[1];
	n[2]=a[2]-b[2];
	n[3]=a[3]-b[3];
	n[4]=a[4]-b[4];
	n[5]=a[5]-b[5];
	n[6]=a[6]-b[6];
	n[7]=a[7]-b[7];
	n[8]=a[8]-b[8];
	n[9]=a[9]-b[9];
	n[10]=a[10]-b[10];
	n[11]=a[11]-b[11];
	n[12]=a[12]-b[12];
	n[13]=a[13]-b[13];
	n[14]=a[14]-b[14];
	n[15]=a[15]-b[15];
	n[16]=a[16]-b[16];
	n[17]=a[17]-b[17];
	n[18]=a[18]-b[18];
	n[19]=a[19]-b[19];
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[19]+=((spint)0x4d0000u)&carry;
	(void)prop(n);
}

//Modular negation
static void modneg(const spint *b,spint *n) {
	spint carry;
	n[0]=(spint)0-b[0];
	n[1]=(spint)0-b[1];
	n[2]=(spint)0-b[2];
	n[3]=(spint)0-b[3];
	n[4]=(spint)0-b[4];
	n[5]=(spint)0-b[5];
	n[6]=(spint)0-b[6];
	n[7]=(spint)0-b[7];
	n[8]=(spint)0-b[8];
	n[9]=(spint)0-b[9];
	n[10]=(spint)0-b[10];
	n[11]=(spint)0-b[11];
	n[12]=(spint)0-b[12];
	n[13]=(spint)0-b[13];
	n[14]=(spint)0-b[14];
	n[15]=(spint)0-b[15];
	n[16]=(spint)0-b[16];
	n[17]=(spint)0-b[17];
	n[18]=(spint)0-b[18];
	n[19]=(spint)0-b[19];
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[19]+=((spint)0x4d0000u)&carry;
	(void)prop(n);
}

// Overflow limit   = 18446744073709551616
// maximum possible = 5765962134242033684
// Modular multiplication, c=a*b mod 2p
static void modmul(const spint *a,const spint *b,spint *c) {
	dpint t=0;
	spint p19=0x268000u;
	spint q=((spint)1<<29u); // q is unsaturated radix 
	spint mask=(spint)(q-(spint)1);
	t+=(dpint)a[0]*b[0]; spint v0=((spint)t & mask); t>>=29;
	t+=(dpint)a[0]*b[1]; t+=(dpint)a[1]*b[0]; spint v1=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[2]; t+=(dpint)a[1]*b[1]; t+=(dpint)a[2]*b[0]; spint v2=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[3]; t+=(dpint)a[1]*b[2]; t+=(dpint)a[2]*b[1]; t+=(dpint)a[3]*b[0]; spint v3=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[4]; t+=(dpint)a[1]*b[3]; t+=(dpint)a[2]*b[2]; t+=(dpint)a[3]*b[1]; t+=(dpint)a[4]*b[0]; spint v4=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[5]; t+=(dpint)a[1]*b[4]; t+=(dpint)a[2]*b[3]; t+=(dpint)a[3]*b[2]; t+=(dpint)a[4]*b[1]; t+=(dpint)a[5]*b[0]; spint v5=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[6]; t+=(dpint)a[1]*b[5]; t+=(dpint)a[2]*b[4]; t+=(dpint)a[3]*b[3]; t+=(dpint)a[4]*b[2]; t+=(dpint)a[5]*b[1]; t+=(dpint)a[6]*b[0]; spint v6=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[7]; t+=(dpint)a[1]*b[6]; t+=(dpint)a[2]*b[5]; t+=(dpint)a[3]*b[4]; t+=(dpint)a[4]*b[3]; t+=(dpint)a[5]*b[2]; t+=(dpint)a[6]*b[1]; t+=(dpint)a[7]*b[0]; spint v7=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[8]; t+=(dpint)a[1]*b[7]; t+=(dpint)a[2]*b[6]; t+=(dpint)a[3]*b[5]; t+=(dpint)a[4]*b[4]; t+=(dpint)a[5]*b[3]; t+=(dpint)a[6]*b[2]; t+=(dpint)a[7]*b[1]; t+=(dpint)a[8]*b[0]; spint v8=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[9]; t+=(dpint)a[1]*b[8]; t+=(dpint)a[2]*b[7]; t+=(dpint)a[3]*b[6]; t+=(dpint)a[4]*b[5]; t+=(dpint)a[5]*b[4]; t+=(dpint)a[6]*b[3]; t+=(dpint)a[7]*b[2]; t+=(dpint)a[8]*b[1]; t+=(dpint)a[9]*b[0]; spint v9=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[10]; t+=(dpint)a[1]*b[9]; t+=(dpint)a[2]*b[8]; t+=(dpint)a[3]*b[7]; t+=(dpint)a[4]*b[6]; t+=(dpint)a[5]*b[5]; t+=(dpint)a[6]*b[4]; t+=(dpint)a[7]*b[3]; t+=(dpint)a[8]*b[2]; t+=(dpint)a[9]*b[1]; t+=(dpint)a[10]*b[0]; spint v10=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[11]; t+=(dpint)a[1]*b[10]; t+=(dpint)a[2]*b[9]; t+=(dpint)a[3]*b[8]; t+=(dpint)a[4]*b[7]; t+=(dpint)a[5]*b[6]; t+=(dpint)a[6]*b[5]; t+=(dpint)a[7]*b[4]; t+=(dpint)a[8]*b[3]; t+=(dpint)a[9]*b[2]; t+=(dpint)a[10]*b[1]; t+=(dpint)a[11]*b[0]; spint v11=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[12]; t+=(dpint)a[1]*b[11]; t+=(dpint)a[2]*b[10]; t+=(dpint)a[3]*b[9]; t+=(dpint)a[4]*b[8]; t+=(dpint)a[5]*b[7]; t+=(dpint)a[6]*b[6]; t+=(dpint)a[7]*b[5]; t+=(dpint)a[8]*b[4]; t+=(dpint)a[9]*b[3]; t+=(dpint)a[10]*b[2]; t+=(dpint)a[11]*b[1]; t+=(dpint)a[12]*b[0]; spint v12=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[13]; t+=(dpint)a[1]*b[12]; t+=(dpint)a[2]*b[11]; t+=(dpint)a[3]*b[10]; t+=(dpint)a[4]*b[9]; t+=(dpint)a[5]*b[8]; t+=(dpint)a[6]*b[7]; t+=(dpint)a[7]*b[6]; t+=(dpint)a[8]*b[5]; t+=(dpint)a[9]*b[4]; t+=(dpint)a[10]*b[3]; t+=(dpint)a[11]*b[2]; t+=(dpint)a[12]*b[1]; t+=(dpint)a[13]*b[0]; spint v13=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[14]; t+=(dpint)a[1]*b[13]; t+=(dpint)a[2]*b[12]; t+=(dpint)a[3]*b[11]; t+=(dpint)a[4]*b[10]; t+=(dpint)a[5]*b[9]; t+=(dpint)a[6]*b[8]; t+=(dpint)a[7]*b[7]; t+=(dpint)a[8]*b[6]; t+=(dpint)a[9]*b[5]; t+=(dpint)a[10]*b[4]; t+=(dpint)a[11]*b[3]; t+=(dpint)a[12]*b[2]; t+=(dpint)a[13]*b[1]; t+=(dpint)a[14]*b[0]; spint v14=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[15]; t+=(dpint)a[1]*b[14]; t+=(dpint)a[2]*b[13]; t+=(dpint)a[3]*b[12]; t+=(dpint)a[4]*b[11]; t+=(dpint)a[5]*b[10]; t+=(dpint)a[6]*b[9]; t+=(dpint)a[7]*b[8]; t+=(dpint)a[8]*b[7]; t+=(dpint)a[9]*b[6]; t+=(dpint)a[10]*b[5]; t+=(dpint)a[11]*b[4]; t+=(dpint)a[12]*b[3]; t+=(dpint)a[13]*b[2]; t+=(dpint)a[14]*b[1]; t+=(dpint)a[15]*b[0]; spint v15=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[16]; t+=(dpint)a[1]*b[15]; t+=(dpint)a[2]*b[14]; t+=(dpint)a[3]*b[13]; t+=(dpint)a[4]*b[12]; t+=(dpint)a[5]*b[11]; t+=(dpint)a[6]*b[10]; t+=(dpint)a[7]*b[9]; t+=(dpint)a[8]*b[8]; t+=(dpint)a[9]*b[7]; t+=(dpint)a[10]*b[6]; t+=(dpint)a[11]*b[5]; t+=(dpint)a[12]*b[4]; t+=(dpint)a[13]*b[3]; t+=(dpint)a[14]*b[2]; t+=(dpint)a[15]*b[1]; t+=(dpint)a[16]*b[0]; spint v16=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[17]; t+=(dpint)a[1]*b[16]; t+=(dpint)a[2]*b[15]; t+=(dpint)a[3]*b[14]; t+=(dpint)a[4]*b[13]; t+=(dpint)a[5]*b[12]; t+=(dpint)a[6]*b[11]; t+=(dpint)a[7]*b[10]; t+=(dpint)a[8]*b[9]; t+=(dpint)a[9]*b[8]; t+=(dpint)a[10]*b[7]; t+=(dpint)a[11]*b[6]; t+=(dpint)a[12]*b[5]; t+=(dpint)a[13]*b[4]; t+=(dpint)a[14]*b[3]; t+=(dpint)a[15]*b[2]; t+=(dpint)a[16]*b[1]; t+=(dpint)a[17]*b[0]; spint v17=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[18]; t+=(dpint)a[1]*b[17]; t+=(dpint)a[2]*b[16]; t+=(dpint)a[3]*b[15]; t+=(dpint)a[4]*b[14]; t+=(dpint)a[5]*b[13]; t+=(dpint)a[6]*b[12]; t+=(dpint)a[7]*b[11]; t+=(dpint)a[8]*b[10]; t+=(dpint)a[9]*b[9]; t+=(dpint)a[10]*b[8]; t+=(dpint)a[11]*b[7]; t+=(dpint)a[12]*b[6]; t+=(dpint)a[13]*b[5]; t+=(dpint)a[14]*b[4]; t+=(dpint)a[15]*b[3]; t+=(dpint)a[16]*b[2]; t+=(dpint)a[17]*b[1]; t+=(dpint)a[18]*b[0]; spint v18=((spint)t & mask);  t>>=29;
	t+=(dpint)a[0]*b[19]; t+=(dpint)a[1]*b[18]; t+=(dpint)a[2]*b[17]; t+=(dpint)a[3]*b[16]; t+=(dpint)a[4]*b[15]; t+=(dpint)a[5]*b[14]; t+=(dpint)a[6]*b[13]; t+=(dpint)a[7]*b[12]; t+=(dpint)a[8]*b[11]; t+=(dpint)a[9]*b[10]; t+=(dpint)a[10]*b[9]; t+=(dpint)a[11]*b[8]; t+=(dpint)a[12]*b[7]; t+=(dpint)a[13]*b[6]; t+=(dpint)a[14]*b[5]; t+=(dpint)a[15]*b[4]; t+=(dpint)a[16]*b[3]; t+=(dpint)a[17]*b[2]; t+=(dpint)a[18]*b[1]; t+=(dpint)a[19]*b[0]; t+=(dpint)v0*(dpint)p19;  spint v19=((spint)t & mask);  t>>=29;
	t+=(dpint)a[1]*b[19]; t+=(dpint)a[2]*b[18]; t+=(dpint)a[3]*b[17]; t+=(dpint)a[4]*b[16]; t+=(dpint)a[5]*b[15]; t+=(dpint)a[6]*b[14]; t+=(dpint)a[7]*b[13]; t+=(dpint)a[8]*b[12]; t+=(dpint)a[9]*b[11]; t+=(dpint)a[10]*b[10]; t+=(dpint)a[11]*b[9]; t+=(dpint)a[12]*b[8]; t+=(dpint)a[13]*b[7]; t+=(dpint)a[14]*b[6]; t+=(dpint)a[15]*b[5]; t+=(dpint)a[16]*b[4]; t+=(dpint)a[17]*b[3]; t+=(dpint)a[18]*b[2]; t+=(dpint)a[19]*b[1]; t+=(dpint)v1*(dpint)p19;  c[0]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[2]*b[19]; t+=(dpint)a[3]*b[18]; t+=(dpint)a[4]*b[17]; t+=(dpint)a[5]*b[16]; t+=(dpint)a[6]*b[15]; t+=(dpint)a[7]*b[14]; t+=(dpint)a[8]*b[13]; t+=(dpint)a[9]*b[12]; t+=(dpint)a[10]*b[11]; t+=(dpint)a[11]*b[10]; t+=(dpint)a[12]*b[9]; t+=(dpint)a[13]*b[8]; t+=(dpint)a[14]*b[7]; t+=(dpint)a[15]*b[6]; t+=(dpint)a[16]*b[5]; t+=(dpint)a[17]*b[4]; t+=(dpint)a[18]*b[3]; t+=(dpint)a[19]*b[2]; t+=(dpint)v2*(dpint)p19;  c[1]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[3]*b[19]; t+=(dpint)a[4]*b[18]; t+=(dpint)a[5]*b[17]; t+=(dpint)a[6]*b[16]; t+=(dpint)a[7]*b[15]; t+=(dpint)a[8]*b[14]; t+=(dpint)a[9]*b[13]; t+=(dpint)a[10]*b[12]; t+=(dpint)a[11]*b[11]; t+=(dpint)a[12]*b[10]; t+=(dpint)a[13]*b[9]; t+=(dpint)a[14]*b[8]; t+=(dpint)a[15]*b[7]; t+=(dpint)a[16]*b[6]; t+=(dpint)a[17]*b[5]; t+=(dpint)a[18]*b[4]; t+=(dpint)a[19]*b[3]; t+=(dpint)v3*(dpint)p19;  c[2]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[4]*b[19]; t+=(dpint)a[5]*b[18]; t+=(dpint)a[6]*b[17]; t+=(dpint)a[7]*b[16]; t+=(dpint)a[8]*b[15]; t+=(dpint)a[9]*b[14]; t+=(dpint)a[10]*b[13]; t+=(dpint)a[11]*b[12]; t+=(dpint)a[12]*b[11]; t+=(dpint)a[13]*b[10]; t+=(dpint)a[14]*b[9]; t+=(dpint)a[15]*b[8]; t+=(dpint)a[16]*b[7]; t+=(dpint)a[17]*b[6]; t+=(dpint)a[18]*b[5]; t+=(dpint)a[19]*b[4]; t+=(dpint)v4*(dpint)p19;  c[3]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[5]*b[19]; t+=(dpint)a[6]*b[18]; t+=(dpint)a[7]*b[17]; t+=(dpint)a[8]*b[16]; t+=(dpint)a[9]*b[15]; t+=(dpint)a[10]*b[14]; t+=(dpint)a[11]*b[13]; t+=(dpint)a[12]*b[12]; t+=(dpint)a[13]*b[11]; t+=(dpint)a[14]*b[10]; t+=(dpint)a[15]*b[9]; t+=(dpint)a[16]*b[8]; t+=(dpint)a[17]*b[7]; t+=(dpint)a[18]*b[6]; t+=(dpint)a[19]*b[5]; t+=(dpint)v5*(dpint)p19;  c[4]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[6]*b[19]; t+=(dpint)a[7]*b[18]; t+=(dpint)a[8]*b[17]; t+=(dpint)a[9]*b[16]; t+=(dpint)a[10]*b[15]; t+=(dpint)a[11]*b[14]; t+=(dpint)a[12]*b[13]; t+=(dpint)a[13]*b[12]; t+=(dpint)a[14]*b[11]; t+=(dpint)a[15]*b[10]; t+=(dpint)a[16]*b[9]; t+=(dpint)a[17]*b[8]; t+=(dpint)a[18]*b[7]; t+=(dpint)a[19]*b[6]; t+=(dpint)v6*(dpint)p19;  c[5]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[7]*b[19]; t+=(dpint)a[8]*b[18]; t+=(dpint)a[9]*b[17]; t+=(dpint)a[10]*b[16]; t+=(dpint)a[11]*b[15]; t+=(dpint)a[12]*b[14]; t+=(dpint)a[13]*b[13]; t+=(dpint)a[14]*b[12]; t+=(dpint)a[15]*b[11]; t+=(dpint)a[16]*b[10]; t+=(dpint)a[17]*b[9]; t+=(dpint)a[18]*b[8]; t+=(dpint)a[19]*b[7]; t+=(dpint)v7*(dpint)p19;  c[6]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[8]*b[19]; t+=(dpint)a[9]*b[18]; t+=(dpint)a[10]*b[17]; t+=(dpint)a[11]*b[16]; t+=(dpint)a[12]*b[15]; t+=(dpint)a[13]*b[14]; t+=(dpint)a[14]*b[13]; t+=(dpint)a[15]*b[12]; t+=(dpint)a[16]*b[11]; t+=(dpint)a[17]*b[10]; t+=(dpint)a[18]*b[9]; t+=(dpint)a[19]*b[8]; t+=(dpint)v8*(dpint)p19;  c[7]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[9]*b[19]; t+=(dpint)a[10]*b[18]; t+=(dpint)a[11]*b[17]; t+=(dpint)a[12]*b[16]; t+=(dpint)a[13]*b[15]; t+=(dpint)a[14]*b[14]; t+=(dpint)a[15]*b[13]; t+=(dpint)a[16]*b[12]; t+=(dpint)a[17]*b[11]; t+=(dpint)a[18]*b[10]; t+=(dpint)a[19]*b[9]; t+=(dpint)v9*(dpint)p19;  c[8]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[10]*b[19]; t+=(dpint)a[11]*b[18]; t+=(dpint)a[12]*b[17]; t+=(dpint)a[13]*b[16]; t+=(dpint)a[14]*b[15]; t+=(dpint)a[15]*b[14]; t+=(dpint)a[16]*b[13]; t+=(dpint)a[17]*b[12]; t+=(dpint)a[18]*b[11]; t+=(dpint)a[19]*b[10]; t+=(dpint)v10*(dpint)p19;  c[9]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[11]*b[19]; t+=(dpint)a[12]*b[18]; t+=(dpint)a[13]*b[17]; t+=(dpint)a[14]*b[16]; t+=(dpint)a[15]*b[15]; t+=(dpint)a[16]*b[14]; t+=(dpint)a[17]*b[13]; t+=(dpint)a[18]*b[12]; t+=(dpint)a[19]*b[11]; t+=(dpint)v11*(dpint)p19;  c[10]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[12]*b[19]; t+=(dpint)a[13]*b[18]; t+=(dpint)a[14]*b[17]; t+=(dpint)a[15]*b[16]; t+=(dpint)a[16]*b[15]; t+=(dpint)a[17]*b[14]; t+=(dpint)a[18]*b[13]; t+=(dpint)a[19]*b[12]; t+=(dpint)v12*(dpint)p19;  c[11]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[13]*b[19]; t+=(dpint)a[14]*b[18]; t+=(dpint)a[15]*b[17]; t+=(dpint)a[16]*b[16]; t+=(dpint)a[17]*b[15]; t+=(dpint)a[18]*b[14]; t+=(dpint)a[19]*b[13]; t+=(dpint)v13*(dpint)p19;  c[12]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[14]*b[19]; t+=(dpint)a[15]*b[18]; t+=(dpint)a[16]*b[17]; t+=(dpint)a[17]*b[16]; t+=(dpint)a[18]*b[15]; t+=(dpint)a[19]*b[14]; t+=(dpint)v14*(dpint)p19;  c[13]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[15]*b[19]; t+=(dpint)a[16]*b[18]; t+=(dpint)a[17]*b[17]; t+=(dpint)a[18]*b[16]; t+=(dpint)a[19]*b[15]; t+=(dpint)v15*(dpint)p19;  c[14]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[16]*b[19]; t+=(dpint)a[17]*b[18]; t+=(dpint)a[18]*b[17]; t+=(dpint)a[19]*b[16]; t+=(dpint)v16*(dpint)p19;  c[15]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[17]*b[19]; t+=(dpint)a[18]*b[18]; t+=(dpint)a[19]*b[17]; t+=(dpint)v17*(dpint)p19;  c[16]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[18]*b[19]; t+=(dpint)a[19]*b[18]; t+=(dpint)v18*(dpint)p19;  c[17]=((spint)t & mask);  t>>=29;
	t+=(dpint)a[19]*b[19]; t+=(dpint)v19*(dpint)p19;  c[18]=((spint)t & mask);  t>>=29;
	c[19] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
static void modsqr(const spint *a,spint *c) {
	dpint tot;
	dpint t=0;
	spint p19=0x268000u;
	spint q=((spint)1<<29u); // q is unsaturated radix 
	spint mask=(spint)(q-(spint)1);
	tot=(dpint)a[0]*a[0]; t=tot; spint v0=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[1]; tot*=2; t+=tot;  spint v1=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[2]; tot*=2; tot+=(dpint)a[1]*a[1]; t+=tot;  spint v2=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[3]; tot+=(dpint)a[1]*a[2]; tot*=2; t+=tot;  spint v3=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[4]; tot+=(dpint)a[1]*a[3]; tot*=2; tot+=(dpint)a[2]*a[2]; t+=tot;  spint v4=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[5]; tot+=(dpint)a[1]*a[4]; tot+=(dpint)a[2]*a[3]; tot*=2; t+=tot;  spint v5=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[6]; tot+=(dpint)a[1]*a[5]; tot+=(dpint)a[2]*a[4]; tot*=2; tot+=(dpint)a[3]*a[3]; t+=tot;  spint v6=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[7]; tot+=(dpint)a[1]*a[6]; tot+=(dpint)a[2]*a[5]; tot+=(dpint)a[3]*a[4]; tot*=2; t+=tot;  spint v7=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[8]; tot+=(dpint)a[1]*a[7]; tot+=(dpint)a[2]*a[6]; tot+=(dpint)a[3]*a[5]; tot*=2; tot+=(dpint)a[4]*a[4]; t+=tot;  spint v8=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[9]; tot+=(dpint)a[1]*a[8]; tot+=(dpint)a[2]*a[7]; tot+=(dpint)a[3]*a[6]; tot+=(dpint)a[4]*a[5]; tot*=2; t+=tot;  spint v9=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[10]; tot+=(dpint)a[1]*a[9]; tot+=(dpint)a[2]*a[8]; tot+=(dpint)a[3]*a[7]; tot+=(dpint)a[4]*a[6]; tot*=2; tot+=(dpint)a[5]*a[5]; t+=tot;  spint v10=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[11]; tot+=(dpint)a[1]*a[10]; tot+=(dpint)a[2]*a[9]; tot+=(dpint)a[3]*a[8]; tot+=(dpint)a[4]*a[7]; tot+=(dpint)a[5]*a[6]; tot*=2; t+=tot;  spint v11=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[12]; tot+=(dpint)a[1]*a[11]; tot+=(dpint)a[2]*a[10]; tot+=(dpint)a[3]*a[9]; tot+=(dpint)a[4]*a[8]; tot+=(dpint)a[5]*a[7]; tot*=2; tot+=(dpint)a[6]*a[6]; t+=tot;  spint v12=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[13]; tot+=(dpint)a[1]*a[12]; tot+=(dpint)a[2]*a[11]; tot+=(dpint)a[3]*a[10]; tot+=(dpint)a[4]*a[9]; tot+=(dpint)a[5]*a[8]; tot+=(dpint)a[6]*a[7]; tot*=2; t+=tot;  spint v13=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[14]; tot+=(dpint)a[1]*a[13]; tot+=(dpint)a[2]*a[12]; tot+=(dpint)a[3]*a[11]; tot+=(dpint)a[4]*a[10]; tot+=(dpint)a[5]*a[9]; tot+=(dpint)a[6]*a[8]; tot*=2; tot+=(dpint)a[7]*a[7]; t+=tot;  spint v14=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[15]; tot+=(dpint)a[1]*a[14]; tot+=(dpint)a[2]*a[13]; tot+=(dpint)a[3]*a[12]; tot+=(dpint)a[4]*a[11]; tot+=(dpint)a[5]*a[10]; tot+=(dpint)a[6]*a[9]; tot+=(dpint)a[7]*a[8]; tot*=2; t+=tot;  spint v15=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[16]; tot+=(dpint)a[1]*a[15]; tot+=(dpint)a[2]*a[14]; tot+=(dpint)a[3]*a[13]; tot+=(dpint)a[4]*a[12]; tot+=(dpint)a[5]*a[11]; tot+=(dpint)a[6]*a[10]; tot+=(dpint)a[7]*a[9]; tot*=2; tot+=(dpint)a[8]*a[8]; t+=tot;  spint v16=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[17]; tot+=(dpint)a[1]*a[16]; tot+=(dpint)a[2]*a[15]; tot+=(dpint)a[3]*a[14]; tot+=(dpint)a[4]*a[13]; tot+=(dpint)a[5]*a[12]; tot+=(dpint)a[6]*a[11]; tot+=(dpint)a[7]*a[10]; tot+=(dpint)a[8]*a[9]; tot*=2; t+=tot;  spint v17=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[18]; tot+=(dpint)a[1]*a[17]; tot+=(dpint)a[2]*a[16]; tot+=(dpint)a[3]*a[15]; tot+=(dpint)a[4]*a[14]; tot+=(dpint)a[5]*a[13]; tot+=(dpint)a[6]*a[12]; tot+=(dpint)a[7]*a[11]; tot+=(dpint)a[8]*a[10]; tot*=2; tot+=(dpint)a[9]*a[9]; t+=tot;  spint v18=((spint)t & mask); t>>=29;
	tot=(dpint)a[0]*a[19]; tot+=(dpint)a[1]*a[18]; tot+=(dpint)a[2]*a[17]; tot+=(dpint)a[3]*a[16]; tot+=(dpint)a[4]*a[15]; tot+=(dpint)a[5]*a[14]; tot+=(dpint)a[6]*a[13]; tot+=(dpint)a[7]*a[12]; tot+=(dpint)a[8]*a[11]; tot+=(dpint)a[9]*a[10]; tot*=2; t+=tot;  t+=(dpint)v0*p19;  spint v19=((spint)t & mask); t>>=29;
	tot=(dpint)a[1]*a[19]; tot+=(dpint)a[2]*a[18]; tot+=(dpint)a[3]*a[17]; tot+=(dpint)a[4]*a[16]; tot+=(dpint)a[5]*a[15]; tot+=(dpint)a[6]*a[14]; tot+=(dpint)a[7]*a[13]; tot+=(dpint)a[8]*a[12]; tot+=(dpint)a[9]*a[11]; tot*=2; tot+=(dpint)a[10]*a[10]; t+=tot;  t+=(dpint)v1*p19;  c[0]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[2]*a[19]; tot+=(dpint)a[3]*a[18]; tot+=(dpint)a[4]*a[17]; tot+=(dpint)a[5]*a[16]; tot+=(dpint)a[6]*a[15]; tot+=(dpint)a[7]*a[14]; tot+=(dpint)a[8]*a[13]; tot+=(dpint)a[9]*a[12]; tot+=(dpint)a[10]*a[11]; tot*=2; t+=tot;  t+=(dpint)v2*p19;  c[1]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[3]*a[19]; tot+=(dpint)a[4]*a[18]; tot+=(dpint)a[5]*a[17]; tot+=(dpint)a[6]*a[16]; tot+=(dpint)a[7]*a[15]; tot+=(dpint)a[8]*a[14]; tot+=(dpint)a[9]*a[13]; tot+=(dpint)a[10]*a[12]; tot*=2; tot+=(dpint)a[11]*a[11]; t+=tot;  t+=(dpint)v3*p19;  c[2]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[4]*a[19]; tot+=(dpint)a[5]*a[18]; tot+=(dpint)a[6]*a[17]; tot+=(dpint)a[7]*a[16]; tot+=(dpint)a[8]*a[15]; tot+=(dpint)a[9]*a[14]; tot+=(dpint)a[10]*a[13]; tot+=(dpint)a[11]*a[12]; tot*=2; t+=tot;  t+=(dpint)v4*p19;  c[3]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[5]*a[19]; tot+=(dpint)a[6]*a[18]; tot+=(dpint)a[7]*a[17]; tot+=(dpint)a[8]*a[16]; tot+=(dpint)a[9]*a[15]; tot+=(dpint)a[10]*a[14]; tot+=(dpint)a[11]*a[13]; tot*=2; tot+=(dpint)a[12]*a[12]; t+=tot;  t+=(dpint)v5*p19;  c[4]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[6]*a[19]; tot+=(dpint)a[7]*a[18]; tot+=(dpint)a[8]*a[17]; tot+=(dpint)a[9]*a[16]; tot+=(dpint)a[10]*a[15]; tot+=(dpint)a[11]*a[14]; tot+=(dpint)a[12]*a[13]; tot*=2; t+=tot;  t+=(dpint)v6*p19;  c[5]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[7]*a[19]; tot+=(dpint)a[8]*a[18]; tot+=(dpint)a[9]*a[17]; tot+=(dpint)a[10]*a[16]; tot+=(dpint)a[11]*a[15]; tot+=(dpint)a[12]*a[14]; tot*=2; tot+=(dpint)a[13]*a[13]; t+=tot;  t+=(dpint)v7*p19;  c[6]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[8]*a[19]; tot+=(dpint)a[9]*a[18]; tot+=(dpint)a[10]*a[17]; tot+=(dpint)a[11]*a[16]; tot+=(dpint)a[12]*a[15]; tot+=(dpint)a[13]*a[14]; tot*=2; t+=tot;  t+=(dpint)v8*p19;  c[7]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[9]*a[19]; tot+=(dpint)a[10]*a[18]; tot+=(dpint)a[11]*a[17]; tot+=(dpint)a[12]*a[16]; tot+=(dpint)a[13]*a[15]; tot*=2; tot+=(dpint)a[14]*a[14]; t+=tot;  t+=(dpint)v9*p19;  c[8]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[10]*a[19]; tot+=(dpint)a[11]*a[18]; tot+=(dpint)a[12]*a[17]; tot+=(dpint)a[13]*a[16]; tot+=(dpint)a[14]*a[15]; tot*=2; t+=tot;  t+=(dpint)v10*p19;  c[9]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[11]*a[19]; tot+=(dpint)a[12]*a[18]; tot+=(dpint)a[13]*a[17]; tot+=(dpint)a[14]*a[16]; tot*=2; tot+=(dpint)a[15]*a[15]; t+=tot;  t+=(dpint)v11*p19;  c[10]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[12]*a[19]; tot+=(dpint)a[13]*a[18]; tot+=(dpint)a[14]*a[17]; tot+=(dpint)a[15]*a[16]; tot*=2; t+=tot;  t+=(dpint)v12*p19;  c[11]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[13]*a[19]; tot+=(dpint)a[14]*a[18]; tot+=(dpint)a[15]*a[17]; tot*=2; tot+=(dpint)a[16]*a[16]; t+=tot;  t+=(dpint)v13*p19;  c[12]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[14]*a[19]; tot+=(dpint)a[15]*a[18]; tot+=(dpint)a[16]*a[17]; tot*=2; t+=tot;  t+=(dpint)v14*p19;  c[13]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[15]*a[19]; tot+=(dpint)a[16]*a[18]; tot*=2; tot+=(dpint)a[17]*a[17]; t+=tot;  t+=(dpint)v15*p19;  c[14]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[16]*a[19]; tot+=(dpint)a[17]*a[18]; tot*=2; t+=tot;  t+=(dpint)v16*p19;  c[15]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[17]*a[19]; tot*=2; tot+=(dpint)a[18]*a[18]; t+=tot;  t+=(dpint)v17*p19;  c[16]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[18]*a[19]; tot*=2; t+=tot;  t+=(dpint)v18*p19;  c[17]=((spint)t & mask);  t>>=29;
	tot=(dpint)a[19]*a[19]; t+=tot;  t+=(dpint)v19*p19;  c[18]=((spint)t & mask);  t>>=29;
	c[19] = (spint)t;
}

//copy
static void modcpy(const spint *a,spint *c) {
	int i;
	for (i=0;i<20;i++) {
		c[i]=a[i];
	}
}

//square n times
static void modnsqr(spint *a,int n) {
	int i;
	for (i=0;i<n;i++) {
		modsqr(a,a);
	}
}

//Calculate progenitor
static void modpro(const spint *w,spint *z) {
	spint x[20];
	spint t0[20];
	spint t1[20];
	spint t2[20];
	spint t3[20];
	spint t4[20];
	spint t5[20];
	modcpy(w,x);
	modsqr(x,t0);
	modsqr(t0,z);
	modmul(t0,z,t1);
	modmul(x,t1,z);
	modsqr(z,t2);
	modmul(z,t2,t2);
	modsqr(t2,t2);
	modmul(z,t2,t2);
	modmul(t0,t2,t0);
	modcpy(t0,t2);
	modnsqr(t2,2);
	modmul(t0,t2,t3);
	modsqr(t3,t4);
	modcpy(t4,t2);
	modnsqr(t2,2);
	modcpy(t2,t5);
	modnsqr(t5,6);
	modmul(t4,t5,t4);
	modnsqr(t4,10);
	modmul(t2,t4,t2);
	modcpy(t2,t4);
	modnsqr(t4,5);
	modmul(t3,t4,t3);
	modnsqr(t3,27);
	modmul(t2,t3,t2);
	modmul(z,t2,z);
	modcpy(z,t2);
	modnsqr(t2,3);
	modcpy(t2,t3);
	modnsqr(t3,59);
	modmul(t2,t3,t2);
	modmul(t1,t2,t1);
	modsqr(t1,t2);
	modmul(t1,t2,t2);
	modnsqr(t2,3);
	modmul(t1,t2,t2);
	modmul(t0,t2,t0);
	modmul(t1,t0,t1);
	modmul(t0,t1,t0);
	modmul(t1,t0,t2);
	modmul(t0,t2,t0);
	modmul(t1,t0,t1);
	modnsqr(t1,127);
	modmul(t0,t1,t1);
	modnsqr(t1,128);
	modmul(t0,t1,t1);
	modnsqr(t1,128);
	modmul(t0,t1,t0);
	modnsqr(t0,59);
	modmul(z,t0,z);
}

//calculate inverse, provide progenitor h if available
static void modinv(const spint *x,const spint *h,spint *z) {
	spint s[20];
	spint t[20];
	if (h==NULL) {
		modpro(x,t);
	} else {
		modcpy(h,t);
	}
	modcpy(x,s);
	modnsqr(t,2);
	modmul(s,t,z);
}

//Convert m to n-residue form, n=nres(m) 
static void nres(const spint *m,spint *n) {
	const spint c[20]={0x1b0403f8u,0x1606a63bu,0xc0d4c77u,0x181a98efu,0x103531deu,0x6a63bdu,0xd4c77bu,0x1a98ef6u,0x3531decu,0x6a63bd8u,0xd4c77b0u,0x1a98ef60u,0x1531dec0u,0xa63bd81u,0x14c77b03u,0x98ef606u,0x131dec0du,0x63bd81au,0xc77b035u,0x23606au};
	modmul(m,c,n);
}

//Convert n back to normal form, m=redc(n) 
static void redc(const spint *n,spint *m) {
	int i;
	spint c[20];
	c[0]=1;
	for (i=1;i<20;i++) {
		c[i]=0;
	}
	modmul(n,c,m);
	(void)modfsb(m);
}

//is unity?
static int modis1(const spint *a) {
	int i;
	spint c[20];
	spint c0;
	spint d=0;
	redc(a,c);
	for (i=1;i<20;i++) {
		d|=c[i];
	}
	c0=(spint)c[0];
	return ((spint)1 & ((d-(spint)1)>>29u) & (((c0^(spint)1)-(spint)1)>>29u));
}

//is zero?
static int modis0(const spint *a) {
	int i;
	spint c[20];
	spint d=0;
	redc(a,c);
	for (i=0;i<20;i++) {
		d|=c[i];
	}
	return ((spint)1 & ((d-(spint)1)>>29u));
}

//set to zero
static void modzer(spint *a) {
	int i;
	for (i=0;i<20;i++) {
		a[i]=0;
	}
}

//set to one
static void modone(spint *a) {
	int i;
	a[0]=1;
	for (i=1;i<20;i++) {
		a[i]=0;
	}
	nres(a,a);
}

//set to integer
static void modint(int x,spint *a) {
	int i;
	a[0]=(spint)x;
	for (i=1;i<20;i++) {
		a[i]=0;
	}
	nres(a,a);
}

// Modular multiplication by an integer, c=a*b mod 2p
// uses special method for trinomials, otherwise Barrett-Dhem reduction
static void modmli(const spint *a,int b,spint *c) {
	spint p19=0x268000u;
	spint mask=((spint)1<<29u)-(spint)1;
	dpint t=0;
	spint q,h,r=0x3531dec0;
	t+=(dpint)a[0]*(dpint)b; c[0]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[1]*(dpint)b; c[1]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[2]*(dpint)b; c[2]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[3]*(dpint)b; c[3]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[4]*(dpint)b; c[4]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[5]*(dpint)b; c[5]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[6]*(dpint)b; c[6]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[7]*(dpint)b; c[7]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[8]*(dpint)b; c[8]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[9]*(dpint)b; c[9]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[10]*(dpint)b; c[10]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[11]*(dpint)b; c[11]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[12]*(dpint)b; c[12]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[13]*(dpint)b; c[13]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[14]*(dpint)b; c[14]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[15]*(dpint)b; c[15]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[16]*(dpint)b; c[16]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[17]*(dpint)b; c[17]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[18]*(dpint)b; c[18]=(spint)t & mask; t=t>>29u;
	t+=(dpint)a[19]*(dpint)b; c[19]=(spint)t;
	
//Barrett-Dhem reduction
	h = (spint)(t>>19u);
	q=(spint)(((dpint)h*(dpint)r)>>32u);
	c[0]+=q;
	c[19]-=q*p19;
}

//Test for quadratic residue 
static int modqr(const spint *h,const spint *x) {
	spint r[20];
	if (h==NULL) {
		modpro(x,r);
		modsqr(r,r);
	} else {
		modsqr(h,r);
	}
	modmul(r,x,r);
	return modis1(r) | modis0(x);
}

//conditional move g to f if d=1
//strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__ ((noinline)) modcmv(int b,const spint *g,volatile spint *f) {
	int i;
	spint c0,c1,s,t,w,aux;
	static spint R=0;
	R+=0x5aa5a55au;
	w=R;
	c0=(~b)&(w+1);
	c1=b+w;
	for (i=0;i<20;i++) {
		s=g[i]; t=f[i];
		f[i] = aux = c0*t+c1*s;
		f[i] = aux - w*(t+s);
	}
}

//conditional swap g and f if d=1
//strongly recommend inlining be disabled using compiler specific syntax
static void __attribute__ ((noinline)) modcsw(int b,volatile spint *g,volatile spint *f) {
	int i;
	spint c0,c1,s,t,w,v,aux;
	static spint R=0;
	R+=0x5aa5a55au;
	w=R;
	c0=(~b)&(w+1);
	c1=b+w;
	for (i=0;i<20;i++) {
		s=g[i]; t=f[i];
		v=w*(t+s);
		f[i] = aux = c0*t+c1*s;
		f[i] = aux - v;
		g[i] = aux = c0*s+c1*t;
		g[i] = aux - v;
	}
}

//Modular square root, provide progenitor h if available, NULL if not
static void modsqrt(const spint *x,const spint *h,spint *r) {
	spint s[20];
	spint y[20];
	if (h==NULL) {
		modpro(x,y);
	} else {
		modcpy(h,y);
	}
	modmul(y,x,s);
	modcpy(s,r);
}

//shift left by less than a word
static void modshl(unsigned int n,spint *a) {
	int i;
	a[19]=((a[19]<<n)) + (a[18]>>(29u-n));
	for (i=18;i>0;i--) {
		a[i]=((a[i]<<n)&(spint)0x1fffffff) + (a[i-1]>>(29u-n));
	}
	a[0]=(a[0]<<n)&(spint)0x1fffffff;
}

//shift right by less than a word. Return shifted out part
static int modshr(unsigned int n,spint *a) {
	int i;
	spint r=a[0]&(((spint)1<<n)-(spint)1);
	for (i=0;i<19;i++) {
		a[i]=(a[i]>>n) + ((a[i+1]<<(29u-n))&(spint)0x1fffffff);
	}
	a[19]=a[19]>>n;
	return r;
}

//divide by 2. Shift right 1 bit (or add p and shift right one bit)
static void modhaf(spint *n) {
	int lsb;
	spint t[20];
	(void)prop(n);
	modcpy(n,t);
	lsb=modshr(1,t);
	n[0]-=(spint)1;
	n[19]+=((spint)0x268000u);
	(void)prop(n);
	modshr(1,n);
	modcmv(1-lsb,t,n);
}

//set a= 2^r
static void mod2r(unsigned int r,spint *a) {
	unsigned int n=r/29u;
	unsigned int m=r%29u;
	modzer(a);
	if (r>=72*8) return;
	a[n]=1; a[n]<<=m;
nres(a,a);
}

//export to byte array
static void modexp(const spint *a,char *b) {
	int i;
	spint c[20];
	redc(a,c);
	for (i=71;i>=0;i--) {
		b[i]=c[0]&(spint)0xff;
		(void)modshr(8,c);
	}
}

//import from byte array
//returns 1 if in range, else 0
static int modimp(const char *b, spint *a) {
	int i,res;
	for (i=0;i<20;i++) {
		a[i]=0;
	}
	for (i=0;i<72;i++) {
		modshl(8,a);
		a[0]+=(spint)(unsigned char)b[i];
	}
	res=(int)modfsb(a);
	nres(a,a);
	return res;
}

//determine sign
static int modsign(const spint *a) {
	spint c[20];
	redc(a,c);
	return c[0]%2;
}

//return true if equal
static int modcmp(const spint *a,const spint *b) {
	spint c[20],d[20];
	int i,eq=1;
	redc(a,c);
	redc(b,d);
	for (i=0;i<20;i++) {
		eq&=(((c[i]^d[i])-1)>>29)&1;
	}
	return eq;
}


/******************************************************************************
MIKE API functions calling generated code above
******************************************************************************/

#include <fp.h>

const digit_t ZERO[NWORDS_FIELD] = {0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0};
const digit_t ONE[NWORDS_FIELD] = {0xd4,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x1e0000};
// Montgomery representation of 2^-1
static const digit_t TWO_INV[NWORDS_FIELD] = {0x6a,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0x0,
0xf0000};
// Montgomery representation of 3^-1
static const digit_t THREE_INV[NWORDS_FIELD] = {0xaaaaaf1,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x15555555,
0xaaaaaaa,
0x16d555};

void
fp_set_small(fp_t *x, const digit_t val)
{
    modint((int)val, *x);
}

void
fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val)
{
    modmli(*a, (int)val, *x);
}

void
fp_set_zero(fp_t *x)
{
    modzer(*x);
}

void
fp_set_one(fp_t *x)
{
    modone(*x);
}

uint32_t
fp_is_equal(const fp_t *a, const fp_t *b)
{
    return -(uint32_t)modcmp(*a, *b);
}

uint32_t
fp_is_zero(const fp_t *a)
{
    return -(uint32_t)modis0(*a);
}

void
fp_copy(fp_t *out, const fp_t *a)
{
    modcpy(*a, *out);
}

void
fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
{
    modcsw((int)(ctl & 0x1), *a, *b);
}

void
fp_add(fp_t *out, const fp_t *a, const fp_t *b)
{
    modadd(*a, *b, *out);
}

void
fp_sub(fp_t *out, const fp_t *a, const fp_t *b)
{
    modsub(*a, *b, *out);
}

void
fp_neg(fp_t *out, const fp_t *a)
{
    modneg(*a, *out);
}

void
fp_sqr(fp_t *out, const fp_t *a)
{
    modsqr(*a, *out);
}

void
fp_mul(fp_t *out, const fp_t *a, const fp_t *b)
{
    modmul(*a, *b, *out);
}

void
fp_inv(fp_t *x)
{
    modinv(*x, NULL, *x);
}

uint32_t
fp_is_square(const fp_t *a)
{
    return -(uint32_t)modqr(NULL, *a);
}

void
fp_sqrt(fp_t *a)
{
    modsqrt(*a, NULL, *a);
}

void
fp_half(fp_t *out, const fp_t *a)
{
    modmul(TWO_INV, *a, *out);
}

void
fp_exp3div4(fp_t *out, const fp_t *a)
{
    modpro(*a, *out);
}

void
fp_div3(fp_t *out, const fp_t *a)
{
    modmul(THREE_INV, *a, *out);
}

void
fp_encode(void *dst, const fp_t *a)
{
    // little-endian canonical encoding
    int i;
    spint c[Nlimbs]= {0};
    redc(*a, c);
    for (i = 0; i < Nbytes; i++) {
        ((char *)dst)[i] = c[0] & (spint)0xff;
        (void)modshr(8, c);
    }
}

uint32_t
fp_decode(fp_t *d, const void *src)
{
    // inverse of fp_encode; returns 0xFFFFFFFF iff the input was canonical
    int i;
    spint res;
    const unsigned char *b = src;
    for (i = 0; i < Nlimbs; i++) {
        (*d)[i] = 0;
    }
    for (i = Nbytes - 1; i >= 0; i--) {
        modshl(8, *d);
        (*d)[0] += (spint)b[i];
    }
    res = (spint)-modfsb(*d);
    nres(*d, *d);
    for (i = 0; i < Nlimbs; i++) {
        (*d)[i] &= res;
    }
    return (uint32_t)res;
}

void
fp_decode_reduce(fp_t *d, const void *src, size_t len)
{
    // Reduce a little-endian byte string of arbitrary length mod p.
    // Radix-256 Horner using only fp_add/fp_set_small: correct for ANY prime,
    // no prime-shape-specific fast reduction required.
    const unsigned char *b = src;
    size_t i;
    fp_set_zero(d);
    for (i = len; i-- > 0;) {
        fp_t t;
        int k;
        for (k = 0; k < 8; k++) {
            fp_add(d, d, d); // d *= 2  (eight times => d *= 256)
        }
        fp_set_small(&t, (digit_t)b[i]);
        fp_add(d, d, &t);
    }
}

