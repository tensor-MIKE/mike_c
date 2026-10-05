
//Automatically generated modular arithmetic C code
//Command line : python external/modarith/monty.py 64 330097605893050801307629086915362912242068350256777637543085287536978780864049186563607177986047
//Python Script by Mike Scott (Technology Innovation Institute, UAE, 2025)

#include <stdio.h>
#include <stdint.h>

#define sspint int64_t
#define spint uint64_t
#define dpint __uint128_t
#define sdpint __int128_t
#define Wordlength 64
#define Nlimbs 6
#define Radix 54
#define Nbits 318
#define Nbytes 40

#define MONTGOMERY
//propagate carries
static inline spint prop(spint *n) {
	int i;
	spint mask=((spint)1<<54u)-(spint)1;
	sspint carry=(sspint)n[0];
	carry>>=54u;
	n[0]&=mask;
	for (i=1;i<5;i++) {
		carry+=(sspint)n[i];
		n[i] = (spint)carry & mask;
		carry>>=54u;
	}
	n[5]+=(spint)carry;
	return -((n[5]>>1)>>62u);
}

//propagate carries and add p if negative, propagate carries again
static inline spint flatten(spint *n) {
	spint carry=prop(n);
	n[0]-=(spint)1u&carry;
	n[5]+=((spint)0x9e4000000000u)&carry;
	(void)prop(n);
	return (carry&1);
}

//Montgomery final subtract
static inline spint modfsb(spint *n) {
	n[0]+=(spint)1u;
	n[5]-=(spint)0x9e4000000000u;
	return flatten(n);
}

//Modular addition - reduce less than 2p
static inline void modadd(const spint *a,const spint *b,spint *n) {
	spint carry;
	n[0]=a[0]+b[0];
	n[1]=a[1]+b[1];
	n[2]=a[2]+b[2];
	n[3]=a[3]+b[3];
	n[4]=a[4]+b[4];
	n[5]=a[5]+b[5];
	n[0]+=(spint)2u;
	n[5]-=(spint)0x13c8000000000u;
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[5]+=((spint)0x13c8000000000u)&carry;
	(void)prop(n);
}

//Modular subtraction - reduce less than 2p
static inline void modsub(const spint *a,const spint *b,spint *n) {
	spint carry;
	n[0]=a[0]-b[0];
	n[1]=a[1]-b[1];
	n[2]=a[2]-b[2];
	n[3]=a[3]-b[3];
	n[4]=a[4]-b[4];
	n[5]=a[5]-b[5];
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[5]+=((spint)0x13c8000000000u)&carry;
	(void)prop(n);
}

//Modular negation
static inline void modneg(const spint *b,spint *n) {
	spint carry;
	n[0]=(spint)0-b[0];
	n[1]=(spint)0-b[1];
	n[2]=(spint)0-b[2];
	n[3]=(spint)0-b[3];
	n[4]=(spint)0-b[4];
	n[5]=(spint)0-b[5];
	carry=prop(n);
	n[0]-=(spint)2u&carry;
	n[5]+=((spint)0x13c8000000000u)&carry;
	(void)prop(n);
}

// Overflow limit   = 340282366920938463463374607431768211456
// maximum possible = 1950245786148920193139679457968134
// Modular multiplication, c=a*b mod 2p
static inline void modmul(const spint *a,const spint *b,spint *c) {
	dpint t=0;
	spint p5=0x9e4000000000u;
	spint q=((spint)1<<54u); // q is unsaturated radix 
	spint mask=(spint)(q-(spint)1);
	t+=(dpint)a[0]*b[0]; spint v0=((spint)t & mask); t>>=54;
	t+=(dpint)a[0]*b[1]; t+=(dpint)a[1]*b[0]; spint v1=((spint)t & mask);  t>>=54;
	t+=(dpint)a[0]*b[2]; t+=(dpint)a[1]*b[1]; t+=(dpint)a[2]*b[0]; spint v2=((spint)t & mask);  t>>=54;
	t+=(dpint)a[0]*b[3]; t+=(dpint)a[1]*b[2]; t+=(dpint)a[2]*b[1]; t+=(dpint)a[3]*b[0]; spint v3=((spint)t & mask);  t>>=54;
	t+=(dpint)a[0]*b[4]; t+=(dpint)a[1]*b[3]; t+=(dpint)a[2]*b[2]; t+=(dpint)a[3]*b[1]; t+=(dpint)a[4]*b[0]; spint v4=((spint)t & mask);  t>>=54;
	t+=(dpint)a[0]*b[5]; t+=(dpint)a[1]*b[4]; t+=(dpint)a[2]*b[3]; t+=(dpint)a[3]*b[2]; t+=(dpint)a[4]*b[1]; t+=(dpint)a[5]*b[0]; t+=(dpint)v0*(dpint)p5;  spint v5=((spint)t & mask);  t>>=54;
	t+=(dpint)a[1]*b[5]; t+=(dpint)a[2]*b[4]; t+=(dpint)a[3]*b[3]; t+=(dpint)a[4]*b[2]; t+=(dpint)a[5]*b[1]; t+=(dpint)v1*(dpint)p5;  c[0]=((spint)t & mask);  t>>=54;
	t+=(dpint)a[2]*b[5]; t+=(dpint)a[3]*b[4]; t+=(dpint)a[4]*b[3]; t+=(dpint)a[5]*b[2]; t+=(dpint)v2*(dpint)p5;  c[1]=((spint)t & mask);  t>>=54;
	t+=(dpint)a[3]*b[5]; t+=(dpint)a[4]*b[4]; t+=(dpint)a[5]*b[3]; t+=(dpint)v3*(dpint)p5;  c[2]=((spint)t & mask);  t>>=54;
	t+=(dpint)a[4]*b[5]; t+=(dpint)a[5]*b[4]; t+=(dpint)v4*(dpint)p5;  c[3]=((spint)t & mask);  t>>=54;
	t+=(dpint)a[5]*b[5]; t+=(dpint)v5*(dpint)p5;  c[4]=((spint)t & mask);  t>>=54;
	c[5] = (spint)t;
}

// Modular squaring, c=a*a  mod 2p
static inline void modsqr(const spint *a,spint *c) {
	dpint tot;
	dpint t=0;
	spint p5=0x9e4000000000u;
	spint q=((spint)1<<54u); // q is unsaturated radix 
	spint mask=(spint)(q-(spint)1);
	tot=(dpint)a[0]*a[0]; t=tot; spint v0=((spint)t & mask); t>>=54;
	tot=(dpint)a[0]*a[1]; tot*=2; t+=tot;  spint v1=((spint)t & mask); t>>=54;
	tot=(dpint)a[0]*a[2]; tot*=2; tot+=(dpint)a[1]*a[1]; t+=tot;  spint v2=((spint)t & mask); t>>=54;
	tot=(dpint)a[0]*a[3]; tot+=(dpint)a[1]*a[2]; tot*=2; t+=tot;  spint v3=((spint)t & mask); t>>=54;
	tot=(dpint)a[0]*a[4]; tot+=(dpint)a[1]*a[3]; tot*=2; tot+=(dpint)a[2]*a[2]; t+=tot;  spint v4=((spint)t & mask); t>>=54;
	tot=(dpint)a[0]*a[5]; tot+=(dpint)a[1]*a[4]; tot+=(dpint)a[2]*a[3]; tot*=2; t+=tot;  t+=(dpint)v0*p5;  spint v5=((spint)t & mask); t>>=54;
	tot=(dpint)a[1]*a[5]; tot+=(dpint)a[2]*a[4]; tot*=2; tot+=(dpint)a[3]*a[3]; t+=tot;  t+=(dpint)v1*p5;  c[0]=((spint)t & mask);  t>>=54;
	tot=(dpint)a[2]*a[5]; tot+=(dpint)a[3]*a[4]; tot*=2; t+=tot;  t+=(dpint)v2*p5;  c[1]=((spint)t & mask);  t>>=54;
	tot=(dpint)a[3]*a[5]; tot*=2; tot+=(dpint)a[4]*a[4]; t+=tot;  t+=(dpint)v3*p5;  c[2]=((spint)t & mask);  t>>=54;
	tot=(dpint)a[4]*a[5]; tot*=2; t+=tot;  t+=(dpint)v4*p5;  c[3]=((spint)t & mask);  t>>=54;
	tot=(dpint)a[5]*a[5]; t+=tot;  t+=(dpint)v5*p5;  c[4]=((spint)t & mask);  t>>=54;
	c[5] = (spint)t;
}

//copy
static inline void modcpy(const spint *a,spint *c) {
	int i;
	for (i=0;i<6;i++) {
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
	spint x[6];
	spint t0[6];
	spint t1[6];
	spint t2[6];
	spint t3[6];
	modcpy(w,x);
	modsqr(x,z);
	modmul(x,z,z);
	modsqr(z,z);
	modmul(x,z,z);
	modsqr(z,t0);
	modmul(x,t0,t1);
	modmul(x,t1,t0);
	modmul(t1,t0,t1);
	modsqr(t1,t2);
	modmul(t0,t2,t0);
	modcpy(t0,t2);
	modnsqr(t2,2);
	modsqr(t2,t3);
	modmul(t2,t3,t2);
	modmul(t3,t2,t3);
	modmul(t2,t3,t2);
	modmul(t3,t2,t3);
	modnsqr(t3,6);
	modmul(t2,t3,t2);
	modcpy(t2,t3);
	modnsqr(t3,12);
	modmul(t2,t3,t2);
	modcpy(t2,t3);
	modnsqr(t3,24);
	modmul(t2,t3,t2);
	modmul(z,t2,z);
	modmul(t0,z,t0);
	modmul(z,t0,z);
	modmul(t1,z,t1);
	modmul(t0,t1,t0);
	modsqr(t0,t2);
	modmul(z,t2,z);
	modsqr(z,t2);
	modmul(z,t2,t2);
	modmul(t0,t2,t0);
	modmul(t1,t0,t1);
	modmul(t0,t1,t0);
	modmul(z,t0,z);
	modmul(t1,z,t1);
	modmul(t0,t1,t0);
	modmul(t1,t0,t1);
	modmul(t0,t1,t0);
	modmul(t1,t0,t2);
	modmul(t0,t2,t0);
	modmul(t1,t0,t1);
	modnsqr(t1,63);
	modmul(t0,t1,t1);
	modnsqr(t1,64);
	modmul(t0,t1,t1);
	modnsqr(t1,64);
	modmul(t0,t1,t0);
	modnsqr(t0,60);
	modmul(z,t0,z);
}

//calculate inverse, provide progenitor h if available
static void modinv(const spint *x,const spint *h,spint *z) {
	spint s[6];
	spint t[6];
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
	const spint c[6]={0x34458f91ff5aceu,0x4a69f57c28713u,0x21774e0dc019e2u,0x2f6b2c1507af1du,0x13d1163e47fcc3u,0x3b1a7d5f0a1cu};
	modmul(m,c,n);
}

//Convert n back to normal form, m=redc(n) 
static void redc(const spint *n,spint *m) {
	int i;
	spint c[6];
	c[0]=1;
	for (i=1;i<6;i++) {
		c[i]=0;
	}
	modmul(n,c,m);
	(void)modfsb(m);
}

//is unity?
static int modis1(const spint *a) {
	int i;
	spint c[6];
	spint c0;
	spint d=0;
	redc(a,c);
	for (i=1;i<6;i++) {
		d|=c[i];
	}
	c0=(spint)c[0];
	return ((spint)1 & ((d-(spint)1)>>54u) & (((c0^(spint)1)-(spint)1)>>54u));
}

//is zero?
static int modis0(const spint *a) {
	int i;
	spint c[6];
	spint d=0;
	redc(a,c);
	for (i=0;i<6;i++) {
		d|=c[i];
	}
	return ((spint)1 & ((d-(spint)1)>>54u));
}

//set to zero
static void modzer(spint *a) {
	int i;
	for (i=0;i<6;i++) {
		a[i]=0;
	}
}

//set to one
static void modone(spint *a) {
	int i;
	a[0]=1;
	for (i=1;i<6;i++) {
		a[i]=0;
	}
	nres(a,a);
}

//set to integer
static void modint(int x,spint *a) {
	int i;
	a[0]=(spint)x;
	for (i=1;i<6;i++) {
		a[i]=0;
	}
	nres(a,a);
}

// Modular multiplication by an integer, c=a*b mod 2p
// uses special method for trinomials, otherwise Barrett-Dhem reduction
static inline void modmli(const spint *a,int b,spint *c) {
	spint p5=0x9e4000000000u;
	spint mask=((spint)1<<54u)-(spint)1;
	dpint t=0;
	spint q,h,r=0x67884a69f57c28;
	t+=(dpint)a[0]*(dpint)b; c[0]=(spint)t & mask; t=t>>54u;
	t+=(dpint)a[1]*(dpint)b; c[1]=(spint)t & mask; t=t>>54u;
	t+=(dpint)a[2]*(dpint)b; c[2]=(spint)t & mask; t=t>>54u;
	t+=(dpint)a[3]*(dpint)b; c[3]=(spint)t & mask; t=t>>54u;
	t+=(dpint)a[4]*(dpint)b; c[4]=(spint)t & mask; t=t>>54u;
	t+=(dpint)a[5]*(dpint)b; c[5]=(spint)t;
	
//Barrett-Dhem reduction
	h = (spint)(t>>38u);
	q=(spint)(((dpint)h*(dpint)r)>>64u);
	c[0]+=q;
	c[5]-=q*p5;
}

//Test for quadratic residue 
static int modqr(const spint *h,const spint *x) {
	spint r[6];
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
	R+=0x3cc3c33c5aa5a55au;
	w=R;
	c0=(~b)&(w+1);
	c1=b+w;
	for (i=0;i<6;i++) {
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
	R+=0x3cc3c33c5aa5a55au;
	w=R;
	c0=(~b)&(w+1);
	c1=b+w;
	for (i=0;i<6;i++) {
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
	spint s[6];
	spint y[6];
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
	a[5]=((a[5]<<n)) + (a[4]>>(54u-n));
	for (i=4;i>0;i--) {
		a[i]=((a[i]<<n)&(spint)0x3fffffffffffff) + (a[i-1]>>(54u-n));
	}
	a[0]=(a[0]<<n)&(spint)0x3fffffffffffff;
}

//shift right by less than a word. Return shifted out part
static int modshr(unsigned int n,spint *a) {
	int i;
	spint r=a[0]&(((spint)1<<n)-(spint)1);
	for (i=0;i<5;i++) {
		a[i]=(a[i]>>n) + ((a[i+1]<<(54u-n))&(spint)0x3fffffffffffff);
	}
	a[5]=a[5]>>n;
	return r;
}

//divide by 2. Shift right 1 bit (or add p and shift right one bit)
static void modhaf(spint *n) {
	int lsb;
	spint t[6];
	(void)prop(n);
	modcpy(n,t);
	lsb=modshr(1,t);
	n[0]-=(spint)1;
	n[5]+=((spint)0x9e4000000000u);
	(void)prop(n);
	modshr(1,n);
	modcmv(1-lsb,t,n);
}

//set a= 2^r
static void mod2r(unsigned int r,spint *a) {
	unsigned int n=r/54u;
	unsigned int m=r%54u;
	modzer(a);
	if (r>=40*8) return;
	a[n]=1; a[n]<<=m;
nres(a,a);
}

//export to byte array
static void modexp(const spint *a,char *b) {
	int i;
	spint c[6];
	redc(a,c);
	for (i=39;i>=0;i--) {
		b[i]=c[0]&(spint)0xff;
		(void)modshr(8,c);
	}
}

//import from byte array
//returns 1 if in range, else 0
static int modimp(const char *b, spint *a) {
	int i,res;
	for (i=0;i<6;i++) {
		a[i]=0;
	}
	for (i=0;i<40;i++) {
		modshl(8,a);
		a[0]+=(spint)(unsigned char)b[i];
	}
	res=(int)modfsb(a);
	nres(a,a);
	return res;
}

//determine sign
static int modsign(const spint *a) {
	spint c[6];
	redc(a,c);
	return c[0]%2;
}

//return true if equal
static int modcmp(const spint *a,const spint *b) {
	spint c[6],d[6];
	int i,eq=1;
	redc(a,c);
	redc(b,d);
	for (i=0;i<6;i++) {
		eq&=(((c[i]^d[i])-1)>>54)&1;
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
0x0};
const digit_t ONE[NWORDS_FIELD] = {0x67,
0x0,
0x0,
0x0,
0x0,
0x544000000000};
// Montgomery representation of 2^-1
static const digit_t TWO_INV[NWORDS_FIELD] = {0x33,
0x0,
0x0,
0x0,
0x0,
0x794000000000};
// Montgomery representation of 3^-1
static const digit_t THREE_INV[NWORDS_FIELD] = {0x15555555555577,
0x15555555555555,
0x15555555555555,
0x15555555555555,
0x15555555555555,
0x859555555555};

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

