#ifndef BINARY_H
#define BINARY_H
#include <cstdint>

//binary
extern uint32_t binact[32];
extern uint8_t binnum[256];
extern uint8_t binaryReverseMatrix16bit[65536][17];
#define BINARY_SIZE 4

#define MaxBinaryLength 3000
#define binaryStrategy 0 // which 0 if we use 32bit and is 1 if we use 64bit 
#define binarySize 32
#define binaryShift 5
#define binaryR 31
#define blockLength 16
#define binary_type unsigned int

//#define binaryArrlength(len) ((len)%32 == 0)?((len)/32):(((len)/32)+1)
#define binaryArrlength(len) (((len)/32)+1)
#define binaryBlockArrlength(len,size) (((len)/size)+1)

// Macros
#define addbinSingle(f,i) (f)|=binact[(i)]
#define removebinSingle(f,i) (f)&=~binact[(i)]
#define checkbinSingle(f,i) (((f)&binact[(i)])?1:0)

#define addbin(f,i)  f[(i)>>binaryShift]|=binact[(i)&binaryR];
#define removebin(f,i) f[(i)>>binaryShift]&=~binact[(i)&binaryR];
#define checkbin(f,i) (((f[(i)>>binaryShift]&binact[(i)&binaryR]))?1:0)
#define check2binOR(f,g,i) ((((f[(i)>>binaryShift]|g[(i)>>binaryShift])&binact[(i)&binaryR]))?1:0)
#define check3binOR(f,g,h,i) ((((f[(i)>>binaryShift]|g[(i)>>binaryShift]|h[(i)>>binaryShift])&binact[(i)&binaryR]))?1:0)
#define check2binAND(f,g,i) ((((f[(i)>>binaryShift]&g[(i)>>binaryShift])&binact[(i)&binaryR]))?1:0)
#define check3binAND(f,g,h,i) ((((f[(i)>>binaryShift]&g[(i)>>binaryShift]&h[(i)>>binaryShift])&binact[(i)&binaryR]))?1:0)

#define binary1st16Bit(f) ((uint16_t) (f))
#define binary2nd16Bit(f) ((f)>>16)

#define binary1st8Bit(f) ((uint8_t) (f))
#define binary2nd8Bit(f) (((uint16_t) (f))>>8)
#define binary3rd8Bit(f) ((uint8_t) ((f)>>16))
#define binary4th8Bit(f) ((f)>>24)

// binary reverse macros
#define binaryItems16(f,block,job) ((((block&1)==0)?(binaryReverseMatrix16bit[binary1st16Bit(f)][job]):(binaryReverseMatrix16bit[binary2nd16Bit(f)][job]))+((job!=16)?16*block:0))
#define binaryNumItems16(f,position) binaryItems16(f,2*position, 16) + binaryItems16(f,2*position+1,16) 
#define batchbinaryNumItems16(f,size,copyarrayi,result) result=0; for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {result += binaryNumItems16(f[copyarrayi],copyarrayi);}
#define convertBintoArr(f_bin,f,size,j,k,tmp,nactbin)	size =  0;  for (k = 0; k<nactbin*2; k++){tmp = binaryItems16(f_bin[(int)k/2],k,16);for (j = 0; j<tmp; j++){f[size++] = binaryItems16(f_bin[(int)k/2],k,j); }}


// binary merge macro


#define BINARYADD(bin,addition,size,copyarrayi) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {bin[copyarrayi] |= addition[copyarrayi];}
#define BINARYREMOVE(bin,exclusion,size,copyarrayi) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {bin[copyarrayi] &=~ exclusion[copyarrayi];}
#define BINARYUNION(uni,bin1,bin2,size,copyarrayi) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {uni[copyarrayi] = bin1[copyarrayi] | bin2[copyarrayi];}
#define BINARYINTERSECT(intersect,bin1,bin2,size,copyarrayi) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {intersect[copyarrayi] = bin1[copyarrayi] & bin2[copyarrayi];}
#define BINARYEXCLUDE_bin2_from_bin1(exc,bin1,bin2,size,copyarrayi) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {exc[copyarrayi] =(bin1[copyarrayi] & (~bin2[copyarrayi]));}
#define BINARYCHECK(bin1,bin2,size, copyarrayi,result) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {if ((bin1[copyarrayi] & bin2[copyarrayi])!=bin2[copyarrayi]) {result = false;break;}} if (copyarrayi == size) result = true;
#define BINARYCMP(bin1,bin2,size, copyarrayi,result) for(copyarrayi = 0; copyarrayi<(size); copyarrayi++) {if (bin1[copyarrayi] != bin2[copyarrayi]) {result = false;break;}} if (copyarrayi == size) result = true;

void initBinary();
void deleteBinary();
void print_binary(uint32_t b);
char* binaryArrtoIntArrtoString(uint32_t * b, int nact, int nactbin, char* s);
void binaryArrtoIntArr(uint32_t * b, int* Arr, int nactbin);
#endif //BINARY_H