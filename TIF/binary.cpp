#include "binary.h"
#include <stdio.h>

uint32_t binact[32];
void initBinaryActBin();
uint8_t copyarrayi;
uint8_t binaryReverseMatrix16bit[65536][17];
uint8_t binnum[256];

void initBinary() {
	initBinaryActBin();
}
void initBinaryActBin() {
	binact[0] = 1;
	for (int i = 1; i < 32; i++) binact[i] = binact[i - 1] << 1;

	int num;
	for (uint32_t i = 0; i < 256; i++) {
		num = 0;
		for (int j = 0; j < 8; j++) { if (checkbinSingle(i, j)) { num++; } }
		binnum[i] = num;
	}

	for (int k = 0; k < 65536; k++)
	{
		int count = 0;
		for (int j = 0; j < 16; j++) { if (checkbinSingle(k, j)) binaryReverseMatrix16bit[k][count++] = j; }
		binaryReverseMatrix16bit[k][16] = count;
	}
}
void deleteBinary() {
	//delete[] binaryReverseMatrix16bit;
}

void print_binary(uint32_t b) {
	printf("\n");
	for (int j = 31; j >= 0; j--) { if (checkbinSingle(b, j)) { printf("1"); } else { printf("0"); } if (j % 8 == 0) { printf(" "); } }
}

uint8_t calc_binary_n_func(uint32_t b) {
	return binnum[(uint8_t)(b)] + binnum[(uint8_t)(b >> 8)] + binnum[(uint8_t)(b >> 16)] + binnum[(uint8_t)(b >> 24)];
}

void binaryArrtoIntArr(uint32_t* b, int* Arr, int nactbin) {
	int size, j, k, tmp;
	convertBintoArr(b, Arr, size, j, k, tmp, nactbin);
}

char* binaryArrtoIntArrtoString(uint32_t* b, int nact, int nactbin, char* s) {
	//char s[500];
	//int size, j, k, tmp;
	//uint8_t* B = new uint8_t[nact];
	//convertBintoArr(b, B, size, j, k, tmp, nactbin);
	//k = sprintf(s, " ");
	//for (j = 0; j < size; j++) {
	//	//k += sprintf_s(s + k, "%d ", B[j]);
	//}
	//delete[] B;
	return s;
}

