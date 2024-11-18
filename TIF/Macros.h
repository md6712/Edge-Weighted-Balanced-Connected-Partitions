#ifndef MACROS_H
#define MACROS_H

// max min macros
#define MAX(a,b)  ((a)>(b)? (a):(b))
#define MIN(a,b)  ((a)<(b)? (a):(b))
#define SQR(x)    ((x)*(x))

#define ROUND3(x) (round( x * 1000.0 ) / 1000.0)

// array macro
#define CHECK_VALUE(k,value) (((k)==(value))?1:0)\

#define _PRINT_LOG_L0(prob,a,b,...) _PrintL0(a,b,__VA_ARGS__);_LogL0(prob,a,b,__VA_ARGS__);
#define _PRINT_LOG_L1(prob,a,b,...) _PrintL1(a,b,__VA_ARGS__);_LogL1(prob,a,b,__VA_ARGS__);
#define _PRINT_LOG_L2(prob,a,b,...) _PrintL2(a,b,__VA_ARGS__);_LogL2(prob,a,b,__VA_ARGS__);
#define _PRINT_LOG_L3(prob,a,b,...) _PrintL3(a,b,__VA_ARGS__);_LogL3(prob,a,b,__VA_ARGS__);
#define _PRINT_LOG_L4(prob,a,b,...) _PrintL4(a,b,__VA_ARGS__);_LogL4(prob,a,b,__VA_ARGS__);
#define _PRINT_LOG_L5(prob,a,b,...) _PrintL5(a,b,__VA_ARGS__);_LogL5(prob,a,b,__VA_ARGS__);

#define fscanfCheck(checkscanf,fp,s,ch,...) {checkscanf  =0; while(checkscanf != 1){checkscanf = fscanf(fp, s, __VA_ARGS__);if (checkscanf!= 1) fscanf(fp, "%c", ch);}}

#endif