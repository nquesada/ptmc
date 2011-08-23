#include"lj_params.h"
/*! \file lj_params.c
   \brief See documentation of lj_params.h - this file contains source code only, with function prototypes in lj_params.h
 */


inline double
lj_sigma2 (int i, int j)
{
  return 1;
}

inline double
lj_epsilon (int i, int j)
{
  return 1;
}

/*
#define sigmaAA2 1
#define sigmaAB2 1.15348
#define sigmaBB2 1  //the one with index m>
#define epsilonAA 1
#define epsilonAB 1.48
#define epsilonBB 1

#define m 1

inline double lj_sigma2(int i,int j){
  if(i<m&&j<m){
    return sigmaAA2;
  }
  else if(i>=m&&j>=m){
    return sigmaBB2;
  }
  else{
    return sigmaAB2;
  }
}

inline double lj_epsilon(int i,int j){
  if(i<m&&j<m){
    return epsilonAA;
  }
  else if(i>=m&&j>=m){
    return epsilonBB;
  }
  else{
    return epsilonAB;
  }
}
*/
