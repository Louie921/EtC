#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MT.h" // メルセンヌツイスタ法

#define rnd() (genrand_real1()) // 乱数生成マクロ[0, 1]

double awgn(double sn);

// Box-Muller法に基づいて，各SNRに対するAWGNを作成する関数
double awgn(double sn)
{
  double nvar;
  double x,y;
  double noise;
  
  nvar= pow( 10.0,(-1.0/10.0*sn)); // 分散の計算
  x=0.0;
  while(x==0.0)	x=rnd();
  y=rnd();
  noise=sqrt(-2.0*nvar*log(x)) * cos(2.0*M_PI*y); // AWGN生成
  
  return noise;
}
