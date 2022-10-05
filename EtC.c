// 暗号化後圧縮

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "MT.h" // メルセンヌツイスタ法
#define rnd() (genrand_real1()) // 乱数生成マクロ[0, 1]
#define FLAMENUM 100 // フレーム数

double awgn(double sn);

typedef struct matrix{
    int weight;
    int column;
    int row;
}matrix;

void shuffle(int array[], int size) {
    int i, j;
    int tmp;
    /* シャッフル対象の末尾を設定 */
    i = size - 1;
    while (i > 0) {
        /* シャッフル対象(0〜i)から位置をランダム決定 */
        j = rand() % (i + 1);
        /* ランダムに決めた位置と
           シャッフル対象の末尾の位置のデータを交換 */
        tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
        /* シャッフル対象の範囲を狭める */
        i--;
    } 
}

/* 数値を昇順にソート */
void sort(int total, int number[]){
    for (int i=0; i<total; ++i) {
        for (int j=i+1; j<total; ++j) {
            if (number[i] > number[j]) {
                int tmp =  number[i];
                number[i] = number[j];
                number[j] = tmp;
            }
        }
    }
} 


void makeP(int array[], matrix mP0, int P0[mP0.column][mP0.weight]){
    int tmp=0;
    shuffle(array,mP0.row);
    for(int i=0;i<mP0.column;i++){
        for(int k=0;k<mP0.weight;k++){
            P0[i][k]=array[tmp];
            tmp++;
        }
        sort(mP0.weight,P0[i]);
    }
}

void makePt(int Pt[], int weight, int row, matrix mH, int H[mH.column][mH.weight]){//転置行列の作成
    // for(int l=0;l<weight;l++){
        int l=0;
        for(int i=0;i<mH.column;i++){
            for(int k=0;k<mH.weight;k++){
                if(H[i][k]==row){//同じなら値を入れて次のループへ
                    Pt[l]=i;
                    l=l+1;
                    break;
                }
                if(H[i][k]>row){//求めたいものよりも大きければ終了し次へ
                    break;
                }
            }
        }
    // }
}

void xor(int M,int X1[M],int X2[M],int X[M]){//排他的論理和 
    for(int m=0;m<M;m++){
        X[m]=X1[m]^X2[m];
    }
}

void mkplain(int M,int X[M]){//平文ブロック生成
    for(int m=0;m<M;m++){
        if(rnd() >= 0.5){
            X[m]= 1;
        } else{
            X[m]= 0;
        }
    }
}

double correlation(int M,int key[M],int encr[M]){
  int n;
  int err=0;

  double p;

  for(n=0;n<M;n++){
    if(key[n] != encr[n]){
      err++;
    }
  }

  p = ((double)err) / ((double)M);

  return p;
}

void compress(matrix mH ,int N,int M,int J,int EX[N][M],int H[mH.column][mH.weight],int s[N-1][J]){
    int m=0;
    for(int n=0;n<N-1;n++){
        for(int j=0;j<J;j++){
            s[n][j]=0;
            for(int d=0;d<mH.weight;d++){
                m=H[j][d];
                s[n][j]=s[n][j]+EX[n][m];
            }
            s[n][j]=s[n][j]%2;
        }

    }
}

int main(void){
    srand((unsigned int) time(NULL));
    init_genrand((unsigned)time(NULL));

    int N=3;//ブロック数
    int M=8;//ブロック長
    // int IV[M];//初期ベクトル
    int X[N][M];//平文ブロック X0初期ベクトル
    int EX[N][M];//暗号ブロック
    int n,m,j,k;

    // 平文ブロック生成
    for(n=0;n<N+1;n++){
        mkplain(M,X[n]);
    }
    for(n=0;n<N+1;n++){
        printf("%d\n",n);
        for(m=0;m<M;m++){
            printf("\t%d\n",X[n][m]);
        }
    }
    // // 初期ベクトル生成
    // mkplain(M,IV);

    // 暗号ブロック生成
    for(m=0;m<M;m++){
        EX[0][m]=X[0][m];
    }
    for(n=1;n<N+1;n++){
        xor(M,EX[n-1],X[n],EX[n]);
    }
    for(n=0;n<N+1;n++){
        printf("%d\n",n);
        for(m=0;m<M;m++){
            printf("\t%d\n",EX[n][m]);
        }
    }

    //Hの作成
    matrix mP0;
    mP0.weight=2;
    int R=0.5;//圧縮率
    mP0.row=M*R;
    mP0.column=mP0.row/2;
    int P0[mP0.column][mP0.weight];
    int C;//部分行列の数
    C=mP0.row/mP0.column;
    matrix mH;
    mH.weight=mP0.weight+1;
    mH.column=mP0.row;
    mH.row=2*mP0.row;
    
    int array[mP0.row];
    for (m=0; m< mP0.row; m++) {
        array[m] = m;
    }
    int H[mH.column][mH.weight];
    int d=mP0.row;
    for(j=0;j<C;j++){
        makeP(array, mP0, P0);
        for(n=j*mP0.column;n<(j+1)*mP0.column;n++){
            for(k=0;k<mP0.weight;k++){
                H[n][k]=P0[n-(j*mP0.column)][k];
            }
            H[n][mP0.weight]=d;
            d=d+1;
        }
    }

    for(n=0;n<mH.column;n++){
        for(k=0;k<mH.weight;k++){
            printf("%d\t",H[n][k]);
        }
        printf("\n");
    }
    //圧縮
    int J=mH.column;
    int s[N-1][J];//圧縮ブロック

    compress(mH,N,M,J,EX,H,s);
    for(n=0;n<N;n++){
        printf("%d\n",n);
        for(j=0;j<J;j++){
            printf("\t%d\n",s[n][j]);
        }
    }

    




    return 0;

}