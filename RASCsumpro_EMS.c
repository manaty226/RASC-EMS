#include<math.h>

int PRINTSWITCH = 1;
LLR **metric;
Complex *Xl;

#ifndef RASC_ADD
int RASCadd(int a, int b)
{

  int c;
  int i,n_a,n_b;
  int tmp_a,tmp_b,tmp_c;


  c = 0;
  for(i = 0; i < 2 * Nb; i++) {
    tmp_a = a % L;
    tmp_b = b % L;
    a /= L;
    b /= L;
    tmp_c = (tmp_a+tmp_b) % L;
    c += tmp_c * (int)pow(L,i);
  }


  return c;
}

int RASCinvMod(int a) {

  int c;
  int i;
  int tmp;

  c = 0;
  for(i = 0; i < 2 * Nb; i++) {
    tmp = a % L;
    a /= L;

    tmp = (-tmp + L) % L;
    c += tmp * (int)pow(L,i);
  }

  return c;
}
#define RASC_ADD
#endif


void Sort(SORTER *sorter)
{

  int i,j,k;
  SORTER temp;


  for(i = 1; i < Nc; i++) {
    k = i;
    while(k > 0 && sorter[k].l.llr > sorter[k-1].l.llr) {
      temp = sorter[k-1];
      sorter[k-1] = sorter[k];
      sorter[k] = temp;
      k--;
    }
  }
}

void llrsorting(LLR *sorter, int size)
{

  int i,j,k;
  LLR temp;

  for(i = 1; i < size; i++) {
    k = i;
    while(k > 0 && (sorter+k)->llr > (sorter+(k-1))->llr) {
      temp = *(sorter+(k-1));
      *(sorter+(k-1)) = *(sorter+k);
      *(sorter+k) = temp;
      k--;
    }
  }
}

void EMS(Message *I, Message U)
{

  int i,j,k;
  int matchflag;
  double gamma;
  SORTER sorter[Nc];
  Message V;
  int index[Nc];
  int endflag = 0;


  /* printf("EMS \n"); */
  for(i = 0; i < Nc; i++) {
    sorter[i].l.llr = I->llr[i].llr + U.llr[0].llr;
    sorter[i].l.symbol = RASCinvMod(RASCadd(I->llr[i].symbol,U.llr[0].symbol));
    sorter[i].index = i;
  }

  for(i = 0; i < Nc; i++) {
    V.llr[i].symbol = -1;
    index[i] = 1;
  }

  k = 0;
  for(i = 0; i < 2 * Nc; i++) {
    matchflag = 0;
    for(j = 0; j < Nc; j++) {
      /* printf("%d %d \n",V.llr[j].symbol,sorter[0].l.symbol); */
      if(V.llr[j].symbol == sorter[0].l.symbol) {
        matchflag = 1;
      }
    }
    if(matchflag != 1) {
      V.llr[k] = sorter[0].l;
      k++;
    }

    if(index[sorter[0].index] < Nc) {
      sorter[0].l.llr = I->llr[sorter[0].index].llr + U.llr[index[sorter[0].index]].llr;
      sorter[0].l.symbol = RASCinvMod(RASCadd(I->llr[sorter[0].index].symbol,U.llr[index[sorter[0].index]].symbol));
      index[sorter[0].index]++;
      Sort(sorter);
    }
    else if(sorter[0].index < Nc-1) {
      sorter[0].l.llr = I->llr[sorter[0].index+1].llr + U.llr[index[sorter[0].index+1]].llr;
      sorter[0].l.symbol = RASCinvMod(RASCadd(I->llr[sorter[0].index+1].symbol,U.llr[index[sorter[0].index+1]].symbol));
      sorter[0].index++;
      index[sorter[0].index]++;
      Sort(sorter);
    }
    else {
      endflag = 1;
      break;
    }

    if(k > Nc-1) {
      endflag = 2;
      break;
    }

  }

  if(endflag == 0 || endflag == 1) {
    gamma = V.llr[k-1].llr - log(STNUM-(k-1));
    while(k-1 < Nc) {
      for(i = 0; i < STNUM; i++) {
        matchflag = 0;
        for(j = 0; j < Nc; j++) {
          if(i == V.llr[j].symbol) {
            matchflag = 1;
            break;
          }
        }
        if(matchflag != 1) {
          V.llr[k].symbol = i;
          V.llr[k].llr = gamma;
          break;
        }
      }
      k++;
    }
  }

  for(i = 0; i < Nc; i++) {
    if(V.llr[i].symbol != -1) {
      I->llr[i] = V.llr[i];
    }
  }

  if(endflag == 0 || endflag == 1) {
    I->gamma = gamma - ALPHA;
  }
  else {
    I->gamma = V.llr[Nc-1].llr - log(STNUM-(Nc-1)) - ALPHA;
  }



}

void CNproc(Message *I, int row, int column, Sump **V)
{
  int i,j,k;

  if(row < 0 || column < 0) {
    printf("行列の番号が不正 @ RowMsg\n");
    exit(1);
  }

  for(i=0; i<CW; i++) {
    if(V[column][i].nmb == row) {
      EMS(I,V[column][i].msg);
    }
  }
}

void CHproc(Message *I, int row,int column,LLR **metric, int metric_sort[][CNUM])
{
  int i,j,k;
  int index[STNUM];
  int match = -1;
  LLR T[2*Nc];


  if(row < 0 || column < 0) {
    printf("行列の番号が不正 @ RowMsg\n");
    exit(1);
  }

  for(i = 0; i < 2*Nc; i++) {
    T[i].llr = -999999;
    T[i].symbol = -1;
  }


  for(i = 0; i < STNUM; i++) {
    index[i] = 0;
  }

  for(j = 0; j < Nc; j++) {
    if((I->llr[j].symbol != -1) && (I->llr[j].symbol < STNUM)) {
      T[j].llr = I->llr[j].llr + metric[I->llr[j].symbol][column].llr;
      T[j].symbol = I->llr[j].symbol;
      // index[I->llr[j].symbol] = 1;
    } else if (I->llr[j].symbol == -1){

    } else {
      printf("error %d \n",T[j].symbol);
    }
  }

  k = 0;
  for(j = 0; j < STNUM && k < Nc; j++) {
      T[Nc+k].llr = I->gamma + metric[metric_sort[j][column]][column].llr;
      T[Nc+k].symbol = metric[metric_sort[j][column]][column].symbol;
      if(T[Nc+k].symbol > STNUM) {
        printf("error \n");
        printf("error symbol %d \n",T[Nc+j].symbol);
      }
      k++;
  }

  llrsorting(T,2*Nc);


  for(i = 0; i < Nc; i++) {
    I->llr[i].llr = T[i].llr - T[Nc-1].llr;
    I->llr[i].symbol = T[i].symbol;
    if(T[i].symbol > STNUM) {
      printf("error %d \n",T[i].symbol);
    }

  }

  I->gamma = T[Nc].llr - T[Nc-1].llr - log(STNUM-Nc) - ALPHA;

}


void VNproc(Message *I, Sump **A, int row, int column) {
  int i,j,k;
  int index[Nc];
  int match=-1;
  LLR T[2*Nc];

  if(row < 0 || column < 0) {
    printf("invalid row number @ RowMsg\n");
    exit(1);
  }

  for(k =0; k < Nc; k++) {
    index[k] = 0;
  }

  for(i = 0; i < 2*Nc; i++) {
    T[i].llr = -999999;
    T[i].symbol = -1;
  }


  for(i=0; i<RW; i++) {
    if(A[row][i].nmb == column) {
      for(j = 0; j < Nc; j++) {
        match = -1;
        for(k = 0; k < Nc; k++) {
          if(I->llr[j].symbol == A[row][i].msg.llr[k].symbol) {
            match = 1;
            T[j].llr = I->llr[j].llr + A[row][i].msg.llr[k].llr;
            T[j].symbol = I->llr[j].symbol;
            index[k] = 1;
            if(T[j].symbol > STNUM) {
              printf("error %d \n",T[j].symbol);
            }
          }
        }
        if(-1 == match) {
          T[j].llr = I->llr[j].llr + A[row][i].msg.gamma;
          T[j].symbol = I->llr[j].symbol;
          if(T[j].symbol > STNUM) {
            printf("error %d \n",T[j].symbol);
          }
        }
      }

      for(j = 0; j < Nc; j++) {
        if(0 == index[j]) {
          T[Nc+j].llr = A[row][i].msg.llr[j].llr + I->gamma;
          T[Nc+j].symbol = A[row][i].msg.llr[j].symbol;
          if(T[Nc+j].symbol > STNUM) {
            printf("error \n");
            printf("error symbol %d \n",T[Nc+j].symbol);
          }
        }
      }
    }
  }


  llrsorting(T,2*Nc);


  for(i = 0; i < Nc; i++) {
    I->llr[i].llr = T[i].llr - T[Nc-1].llr;
    I->llr[i].symbol = T[i].symbol;
    if(T[i].symbol > STNUM) {
      printf("error %d \n",T[i].symbol);
    }

  }

  I->gamma = T[Nc].llr - T[Nc-1].llr - log(STNUM-Nc) - ALPHA;

}

void VNproc_hidden(Message *I, Sump **A, int row, int column) {
  int i,j,k;
  int index[Nc];
  int match=-1;
  LLR T[2*Nc];
  int last;
  int symidx[CONSTSIZE];

  if(row < 0 || column < 0) {
    printf("行列の番号が不正 @ RowMsg\n");
    exit(1);
  }

  for(k =0; k < Nc; k++) {
    index[k] = 0;
  }
  for(i = 0; i < CONSTSIZE; i++) {
    symidx[i] = 0;
  }

  for(i = 0; i < 2*Nc; i++) {
    T[i].llr = -999999;
    T[i].symbol = -1;
  }


  for(i=0; i<RW; i++) {
    if(A[row][i].nmb == column) {
      for(j = 0; j < Nc; j++) {
        match = -1;
        for(k = 0; k < Nc; k++) {
          if(I->llr[j].symbol == A[row][i].msg.llr[k].symbol) {
            match = 1;
            T[j].llr = I->llr[j].llr + A[row][i].msg.llr[k].llr;
            T[j].symbol = I->llr[j].symbol;
            index[k] = 1;
            if(T[j].symbol > STNUM) {
              printf("error %d \n",T[j].symbol);
            }
          }
        }
        if(-1 == match) {
          T[j].llr = I->llr[j].llr + A[row][i].msg.gamma;
          T[j].symbol = I->llr[j].symbol;
          if(T[j].symbol > STNUM) {
            printf("error %d \n",T[j].symbol);
          }
        }
        if(T[j].symbol > CONSTSIZE-1) {
          T[j].llr = -99;
        }
      }

      for(j = 0; j < Nc; j++) {
        if(0 == index[j]) {
          T[Nc+j].llr = A[row][i].msg.llr[j].llr + I->gamma;
          T[Nc+j].symbol = A[row][i].msg.llr[j].symbol;
          if(T[Nc+j].symbol > STNUM) {
            printf("error \n");
            printf("error symbol %d \n",T[Nc+j].symbol);
          }
          if(T[Nc+j].symbol > CONSTSIZE-1) {
            T[Nc+j].llr = -99;
          }
        }
      }
    }
  }

  llrsorting(T,2*Nc);


  for(i = 0; i < Nc; i++) {
    if(T[i].symbol > CONSTSIZE) {
      last = i;
      break;
    }
    else {
      symidx[T[i].symbol] = 1;
    }
  }

  for(i = last; i < Nc; i++) {
    if(last != 0) {
      for(j = 0; j < CONSTSIZE + 1; j++) {
        if( (j < CONSTSIZE) && (symidx[j] == 0)) {
          symidx[j] = 1;
          break;
        }
      }
      if(j < CONSTSIZE) {
        T[i].llr = T[last-1].llr - log(STNUM-last) - ALPHA;
        T[i].symbol = j;
      }
      else {
        break;
      }
    }
    else {
      for(j = 0; j < CONSTSIZE+1; j++) {
        if(j < CONSTSIZE && symidx[j] == 0) {
          symidx[j] = 1;
          break;
        }
      }
      if(j < CONSTSIZE) {
        T[i].llr = - log(STNUM-last) - ALPHA;
        T[i].symbol = j;
      }
      else {
        break;
      }
    }
  }


  for(i = 0; i < Nc; i++) {
    if(last != 0) {
      I->llr[i].llr = T[i].llr - T[last-1].llr;
    }
    else {
      I->llr[i].llr = T[i].llr;
    }
    I->llr[i].symbol = T[i].symbol;
    if(T[i].symbol > STNUM) {
      printf("error VN hidden= %d last = %d, i = %d \n",T[i].symbol,last,i);
    }

  }

  I->gamma = -999;

}

void VNmessageset(Message *I, int row, int column, Sump **A)
{

  int i,j,k;

  if(row < 0 || column < 0) {
    printf("行列の番号が不正 @ ColMsg\n");
    exit(1);
  }

  for(i=0; i<RW; i++) {
    if(A[row][i].nmb == column) {
      *I = A[row][i].msg;
    }
  }
}


void CNmessageset(Message *I, int row, int column, Sump **B)
{
  int i,j,k;

  if(row < 0 || column < 0) {
    printf("行列の番号が不正 @ ColMsg\n");
    exit(1);
  }

  for(i=0; i<CW; i++) {
    if(B[column][i].nmb == row) {
      *I = B[column][i].msg;
    }
  }

}



void RASCMod(Complex X[]) {
  int i,j,k;

  for(k=0; k<Nb; k++) {
    X[k].re = fmod(fmod(X[k].re,L)+L,L);
    X[k].im = fmod(fmod(X[k].im,L)+L,L);
  }

}


void RASCDem(Complex Y[],LLR **metric,Complex s[][STNUM],double noise, Complex *replica, int metric_sort[][CNUM]) {

  int i,j,k;
  int met_idx;
  Complex temp;
  double metzero,sigma,mettemp;
  LLR llrtemp;
  double max,sum;
  double met_sorttemp[Nc];

  sigma = 2*pow(noise,2);


  for(i=SNUM; i<CNUM; i++) {
    max = 0;
    sum = 0;
    for(j = 0; j < STNUM; j++) {
      metric[j][i].llr = -999;
      metric[j][i].symbol = -1;
      if(j < Nc) {
        met_sorttemp[j] = -999;
      }
    }
    for(j=0; j<STNUM; j++) {
      temp = csub(Y[i],replica[j]);
      mettemp = -(pow(temp.re,2)+pow(temp.im,2))/(sigma);

      if(j == 0)
      metzero = mettemp;

      mettemp -= metzero;
      metric[j][i].llr = mettemp;
      metric[j][i].symbol = j;

      if(met_sorttemp[Nc-1] < mettemp) {
        met_sorttemp[Nc-1] = mettemp;
        metric_sort[Nc-1][i] = j;
        k = Nc-1;
        while(k > 0 && met_sorttemp[k] > met_sorttemp[k-1]) {
          mettemp = met_sorttemp[k-1];
          met_idx = metric_sort[k-1][i];

          met_sorttemp[k-1] = met_sorttemp[k];
          metric_sort[k-1][i] = metric_sort[k][i];

          met_sorttemp[k] = mettemp;
          metric_sort[k][i] = met_idx;
          k -= 1;
        }
      }
    }

  }

}

//
// void InvLLR(double prob[][CNUM],double llr[][CNUM]) {
//   int i,j;
//   double sum;
//
//   for(j=0; j<CNUM; j++) {
//     sum = 0.;
//     for(i=0; i<STNUM; i++) {
//       prob[i][j] = exp(llr[i][j]);
//       sum += prob[i][j];
//     }
//     for(i=0; i<STNUM; i++) {
//       prob[i][j] /= sum;
//     }
//   }
//
// }

// void Permutate(double org[], int perm[]) {
//   int i,j,k;
//   double copy[STNUM];
//
//   for(i=0; i<STNUM; i++) {
//     copy[i] = org[i];
//   }
//
//   for(i=0; i<STNUM; i++) {
//     org[perm[i]] = copy[i];
//   }
//
// }
//
// void Depermutate(double org[], int perm[]) {
//   int i,j,k;
//   double copy[STNUM];
//
//   for(i=0; i<STNUM; i++) {
//     copy[i] = org[i];
//   }
//
//   for(i=0; i<STNUM; i++) {
//     org[i] = copy[perm[i]];
//   }
//
// }

void RASCsumpro(Complex Y[], Complex R[], Complex F[], double noise, int inter[],Sump **A, Sump **B, Complex *replica) {
  int i,j,k,n,u,t,s,decnum,symnum,CNcnt,VNcnt;
  int matchcount,purm[STNUM],depurm[STNUM];
  int metric_sort[Nc][CNUM];
  Complex conj,fbtemp;
  Complex bs[STNUM][Nb],cs[Nb];
  Complex state[Nb][STNUM];
  double sum,max;
  int maxnum,count;
  Message intmdt,out;

  metric = (LLR **)calloc(STNUM, sizeof(LLR *));
  metric[0] = calloc((STNUM*CNUM), sizeof(LLR));
  for(i=1; i<STNUM; i++) metric[i] = metric[0]+i*CNUM;


  conj.re = 0;
  conj.im = 1;

  /* 状態の初期化 */
  for(i=0; i<STNUM; i++) {
    n = i;
    for(j=0; j<Nb; j++) {
      state[j][i].re = n%L;
      n /= L;
      state[j][i].im = n%L;
      n /= L;
    }
  }

  for(i=0; i<STNUM; i++) {
    /* gx_i-1の計算 */
    for(u=0; u<Nb; u++) {
      bs[i][u].re = 0;
      bs[i][u].im = 0;
      for(t=0; t<Nb; t++) {
        for(s=0; s<Nb; s++) {
          if((s+t)%Nb == u) {
            if((s+t) == u)
            bs[i][u] = cadd(bs[i][u],cmul(state[s][i],F[t]));
            else
            bs[i][u] = cadd(bs[i][u],cmul(conj,cmul(state[s][i],F[t])));
          }
        }
      }
    }
    for(j=0; j<STNUM; j++) {
      RASCMod(bs[i]);
      matchcount = 0;
      for(u=0; u<Nb; u++) {
        matchcount += pow(L,2*u)*bs[i][u].re;
        matchcount += pow(L,2*u+1)*bs[i][u].im;
      }
      if(matchcount == j) {
        purm[i] = j;
        depurm[j] = i;
      }
    }
  }

// Calculate channel LLRs
  RASCDem(Y,metric,state,noise,replica,metric_sort);

  /* Initialize messages */
  for(i=0; i<CNUM; i++) {
    for(j=0; j<CW; j++) {
      for(k=0; k<Nc; k++) {
        B[i][j].msg.llr[k].llr = -9999999;
        B[i][j].msg.llr[k].symbol = -1;
        B[i][j].msg.gamma = 0;
      }
    }
    if(i < NK) {
      for(j=0; j<RW; j++) {
        if(A[i][j].nmb != -1) {
          for(k=0; k<Nc; k++) {
            A[i][j].msg.llr[k].llr = -999999;
            A[i][j].msg.llr[k].symbol = -1;
          }
          A[i][j].msg.gamma = 0;
        }
      }
    }
  }

  count = 0;
  for(decnum=0; decnum<DECNUM; decnum++) {
    /* Variable node update for observation nodes */
    for(i=SNUM; i<CNUM; i++) {
      for(j=0; j<CW; j++) {
        if(B[i][j].nmb != -1) {
          VNcnt = 0;
          for(k=0; k<Nc; k++) {
            intmdt.llr[k].llr = -9999999;
            intmdt.llr[k].symbol = -1;
            intmdt.gamma = 0;
          }

          for(k=0; k<CW; k++) {
            if(B[i][k].nmb != B[i][j].nmb && B[i][k].nmb != -1) {
              if(VNcnt == 0) {
                VNmessageset(&intmdt,B[i][k].nmb,i,A);
                CHproc(&intmdt,B[i][k].nmb,i,metric,metric_sort);
                VNcnt = 1;
              }
            }
          }
          if(i == CNUM-1) {
            CHproc(&intmdt,B[i][k].nmb,i,metric,metric_sort);
          }
          B[i][j].msg = intmdt;
        }
      }
    }


// permutation
    for(i=SNUM; i<CNUM; i++) {
      for(j=0; j<CW; j++) {
        if(B[i][j].pflag == 1) {
          for(k = 0; k < Nc; k++) {
            B[i][j].msg.llr[k].symbol = purm[B[i][j].msg.llr[k].symbol];
          }
        }
      }
    }


// Check node update for hidden nodes
    for(i=0; i<NK; i++) {
      for(j=0; j<RW; j++) {
        if(A[i][j].nmb != -1 && A[i][j].nmb < SNUM) {
          CNcnt = 0;
          for(k=0; k<RW; k++) {
            if(A[i][k].nmb != A[i][j].nmb && A[i][k].nmb != -1) {
              if(CNcnt == 0) {
                CNmessageset(&intmdt,i,A[i][k].nmb,B);
                CNcnt = 1;
              }
              else {
                CNproc(&intmdt,i,A[i][k].nmb,B);
              }
            }
          }
          if(0 == i) {
            for(k = 0; k < Nc; k++) {
              intmdt.llr[k].symbol = RASCinvMod(intmdt.llr[k].symbol);
            }
            A[i][j].msg = intmdt;
          }
          else {
            A[i][j].msg = intmdt;
          }
        }
      }
    }


// Variable node update for hidden nodes
    for(i=0; i<SNUM; i++) {
      for(j=0; j<CW; j++) {
        if(B[i][j].nmb != -1) {
          for(k=0; k<Nc; k++) {
            intmdt.llr[k].llr = -999;
            intmdt.llr[k].symbol = -1;
            intmdt.gamma = 0;
          }

          for(k=0; k<CW; k++) {
            if(B[i][k].nmb != B[i][j].nmb && B[i][k].nmb != -1) {
              #ifndef FULLRATE
              VNproc_hidden(&intmdt,A,B[i][k].nmb,i);
              #else
              VNproc(&intmdt,A,B[i][k].nmb,i);
              #endif
            }
          }
          B[i][j].msg = intmdt;
        }
      }

    }

// Check node update for observation nodes
    for(i=0; i<NK; i++) {
      for(j=0; j<RW; j++) {
        if(A[i][j].nmb != -1 && A[i][j].nmb > SNUM-1) {
          CNcnt = 0;
          for(k=0; k<RW; k++) {
            if(A[i][k].nmb != A[i][j].nmb && A[i][k].nmb != -1) {
              if(CNcnt == 0) {
                CNmessageset(&intmdt,i,A[i][k].nmb,B);
                CNcnt = 1;
              }
              else {
                CNproc(&intmdt,i,A[i][k].nmb,B);
              }
            }
          }
          if(0 == i) {
            for(k = 0; k < Nc; k++) {
              intmdt.llr[k].symbol = RASCinvMod(intmdt.llr[k].symbol);
            }
            A[i][j].msg = intmdt;
          }
          else {
            A[i][j].msg = intmdt;
          }


        }
      }
    }

// depermutation
    for(i=0; i<NK; i++) {
      for(j=0; j<RW; j++) {
        if(A[i][j].nmb != -1 && A[i][j].nmb > SNUM-1) {
          if(A[i][j].pflag == 1) {
            for(k = 0; k < Nc ; k++) {
              A[i][j].msg.llr[k].symbol = depurm[A[i][j].msg.llr[k].symbol];
            }
          }
        }
      }
    }

// Tentative decision
    for(i=0; i<SNUM; i++) {
      for(k=0; k<Nc; k++) {
        intmdt.llr[k].llr = -9999999;
        intmdt.llr[k].symbol = -1;
        intmdt.gamma = 0;
      }
      for(j=0; j<CW; j++) {
        if(B[i][j].nmb != -1) {
          VNproc(&intmdt,A,B[i][j].nmb,i);

        }
      }
      R[i] = cadd(state[0][intmdt.llr[0].symbol],cmul(state[1][intmdt.llr[0].symbol],Polar(1,(double)1./(2.*Nb))));
    }
  }

  free(metric[0]);
  free(metric);
}
