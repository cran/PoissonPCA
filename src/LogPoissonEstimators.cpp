/*

   Copyright 2018 Toby Kenney

   This file is part of the R package PoissonPCA.

   PoissonPCA is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   PoissonPCA is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with PoissonPCA.  If not, see <http://www.gnu.org/licenses/>.

 */

#include <R.h>

#define NUMTERMS 9
#define NUMLOWVALS 20
  
long double ankdiag[NUMTERMS];
long double *lowvals;
long double* lowgvalues;

unsigned int numlowgvals;

class vallistelt{
  int ref;
  long double akd[NUMTERMS];
  unsigned int numlowgvals;
  long double* lowgvalues;
  long double* lowvals;
  vallistelt *next;
  friend class vallist;
  unsigned int checkvalue(int refno);
};

unsigned int vallistelt::checkvalue(int refno){
  //Checks if reference number is already used.
  if(ref==refno){
    return 1;
  }else if(next==NULL){
    return 0;
  }else{
    return next->checkvalue(refno);
  };
};

class vallist{
  static int lastref;
  static int currentref;
  static vallistelt *start;
public:
  static int insert(long double ak[],long double lwEC[] ,unsigned int nlg,long double *lwg);
  static void insert(long double ak[],long double lwEC[] ,unsigned int nlg,long double *lwg,int refno);
  static void lookupvalues(int ref);
};


int vallist::lastref=0;
int vallist::currentref=0;
vallistelt *vallist::start=NULL;

int vallist::insert(long double ak[],long double * lwEC ,unsigned int nlg,long double *lwg){
  vallistelt *tmp=new vallistelt;
  tmp->next=start;
  start=tmp;
  start->ref=++lastref;
  for(unsigned int i=0;i<NUMTERMS;i++){
    start->akd[i]=ak[i];
  };
  start->numlowgvals=nlg;
  start->lowgvalues=lwg;
  start->lowvals=lwEC;
  return lastref;
};


void vallist::insert(long double ak[],long double *lwEC ,unsigned int nlg,long double *lwg,int refno){
  if(start&&start->checkvalue(refno)){
    throw("Reference already defined.\n");
  };
  vallistelt *tmp=new vallistelt;
  tmp->next=start;
  start=tmp;
  start->ref=refno;
  if(refno>lastref){
    lastref=refno;
  };
  for(unsigned int i=0;i<NUMTERMS;i++){
    start->akd[i]=ak[i];
  };
  start->numlowgvals=nlg;
  start->lowgvalues=lwg;
  start->lowvals=lwEC;
};


void vallist::lookupvalues(int ref){
  if(ref!=currentref){
    //if values currently loaded then nothing to be done.
    vallistelt *current=start;
    for(;current&&current->ref!=ref;current=current->next);
    if(current){
      for(unsigned int i=0;i<NUMTERMS;i++){
	ankdiag[i]=current->akd[i];
      };
      lowgvalues=current->lowgvalues;
      lowvals=current->lowvals;
      numlowgvals=current->numlowgvals;
      currentref=ref;
    }else{
      throw "Undefined Reference";
    };
  };
};

void setuplowgvalues(long double a,unsigned int N,unsigned int nlg=7){
  numlowgvals=nlg;
  lowgvalues=new long double[numlowgvals];
  long double start=log(a);
  for(unsigned int i=N;i>0;i--){
    start-=1/((long double) i);
  };
  for(unsigned int i=0;i<numlowgvals;i++){
    long double term=-(N*i/a);
    *(lowgvalues+i)=start-term;
    for(unsigned int j=1;term!=0&&j<N;j++){
      term*=(-1)*((int)i-(int)j)*((N-j)/a)*j/((long double)(j+1)*(j+1));
      *(lowgvalues+i)-=term;
    };
  };
};

inline long double getlogg(unsigned int x){
    if(x<numlowgvals){
      if(lowgvalues==NULL){
	throw("g values table not set up.");
      };
      return *(lowgvalues+x);
    }else{
      return log(x);      
    };
};



void setuplowECVar(long double a,unsigned int N){
  lowvals=new long double[numlowgvals];
  long double la=log(a);
  long double lasq=la*la;
  long double sm1n=0;
  for(unsigned int i=N;i>0;i--){
    sm1n+=1/((long double )i);
  };
  long double *Nchkk=new long double[N];//N choose k divided by k
  *Nchkk=N;
  for(unsigned int i=1;i<N;i++){
    *(Nchkk+i)=((*(Nchkk+i-1)*(N-i))*i)/((i+1)*(i+1));
  };
  long double *Nchkl=new long double[2*N-1];//sum_{k+l=m}Nchkk*Nchll
  for(unsigned int i=0;i<2*N-1;i++){
    *(Nchkl+i)=0;
  };
  for(unsigned int i=0;i<N;i++){
    for(unsigned int j=0;j<N;j++){
      *(Nchkl+i+j)+=*(Nchkk+i)*(*(Nchkk+j));
    };
  };

  for(unsigned int i=0;i<numlowgvals;i++){    
    lowvals[i]=lasq;
    long double sum=sm1n;
    long double apow=-1/a;
    unsigned int cumprod=i;
    for(unsigned int j=1;j<=i&&j<=N;j++){
      sum+=apow*(*(Nchkk+j-1))*cumprod;
      cumprod*=i-j;
      apow/=-a;
    };
    lowvals[i]-=2*(la-sm1n)*sum;
    lowvals[i]-=sm1n*sm1n;
    apow=1/(a*a);
    cumprod=i*(i-1);
    for(unsigned int j=0;j<2*N-1;j++){
      lowvals[i]+=apow*(*(Nchkl+j))*cumprod;
      cumprod*=(i-j-2);
      apow/=-a;
    };
    long double gx=getlogg(i);
    lowvals[i]=gx*gx-lowvals[i];
  };

  for(unsigned int i=numlowgvals;i<numlowgvals;i++){
    lowvals[i]=0;
    long double cpr=1;
    for(unsigned int j=0;j<NUMTERMS;j++){
      cpr/=(i+j+1);
      lowvals[i]+=cpr*ankdiag[j];
    };
  };
};



inline long double getaaval(unsigned int n){
  if(n<=numlowgvals){
    //lookup value in table
    return *(lowvals+n);
  };
  long double lnt=log(((long double)n)/2);
  return lnt*lnt-(lnt+1)/n-(n+1)*(lnt/2+((long double)5)/12)/(2*((long double)n)*n*n);
};




extern "C"{

  void log_g(double *answers,int *n,const int *xpointer,const int *refpointer){
    vallist::lookupvalues(*refpointer);
    for(unsigned int i=0;i<(unsigned) (*n);i++){
      if((unsigned int)(*(xpointer+i))<numlowgvals){
	*(answers+i)=*(lowgvalues+*(xpointer+i));
      }else{
	*(answers+i)=log(*(xpointer+i));      
      };
    }
  };

  
  void log_ECVar(double *answers,const int *npointer,const double *apointer,const int *Npointer,const int *x){
    unsigned int n=*npointer;
    if(lowgvalues==NULL){
      setuplowgvalues(*apointer,*Npointer);
    };      
    long double totgvalsq=0;
    long double totmsq=0;
    for(unsigned int i=0;i<n;i++){            
      totgvalsq+=getlogg(*(x+i))*(getlogg(*(x+i)));
      long double numterms=0;
      long double msqi=0;
      for(unsigned int j=0;j<n;j++){
	long double logweight=0;
	int sum=*(x+i)+(*(x+j));
	if(*(x+i)&&(*(x+j))){	  
	  logweight=sum*log(sum)-*(x+i)*log(*(x+i)*2)-*(x+j)*log(*(x+j)*2);
	}else{
	  logweight=-sum*log(2);
	};
	long double weight=exp(logweight);
	numterms+=weight;
	msqi+=weight*getaaval(sum);
      };
      msqi/=numterms;
      totmsq+=msqi;
    };
    *answers=(totgvalsq-totmsq)/n;
  };


  
  void precalculate(int *ref, const double *a,const int *N,const int *nlg){
    //first calculate low values of the estimator from first principles.

    setuplowgvalues(*a,*N,*nlg);  
    
  //For larger values use g(X)=log(X) and Taylor series expansion.
    //stirling[i][j] is the number of ways to divide i+1 things into
    //j+1 partitions all of size at least 2.
    long unsigned int stirling[NUMTERMS+1][(NUMTERMS+1)/2];
    for(unsigned int i=1;i<NUMTERMS+1;i++){
      stirling[i][0]=1;
    };
    for(unsigned int i=0;i<(NUMTERMS+1)/2;i++){
      stirling[0][i]=0;
    };
    for(unsigned int i=1;i<(NUMTERMS+1);i++){
      for(unsigned int j=1;j<(NUMTERMS+1)/2;j++){
	stirling[i][j]=stirling[i-1][j]*(j+1);
	if(i>1){
	  stirling[i][j]+=i*stirling[i-2][j-1];
	};
      };
    };
    long double stirlingscaled[NUMTERMS+1][(NUMTERMS+1)/2];
    for(unsigned int i=0;i<(NUMTERMS+1);i++){
      for(unsigned int j=0;j<(NUMTERMS+1)/2;j++){
	stirlingscaled[i][j]=((long double)stirling[i][j])/(i+1);
      };
    };
    
    long double coeff[NUMTERMS+1][(NUMTERMS+1)/2];

    for(unsigned int i=0;i<(NUMTERMS+1);i++){
      for(unsigned int j=0;j<(NUMTERMS+1)/2;j++){
	coeff[i][j]=0;
	for(unsigned int k=0;k<i;k++){
	  for(unsigned int l=0;l<j;l++){
	    coeff[i][j]+=stirlingscaled[k][l]*stirlingscaled[i-1-k][j-1-l];
	  };
	};
      };
    };
    
    long double  scale[NUMTERMS+1];
    for(unsigned int i=0;i<NUMTERMS;i++){
      ankdiag[i]=0;
      if(i){
	scale[i]=(scale[i-1]*(i+1)+2/((long double)i+1))/(i+2);
      }else{
	scale[i]=1;//1;//0??
      };
      int n1p=1;
      for(unsigned int j=0;j<(NUMTERMS+1)/2&&j<i;j++){
	ankdiag[i]+=n1p*(scale[i-j]*stirling[i-j][j]-coeff[i-j][j]);
	n1p*=-1;
      };
    }; 
    for(unsigned int i=0;i<NUMTERMS-1;i++){
      ankdiag[i]=ankdiag[i+1];
    };
    scale[NUMTERMS]=(scale[NUMTERMS-1]*(NUMTERMS+1)+2/((long double)NUMTERMS+1))/(NUMTERMS+2);
    ankdiag[NUMTERMS-1]=0;
    int n1p=1;
    for(unsigned int j=0;j<(NUMTERMS+1)/2;j++){
      ankdiag[NUMTERMS-1]+=n1p*(scale[NUMTERMS-j]*stirling[NUMTERMS-j][j]-coeff[NUMTERMS-j][j]);
            n1p*=-1;
    };

    setuplowECVar(*a,*N);
    if(*ref==-1){
      *ref=vallist::insert(ankdiag,lowvals , numlowgvals,lowgvalues);
    }else{
      vallist::insert(ankdiag,lowvals , numlowgvals,lowgvalues,*ref);
    };
  };
  
  void log_ECVar_unbiassed(double *answers,const int *npointer,const int *ref,const int *x, const double *a,const int *N,const int *nlg){
    unsigned int n=*npointer;
    try{
      vallist::lookupvalues(*ref);
    }catch(...){
      //Reference not found, probably deleted on closing R
      int refno=*ref;
      precalculate(&refno,a,N,nlg);
    };
    *answers=0;
    for(unsigned int i=0;i<n;i++){
      if((unsigned int)(*(x+i))<numlowgvals){
	*answers+=lowvals[*(x+i)];
      }else{
	long double cpr=1;
	for(unsigned int j=0;j<NUMTERMS;j++){
	  cpr/=(*(x+i)+j+1);
	  *answers+=cpr*ankdiag[j];
	};
      };
    };
    *answers/=n;
  };

};




