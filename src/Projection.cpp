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

//#include <iostream>
//using namespace std;


#include <R.h>

void solve(unsigned int p,double *answer,const double *M,const double *lambda,double midd,const double *x){
  //solves (Lambda+M)a=x
  //p is the dimension
  //answer is the vector into which the answer is written
  //lambda is a diagonal matrix
  //M is a p*p matrix (column-major format)
  //x is a length p vector

  //Calculates an approximate solution based on lambda+midrange of d.

  double *temp=new double[p];
  for(unsigned int i=0;i<p;i++){
    *(answer+i)=*(temp+i)=*(x+i)/(midd+*(lambda+i));
  };
  delete[] temp;
//  for(unsigned int i=0;i<p;i++){
//    for(unsigned int j=0;j<p;j++){
//      *(answer+i)-=((*(M+i*p+j)-midd)/(midd+*(lambda+j)))*(*(temp+j));
//    };
//  };
};

extern "C"{

  void get_scores_log(double *answers,int *p,const int *x,const double *V, const double *d,const int *k,const double *mu,const double *M){
    //p is the dimension.
    //xpointer is the data - a vector of length p
    //V is the p*p matrix of eigenvectors of the estimated variance
    //d is the vector of eigenvalues of the estimated variance
    //k is the number of components
    //mu is the esimated mean of lambda
    //M=VDV^T

    //    cout<<"p="<<*p<<"\n";

    long double movesize=1;

    long double mind=*d;
    long double maxd=*d;
    for(unsigned int i=0;i<(unsigned int)(*p);i++){
      if(*(d+i)>maxd){
	maxd=*(d+i);
      }else if(*(d+i)<mind){
	mind=*(d+i);
      };
    };
    long double midd=(mind+maxd)/2;

    //    cout<<"midd="<<midd<<"p="<<*p<<"\n";
    
    double *lambda=new double[(*p)];
    double *llmm=new  double[(*p)];
    double *llmmchange=new  double[(*p)];
    double *s=new  double[(*p)];
    double *t=new  double[(*p)];
    double *step=new double[(*p)];
    for(unsigned int i=0;i<(unsigned int)(*p);i++){
      *(lambda+i)=*(x+i);
      if(*(x+i)==0){
	*(lambda+i)=0.05;
      };
    };
    
    
#define MAXITER 2000
#define TOLERANCE 1e-16
#define NUMEIGENAPPROX 1

    unsigned int calculatepropers= 1;

    //Initial value lambda=x
    
    
    for(unsigned int i=0;i<(unsigned int)(*p);i++){
      *(s+i)=0;
    };

    for(unsigned int numiter=0;movesize>TOLERANCE&&numiter<MAXITER;numiter++){
      for(unsigned int i=0;i<(unsigned int)(*p);i++){
	//		cout<<*(lambda+i)<<"  ";
	long double temp=log(*(lambda+i))-*(mu+i);
	if(numiter){
	  *(llmmchange+i)=temp-*(llmm+i);
	}else{
	  *(llmmchange+i)=temp;
	};
	*(llmm+i)=temp;
      };

      //            cout<<"\n";
      //llmm = log(lambda)-mu
      //llmmchange= log(lambda)-log(lambda_old)
      if(numiter%calculatepropers){
	//Use only first few eigenvectors of V.
	for(unsigned int i=0;i<NUMEIGENAPPROX;i++){
	  double vvalue=0;
	  for(unsigned int j=0;j<(unsigned int)(*p);j++){
	    vvalue+=*(llmmchange+j)*(*(V+i*(*p)+j));
	  };
	  vvalue*=*(d+i);
	  for(unsigned int j=0;j<(unsigned int)(*p);j++){
	    *(s+j)+=*(V+i*(*p)+j)*vvalue;
	  };
	};
      }else{
	//	cout<<"Calculating full Hessian matrix.\n";
	for(unsigned int j=0;j<(unsigned int)(*p);j++){
	  *(s+j)=0;
	};	
	//recalculate s fully every 10 iterations
	for(unsigned int i=0;i<(unsigned int)(*p);i++){      
	  for(unsigned int j=0;j<(unsigned int)(*p);j++){
	    //*(s+j)+=*(M+i*(*p)+j)*(*(llmmchange+i));
	    *(s+j)+=*(M+i*(*p)+j)*(*(llmm+i));
	  };
	};
      };
      //s=M llmm
      //    cout<<"t=";
      for(unsigned int i=0;i<(unsigned int)(*p);i++){
	*(t+i)=*(x+i)-*(lambda+i)-*(s+i);
	//		cout<<*(t+i)<<" ";
      };
      //   cout<<"\n";
      
      solve((*p),step,M,lambda,midd,t);

      //Now a_new=a + V^T step
      //log(lambda_new)=log(lambda)+step
      //lambda_new=lambda exp(step)
      
      movesize=0;
      for(unsigned int i=0;i<(unsigned int)(*p);i++){
	//	*(step+i)-=1;
	if((1+*(step+i))<0.5){
	  //too close to zero. optimise with other values fixed.
	  double c=-*(x+i);
	  for(unsigned int j=0;j<(unsigned int)(*p);j++){
	    if(j!=i){
	      if(*(llmm+j)<-2){
		c-=*(M+i*(*p)+j)*2;		
	      }else if(*(llmm+j)>2){
		c+=*(M+i*(*p)+j)*2;		
	      }else{
		c+=*(M+i*(*p)+j)*(*(llmm+j));
	      };
	    };
	  };
	  //solve M log(1+step)+c+lambda_i(1+step)=0;
	  double s=0.01;
	  //This loop should converge, but in case it doesn't, we add a maximum number of iterations.
	  unsigned int onevariableiterations=0;
	  for(double val=1;val*val>1e-10&&onevariableiterations<10000;){
	    onevariableiterations++;
	    val=log(s)*(*(M+i*(*p)+i))-c+*(lambda+i)*s;
	    double dval=*(M+i*(*p)+i)/s+*(lambda+i);
	    if(val/(dval*s)>0.5){
	      s/=2;
	    }else{
	      s-=val/dval;
	    };
	  };
	  *(step+i)=s-1;
	  //	  	  cout<<"i="<<i<<"s="<<s<<"\n";
	  	  if(*(step+i)*(*(step+i))>0.99){
		    //	    cout<<numiter<<"  "<<i<<"  "<<*(step+i)<<"  "<<*(M+i*(*p)+i)<<"  "<<c<<"  "<<*(lambda+i)<<"\n";
	  	  };
	};
	if(*(step+i)<-0.99999){
	  //	  cout<<numiter<<"  "<<i<<"  "<<*(step+i)<<"\n";
	  *(step+i)=-0.999;//limit stepsize
	};
	*(lambda+i)+=*(lambda+i)*(*(step+i));
	movesize+=*(step+i)*(*(step+i));
      };
      //      cout<<"movesize="<<movesize<<"\n";
      if((numiter==MAXITER-2)&&calculatepropers>1){
	calculatepropers=1;
      };
            *(answers+(*p)+(*k))=numiter;//convergence information.
	    *(answers+(*k)+(*p)+1)=movesize;
      if(movesize<TOLERANCE||numiter==MAXITER-1){
	if((numiter<MAXITER-1)&&(numiter%calculatepropers)){
	  calculatepropers=1;
	  movesize=TOLERANCE+1;
	}else{
	  *(answers+(*p)+(*k))=numiter;//convergence information.
	  *(answers+(*k)+(*p)+1)=movesize;
	};
      };
    };

    long double sumt=0;
    for(unsigned int i=0;i<(unsigned int)(*p);i++){
      sumt+=*(t+i)*(*(t+i));
    };
    for(unsigned int i=0;i<(unsigned int)(*k);i++){
      *(answers+i)=0;
      for(unsigned int j=0;j<(unsigned int)(*p);j++){
	*(answers+i)+=*(V+i*(*p)+j)*(*(llmm+j));
      };
    };
    
    for(unsigned int i=0;i<(unsigned int)(*p);i++){
      *(answers+i+(*k))=*(lambda+i);
    };
    
    delete[] lambda;
    delete[] llmm;
    delete[] llmmchange;
    delete[] s;
    delete[] t;
    delete[] step;
  };
 
};
