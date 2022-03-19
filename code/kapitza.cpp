#include <iostream>
#include <cmath>
#include<iomanip>
#include <fstream>


#define STAGE 1
#define METHOD 2

using namespace std;
void rhs(double t, double *X, double *R);
void RK4Step(double t, double *Y, void (*RHS_Func)(double, double *, double *), double h, int neq);
void I_trap_kapitza(double t, double *Y, double h);
void forward_time(double t, double *y, void (*RHS_Func)(double, double *, double *), 
  double dt, int neq, int n_periods, bool verbose);
double Bisection(double (*F)(double *Y, double x, double h, double t),
  double *Y,double a, double b, double tol, double h, double t);
double l=10, b=0.5, g=1, omega=2*sqrt(2*100)/b;


int main(){
  double t=0.0;
  double dt;

  double y[2]={-M_PI*0.80, 0.0}; 
  int N_periods;
  #if STAGE==1     
  N_periods=300;
  dt=0.01;
  bool verbose=true;
  #endif
  #if STAGE==2
  N_periods=3000;
  dt=0.001;
  bool verbose=true;
  #endif
  #if STAGE==1 || STAGE==2
  forward_time(t,y,rhs,dt,2,N_periods, verbose);
  #endif
  #if STAGE==3
  N_periods=3000;
  bool verbose=false;
  g=10;
  l=10;
  for(b=1;b>0.15;b=b/sqrt(2)){
    for(dt=0.0001;dt<1;dt=10*dt){
      t=0.0;
      forward_time(t,y,rhs,dt,2,N_periods,verbose);
    }
  }
  #endif
  return 0;
}



void forward_time(double t, double *y, void (*RHS_Func)(double, double *, double *),
    double dt, int neq, int n_periods, bool verbose){
  #if STAGE==1
  ofstream fdata;
  fdata.open("kapitza.dat");
  int cycle=0;
  double vaux=0.0;
  for(int i=0;(cycle<2*n_periods);i++){
    vaux=y[1];
    #if METHOD==1
    RK4Step(t, y, rhs, dt, 2);
    #endif
    #if METHOD==2
    I_trap_kapitza(t,y,dt);
    #endif
    t+=dt;
    if (verbose==true){
      cout <<t << "    " <<cycle <<endl;
    }
    fdata <<t << "    " << 1 << "    "<<0.0 <<" " <<-b*cos(omega*t) <<endl;
    fdata <<t << "    " << 2 << "    "<< l*sin(y[0]) <<" " 
      <<-l*cos(y[0])-b*cos(omega*t)  <<endl;
    if(vaux*y[1]<=0.0){
      cycle++;
    }
  }
  fdata <<endl <<endl;
  fdata.close();
  #endif
  #if STAGE==2 || STAGE==3
  double tol=0.001, vaux, cycle;
  double low_omega=1.;
  double upper_omega=81;
  bool stable;
  double a,d;
  while(fabs(upper_omega-low_omega)>tol){
    stable=false;
    omega=0.5*(low_omega+upper_omega);
    t=0.0;
    y[0]=M_PI*(-0.96);
    y[1]=0.; 
    cycle=0.0;
    for(int i=0;(cycle<2*n_periods)&&(t<300);i++){
      vaux=y[1];
      #if METHOD==1
      RK4Step(t, y, rhs, dt, 2);
      #endif
      #if METHOD==2
      I_trap_kapitza(t,y,dt);
      #endif
      t+=dt;
      //cout <<t << "    " <<cycle <<endl;
      if(vaux*y[1]<=0.0){
        cycle++;
      }
    }
    if (fabs(y[0]+M_PI)<0.041*M_PI){
      stable=true;
    }
    if(stable==true){
      upper_omega=omega;
      if(verbose==true){
        cout <<omega <<": Stable"  <<endl;
      }
    }
    if(stable==false){
      low_omega=omega;
      if(verbose==true){
        cout <<omega <<": Unstable" <<endl;
      }
    }
    if(verbose==true){
      cout <<"Lower:" <<low_omega <<" Upper:" <<upper_omega <<endl;
    }
  }
  a=0.5*(low_omega+upper_omega);
  d=sqrt(2*g*l)/b;
  cout <<"Critical Omega" <<a  <<"   Expected Omega:" <<d 
    <<"   Difference: " <<(a-d) <<"   dt= " <<dt <<endl;
  #endif
}


void rhs(double t, double *X, double *R){
  double mu=0.0;
  R[0]=X[1];
  R[1]=-g/l*sin(X[0])-b/l*omega*omega*cos(omega*t)*sin(X[0])-mu*X[1];
}

void RK4Step(double t, double *Y, void (*RHS_Func)
    (double, double *, double *), double h, int neq){
  double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];
  int i;    
  RHS_Func(t,Y,k1);
  for(i=0; i<neq; i++){
      Y1[i] = Y[i]+0.5*h*k1[i];
  }    
  RHS_Func(t+0.5*h,Y1,k2);
  for(i=0; i<neq; i++){
      Y1[i] = Y[i]+0.5*h*k2[i];
  }    
  RHS_Func(t+0.5*h,Y1,k3);
  for(i=0; i<neq; i++){
      Y1[i] = Y[i]+h*k3[i];
  }    
  RHS_Func(t+h,Y1,k4);
  for(i=0; i<neq; i++){
      Y[i] += h/6.0*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
  }
}


double func_for_trap(double *Y, double x, double h, double t){
  return x-Y[0]-h*Y[1]+h*h/4*((b/l*omega*omega*cos(omega*(t+h))+g/l)*sin(x)
    +(b/l*omega*omega*cos(omega*(t))+g/l)*sin(Y[0]));
}

//Had to rewrite bisection in order to pass more arguments
void I_trap_kapitza(double t, double *Y, double h){
  double Y1[2];
  Y1[0]=Bisection(func_for_trap,Y, Y[0]-1,Y[0]+1,1e-7 ,h,t);
  Y1[1]=Y[1]-h/2*((b/l*omega*omega*cos(omega*t)+g/l)*sin(Y[0])+
    (b/l*omega*omega*cos(omega*(t+h))+g/l)*sin(Y1[0]));
  Y[0]=Y1[0];
  Y[1]=Y1[1];
}

double Bisection(double (*F)(double *Y, double x, double h, double t),
    double *Y,double a, double b, double tol, double h, double t){
  double fa=F(Y,a,h,t);
  double fb=F(Y,b,h,t);
  double xm;
  double fm;
  double dx=fabs(b-a);
  for(int k=1;fabs(a-b)>tol; k++){
    xm=(a+b)*0.5;
    fm=F(Y,xm,h,t);
    //cout <<"Bisection:  k=" <<k <<"; [a,b]= [" <<a <<", " <<b <<"]" <<"; xm= "
      //<<xm <<endl;
    if(fm*fa>0){
      a=xm;
      fa=fm;
    }
    else if(fm*fa<0){
      b=xm;
      fb=fm;
    }
    else if(fm*fa==0){
      return xm;
    }   
  }
  return xm;
}

