#include <iostream>
#include <cmath>
#include<iomanip>
#include <fstream>

using namespace std;

void rhs(double t, double *X, double *R);
double func_for_Euler(double *Y, double x, double h, double t);
double func_for_Euler_old(double *Y, double x, double h, double t);
void I_Euler_kapitza(double t, double *Y, double h);
void I_Euler_kapitza_old(double t, double *Y, double h);
void RK2Step(double t, double *Y, void (*RHS_Func)(double, double *, double *), double h, int neq);
void RK4Step(double t, double *Y, void (*RHS_Func)
    (double, double *, double *), double h, int neq);
double Bisection(double (*F)(double *Y, double x, double h, double t),
    double *Y,double a, double b, double tol, double h, double t);
double g=10, l=10;

int main(){
  double dt=0.1,t=0.0;
  double y[2]={M_PI*0.25, 0.0};     
  double E=0.0;
  //g,m,L=1
  int N_periods=500;
  ofstream fdata;
  fdata.open("pendulum.dat");
  double vaux, cycle=0.0;
  for(int i=0;cycle<2*N_periods;i++){
    vaux=y[1];
    I_Euler_kapitza(t,y,dt);
    t+=dt;
    fdata << t << "   " <<y[0] <<" " <<y[1] <<" " <<0.5*y[1]*y[1]+1-cos(y[0])  << endl;
    if(vaux*y[1]<=0){
      cycle++;
    }
  }
  cycle=0;
  t=0.0;
  fdata <<endl <<endl;

  cycle=0;
  t=0;
  y[0]=M_PI*0.25;
  y[1]=0.0;
  for(int i=0;cycle<2*N_periods;i++){
    vaux=y[1];
    I_Euler_kapitza_old(t,y,dt);
    t+=dt;
    fdata << t << "   " <<y[0] <<" " <<y[1] <<" " <<0.5*y[1]*y[1]+1-cos(y[0])  << endl;
    if(vaux*y[1]<=0){
      cycle++;
    }
  }
  cycle=0;
  t=0.0;
  fdata <<endl <<endl;

  y[0]=M_PI*0.25;
  y[1]=0.0;
  for(int i=0;cycle<2*N_periods;i++){
    vaux=y[1];
    RK2Step(t, y, rhs, dt, 2);
    t+=dt;
    fdata <<t << "    " << y[0] <<" " <<y[1] <<" " <<0.5*y[1]*y[1]+1-cos(y[0]) <<endl;
    if(vaux*y[1]<=0.0){
      cycle++;
    }
  }
  fdata <<endl <<endl;

  t=0.0;
  cycle=0;
  y[0]=M_PI*0.25;
  y[1]=0.0;
  for(int i=0;cycle<2*N_periods;i++){
    vaux=y[1];
    RK4Step(t, y, rhs, dt, 2);
    t+=dt;
    fdata <<t << "    " << y[0] <<" " <<y[1] <<" " <<0.5*y[1]*y[1]+1-cos(y[0]) <<endl;
    if(vaux*y[1]<=0.0){
      cycle++;
    }
  }
  fdata <<endl <<endl;

  fdata.close();
  return 0;
}


double func_for_Euler(double *Y, double x, double h, double t){
  return x-Y[0]-h*Y[1]+h*h/4*(g/l)*(sin(x)+sin(Y[0]));
}

double func_for_Euler_old(double *Y, double x, double h, double t){
  return x-Y[0]-h*Y[1]+h*h*(g/l)*(sin(x));
}

//Had to rewrite bisection in order to pass more arguments
void I_Euler_kapitza(double t, double *Y, double h){
  double Y1[2];
  Y1[0]=Bisection(func_for_Euler,Y, Y[0]-1,Y[0]+1,1e-7 ,h,t);
  Y1[1]=Y[1]-h/2*(g/l)*(sin(Y1[0])+sin(Y[0]));
  Y[0]=Y1[0];
  Y[1]=Y1[1];
}

void I_Euler_kapitza_old(double t, double *Y, double h){
  double Y1[2];
  Y1[0]=Bisection(func_for_Euler_old,Y, Y[0]-1,Y[0]+1,1e-7 ,h,t);
  Y1[1]=Y[1]-h*(g/l)*sin(Y1[0]);
  Y[0]=Y1[0];
  Y[1]=Y1[1];
}


void rhs(double t, double *X, double *R){
  double mu=0.0;
  R[0]=X[1];
  R[1]=-g/l*sin(X[0]);
}

void RK2Step(double t, double *Y, void (*RHS_Func)(double, double *, double *), double h, int neq){
    double Y1[neq], k1[neq], k2[neq];
    int i;   
    RHS_Func(t,Y,k1);
    for(i=0; i<neq; i++){
        Y1[i] = Y[i]+0.5*h*k1[i];
    }    
    RHS_Func(t+0.5*h,Y1,k2);    
    for(i=0; i<neq; i++){
        Y[i] += h*k2[i];
    }
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



/* plot "pendulum.dat" using 2:3 index 1 title "RK4"
gnuplot> replot "pendulum.dat" using 2:3 index 0 title "Position Verlet"
gnuplot> set label 1 "dt=0.1, Nperiodi=500, g/L=1" at 1.2,1.2
gnuplot> replot
gnuplot> unset label 1
gnuplot> replot
gnuplot> set label 1 "dt=0.1, Nperiodi=500, g/L=1" at 1,1.2
gnuplot> replot
gnuplot> set label 2 "theta0=0.25*pi, dtheta/dt0=0" at 1,1.1
gnuplot> replot
gnuplot> set label 3 "Come atteso, Verlet preserva l'energia " at 1,0.5
gnuplot> replot
gnuplot> unset label 3
gnuplot> set label 3 "Come atteso, Verlet preserva l'energia " at 1,-0.5
gnuplot> replot
gnuplot> set label 4 "e questo comporta che la curva resti " at 1,-0.6
gnuplot> set label 5 "chiusa nello spazio delle fasi " at 1,-0.7
gnuplot> replot
gnuplot> set title "Caretti Pendulum"
gnuplot> set xlabel "theta"
gnuplot> set ylabel "dtheta/dt"
gnuplot> replot*/
