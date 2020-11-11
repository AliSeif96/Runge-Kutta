/************************************************************************************************/
/*** Topic: Runge-Kutta 4th Order Method to Solve Differential Equation          Ali-Seif     ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 11/9/2020                                                                          ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/

#include <iostream>
using namespace std;

//________________________Differential Equation_______________________//

double diffOfy(double t, double x) {
   return   (x-(t*t)+1);                                              //function x'=x-t^2+1
}
//_______________________Runge-Kutta calculation______________________//

double rk4thOrder(double t0, double x0, double t, double h) {

    double  k1, k2, k3, k4;
    double  m=      (t - t0)/h;                                       //calculate number of iterations
    int     n=      m;
    double  x=      x0;                                               //initially x is f(x0)

    for(int  i=1; i<=n; i++) {
            k1=     h*diffOfy(t0, x);
            k2=     h*diffOfy((t0+h/2), (x+k1/2));
            k3=     h*diffOfy((t0+h/2), (x+k2/2));
            k4=     h*diffOfy((t0+h), (x+k3));

            x+=     double((1.0/6.0)*(k1+2*k2+2*k3+k4));              //update x using del x
            t0+=    h;                                                //update t0 by h
   }
   return   x;                                                        //f(x) value
}
//____________________The principle of the program____________________//

int main() {

    double  t0, x0, t, h;
    cout << "Enter t0(start point): "; cin >> t0;                     //for example t0=0
    cout << "Enter t (The final amount required): "; cin >> t;        //for example t=2
    cout << "Enter x0(Starting point value): "; cin  >> x0;           //for example t=0.5
    cout << "Enter h (Time step): "; cin >> h;                        //for example h=0.02
    cout << "Answer of differential equation: " << rk4thOrder(t0, x0, t, h);
    return 0;
}
