/************************************************************************************************/
/*** Topic: Hodgkin–Huxley model               Ali-Seif                                       ***/
/*** Version Release 17.12 rev 11256                                                          ***/
/*** Date: 1397/3/1                                                                           ***/
/*** Code implemented in CodeBlocks C++ compiler (v. 17.12),                                  ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/************************************************************************************************/
#include <iostream>
#include <math.h>

using namespace std;


int main()
{
//__________Define the parameters____________/

    int   step=         100;
    float h=           0.01;
    float Iapp =        0.75;

//=============First neurons===========

    float v1[step] =    { 0.0 };
    float h1[step] =    { 0.0 };
    float n1[step] =    { 0.0 };

    float INa1=         0.0;
    float IK1=          0.0;
    float Il1=          0.0;

//_____________Fixed values___________________/

    float C_m =         1.0;        //microFarad

    int   E_Na=         55;         //miliVolt
    int   E_K=          -90;        //miliVolt
    int   E_l=          -65;        //miliVolt
    int   g_Na=         35;         //mScm-2
    int   g_K=          9;          //mScm-2
    float g_l=          0.1;        //mScm-2

    int   phi=          5;

//______________Initial values________________/

          v1[0]=        -63.0;

//=============First neurons===========

    float alpha_m_1=    -0.1*(v1[0]+35)/(exp(-0.1*(v1[0]+35))-1);
    float alpha_h_1=    0.07*exp(-(v1[0]+58)/20);
    float alpha_n_1=    -0.01*(v1[0]+34)/(exp(-0.1*(v1[0]+34))-1);

    float beta_m_1=     4*exp(-(v1[0]+60)/18);
    float beta_h_1=     1/(exp(-0.1*(v1[0]+28))+1);
    float beta_n_1=     0.125*exp(-(v1[0]+44)/80);


    float m1 =          alpha_m_1/(alpha_m_1+beta_m_1);
          n1[0] =       alpha_n_1/(alpha_n_1+beta_n_1);
          h1[0] =       alpha_h_1/(alpha_h_1+beta_h_1);


    //cout <<'\t'<<m1<<'\t'<<n1[0]<<'\t'<<h1[0]<<endl;


//_______________Start changes________________/

    for (int i=1 ; i<=step ;i++){




        n1[i] =(n1[i-1] + h*phi*((alpha_n_1*(1-n1[i-1]))-beta_n_1*n1[i-1]));
        h1[i] =(h1[i-1] + h*phi*((alpha_h_1*(1-h1[i-1]))-beta_h_1*h1[i-1]));

        cout <<i<<'\t'<<m1<<'\t'<<n1[i]<<'\t'<<h1[i]<<endl;



        INa1 = g_Na*h1[i-1]*pow(m1,3)*(v1[i-1]-E_Na);
        IK1 = g_K*pow(n1[i-1],4)*(v1[i-1]-E_K);
        Il1=g_l*(v1[i-1]-E_l);



        v1[i]=(v1[i-1]+(h)*((1/C_m)*(Iapp-(INa1+IK1+Il1))));


        alpha_m_1=    -0.1*(v1[i]+35)/(exp(-0.1*(v1[i]+35))-1);
        alpha_h_1=    0.07*exp(-(v1[i]+58)/20);
        alpha_n_1=    -0.01*(v1[i]+34)/(exp(-0.1*(v1[i]+34))-1);

        beta_m_1=     4*exp(-(v1[i]+60)/18);
        beta_h_1=     1/(exp(-0.1*(v1[i]+28))+1);
        beta_n_1=     0.125*exp(-(v1[i]+44)/80);

        m1 =          alpha_m_1/(alpha_m_1+beta_m_1);



        //cout<< i <<'\t' << v1[i]<<endl;
    }


/*
    //cout<<v[0];
    cout<<"DON ;D";

    ofstream temp("temp.txt", ios::out | ios::trunc);              //یه فایل بساز
    for(int i=0; i<=step; i++)
    temp << i << '\t' <<v1[i]<< endl;     //چاپ کن هر کدومو داخل فایل
    temp.close();
    cout << "\nFinish" << endl;

*/
    return 0;
}
