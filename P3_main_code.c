/*
RK2 method for plotting Phase Diagrams of Systems of Differential Equations.
Author: Tomas Coronado Gonzalez
Date: Augost 2020

IVP:    {  y1'(t) = -2*y1(t)*y2(t)
        {  y2'(t) = y2(t)*y2(t)+y1(t)-1
*/

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define T 1e6

double f1 (double w1_fun, double w2_fun, double t_fun_w);        // y1'(t) = -2*y1(t)*y2(t)       --> dy1/dt = f1(y1,y2,t) = -2*y1(t)*y2(t)       (ODE expressed in its normal form)
double f2 (double w1_fun, double w2_fun, double t_fun_w);        // y2'(t) = y2(t)*y2(t)+y1(t)-1  --> dy2/dt = f2(y1,y2,t) = y2(t)*y2(t)+y1(t)-1  (ODE expressed in its normal form)                                                 

int main()
{
    FILE *pf;
    
    int n=100;                                                   // Number of steps of the interval
    int j,k;

    int M1=0,M2=0;                                               // M1 is the number of curves plotted on the Phase Diagram, or the number of initial points
    
    int p,q,m;
    double s1=6,s1_bis=4,r1=0.5,r1_bis=0.25;
    double s2=0.5,r2=0.1;
    double s3=0.5,r3=0.1;
    double s4=0.5,r4=0.1;

    // Initial points around all the domain:
    double m1=-s1,m2=-s1;
    for (p=0;m1<=s1;p++)   
    {
        for (q=0;m2<=s1;q++)   
        {
            M1++;     
            m2=m2+r1;   
        } 
        m1=m1+r1;
        m2=-s1;
    }

    m1=-s1_bis,m2=-s1_bis;
    for (p=0;m1<=s1_bis;p++)   
    {
        for (q=0;m2<=s1_bis;q++)   
        {
            M1++;     
            m2=m2+r1_bis;   
        } 
        m1=m1+r1_bis;
        m2=-s1_bis;
    }

    // Initial points around the point (1,0)
    m1=1-s2,m2=-s2;
    for (p=0;m1<=1+s2;p++)   
    {
        for (q=0;m2<=s2;q++)   
        {
            M1++;     
            m2=m2+r2;   
        } 
        m1=m1+r2;
        m2=-s2;
    }

    // Initial points around the point (0,1)
    m1=-s3,m2=1-s3;
    for (p=0;m1<=s3;p++)   
    {
        for (q=0;m2<=1+s3;q++)   
        {
            M1++;     
            m2=m2+r3;   
        } 
        m1=m1+r3;
        m2=1-s3;
    }

    // Initial points around the point (0,-1)
    m1=-s4,m2=-1-s4;
    for (p=0;m1<=s4;p++)   
    {
        for (q=0;m2<=-1+s4;q++)   
        {
            M1++;     
            m2=m2+r4;   
        } 
        m1=m1+r4;
        m2=-1-s4;
    }

    printf("%d\n",M1);
    double a = 0;                                                // Begining of the interval
    double b = 0.80;                                             // End of the interval             
    double h = (b-a)/n;                                          // Lengh of a subinterval
    
    double t[n+1];                                               // Time array
    t[0] = a;

    // Aproximate solution for RK2 method:
    double w_initial_value_RK2[2*M1], w_RK2[2][n+1], k1, k2, l1, l2;
                                                                 // w_initial_value_RK2[2*M1] is the array with the coordinates x and y of each initial point for each curve                       
                                                                 // w_RK2[2][n+1] is the array where we are going to save the solution for each curve; we are going to overwrite this array, for each curve, in a for loop 
    // Initial points around all the domain:
    m1=1-s1,m2=-s1;
    for (p=0;m1<=1+s1;p++)   
    {
        for (q=0;m2<=s1;q++)   
        {
            w_initial_value_RK2[M2]   = m1;
            w_initial_value_RK2[M2+1] = m2;
            
            M2=M2+2;
            m2=m2+r1;   
        } 
        m1=m1+r1;
        m2=-s1;
    }

    m1=1-s1_bis,m2=-s1_bis;
    for (p=0;m1<=1+s1_bis;p++)   
    {
        for (q=0;m2<=s1_bis;q++)   
        {
            w_initial_value_RK2[M2]   = m1;
            w_initial_value_RK2[M2+1] = m2;
            
            M2=M2+2;
            m2=m2+r1_bis;   
        } 
        m1=m1+r1_bis;
        m2=-s1_bis;
    }

    // Initial points around the point (1,0)
    m1=1-s2,m2=-s2;
    for (p=0;m1<=1+s2;p++)   
    {
        for (q=0;m2<=s2;q++)   
        {
            w_initial_value_RK2[M2]   = m1;
            w_initial_value_RK2[M2+1] = m2;
            
            M2=M2+2;
            m2=m2+r2;   
        } 
        m1=m1+r2;
        m2=-s2;
    }

    // Initial points around the point (0,1)
    m1=-s3,m2=1-s3;
    for (p=0;m1<=s3;p++)   
    {
        for (q=0;m2<=1+s3;q++)   
        {
            w_initial_value_RK2[M2]   = m1;
            w_initial_value_RK2[M2+1] = m2;
            
            M2=M2+2;
            m2=m2+r3;   
        } 
        m1=m1+r3;
        m2=1-s3;
    }

    // Initial points around the point (0,-1)
    m1=-s4,m2=-1-s4;
    for (p=0;m1<=s4;p++)   
    {
        for (q=0;m2<=-1+s4;q++)   
        {
            w_initial_value_RK2[M2]   = m1;
            w_initial_value_RK2[M2+1] = m2;
            
            M2=M2+2;
            m2=m2+r4;   
        } 
        m1=m1+r4;
        m2=-1-s4;
    }

    M2=M2-1;
    
    printf("Creating file...\n"); 

    // Open file:

    pf = fopen("output_P3.csv", "w");
    if (pf == NULL)
    {
        printf("Error: the file cannot be opened\n");
        exit(1);
    }
    else
    { 
        // Writting on the file:

        for (j=1;j<n+1;j++)
        {
            t[j] = t[j-1] + h;
            fprintf(pf,"%f,", t[j-1]);
        }
        fprintf(pf,"%f\n",t[n]);

        for (m=0;m<M2;m=m+2)        
        {
            w_RK2[0][0] = w_initial_value_RK2[m];
            w_RK2[1][0] = w_initial_value_RK2[m+1];

            // Runge-Kutta with two stages (RK2):
            for (k=1;k<n+1;k++)                             
            {
                k1 = f1(w_RK2[0][k-1],w_RK2[1][k-1],t[k-1]);
                l1 = f2(w_RK2[0][k-1],w_RK2[1][k-1],t[k-1]);

                k2 = f1(w_RK2[0][k-1]+h*k1,w_RK2[1][k-1]+h*l1,t[k-1]);
                l2 = f2(w_RK2[0][k-1]+h*k1,w_RK2[1][k-1]+h*l1,t[k-1]);

                w_RK2[0][k]=w_RK2[0][k-1] + h*(0.5*k1+0.5*k2);
                w_RK2[1][k]=w_RK2[1][k-1] + h*(0.5*l1+0.5*l2);

                fprintf(pf,"%f,", w_RK2[0][k-1]);
            } 
            fprintf(pf,"%f\n", w_RK2[0][n]);

            for (k=1;k<n+1;k++)                             
            {
                fprintf(pf,"%f,", w_RK2[1][k-1]);
            } 
            fprintf(pf,"%f\n", w_RK2[1][n]);
        }
        fprintf(pf,"%d\n", M1);

        // Close file:
        fclose(pf);
        printf("File with the solution (created and closed)\n");
    }
    if (ferror(pf))
    {
        printf("Error when creating the file\n");
        clearerr(pf);
    }
    
    return 0;
}

double f1 (double w1_fun, double w2_fun, double t_fun_w) {       // y1'(t) = -y2(t)
    double funcion_w;

    if (w1_fun > T) w1_fun = T;
    if (w1_fun < -T) w1_fun = -T;
    if (w2_fun > T) w2_fun = T;
    if (w2_fun < -T) w2_fun = -T;

    funcion_w = -2*w1_fun*w2_fun;

    return(funcion_w);
}

double f2 (double w1_fun, double w2_fun, double t_fun_w) {       // y2'(t) = y1(t)
    double funcion_w;

    if (w1_fun > T) w1_fun = T;
    if (w1_fun < -T) w1_fun = -T;
    if (w2_fun > T) w2_fun = T;
    if (w2_fun < -T) w2_fun = -T;

    funcion_w = w2_fun*w2_fun + w1_fun - 1;

    return(funcion_w);
}

