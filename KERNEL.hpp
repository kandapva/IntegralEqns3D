/*
This header file provides acces to the entry of the 
system of linear equation from the discretized integral equation 
                Ax = b
-> Interactions
        Initialize using the object points_dt 
        getMatrixEntry(i,j) - i,j th entry of the matrix A.
        getKii() - provides the diagonal entry of the matrix
*/
#ifndef KERNEL_HPP
#define KERNEL_HPP

#include "points3D.hpp"
#include "Integral3D.hpp"

// This class holds all the interaction information (unsorted)
// The discretization of the integral equation is done as done in 
// https://github.com/klho/FLAM/blob/master/hifie/test/ie_cube1.m
class Interactions
{
points_dt *src;
double Kii;
double h,h3;
size_t N,n;
public:
Interactions(points_dt*& src_)
{
        this->src = src_;
        this->N = src->len();
        n = cbrt(N);
        h = 1.0/n;
        h3 = h*h*h;
        double *a,*b;
        a = new double[3];
        b = new double[3];
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        b[0] = h*0.5;
        b[1] = h*0.5;
        b[2] = h*0.5;
        Kii = triple_integral(a,b);
}
double getMatrixEntry(int i, int j)
{
        double output;
        if(i == j)
                output = Kii + 1; // First Kind integral equation to change it to second kind remove 1
        else
                output = h3 * KernelFunc(src->gridPoints[i], src->gridPoints[j]);
        return output;
}
double getKii()
{
        return Kii;
}
~Interactions(){
        //KD357
}
};

#endif
