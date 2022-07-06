#ifndef points3D_HPP
#define points3D_HPP

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include "Integral3D.hpp"

// Laplacian Kernel single layer in 2D
double Ker(const double r)
{
        return 1.0/(4.0*PI*r);
}


class pts3D
{
public:
double x,y,z;
pts3D(double x_,double y_,double z_)
{
        this->x = x_;
        this->y = y_;
        this->z = z_;
}

friend double KernelFunc(pts3D const &src,pts3D const &dst);
void print_point()
{
        std::cout << "(" << x << ","<< y << "," << z << ")" << std::endl;
}
~pts3D(){
}
};
double KernelFunc(pts3D const &src,pts3D const &dst)
{
        double output;
        double r = (src.x-dst.x)*(src.x-dst.x) + (src.y-dst.y)*(src.y-dst.y) + (src.z-dst.z)*(src.z-dst.z) ;
        // Select the Kernel function
        if(r == double(0))
                output = 0.0;
        else
        {
                r = sqrt(r);
                output = Ker(r);
        }
        return output;
}

class points_dt
{
int N;
Eigen::VectorXd Xdir,Ydir,Zdir;
public:
std::vector<pts3D> gridPoints;     // location of particles in the domain
points_dt(int nPoints, double xc, double yc,double zc,double L)
{
        //std::cout << "Uniform Nodes"<<std::endl;
        this->Xdir = unif_nodes(xc - L,xc + L,nPoints);
        this->Ydir = unif_nodes(yc - L,yc + L,nPoints);
        this->Zdir = unif_nodes(zc - L,zc + L,nPoints);
        std::cout << Xdir << std::endl;
        this->N = nPoints * nPoints * nPoints;
        for(size_t i=0; i<nPoints; i++)
                for(size_t j=0; j<nPoints; j++)
                {
                    for(size_t k=0; k<nPoints; k++)
                    {
                        pts3D temp(Xdir[i],Ydir[j],Zdir[k]);
                        gridPoints.push_back(temp);
                    }
                }
}
Eigen::VectorXd unif_nodes(double a,double b,int n)
{
        std::vector<double> X(n);
        Eigen::VectorXd R;
        double delta = (b-a)/double(n-1);
        for(int k=0; k<n; k++)
                X[k] = a + delta * k;
        R = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(X.data(), X.size());
        return R;
}
void print_points()
{
        for(size_t i=0; i<N; i++)
                gridPoints[i].print_point();

}

int len()
{
        return N;
}
~points_dt(){
}
};


#endif
