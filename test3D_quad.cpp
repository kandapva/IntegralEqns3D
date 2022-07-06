#include <iostream>
#include "Integral3D.hpp"
#include "KERNEL.hpp"
#include "points3D.hpp" 

using namespace std;

int main()
{

points_dt* src;
src = new points_dt(2,0.75,0.75,0.75,0.25);
// Creates a uniform nodes in 3D with center (0.75,.075,.075) and half length 0.25 with 2 points
Interactions* tst3D;
tst3D = new Interactions(src);


double mat_t;
    for (int j=0; j < 8; ++j)
        {
                for (int k=0; k < 8; ++k)
                {
                        mat_t = tst3D->getMatrixEntry(j , k);
                        cout << mat_t << "\t";
                }
                cout << endl;
        }
    return 0;
}