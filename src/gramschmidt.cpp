#include "matrix.h"
#include <math.h>

void gramschmidt(const Matrix<float>& a, Matrix<float>& b)	//b is the orthogonalized matrix
{
    int noofbasis = a.getRows();
    int ineachbasis = a.getCols();
    
    float dot = 0;
    for(int i=0; i<noofbasis; ++i)
    {
        for(int j=0; j<ineachbasis; j++)
        {
            b[i][j] = a[i][j];
        }
        for(int k=i+1; k<noofbasis; k++)
        {
            dot = 0;
            for(int j=0; j<ineachbasis; j++)
            {
                dot+=(a[k][j]*b[i][j]);
            }
            for(int j=0; j<ineachbasis; j++)
            {
                a[k][j]-=(dot*b[i][j]);
            }
        }
    }
}
