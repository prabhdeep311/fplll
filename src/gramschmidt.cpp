#include "matrix.h"
#include <vector>
#include <math.h>

float calcmagnitude(vector<float> v1)
{
    float f = 0;
    for(int i = 0; i<v1.size(); i++){
        f += (v1[i]*v1[i]);
    }
    return sqrtf(f);
}

void gramschmidt(const Matrix<float>& a, Matrix<float>& b)	//b is the orthogonalized matrix
{
    int noOfBasis = a.getCols();
    int inEachBasis = a.getRows();
    vector<float> v;
    
    float mag = 0, dot = 0;
    for (int i = 0; i < noOfBasis; ++i){
        v.clear();
	for (int j = 0; j < inEachBasis; j++){
            v.push_back(a[j][i]);
        }
        mag = calcmagnitude(v);
        for (int j = 0; j < inEachBasis; j++){
            b[j][i] = a[j][i]/mag;
        }
        for(int k = i+1; k < noOfBasis; k++){
            dot = 0;
            for(int j = 0; j < inEachBasis; j++){
                dot += (a[j][k]*b[j][i]);
            }
            for(int j = 0; j < inEachBasis; j++){
                a[j][k] -= (dot*b[j][i]);
            }
        }
    }
}