#ifndef _PCO_PCA_H_FILE_
#define _PCO_PCA_H_FILE_

#include <math.h>
# define PCA_EPS    0.000001   /* 计算精度 */
# define PCA_Iteration    100	 /* 迭代次数 */

int Tri_Symmetry_Diagonal_Eigenvector( int n, double b[], double c[], double q[], double eps, int l );
void Householder_Tri_Symetry_Diagonal( double a[], int n, double q[], double b[], double c[] );
int SymmetricRealMatrix_Eigen( double CovMatrix[], int n, double Eigen[], double EigenVector[] );
int PCA( double * data, int row, int col, double * PC, double * V );

#endif // !_PCO_PCA_H_FILE_