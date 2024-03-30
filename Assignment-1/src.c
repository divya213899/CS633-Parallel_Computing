#include<stdio.h>
#include<stdlib.h>
#include "mpi.h"
#include<math.h>
#include<string.h>

//Compute step value at a point for 5-point Stencil
void computeSmol(double** temp_matrix,double** data, int N, double* up, double* down, double* left, double* right) {
    int denom_corner=5;
 // Compute the stencil for boundary points
    if(up==NULL && left == NULL){
        denom_corner= 3;
    }
   else if(up ==NULL || left == NULL){
        denom_corner= 4;
    }
   else if(down==NULL && right == NULL){
         denom_corner= 3;
    }
    else  if(down==NULL || right == NULL){
        denom_corner = 4;
    }

    // i=0 & j=0
    temp_matrix[0][0] = ((up != NULL ? up[0] : 0) + (left != NULL ? left[0] : 0) + data[1][0] + data[0][1] + data[0][0]) / denom_corner;

    // i=0 & j = N-1
    temp_matrix[0][N - 1] = ((up != NULL ? up[N-1] : 0) + (right != NULL ? right[0] : 0) + data[1][N - 1] + data[0][N - 2] + data[0][N - 1]) / denom_corner;

    // i = N-1, j =0
    temp_matrix[N - 1][0] = ((down != NULL ? down[0] : 0) + (left != NULL ? left[N-1] : 0) + data[N - 2][0] + data[N - 1][1] + data[N - 1][0]) / denom_corner;

    // (i=N-1 and j=N-1)
    temp_matrix[N - 1][N - 1] = ((down != NULL ? down[N-1] : 0) + (right != NULL ? right[N-1] : 0) + data[N - 2][N - 1] + data[N - 1][N - 2] + data[N - 1][N - 1]) / denom_corner;

    int denom_b = 5;

    // Boundary rows or columns
    // These data points must be surrounded by at least 3 neighboring points.
    if (up == NULL || down == NULL || left == NULL || right == NULL) {
        denom_b = 4;
    }
    
    for (int i = 1; i < N - 1; i++) {
        temp_matrix[0][i] = (data[0][i + 1] + data[0][i - 1] + data[1][i] + (up != NULL ? up[i] : 0) + data[0][i]) / denom_b;
        temp_matrix[N - 1][i] = (data[N - 1][i + 1] + data[N - 1][i - 1] + data[N - 2][i] + (down != NULL ? down[i] : 0) + data[N - 1][i]) / denom_b;
        temp_matrix[i][0] = (data[i + 1][0] + data[i - 1][0] + data[i][1] + (left != NULL ? left[i] : 0) + data[i][0]) / denom_b;
        temp_matrix[i][N - 1] = (data[i + 1][N - 1] + data[i - 1][N - 1] + data[i][N - 2] + (right != NULL ? right[i] : 0) + data[i][N - 1]) / denom_b;
    }

    // Compute the stencil for the interior points
    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            // Calculate the new value using the 5-point stencil formula
            temp_matrix[i][j] = (data[i - 1][j] + data[i + 1][j] +
                                 data[i][j - 1] + data[i][j + 1] +
                                 data[i][j]) / 5.0;
        }
    }
    
   
    //copying from temporary matrix back to original matrix
    for(int i = 0; i<N; i++){
    
        for(int j = 0; j<N; j++)
        {
            data[i][j] = temp_matrix[i][j];
        }
    }

}


//Compute step value at a point for 9-point Stencil

void computeBig(double** temp_matrix,double** dat,int n,double up[2][n], double down[2][n], double left[2][n], double right[2][n])
{  
    
    //for non halo elements
    for(int i = 2; i<n-2; i++)
    {
        for(int j = 2; j<n-2; j++)
        {
            temp_matrix[i][j] = (dat[i][j-2] + dat[i][j-1] + dat[i-2][j] + dat[i-1][j] + dat[i][j] + dat[i][j+1] + dat[i][j+2] + dat[i+1][j] + dat[i+2][j])/9;
        }
    }

    for(int j = 2; j<n-2; j++)
    {
        if(!up[0][0])
        {
            temp_matrix[0][j] = (dat[0][j-2] + dat[0][j-1] + dat[0][j] + dat[0][j+1] + dat[0][j+2] + dat[1][j] + dat[2][j])/7;
            temp_matrix[1][j] = (dat[1][j-2] + dat[1][j-1] + dat[0][j] + dat[1][j] + dat[1][j+1] + dat[1][j+2] + dat[2][j] + dat[3][j])/8;
        }else{
            temp_matrix[0][j] = (dat[0][j-2] + dat[0][j-1] + up[0][j] + up[1][j] + dat[0][j] + dat[0][j+1] + dat[0][j+2] + dat[1][j] + dat[2][j])/9;
            temp_matrix[1][j] = (dat[1][j-2] + dat[1][j-1] + up[1][j] + dat[0][j] + dat[1][j] + dat[1][j+1] + dat[1][j+2] + dat[2][j] + dat[3][j])/9;
        }

        if(!down[0][0])
        {
            temp_matrix[n-2][j] = (dat[n-2][j-2] + dat[n-2][j-1] + dat[n-3][j] + dat[n-4][j] + dat[n-2][j] + dat[n-2][j+1] + dat[n-2][j+2] + dat[n-1][j])/8;
            temp_matrix[n-1][j] = (dat[n-1][j-2] + dat[n-1][j-1] + dat[n-2][j] + dat[n-3][j] + dat[n-1][j] + dat[n-1][j+1] + dat[n-1][j+2])/7;
        }else{
            temp_matrix[n-2][j] = (dat[n-2][j-2] + dat[n-2][j-1] + dat[n-3][j] + dat[n-4][j] + dat[n-2][j] + dat[n-2][j+1] + dat[n-2][j+2] + dat[n-1][j] + down[0][j])/9;
            temp_matrix[n-1][j] = (dat[n-1][j-2] + dat[n-1][j-1] + dat[n-2][j] + dat[n-3][j] + dat[n-1][j] + dat[n-1][j+1] + dat[n-1][j+2] + down[0][j] + down[1][j])/9;
        }
    }
    for(int i = 2; i<n-2; i++)
    {
        if(!left[0][0])
        {
            temp_matrix[i][0] = (dat[i-1][0] + dat[i-2][0] + dat[i][0] + dat[i][1] + dat[i][2] + dat[i+1][0] + dat[i+2][0])/7;
            temp_matrix[i][1] = (dat[i][0] + dat[i-1][1] + dat[i-2][1] + dat[i][1] + dat[i][2] + dat[i][3] + dat[i+1][1] + dat[i+2][1])/8;
        }
        else{
            temp_matrix[i][0] = (left[0][i] + left[1][i] + dat[i-1][0] + dat[i-2][0] + dat[i][0] + dat[i][1] + dat[i][2] + dat[i+1][0] + dat[i+2][0])/9;
            temp_matrix[i][1] = (left[1][i] + dat[i][0] + dat[i-1][1] + dat[i-2][1] + dat[i][1] + dat[i][2] + dat[i][3] + dat[i+1][1] + dat[i+2][1])/9;
        }

        if(!right[0][0])
        {
            temp_matrix[i][n-2] = (dat[i][n-3] + dat[i][n-4] + dat[i-2][n-2] + dat[i-1][n-2] + dat[i][n-2] + dat[i][n-1] + dat[i+1][n-2] + dat[i+2][n-2])/8;
            temp_matrix[i][n-1] = (dat[i][n-2] + dat[i][n-3] + dat[i+1][n-1] + dat[i+2][n-1] + dat[i][n-1] + dat[i-1][n-1] + dat[i-2][n-1])/7;
        }
        else{
            temp_matrix[i][n-2] = (dat[i][n-3] + dat[i][n-4] + dat[i-2][n-2] + dat[i-1][n-2] + dat[i][n-2] + dat[i][n-1] + dat[i+1][n-2] + dat[i+2][n-2] + right[0][i])/9;
            temp_matrix[i][n-1] = (dat[i][n-2] + dat[i][n-3] + dat[i+1][n-1] + dat[i+2][n-1] + dat[i][n-1] + dat[i-1][n-1] + dat[i-2][n-1] + right[0][i] + right[1][i])/9;
  
        }
    }
//solving for Top-left corner
    if(!up[0][0])
    {
        
        if(!left[0][0])
        {
            temp_matrix[0][0] = (dat[0][0] + dat[0][1] + dat[0][2] + dat[1][0] + dat[2][0])/5;
            temp_matrix[0][1] = (dat[0][0] + dat[0][1] + dat[0][2] + dat[0][3] + dat[1][1] + dat[2][1])/6;
            temp_matrix[1][0] = (dat[0][0] + dat[1][0] + dat[1][1] + dat[1][2] + dat[2][0] + dat[3][0])/6;
            temp_matrix[1][1] = (dat[1][0] + dat[1][1] + dat[0][1] + dat[1][2] + dat[1][3] + dat[2][1] + dat[3][1])/7;
        }
        else {
            temp_matrix[0][0] = (left[0][0] + left[1][0] + dat[0][0] + dat[0][1] + dat[0][2] + dat[1][0] + dat[2][0])/7;
            temp_matrix[0][1] = (left[1][0] + dat[0][0] + dat[0][1] + dat[0][2] + dat[0][3] + dat[1][1] + dat[2][1])/7;
            temp_matrix[1][0] = (left[0][1] + left[1][1] + dat[0][0] + dat[1][0] + dat[1][1] + dat[1][2] + dat[2][0] + dat[3][0])/8;
            temp_matrix[1][1] = (left[1][1] + dat[1][0] + dat[1][1] + dat[0][1] + dat[1][2] + dat[1][3] + dat[2][1] + dat[3][1])/8;
        }

    }
    else{
        if(!left[0][0])
        {
            temp_matrix[0][0] = (up[0][0] + up[1][0] + dat[0][0] + dat[0][1] + dat[0][2] + dat[1][0] + dat[2][0])/7;
            temp_matrix[0][1] = (up[0][1] + up[1][1] + dat[0][0] + dat[0][1] + dat[0][2] + dat[0][3] + dat[1][1] + dat[2][1])/8;
            temp_matrix[1][0] = (up[1][0] + dat[0][0] + dat[1][0] + dat[1][1] + dat[1][2] + dat[2][0] + dat[3][0])/7;
            temp_matrix[1][1] = (up[1][1] + dat[1][0] + dat[1][1] + dat[0][1] + dat[1][2] + dat[1][3] + dat[2][1] + dat[3][1])/8;
        }
        else {
            temp_matrix[0][0] = (up[0][0] + up[1][0] + left[0][0] + left[1][0] + dat[0][0] + dat[0][1] + dat[0][2] + dat[1][0] + dat[2][0])/9;
            temp_matrix[0][1] = (up[0][1] + up[1][1] + left[1][0] + dat[0][0] + dat[0][1] + dat[0][2] + dat[0][3] + dat[1][1] + dat[2][1])/9;
            temp_matrix[1][0] = (up[1][0] + left[0][1] + left[1][1] + dat[0][0] + dat[1][0] + dat[1][1] + dat[1][2] + dat[2][0] + dat[3][0])/9;
            temp_matrix[1][1] = (up[1][1] + left[1][1] + dat[1][0] + dat[1][1] + dat[0][1] + dat[1][2] + dat[1][3] + dat[2][1] + dat[3][1])/9;
        }
    }
     // solving for Top-right corner

    if(!up[0][0])
    {
        if(!right[0][0])
        {
            temp_matrix[0][n-2] = (dat[0][n-2-2] + dat[0][n-2-1] + dat[0][n-2] + dat[0][n-1] + dat[1][n-2] + dat[2][n-2])/6;
            temp_matrix[0][n-1] = (dat[0][n-1-2] + dat[0][n-1-1] + dat[0][n-1] + dat[1][n-1] + dat[2][n-1])/5;
            temp_matrix[1][n-2] = (dat[1][n-2-2] + dat[1][n-2-1] + dat[1][n-2] + dat[1][n-1] + dat[0][n-2] + dat[2][n-2] + dat[3][n-2])/7;
            temp_matrix[1][n-1] = (dat[1][n-1-2] + dat[1][n-1-1] + dat[1][n-1] + dat[2][n-1] + dat[3][n-1] + dat[0][n-1])/6;
        }
        else{
            temp_matrix[0][n-2] = (right[0][0] + dat[0][n-2-2] + dat[0][n-2-1] + dat[0][n-2] + dat[0][n-1] + dat[1][n-2] + dat[2][n-2])/7;
            temp_matrix[0][n-1] = (right[0][0] + right[1][0] + dat[0][n-1-2] + dat[0][n-1-1] + dat[0][n-1] + dat[1][n-1] + dat[2][n-1])/7;
            temp_matrix[1][n-2] = (right[0][1] + dat[1][n-2-2] + dat[1][n-2-1] + dat[1][n-2] + dat[1][n-1] + dat[0][n-2] + dat[2][n-2] + dat[3][n-2])/8;
            temp_matrix[1][n-1] = (right[0][1] + right[1][1] + dat[1][n-1-2] + dat[1][n-1-1] + dat[1][n-1] + dat[2][n-1] + dat[3][n-1] + dat[0][n-1])/8;
        }
    }
    else{
        if(!right[0][0])
        {
            temp_matrix[0][n-2] = (up[0][n-2] + up[1][n-2] + dat[0][n-2-2] + dat[0][n-2-1] + dat[0][n-2] + dat[0][n-1] + dat[1][n-2] + dat[2][n-2])/8;
            temp_matrix[0][n-1] = (up[0][n-1] + up[1][n-1] + dat[0][n-1-2] + dat[0][n-1-1] + dat[0][n-1] + dat[1][n-1] + dat[2][n-1])/7;
            temp_matrix[1][n-2] = (up[1][n-2] + dat[1][n-2-2] + dat[1][n-2-1] + dat[1][n-2] + dat[1][n-1] + dat[0][n-2] + dat[2][n-2] + dat[3][n-2])/8;
            temp_matrix[1][n-1] = (up[1][n-1] + dat[1][n-1-2] + dat[1][n-1-1] + dat[1][n-1] + dat[2][n-1] + dat[3][n-1] + dat[0][n-1])/7;
        }
        else{
            temp_matrix[0][n-2] = (up[0][n-2] + up[1][n-2] + right[0][0] + dat[0][n-2-2] + dat[0][n-2-1] + dat[0][n-2] + dat[0][n-1] + dat[1][n-2] + dat[2][n-2])/9;
            temp_matrix[0][n-1] = (up[0][n-1] + up[1][n-1] + right[0][0] + right[1][0] + dat[0][n-1-2] + dat[0][n-1-1] + dat[0][n-1] + dat[1][n-1] + dat[2][n-1])/9;
            temp_matrix[1][n-2] = (up[1][n-2] + right[0][1] + dat[1][n-2-2] + dat[1][n-2-1] + dat[1][n-2] + dat[1][n-1] + dat[0][n-2] + dat[2][n-2] + dat[3][n-2])/9;
            temp_matrix[1][n-1] = (up[1][n-1] + right[0][1] + right[1][1] + dat[1][n-1-2] + dat[1][n-1-1] + dat[1][n-1] + dat[2][n-1] + dat[3][n-1] + dat[0][n-1])/9;
        }
    }
    
/// solving for bottom left
    if(!down[0][0])
    {
        if(!left[0][0])
        {
            temp_matrix[n-2][0] = (dat[n-3][0] + dat[n-4][0] + dat[n-2][0] + dat[n-2][0+1] + dat[n-2][0+2] + dat[n-1][0])/6;
            temp_matrix[n-2][1] = (dat[n-3][1] + dat[n-4][1] + dat[n-2][1] + dat[n-2][1+1] + dat[n-2][1+2] + dat[n-1][1] + dat[n-2][0])/7;
            temp_matrix[n-1][0] = (dat[n-2][0] + dat[n-3][0] + dat[n-1][0] + dat[n-1][0+1] + dat[n-1][0+2])/5;
            temp_matrix[n-1][1] = (dat[n-2][1] + dat[n-3][1] + dat[n-1][1] + dat[n-1][1+1] + dat[n-1][1+2] + dat[n-1][0])/6;
        }
        else{
            temp_matrix[n-2][0] = (left[0][n-2] + left[1][n-2] + dat[n-3][0] + dat[n-4][0] + dat[n-2][0] + dat[n-2][0+1] + dat[n-2][0+2] + dat[n-1][0])/8;
            temp_matrix[n-2][1] = (left[1][n-2] + dat[n-3][1] + dat[n-4][1] + dat[n-2][1] + dat[n-2][1+1] + dat[n-2][1+2] + dat[n-1][1] + dat[n-2][0])/8;
            temp_matrix[n-1][0] = (left[0][n-1] + left[1][n-1] + dat[n-2][0] + dat[n-3][0] + dat[n-1][0] + dat[n-1][0+1] + dat[n-1][0+2])/7;
            temp_matrix[n-1][1] = (left[1][n-1] + dat[n-2][1] + dat[n-3][1] + dat[n-1][1] + dat[n-1][1+1] + dat[n-1][1+2] + dat[n-1][0])/7;
        }
    }
    else{
        if(!left[0][0])
        {
            temp_matrix[n-2][0] = (down[0][0] + dat[n-3][0] + dat[n-4][0] + dat[n-2][0] + dat[n-2][0+1] + dat[n-2][0+2] + dat[n-1][0])/7;
            temp_matrix[n-2][1] = (down[0][1] + dat[n-3][1] + dat[n-4][1] + dat[n-2][1] + dat[n-2][1+1] + dat[n-2][1+2] + dat[n-1][1] + dat[n-2][0])/8;
            temp_matrix[n-1][0] = (down[0][0] + down[1][0] + dat[n-2][0] + dat[n-3][0] + dat[n-1][0] + dat[n-1][0+1] + dat[n-1][0+2])/7;
            temp_matrix[n-1][1] = (down[0][1] + down[1][1] + dat[n-2][1] + dat[n-3][1] + dat[n-1][1] + dat[n-1][1+1] + dat[n-1][1+2] + dat[n-1][0])/8;
        }
        else{
            temp_matrix[n-2][0] = (down[0][0] + left[0][n-2] + left[1][n-2] + dat[n-3][0] + dat[n-4][0] + dat[n-2][0] + dat[n-2][0+1] + dat[n-2][0+2] + dat[n-1][0])/9;
            temp_matrix[n-2][1] = (down[0][1] + left[1][n-2] + dat[n-3][1] + dat[n-4][1] + dat[n-2][1] + dat[n-2][1+1] + dat[n-2][1+2] + dat[n-1][1] + dat[n-2][0])/9;
            temp_matrix[n-1][0] = (down[0][0] + down[1][0] + left[0][n-1] + left[1][n-1] + dat[n-2][0] + dat[n-3][0] + dat[n-1][0] + dat[n-1][0+1] + dat[n-1][0+2])/9;
            temp_matrix[n-1][1] = (down[0][1] + down[1][1] + left[1][n-1] + dat[n-2][1] + dat[n-3][1] + dat[n-1][1] + dat[n-1][1+1] + dat[n-1][1+2] + dat[n-1][0])/9;
        }
    }
// solving for bottom right
    if(!down[0][0])
    {
        if(!right[0][0])
        {
            temp_matrix[n-2][n-2] = (dat[n-2][n-2-2] + dat[n-2][n-2-1] + dat[n-2][n-2] + dat[n-2-1][n-2] + dat[n-2-2][n-2] + dat[n-2][n-1] + dat[n-1][n-2])/7;
            temp_matrix[n-2][n-1] = (dat[n-2][n-1-2] + dat[n-2][n-1-1] + dat[n-2][n-1] + dat[n-2-1][n-1] + dat[n-2-2][n-1] + dat[n-1][n-1])/6;
            temp_matrix[n-1][n-2] = (dat[n-1][n-2-2] + dat[n-1][n-2-1] + dat[n-1][n-2] + dat[n-1-1][n-2] + dat[n-1-2][n-2] + dat[n-1][n-1])/6;
            temp_matrix[n-1][n-1] = (dat[n-1][n-1-2] + dat[n-1][n-1-1] + dat[n-1][n-1] + dat[n-1-1][n-1] + dat[n-1-2][n-1])/5;
        }
        else{
            temp_matrix[n-2][n-2] = (right[0][n-2] + dat[n-2][n-2-2] + dat[n-2][n-2-1] + dat[n-2][n-2] + dat[n-2-1][n-2] + dat[n-2-2][n-2] + dat[n-2][n-1] + dat[n-1][n-2])/8;
            temp_matrix[n-2][n-1] = (right[0][n-2] + right[1][n-2] + dat[n-2][n-1-2] + dat[n-2][n-1-1] + dat[n-2][n-1] + dat[n-2-1][n-1] + dat[n-2-2][n-1] + dat[n-1][n-1])/8;
            temp_matrix[n-1][n-2] = (right[0][n-1] + dat[n-1][n-2-2] + dat[n-1][n-2-1] + dat[n-1][n-2] + dat[n-1-1][n-2] + dat[n-1-2][n-2] + dat[n-1][n-1])/7;
            temp_matrix[n-1][n-1] = (right[0][n-1] + right[1][n-1] + dat[n-1][n-1-2] + dat[n-1][n-1-1] + dat[n-1][n-1] + dat[n-1-1][n-1] + dat[n-1-2][n-1])/7;
        }
    }
    else{
        if(!right[0][0])
        {
            temp_matrix[n-2][n-2] = (down[0][n-2] + dat[n-2][n-2-2] + dat[n-2][n-2-1] + dat[n-2][n-2] + dat[n-2-1][n-2] + dat[n-2-2][n-2] + dat[n-2][n-1] + dat[n-1][n-2])/8;
            temp_matrix[n-2][n-1] = (down[0][n-1] + dat[n-2][n-1-2] + dat[n-2][n-1-1] + dat[n-2][n-1] + dat[n-2-1][n-1] + dat[n-2-2][n-1] + dat[n-1][n-1])/7;
            temp_matrix[n-1][n-2] = (down[0][n-2] + down[1][n-2] + dat[n-1][n-2-2] + dat[n-1][n-2-1] + dat[n-1][n-2] + dat[n-1-1][n-2] + dat[n-1-2][n-2] + dat[n-1][n-1])/8;
            temp_matrix[n-1][n-1] = (down[0][n-1] + down[1][n-1] + dat[n-1][n-1-2] + dat[n-1][n-1-1] + dat[n-1][n-1] + dat[n-1-1][n-1] + dat[n-1-2][n-1])/7;
        }
        else{
            temp_matrix[n-2][n-2] = (right[0][n-2] + down[0][n-2] + dat[n-2][n-2-2] + dat[n-2][n-2-1] + dat[n-2][n-2] + dat[n-2-1][n-2] + dat[n-2-2][n-2] + dat[n-2][n-1] + dat[n-1][n-2])/9;
            temp_matrix[n-2][n-1] = (right[0][n-2] + right[1][n-2] + down[0][n-1] + dat[n-2][n-1-2] + dat[n-2][n-1-1] + dat[n-2][n-1] + dat[n-2-1][n-1] + dat[n-2-2][n-1] + dat[n-1][n-1])/9;
            temp_matrix[n-1][n-2] = (right[0][n-1] + down[0][n-2] + down[1][n-2] + dat[n-1][n-2-2] + dat[n-1][n-2-1] + dat[n-1][n-2] + dat[n-1-1][n-2] + dat[n-1-2][n-2] + dat[n-1][n-1])/9;
            temp_matrix[n-1][n-1] = (right[0][n-1] + right[1][n-1] + down[0][n-1] + down[1][n-1] + dat[n-1][n-1-2] + dat[n-1][n-1-1] + dat[n-1][n-1] + dat[n-1-1][n-1] + dat[n-1-2][n-1])/9;
        }
    }


    //copying from temporary matrix back to original matrix
    for(int i = 0; i<n; i++)
    {
        for(int j = 0; j<n; j++)
        {
            dat[i][j] = temp_matrix[i][j];
        }
    }
    
    return;
}
int main(int argc, char* argv[])
{
    MPI_Init(&argc,&argv);

    int myrank, size, Px, N, seed, steps, stencil;

    MPI_Status status;
    MPI_Request request;
    double maxTime;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    Px = atoi(argv[1]); 
    int n = atoi(argv[2]);
    N = sqrt((double)n);
    steps = atoi(argv[3]);
    seed = atoi(argv[4]);
    stencil = atoi(argv[5]); 
    double **data = (double **)malloc(N * sizeof(double *));
    for(int i = 0; i<N; i++)
    {
        data[i] = (double *)malloc(N * sizeof(double));

        for(int j = 0; j<N; j++)
        {
            srand(seed*(myrank + 10));
            data[i][j] = abs(rand() + (i*rand() + j*myrank))/100;
        }
    }
    double** temp_mat= (double**)malloc(N * sizeof(double*));
    
    int lflag=0,rflag=0,uflag=0,dflag=0;
    if(myrank>Px-1)
    {
        uflag=1;
    }
    if(myrank<Px*(size/Px - 1))
    {
        dflag=1;
    }
    if(myrank%Px)
    {
        lflag=1;
    }
    if((myrank+1)%Px)
    {
        rflag=1;
    }
    double** temp_matrix = (double**)malloc(N* sizeof(double*));
    for (int i = 0; i < N; i++) {
        temp_matrix[i] = (double*)malloc(N * sizeof(double));
        if (temp_matrix[i] == NULL) {
            printf("Memory allocation failed\n");
            exit(1);
        }
    }
    if(stencil == 5)
    {
        double left[N],right[N],up[N],down[N];
        double send_left[N],send_right[N],send_up[N],send_down[N];
        double recv_left[N],recv_right[N],recv_up[N],recv_down[N];
        int t=0;
        double stime = MPI_Wtime();

        while(t<steps)
        {
            int position = 0;
            if(uflag) // there is a process above
            {
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[0] + i, 1, MPI_DOUBLE, send_up, N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_up,position,MPI_PACKED,myrank-Px,myrank-Px,MPI_COMM_WORLD,&request);
            }
            if(dflag)
            {
                position = 0;
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[N-1] + i, 1, MPI_DOUBLE, send_down, N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_down, position, MPI_PACKED, myrank+Px, myrank+Px, MPI_COMM_WORLD,&request);
            }
            if(lflag)
            {
                position = 0;
                for(int i = 0; i < N; i++)
                {
                    MPI_Pack(data[i] + 0, 1, MPI_DOUBLE, send_left, N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_left, position, MPI_PACKED, myrank-1, myrank-1, MPI_COMM_WORLD,&request);
            }
            if(rflag)
            {
                position = 0;
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[i] + N-1, 1, MPI_DOUBLE, send_right, N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_right, position, MPI_PACKED, myrank+1, myrank+1, MPI_COMM_WORLD,&request);
            }
            // recieving data from neighbouring processors
            if(uflag)
            {
                position=0;
                MPI_Recv(recv_up,N*sizeof(double),MPI_PACKED,myrank-Px,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_up,N*sizeof(double),&position,up+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(dflag)
            {
                position=0;
                MPI_Recv(recv_down,N*sizeof(double),MPI_PACKED,myrank+Px,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_down,N*sizeof(double),&position,down+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(lflag)
            {
                position=0;
                MPI_Recv(recv_left,N*sizeof(double),MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_left,N*sizeof(double),&position,left+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(rflag)
            {
                position=0;
                MPI_Recv(recv_right,N*sizeof(double),MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_right,N*sizeof(double),&position,right+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }

            computeSmol(temp_matrix,data, N, up, down, left, right);

            //increasing step time
            t++;            
        }
        double time = MPI_Wtime() - stime;
        MPI_Reduce(&time, &maxTime, 1*sizeof(double), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(myrank == 0)
        {
            printf("%lf\n", maxTime);
        }
    }
    else{
        double send_left[2*N],send_right[2*N],send_up[2*N],send_down[2*N];
        double recv_left[2*N],recv_right[2*N],recv_up[2*N],recv_down[2*N];
        double left[2][N],right[2][N],up[2][N],down[2][N];
        int t = 0; 
        double stime = MPI_Wtime();

        while(t<steps)
        {
            int position = 0;
            if(uflag)
            {
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[0] + i, 1, MPI_DOUBLE, send_up, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[1] + i, 1, MPI_DOUBLE, send_up, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_up,position,MPI_PACKED,myrank-Px,myrank-Px,MPI_COMM_WORLD,&request);
            }
            if(dflag)
            {
                position = 0;
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[N-2] + i, 1, MPI_DOUBLE, send_down, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                for(int i = 0; i<N; i++)
                {
                    MPI_Pack(data[N-1] + i, 1, MPI_DOUBLE, send_down, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_down, position, MPI_PACKED, myrank+Px, myrank+Px, MPI_COMM_WORLD,&request);
            }
            if(lflag)
            {
                position = 0;
                for(int i = 0; i < N; i++)
                {
                    MPI_Pack(data[i] + 0, 1, MPI_DOUBLE, send_left, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                for(int i = 0; i < N; i++)
                {
                    MPI_Pack(data[i] + 1, 1, MPI_DOUBLE, send_left, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_left, position, MPI_PACKED, myrank-1, myrank-1, MPI_COMM_WORLD,&request);
            }
            if(rflag)
            {
                position = 0;
                for(int i = 0; i < N; i++)
                {
                    MPI_Pack(data[i] + N-2, 1, MPI_DOUBLE, send_right, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                for(int i = 0; i < N; i++)
                {
                    MPI_Pack(data[i] + N-1, 1, MPI_DOUBLE, send_right, 2*N*sizeof(double), &position, MPI_COMM_WORLD);
                }
                MPI_Isend(send_right, position, MPI_PACKED, myrank+1, myrank+1, MPI_COMM_WORLD,&request);
            }
            //recieving data from neighbouring processors
            if(uflag)
            {
                position = 0;
                MPI_Recv(recv_up, 2*N*sizeof(double), MPI_PACKED, myrank-Px, myrank, MPI_COMM_WORLD, &status);
                for(int i = 0; i<N; i++)
                {
                    MPI_Unpack(recv_up, 2*N*sizeof(double), &position, up[0]+i, 1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
                for(int i = 0; i<N; i++)
                {
                    MPI_Unpack(recv_up, 2*N*sizeof(double), &position, up[1]+i, 1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(dflag)
            {
                position=0;
                MPI_Recv(recv_down,2*N*sizeof(double),MPI_PACKED,myrank+Px,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_down,2*N*sizeof(double),&position,down[0]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_down,2*N*sizeof(double),&position,down[1]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(lflag)
            {
                position = 0;
                MPI_Recv(recv_left,2*N*sizeof(double),MPI_PACKED,myrank-1,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_left,2*N*sizeof(double),&position,left[0]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_left,2*N*sizeof(double),&position,left[1]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            if(rflag)
            {
                position=0;
                MPI_Recv(recv_right,2*N*sizeof(double),MPI_PACKED,myrank+1,myrank,MPI_COMM_WORLD,&status);
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_right,2*N*sizeof(double),&position,right[0]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
                for(int i=0;i<N;i++)
                {
                    MPI_Unpack(recv_right,2*N*sizeof(double),&position,right[1]+i,1,MPI_DOUBLE,MPI_COMM_WORLD);
                }
            }
            computeBig(temp_matrix,data,N, up, down, left, right);

            //increasing step time
            t++;
        }
        double time = MPI_Wtime() - stime;
        MPI_Reduce(&time, &maxTime,1*sizeof(double), MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if(myrank == 0)
        {
            printf("%lf\n", maxTime);
        }
    }

    MPI_Finalize();
    return 0;
}
