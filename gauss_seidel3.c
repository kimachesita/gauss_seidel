/* Gauss Seidel Implementation 
    Solving Ax = b given
    A, b, Xo,tol , max_iteration
*/

#include<stdio.h>
#include<math.h>

#define MAX_SIZE 10


int main(){

    // Initialize Variables
    int n;                              // number of unknowns and equations
    float A[MAX_SIZE][MAX_SIZE];        // array containing A
    float b[MAX_SIZE];                  // array containing b
    float x[MAX_SIZE];                  // solution
    float xo[MAX_SIZE];                 // guess for x
    float tol;                          // tolerance
    int max_iter;                       // max_iter
    
    int i,j,k;
    float sigma,norm;
    
    
    // Ask for number of unknowns
    printf("\nEnter the value of n : \n");
    scanf("%d",&n);
    
    // Ask for coefficient matrix
    printf("\nEnter the coefficients row wise : \n");
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            scanf("%f",&A[i][j]);
        }
    } 
    
    // Ask for constant matrix
    printf("\nEnter the right hand side constants : \n");
    for(i=0;i<n;i++) {
        scanf("%f",&b[i]);
    }
    
    // Ask for guess solution
    printf("\nEnter initial/guess solution : \n");
    for(i=0;i<n;i++) {
        scanf("%f",&x[i]);
    }
    
    // Ask for tolerance
    printf("\nEnter tolerance : \n");
    scanf("%f",&tol);
    
    // Ask for max_iter
    printf("\nEnter maximum iteration : \n");
    scanf("%d",&max_iter); 
       
   
    //Check Inputs
    printf("\nSolving Linear Equation with the ff. parameters:\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            printf("%fx%d ",A[i][j],j);
        }
        printf("= ");
        printf("%f\n",b[i]);
    }
    
    printf("tol: %f\n",tol);
    printf("max_iter %d\n",max_iter);
    printf("Initial Val:");
    for(i=0;i<n;i++){
        printf("%f ",x[i]);
    }
    printf("\n");

   
    //begin iterate process
    k = 1;
    while(k <= max_iter){
        for(i=0;i<n;i++){
            xo[i] = x[i];
        }
        for(i=0;i<n;i++){
            sigma = 0;
            for(j=0;j <= i-1; j++){
                sigma = sigma + A[i][j]*x[j];
            }
            for(j = i + 1; j < n;j++){
                sigma = sigma + A[i][j] * xo[j];
            }
            x[i] = (1/A[i][i])*(b[i] - sigma );         
        }
        
        //compute for norm
        sigma = 0;
        for(i=0;i<n;i++){
            sigma = sigma + pow((xo[i]-x[i]),2);
        }
        norm = sqrt(sigma);
        
        if( norm <= tol){
            break;
        }
        k++;
        
    }
    
    //print solution
    printf("\nSolution of the System is:\n");
    for(i=0;i<n;i++){
        printf("\nx[%d]: %f",i,x[i]);
    }
    printf("\nIterations: %d",k-1);
    printf("\nNorm: %f", norm);
    
    return 0;
    
}