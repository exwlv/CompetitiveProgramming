#include <iostream>
#include <cmath>

typedef double dfunction(double[3]);

// f0, f1 and f2 implement the equations
// x0+x1+x2=0.0
// x0^2+x1^2+x2^2-6.0=0.0
// x0^3+x1^3+x2^3-6.0=0.0
double f0(double x[3]) { return x[0] + x[1] + x[2]; }
double f1(double x[3]) { return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] -6.0; }
double f2(double x[3]) { return x[0] * x[0] * x[0] + x[1] * x[1] * x[1] + x[2] * x[2] * x[2] - 6.0; }

// dfij is the partial derivative of fi with respect to x[j]
double df00(double x[3]) { return 1.0; }
double df01(double x[3]) { return 1.0; }
double df02(double x[3]) { return 1.0; }
double df10(double x[3]) { return 2.0 * x[0]; }
double df11(double x[3]) { return 2.0 * x[1]; }
double df12(double x[3]) { return 2.0 * x[2]; }
double df20(double x[3]) { return 3.0 * x[0] * x[0]; }
double df21(double x[3]) { return 3.0 * x[1] * x[1]; }
double df22(double x[3]) { return 3.0 * x[2] * x[2]; }

// MatrixMul multiplies a 3x3 matrix with a 3 dimensional vector
// and stores the result in ans
void MatrixMul(double A[3][3], double x[3], double ans[3]) {

    ans[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2];
    ans[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2];
    ans[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2];

}

// CalculateInverse inverts a 3x3 matrix using the adjoint matrix
// the return value indicates if the matrix could be inverted
int CalculateInverse(double A[3][3], double Ainv[3][3]) {

    double det=(  A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0]
                 +A[0][2] * A[1][0] * A[2][1] - A[0][0] * A[1][2] * A[2][1]
                 -A[0][1] * A[1][0] * A[2][2] - A[0][2] * A[1][1] * A[2][0] );

    if(det == 0.0)
        return 0;

    Ainv[0][0] =  (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / det;
    Ainv[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / det;
    Ainv[0][2] =  (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / det;
    Ainv[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / det;
    Ainv[1][1] =  (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / det;
    Ainv[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / det;
    Ainv[2][0] =  (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / det;
    Ainv[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]) / det;
    Ainv[2][2] =  (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / det;

    return 1;

}

int Newton(double lambda[3], dfunction *f[3], dfunction *df[3][3], double eps, int maxiter) {

    double temp1[3], temp2[3], m[3][3], minv[3][3];

    int i=0;

    while(((fabs(f[0](lambda)) > eps) || (fabs(f[1](lambda)) > eps) || (fabs(f[2](lambda)) > eps)) && i < maxiter) {

        // calculate the partial derivatives
        m[0][0] = df[0][0](lambda);
        m[0][1] = df[0][1](lambda);
        m[0][2] = df[0][2](lambda);

        m[1][0] = df[1][0](lambda);
        m[1][1] = df[1][1](lambda);
        m[1][2] = df[1][2](lambda);

        m[2][0] = df[2][0](lambda);
        m[2][1] = df[2][1](lambda);
        m[2][2] = df[2][2](lambda);

        // evaluate the functions
        temp1[0] = f[0](lambda);
        temp1[1] = f[1](lambda);
        temp1[2] = f[2](lambda);
        
        if(!CalculateInverse(m, minv))
            return 0;

        MatrixMul(minv, temp1, temp2);

        lambda[0] -= temp2[0];
        lambda[1] -= temp2[1];
        lambda[2] -= temp2[2];

        i++;

    }
    
    if (i >= maxiter)
        return 0;
    
    return 1;

}

int main() {

    dfunction *f[3], *df[3][3];

    // for every one solution found there are two more
    // we choose the initial values to satisfy the first equation
    double lambda[3] = { 0, 1, -1 };

    f[0] = f0;
    f[1] = f1;
    f[2] = f2;

    df[0][0] = df00;
    df[0][1] = df01;
    df[0][2] = df02;

    df[1][0] = df10;
    df[1][1] = df11;
    df[1][2] = df12;

    df[2][0] = df20;
    df[2][1] = df21;
    df[2][2] = df22;

    if(Newton(lambda, f, df, 1e-50, 100)) {

        std::cout << "Lambda 1 : " << lambda[0] << std::endl;
        std::cout << "Lambda 2 : " << lambda[1] << std::endl;
        std::cout << "Lambda 3 : " << lambda[2] << std::endl;

    } else
        std::cout<<"Could not find a solution " << "with the required accuracy." << std::endl;

    return 0;

}

