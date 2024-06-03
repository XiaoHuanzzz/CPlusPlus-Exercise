#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <cblas.h>
#include <lapacke.h>
#include <complex.h>

#define PI (4*atan(1.0))

using namespace std;

class ABCG {
private:
    vector<lapack_complex_double> Hamilton;
    vector<lapack_complex_double> Psi;
    vector<lapack_complex_double> TPreconditioner;
    vector<lapack_complex_double> eigenValues;
    vector<lapack_complex_double> eigenVectors;
    int n;

public:
    ABCG(vector<vector<complex<double>>> H, vector<vector<complex<double>>> initialPsi, vector<vector<complex<double>>> TPre) {
        n = H.size();
        Hamilton.resize(n * n);
        Psi.resize(n * n);
        TPreconditioner.resize(n * n);
        eigenValues.resize(n);
        eigenVectors.resize(n * n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Hamilton[i * n + j] = lapack_make_complex_double(H[i][j].real(), H[i][j].imag());
                Psi[i * n + j] = lapack_make_complex_double(initialPsi[i][j].real(), initialPsi[i][j].imag());
                TPreconditioner[i * n + j] = lapack_make_complex_double(TPre[i][j].real(), TPre[i][j].imag());
            }
        }
    }

    ~ABCG() {
        // 自动管理内存的vector不需要手动delete
    }
    // calculate alpha*A - beta*B
    std::vector<lapack_complex_double> matrix_linear(const std::vector<lapack_complex_double>& A, const std::vector<lapack_complex_double>& B, lapack_complex_double alpha, lapack_complex_double beta, int rows, int cols) {
        if (A.size() != B.size() || A.size() != rows * cols) {
            throw std::invalid_argument("Matrices dimensions do not match");
        }

        std::vector<lapack_complex_double> result(A.size());

        for (size_t i = 0; i < A.size(); ++i) {
            result[i] = A[i] - B[i];
        }

        return result;
    }

    void process() {
        vector<lapack_complex_double> matrix(n * n);
        vector<lapack_complex_double> psiH(n);
        lapack_complex_double zero = lapack_make_complex_double(0.0, 0.0);
        lapack_complex_double one = lapack_make_complex_double(1.0, 0.0);
        vector<lapack_complex_double> zeros(n*n, lapack_make_complex_double(0.0, 0.0));

        // Line 1  
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                cblas_zgemv(CblasRowMajor, CblasNoTrans, n, n, &one, Hamilton.data(), n, &Psi[j * n], 1, &zero, psiH.data(), 1); // psiH = Hamilton * Psi
                cblas_zdotc_sub(n, &Psi[i * n], 1, psiH.data(), 1, &matrix[i * n + j]); // matrix = psiH * Psi
            }
        }

        // diagonalize
        vector<lapack_complex_double> eigenValuesTmp(n);
        vector<lapack_complex_double> eigenVectorsTmp(n * n);
        vector<lapack_complex_double> leftEigenVectorsTmp(n * n);

        int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, 
            matrix.data(), n, 
            eigenValues.data(),
            leftEigenVectorsTmp.data(), n,
            Psi.data(), n); // [eigenvalues, Psi] = diagonalize(matrix)

        if (info > 0) {
            cerr << "The algorithm failed to compute eigenvalues." << endl;
            exit(1);
        } 

        // Line 2
        vector<lapack_complex_double> Phi(n * n);
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Hamilton.data(), n, Psi.data(), n, &zero, Phi.data(), n); // Phi = Hamilton * Psi

        // Line 3
        vector<lapack_complex_double> lambda(n);
        vector<lapack_complex_double> P(n * n, lapack_make_complex_double(0, 0));

        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                // Line 4
                cblas_zdotc_sub(n, &Psi[i * n], 1, &Phi[i * n], 1, &lambda[i]); // lambda = (Phi, Psi)
            }

            // Line 5
            vector<lapack_complex_double> diagLambda(n * n, lapack_make_complex_double(0, 0));
            for (int i = 0; i < n; i++) {
                diagLambda[i * n + i] = lambda[i]; // diagLambda = diag(lambda)
            }

            vector<lapack_complex_double> R(n * n);
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Psi.data(), n, diagLambda.data(), n, &zero, R.data(), n);  
            R = matrix_linear(Phi, R, one, -one, n, n); // R = Phi - Psi * diag(lambda)

            // Line 6-7
            if (k == 0) {
                P = matrix_linear(zeros, P, one, -one, n, n); // P = -R; 
            } 
            // Line 8-10
            else {
              for(int i = 0; i < n; i++){
                lapack_complex_double beta1, beta0;
                cblas_zdotc_sub(n, &R[i*n], 1, &R[i*n], 1, &beta1); // beta1 = (R_i, R_i)
                cblas_zdotc_sub(n, &P[i*n], 1, &P[i*n], 1, &beta0); // beta0 = (P_i, P_i)
                P = matrix_linear(zeros, P, one, -beta1/beta0, n, n); // P = -beta1/beta0 * P;
                P = matrix_linear(R, P, one, one, n, n); // P = (R - beta1/beta0 * P)
                cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, TPreconditioner.data(), n, P.data(), n, &zero, P.data(), n); // P = -TPre * (R - beta1/beta0 * P) 
              }
            }

            // Line 11
            vector<lapack_complex_double> PsiP(n*n);
            cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, &one, Psi.data(), n, P.data(), n, &zero, PsiP.data(), n); // PsiP = Psi^T * P
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Psi.data(), n, PsiP.data(), n, &zero, PsiP.data(), n); // PsiP = Psi * PsiP
            P = matrix_linear(P, PsiP, one, -one, n, n); // P = P - Psi * PsiP

            // Line 12
            vector<lapack_complex_double> Theta(n*n);
            cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, &one, Hamilton.data(), n, P.data(), n, &zero, Theta.data(), n); // Theta = Hamilton * P

            // Line 13
            vector<lapack_complex_double> theta(n);
            for(int i = 0; i < n; i++){
              // Line 14
              lapack_complex_double PTheta, PhiTheta;
              cblas_zdotc_sub(n, &P[i*n], 1, &Theta[i*n], 1, &PTheta);
              cblas_zdotc_sub(n, &Phi[i*n], 1, &Theta[i*n], 1, &PhiTheta);
              if(creal(lambda[i]) - creal(PTheta) > 0){
                theta[i] = 2*PhiTheta / (lambda[i] - PTheta);
                theta[i] = atan(creal(theta[i]))/2 + PI/2;
              }
              else{
                theta[i] = 2*PhiTheta / (lambda[i] - PTheta);
                theta[i] = atan(creal(theta[i]))/2;
              } // Compute the optimal theta
              // Line 15-16
              for(int row = 0; row < n; row++){
                for(int col = 0; col < n; col++){
                  Psi[row*n + col] *= cos(creal(theta[i])) ;
                  Psi[row*n + col] += sin(creal(theta[i])) * P[row*n + col]; // Psi_i = Psi_i * cos(theta_i) + P_i * sin(theta_i)
                  Phi[row*n + col] *= cos(creal(theta[i]));
                  Phi[row*n + col] += sin(creal(theta[i])) * Theta[row*n + col]; // Phi_i = Phi_i * cos(theta_i) + Theta_i * sin(theta_i)
                }
              }
            }
            // Line 18
            cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, &one, Psi.data(), n, Psi.data(), n, &zero, R.data(), n); // R = Psi^T * Psi;
            LAPACKE_zpotrf(LAPACK_ROW_MAJOR, 'L', n, R.data(), n); // R = cholesky(R)
            // Line 19-20
            int ipiv[n];
            info = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, R.data(), n, ipiv);
            info = LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, R.data(), n, ipiv); // R = R^-1
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Psi.data(), n, R.data(), n, &zero, Psi.data(), n); // Psi = Psi * R^-1
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, &one, Phi.data(), n, R.data(), n, &zero, Phi.data(), n); // Phi = Phi * R^-1
        }
        // Line 21
        cblas_zgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, &one, Phi.data(), n, Phi.data(), n, &zero, Phi.data(), n); // Phi = Phi^T * Phi
        info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'V', 'V', n, 
            Phi.data(), n, 
            eigenValues.data(),
            leftEigenVectorsTmp.data(), n,
            eigenVectors.data(), n); // [eigenvalues, eigenVectors] = diagonalize(Phi^T * Phi)

        if (info > 0) {
            cerr << "The algorithm failed to compute eigenvalues." << endl;
            exit(1);
        } 
    }
};

int main(int argc, char* argv[]) {
    vector<vector<complex<double>>> H = {{{1,0},{2,1},{3,2}},{{2,3},{4,1},{2,9}},{{2,1},{4,2},{9,3}}}; // Hamiltonian matrix
    vector<vector<complex<double>>> initialPsi = {{{2,0},{2,2},{2,1}},{{1,0},{1,1},{1,2}},{{3,1},{3,2},{3,3}}}; // every row is a vector 
    vector<vector<complex<double>>> TPre = H; // Teter preconditioner
    ABCG solver(H, initialPsi, TPre);
    solver.process();

    return 0;
}

