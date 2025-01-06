#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <chrono>
using namespace std;

double ***allocate_3d_array(int n) {
    double ***arr = new double **[n];
    for (int i = 0; i < n; ++i) {
        arr[i] = new double *[n];
        for (int j = 0; j < n; ++j) {
            arr[i][j] = new double[n];
        }
    }
    return arr;
}

long long currentTimeInMicroseconds() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    return duration.count();
}

int main(int argc, char *argv[]) {
    //constants
    int n, m, k; // eigen solutions
    n=1;m=1;k=1;
    int N=4; // dimension
    double L=1;
    double dx=L/(N-1);
    double Pi=  4*atan(1);
    //vector solutions 
    double ***phi_o =allocate_3d_array(N);
    double ***phi_n=allocate_3d_array(N);
    double ***RHS=allocate_3d_array(N);
    // double ***error=allocate_3d_array(N);
    // soultios

    double global_error = 20;
    int it = 1;

    for (int x=0;x<N;x++){
        for (int y=0;y<N;y++){
            for (int z=0; z<N;z++) {
                phi_o[x][y][z]=0;
                phi_n[x][y][z]=0;
                // RHS[i]= -(n*n+m*m+k*k)*(Pi*Pi)*sin(n*Pi*(i%N))*cos(m*Pi*(i-i%N)/N)*sin(k*Pi*(i-i%(N*N))/(N*N));
                RHS[x][y][z] = -(n*n+m*m+k*k)*(Pi*Pi)*sin(n*Pi*x*dx)*cos(m*Pi*y*dx)*sin(k*Pi*z*dx);
                // cout <<RHS[x][y][z]<<" ";
                if ((y==0)||(y==N-1)){ 
                    phi_n[x][y][z] = RHS[x][y][z]/(-3*Pi*Pi); 
                    phi_o[x][y][z] = RHS[x][y][z]/(-3*Pi*Pi);                         
                }
            }
            // cout<<endl;
        }   
        // cout<<endl;
    }

    long long t1= currentTimeInMicroseconds();

    while (global_error>0.0000001){
        global_error = 0.0;
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N-1 ; ++j) {
                for (int k = 1; k < N - 1; ++k) {
                    if ((i + j + k) % 2 == 0) { // Red points
                        phi_n[i][j][k] = (phi_o[i + 1][j][k] + phi_o[i - 1][j][k] +
                                            phi_o[i][j + 1][k] + phi_o[i][j - 1][k] +
                                            phi_o[i][j][k + 1] + phi_o[i][j][k - 1] -
                                            dx * dx * RHS[i][j][k]) / 6;   
                }
            }
            }
        }
        for (int i = 1; i < N - 1; ++i) {
            for (int j = 1; j < N-1 ; ++j) {
                for (int k = 1; k < N - 1; ++k) {
                    if ((i + j + k) % 2 == 1) { // Red points
                        phi_n[i][j][k] = (phi_o[i + 1][j][k] + phi_o[i - 1][j][k] +
                                            phi_o[i][j + 1][k] + phi_o[i][j - 1][k] +
                                            phi_o[i][j][k + 1] + phi_o[i][j][k - 1] -
                                            dx * dx * RHS[i][j][k]) / 6.0;
                            
                        }
                }
            }
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    global_error += fabs(phi_o[i][j][k] - phi_n[i][j][k]);
                    phi_o[i][j][k] = phi_n[i][j][k];
                }
            }
        }
        // for (int i=0;i<N*N*N-1;i++){
        //     if (i%2==0){
        //         phi_n[i] = (phi[i+1][j][k] + phi[i-1][j][k] +
        //             phi[i][j+1][k] + phi[i][j-1][k] +
        //             phi[i][j][k+1] + phi[i][j][k-1] -
        //             dx * dx * f[i][j][k]) / 6.0;
        //     }
        // }

        // cout << "Global error = " << global_error << endl;
        it++;
    }

    // for(int i=0; i<N;i++){
    //     cout<<"Z= "<<i*dx<<endl;
    //     for(int j=0; j<N; j++){
    //         for (int k = 0; k < N ; k++) {
    //             cout << phi_n[k][j][i] << " " ;
    //             // cout<<RHS[i][j][k]<<" ";                
    //         }
    //         cout << endl ; 
    //     }
    //     cout << endl ;
    // }

    long long t2 = currentTimeInMicroseconds()-t1;
    cout << "Baseline with N = " << N << endl;
    cout << "Number of iterations: " << it << endl;
    cout << "Global error = " << global_error << endl;
    cout << "Time in microseconds: " << t2 << endl;
    return 0;
}