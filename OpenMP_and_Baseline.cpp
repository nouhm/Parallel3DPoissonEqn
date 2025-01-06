#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <chrono>
using namespace std;


long long currentTimeInMicroseconds() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    return duration.count();
}

int main(int argc, char *argv[]) {
    //constants
    int n, m, k; // eigen solutions
    n=2;m=2;k=2;
    int N=80; // dimension
    double L=1;
    double dx=L/(N-1);
    double Pi=  4*atan(1);
    //vector solutions 
    double *phi_o = new double [N*N*N];
    double *phi_n = new double [N*N*N];
    double *RHS = new double [N*N*N];

    for (int i=0;i<N*N*N;i++)
    {
        RHS[i]= -(n*n+m*m+k*k)*(Pi*Pi)*sin(n*Pi*(i%N)*dx)*cos(m*Pi*(((i-i%N)/N)%N)*dx)*sin(k*Pi*((i-i%(N*N))/(N*N))*dx); 
        // F is not zer at y boundaries
        if(((((i-i%N)/N)%N)==0)||((((i-i%N)/N)%N))==(N-1)){
            phi_o[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
            phi_n[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
        }
        else{
            phi_o[i]=0;
            phi_n[i]=0;
        }
        // cout <<endl<<RHS[i];
    }

    // test for RHS / Phi
    // for (int i=0;i<N;i++){
    //     cout<<"Z= "<<i*dx<<" "<<endl;
    //     for (int j=0;j<N;j++){
    //         for (int k=0;k<N;k++){
    //             cout<<phi_n[i*N*N+j*N+k]<<" ";
    //         }
    //         cout<<endl;
    //     }
    //     cout<<endl<<endl;
    // }

    // soultios

    double global_error = 20;
    int it = 1;
    // omp_get_wtime ()
    long long t1 = currentTimeInMicroseconds();
    while (global_error>0.0000001){
        global_error=0;
        // #pragma omp parallel for reduction(+:global_error) collapse(3)
        for (int i=1;i<N-1;i++){
            for (int j=1;j<N-1;j++){
                for (int k=1;k<N-1;k++){
                    if ((i+j+k)%2==0){
                                            phi_n[i*N*N+j*N+k] = (phi_o[i*N*N+j*N+k + 1] + phi_o[i*N*N+j*N+k - 1] +
                                            phi_o[i*N*N+(j + 1)*N+k] + phi_o[i*N*N+(j - 1)*N+k] +
                                            phi_o[(i+1)*N*N+j*N+k] + phi_o[(i-1)*N*N+j*N+k] -
                                            dx * dx * RHS[i*N*N+j*N+k]) / 6;
                                            global_error+=fabs(phi_n[i*N*N+j*N+k]-phi_o[i*N*N+j*N+k]);
                    }
                }
            }
        }

        // #pragma omp parallel for reduction(+:global_error) collapse(3)
        for (int i=1;i<N-1;i++){
            for (int j=1;j<N-1;j++){
                for (int k=1;k<N-1;k++){
                    if ((i+j+k)%2==1){
                                            phi_n[i*N*N+j*N+k] = (phi_o[i*N*N+j*N+k + 1] + phi_o[i*N*N+j*N+k - 1] +
                                            phi_o[i*N*N+(j + 1)*N+k] + phi_o[i*N*N+(j - 1)*N+k] +
                                            phi_o[(i+1)*N*N+j*N+k] + phi_o[(i-1)*N*N+j*N+k] -
                                            dx * dx * RHS[i*N*N+j*N+k]) / 6;
                                            global_error+=fabs(phi_n[i*N*N+j*N+k]-phi_o[i*N*N+j*N+k]);

                    }
                }
            }
        }
        //    #pragma omp for
        //  for (int i=1;i<N-1;i++){
        //         for (int j=1;j<N-1;j++){
        //             for (int k=1;k<N-1;k++){
        //                 global_error+=fabs(phi_n[i*N*N+j*N+k]-phi_o[i*N*N+j*N+k]);
        //             }
        //         }
        //     }  
        it++;
        double *temp = phi_o;
        phi_o = phi_n;
        phi_n = temp;
        // cout<<" Error is: "<< global_error<<endl;
        if ((it%50)==0)
        cout<<global_error<<endl;
    }
        long long t2 = currentTimeInMicroseconds() - t1;
            // cout<<"finel error is: "<<global_error<<endl;
            // cout<<"iterations: "<<it<<endl;
            // cout<<"Time taken in us: "<<t2<<endl;

        // test for RHS / Phi
        // for (int i=0;i<N;i++){
        //     cout<<"Z= "<<i*dx<<" "<<endl;
        //     for (int j=0;j<N;j++){
        //         for (int k=0;k<N;k++){
        //             cout<<phi_n[i*N*N+j*N+k]<<" ";
        //         }
        //         cout<<endl;
        //     }
        //     cout<<endl<<endl;
        // }


    cout << "Baseline with N = " << N << endl;
    cout << "Number of iterations: " << it << endl;
    cout << "Global error = " << global_error << endl;
    cout << "Time in microseconds: " << t2 << endl;
    return 0;
}