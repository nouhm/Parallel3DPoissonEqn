

#include <iostream>
#include <cmath>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <ctime>
#include <mpi.h>
#include <chrono>

using namespace std;

long long currentTimeInMicroseconds() {
    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    return duration.count();
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // We need exactly 8 ranks
    if (size != 8) {
        if (rank == 0) {
            cerr << "This program requires exactly 8 MPI ranks." << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Problem parameters
    int n = 2, m = 2, k = 2; 
    int N = 5; // global dimension
    double L = 1.0;
    double dx = L / (N - 1);
    double Pi = 4 * atan(1.0);

    // Allocate arrays
    double *phi_o = new double[N * N * N];
    double *phi_n = new double[N * N * N];
    double *RHS   = new double[N * N * N];

    int local_N = N / size; // Decomposition along y-axis
    int starty  = rank * local_N;
    int endy    = (rank + 1) * local_N;

    // Fill RHS and initial guess
    // i-th index corresponds to: 
    // i = k + N*(j + N*i_z) with i_z = i/N^2, j = (i/N)%N, k = i%N
    for (int i = 0; i < N * N * N; i++) {
        int kx = i % N;
        int jy = (i / N) % N;
        int iz = i / (N * N);

        double x = kx * dx;
        double y = jy * dx;
        double z = iz * dx;

        // Given PDE: RHS as per the original code
        RHS[i] = -(n * n + m * m + k * k) * (Pi * Pi)
                 * sin(n * Pi * x) * cos(m * Pi * y) * sin(k * Pi * z);

        // Apply Dirichlet boundary conditions at y=0 and y=L
        if (jy == 0 || jy == N - 1) {
            phi_o[i] = RHS[i] / (-(m * m + n * n + k * k) * Pi * Pi);
            phi_n[i] = phi_o[i];
        } else {
            // interior points
            phi_o[i] = 0.0;
            phi_n[i] = 0.0;
        }
    }

    // Broadcast initial data (if needed). Here, data is already consistent across all ranks.
    MPI_Bcast(RHS, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi_o, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi_n, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Setup neighbors
    int up   = (rank == 0) ? MPI_PROC_NULL : rank - 1;
    int down = (rank == size - 1) ? MPI_PROC_NULL : rank + 1;

    double *ysend_down = new double[N * N]; 
    double *yrecv_down = new double[N * N];
    double *ysend_up   = new double[N * N]; 
    double *yrecv_up   = new double[N * N];

    double global_error = 1e10;
    double tolerance = 1e-7;
    int it = 1;

    long long t1 = currentTimeInMicroseconds();

    while (global_error > tolerance) {
        it++;
        double error = 0.0;

        // Halo exchange along y-direction
        // Prepare data to send: we send boundary planes of phi_o at starty and endy-1
        // Downward neighbor exchange (to 'down')
        int reqCount = 0;
        MPI_Request request[4];
        MPI_Status status[4];

        if (down != MPI_PROC_NULL) {
            for (int iz = 0; iz < N; iz++) {
                for (int kx = 0; kx < N; kx++) {
                    ysend_down[iz * N + kx] = phi_o[(endy - 1) * N + N * N * iz + kx];
                }
            }
            MPI_Irecv(yrecv_down, N * N, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, &request[reqCount++]);
            MPI_Isend(ysend_down, N * N, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &request[reqCount++]);
        }

        // Upward neighbor exchange (to 'up')
        if (up != MPI_PROC_NULL) {
            for (int iz = 0; iz < N; iz++) {
                for (int kx = 0; kx < N; kx++) {
                    ysend_up[iz * N + kx] = phi_o[starty * N + N * N * iz + kx];
                }
            }
            MPI_Irecv(yrecv_up, N * N, MPI_DOUBLE, up, 1, MPI_COMM_WORLD, &request[reqCount++]);
            MPI_Isend(ysend_up, N * N, MPI_DOUBLE, up, 0, MPI_COMM_WORLD, &request[reqCount++]);
        }

        // Wait for all communications to finish before using yrecv arrays
        MPI_Waitall(reqCount, request, status);

        // Now copy received data into halo regions
        if (down != MPI_PROC_NULL) {
            for (int iz = 0; iz < N; iz++) {
                for (int kx = 0; kx < N; kx++) {
                    phi_o[endy * N + N * N * iz + kx] = yrecv_down[iz * N + kx];
                }
            }
        }
        if (up != MPI_PROC_NULL) {
            for (int iz = 0; iz < N; iz++) {
                for (int kx = 0; kx < N; kx++) {
                    phi_o[(starty - 1) * N + N * N * iz + kx] = yrecv_up[iz * N + kx];
                }
            }
        }

        // Jacobi iteration
        // Red-black or standard Jacobi approach:
        // Update interior points owned by this rank
        for (int iz = 1; iz < N - 1; iz++) {
            for (int jy = max(1, starty); jy < min(endy, N - 1); jy++) {
                for (int kx = 1; kx < N - 1; kx++) {
                    // Skip boundaries in y since they are fixed
                    // if (jy == 0 || jy == N - 1) continue;
                    if((iz+jy+kx)%2==0){
                    // Jacobi update
                    double new_val = (phi_o[iz * N * N + jy * N + (kx + 1)] +
                                      phi_o[iz * N * N + jy * N + (kx - 1)] +
                                      phi_o[iz * N * N + (jy + 1) * N + kx] +
                                      phi_o[iz * N * N + (jy - 1) * N + kx] +
                                      phi_o[(iz + 1) * N * N + jy * N + kx] +
                                      phi_o[(iz - 1) * N * N + jy * N + kx] -
                                      dx * dx * RHS[iz * N * N + jy * N + kx]) / 6.0;

                    error += fabs(new_val - phi_o[iz * N * N + jy * N + kx]);
                    phi_n[iz * N * N + jy * N + kx] = new_val;
                    }
                }
            }
        }

        for (int iz = 1; iz < N - 1; iz++) {    
            for (int jy = max(1, starty); jy < min(endy, N - 1); jy++) {
                for (int kx = 1; kx < N - 1; kx++) {
                    // Skip boundaries in y since they are fixed
                    // if (jy == 0 || jy == N - 1) continue;
                    if((iz+jy+kx)%2==1){
                    // Jacobi update
                    double new_val = (phi_o[iz * N * N + jy * N + (kx + 1)] +
                                      phi_o[iz * N * N + jy * N + (kx - 1)] +
                                      phi_o[iz * N * N + (jy + 1) * N + kx] +
                                      phi_o[iz * N * N + (jy - 1) * N + kx] +
                                      phi_o[(iz + 1) * N * N + jy * N + kx] +
                                      phi_o[(iz - 1) * N * N + jy * N + kx] -
                                      dx * dx * RHS[iz * N * N + jy * N + kx]) / 6.0;

                    error += fabs(new_val - phi_o[iz * N * N + jy * N + kx]);
                    phi_n[iz * N * N + jy * N + kx] = new_val;
                    }
                }
            }
        }

        // Re-apply boundary conditions at y=0 and y=N-1
        // (Keeps them fixed so that the solution can converge)
        for (int iz = 0; iz < N; iz++) {
            for (int kx = 0; kx < N; kx++) {
                // y=0 boundary
                int i_bottom = iz * N * N + 0 * N + kx;
                phi_n[i_bottom] = RHS[i_bottom] / (-(m*m+n*n+k*k)*Pi*Pi);

                // y=N-1 boundary
                int i_top = iz * N * N + (N - 1)*N + kx;
                phi_n[i_top] = RHS[i_top] / (-(m*m+n*n+k*k)*Pi*Pi);
            }
        }

        // Compute global error
        double local_error = error;
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Swap phi_o and phi_n for next iteration
        // Instead of MPI_Allreduce for solution arrays, we just swap locally.
        double *temp = phi_o;
        phi_o = phi_n;
        phi_n = temp;
    }

    long long elapsed = currentTimeInMicroseconds() - t1;

    if (rank == 0) {
        cout << "MPI with N = " << N << endl;
        cout << "Number of iterations: " << it << endl;
        cout << "Global error = " << global_error << endl;
        cout << "Time in microseconds: " << elapsed << endl;
    }

    // Clean up
   
    //     MPI_Allreduce(phi_n, phi_o, N * N * N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //         for (int i=0;i<N*N*N;i++){
    //             if(((((i-i%N)/N)%N)==0)||((((i-i%N)/N)%N))==(N-1)){
    //             phi_o[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
    //             phi_n[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
    //             }
    //             else {
    //                 phi_o[i]=phi_o[i]/3;
    //             }
    //         }


    // if(rank==3){
    //     for (int i=0;i<N;i++){
    //         // cout<< " rank Id "<< rank <<endl;
    //         cout<<"Z= "<<i*dx<<" "<<endl;
    //         for (int j=0;j<N;j++){
    //             for (int k=0;k<N;k++){
    //                 cout<<phi_o[i*N*N+j*N+k]<<" ";
    //             }
    //             cout<<endl;
    //         }
    //         cout<<endl<<endl;
    //     }
    // }



    delete[] phi_o;
    delete[] phi_n;
    delete[] RHS;
    delete[] ysend_down;
    delete[] yrecv_down;
    delete[] ysend_up;
    delete[] yrecv_up;

    MPI_Finalize();
    return 0;
}
