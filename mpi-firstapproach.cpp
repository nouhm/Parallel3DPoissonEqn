#include <mpi.h>
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
    int rankID = 0, num_ranks = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);

    //constants
    int n, m, k; // eigen solutions
    n=1;m=1;k=1;
    int N=50; // dimension
    int local_N=N/2;
    double L=1;
    double dx=L/(N-1);
    double Pi=  4*atan(1);

    double *phi_o = new double [N*N*N];
    double *phi_n = new double [N*N*N];
    double *RHS = new double [N*N*N];

    double error = 20;
    double global_error = 20;
    int it = 1;

    for (int x=0;x<N;x++){
        for (int y=0;y<N;y++){
            for (int z=0; z<N;z++) {
                phi_o[x*N*N + y*N + z]=0;
                phi_n[x*N*N + y*N + z]=0;
                RHS[x*N*N + y*N + z] = -(n*n+m*m+k*k)*(Pi*Pi)*sin(n*Pi*x*dx)*cos(m*Pi*y*dx)*sin(k*Pi*z*dx);
                if ((y==0)||(y==N-1)){ 
                    phi_n[x*N*N + y*N + z] = RHS[x*N*N + y*N + z]/(-(m*m+n*n+k*k)*Pi*Pi); 
                    phi_o[x*N*N + y*N + z] = RHS[x*N*N + y*N + z]/(-(m*m+n*n+k*k)*Pi*Pi);                         
                }
            }
        }   
    }

    MPI_Bcast(RHS, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi_o, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(phi_n, N * N * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Define neighbor processors (north, south, east, west, in, out)
    int north = (rankID >= 2 && rankID != 4 && rankID != 5) ? rankID - 2 : MPI_PROC_NULL ; // P 2,3,6,7 
    int south = (rankID < 6 && rankID != 2 && rankID != 3) ? rankID + 2 : MPI_PROC_NULL ; // P 0,1,4,5  
    int west = (rankID % 2 == 0) ? MPI_PROC_NULL : rankID - 1; // P 1,3,5,7
    int east = (rankID % 2 == 1) ? MPI_PROC_NULL : rankID + 1; // P 0,2,4,6
    int in = (rankID >= 4) ? MPI_PROC_NULL : rankID + 4; // P 0,1,2,3
    int out = (rankID < 4) ? MPI_PROC_NULL : rankID - 4; // P 4,5,6,7

    int startx = (rankID < 6 && rankID != 2 && rankID != 3) ? 0 : local_N; // columns
    int endx = (rankID < 6 && rankID != 2 && rankID != 3) ? local_N : N;
    int starty = (rankID % 2 == 1) ? 0 : local_N; // rows
    int endy = (rankID % 2 == 1) ? local_N : N;
    int startz = (rankID < 4) ? 0 : local_N; // depth
    int endz = (rankID < 4) ? local_N : N;

    MPI_Barrier(MPI_COMM_WORLD);

    long long t1= currentTimeInMicroseconds();

    while (global_error>0.0000001){

        error = 0;

        // Communication of boundary values
        double *send_hplane_1 = new double[local_N*local_N]; // S
        double *recv_hplane_1 = new double[local_N*local_N]; // N
        double *send_vplane_1 = new double[local_N*local_N]; // E
        double *recv_vplane_1 = new double[local_N*local_N]; // W
        double *send_rplane_1 = new double[local_N*local_N]; // I
        double *recv_rplane_1 = new double[local_N*local_N]; // O

        double *send_hplane_2 = new double[local_N*local_N]; // N
        double *recv_hplane_2 = new double[local_N*local_N]; // S
        double *send_vplane_2 = new double[local_N*local_N]; // W
        double *recv_vplane_2 = new double[local_N*local_N]; // E
        double *send_rplane_2 = new double[local_N*local_N]; // O
        double *recv_rplane_2 = new double[local_N*local_N]; // I

        for (int i = 0; i<local_N*local_N ; i++) {
            send_vplane_1[i] = 0;
            recv_vplane_1[i] = 0;
            send_hplane_1[i] = 0;
            recv_hplane_1[i] = 0;
            send_rplane_1[i] = 0;
            recv_rplane_1[i] = 0;

            send_vplane_2[i] = 0;
            recv_vplane_2[i] = 0;
            send_hplane_2[i] = 0;
            recv_hplane_2[i] = 0;
            send_rplane_2[i] = 0;
            recv_rplane_2[i] = 0;
        }

        // Communication for south boundary - down
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_hplane_1[i*local_N+k] = phi_o[(local_N-1)*N*N + (starty+i)*N + startz+k]; 
            }
        }
        MPI_Sendrecv(send_hplane_1, local_N*local_N, MPI_DOUBLE, south, 0,
                         recv_hplane_1, local_N*local_N, MPI_DOUBLE, north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (north != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(local_N)*N*N + (starty+i)*N + startz+k] = recv_hplane_1[i*local_N+k];  
                }
            }
        }

        // Communication for north boundary - up
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_hplane_2[i*local_N+k] = phi_o[(local_N)*N*N + (starty+i)*N + startz+k]; 
            }
        }
        MPI_Sendrecv(send_hplane_2, local_N*local_N, MPI_DOUBLE, north, 1,
                         recv_hplane_2, local_N*local_N, MPI_DOUBLE, south, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (south != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(local_N-1)*N*N + (starty+i)*N + startz+k] = recv_hplane_2[i*local_N+k];  
                }
            }
        }

        // Communication for east boundary - right
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_vplane_1[i*local_N+k] = phi_o[(i+startx)*N*N + (local_N-1)*N + k+startz]; 
            }
        }
        MPI_Sendrecv(send_vplane_1, local_N*local_N, MPI_DOUBLE, east, 2,
                         recv_vplane_1, local_N*local_N, MPI_DOUBLE, west, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (west != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(i+startx)*N*N + local_N*N + k+startz] = recv_vplane_1[i*local_N+k]; 
                }
            }
        }


        // Communication for west boundary - left
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_vplane_2[i*local_N+k] = phi_o[(i+startx)*N*N + (local_N)*N + k+startz]; 
            }
        }
        MPI_Sendrecv(send_vplane_2, local_N*local_N, MPI_DOUBLE, west, 3,
                         recv_vplane_2, local_N*local_N, MPI_DOUBLE, east, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (east != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(i+startx)*N*N + (local_N-1)*N + k+startz] = recv_vplane_2[i*local_N+k]; 
                }
            }
        }


        // Communication for in boundary 
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_rplane_1[i*local_N+k] = phi_o[(i+startx)*N*N + (k+starty)*N + local_N-1]; 
            }
        }
        MPI_Sendrecv(send_rplane_1, local_N*local_N, MPI_DOUBLE, in, 4,
                         recv_rplane_1, local_N*local_N, MPI_DOUBLE, out, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (out != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(i+startx)*N*N + (k+starty)*N + local_N] = recv_rplane_1[i*local_N+k]; 
                }
            }
        }


        // Communication for out boundary 
        for (int i = 0; i < local_N; i++) {
            for (int k = 0; k < local_N; k++) {
                send_rplane_2[i*local_N+k] = phi_o[(i+startx)*N*N + (k+starty)*N + local_N]; 
            }
        }
        MPI_Sendrecv(send_rplane_2, local_N*local_N, MPI_DOUBLE, out, 5,
                         recv_rplane_2, local_N*local_N, MPI_DOUBLE, in, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (in != MPI_PROC_NULL) {
            for (int i = 0; i < local_N; i++) {
                for (int k = 0; k < local_N; k++) {
                    phi_o[(i+startx)*N*N + (k+starty)*N + local_N-1] = recv_rplane_2[i*local_N+k]; 
                }
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        delete[] send_hplane_1;
        delete[] recv_hplane_1;
        delete[] send_vplane_1;
        delete[] recv_vplane_1;
        delete[] send_rplane_1;
        delete[] recv_rplane_1;
        delete[] send_hplane_2;
        delete[] recv_hplane_2;
        delete[] send_vplane_2;
        delete[] recv_vplane_2;
        delete[] send_rplane_2;
        delete[] recv_rplane_2;

        for (int i = max(1,startx); i < min(endx,N-1); ++i) {
            for (int j = max(1,starty); j < min(endy,N-1) ; ++j) {
                for (int k = max(1,startz); k < min(endz,N-1); ++k) {
                    if ((i + j + k) % 2 == 0) { // Red points
                        phi_n[i*N*N + j*N + k] = (phi_o[(i + 1)*N*N + j*N + k] + phi_o[(i - 1)*N*N + j*N + k] +
                                                phi_o[i*N*N + (j + 1)*N + k] + phi_o[i*N*N + (j - 1)*N + k] +
                                                phi_o[i*N*N + j*N + k + 1] + phi_o[i*N*N + j*N + k - 1] -
                                                dx * dx * RHS[i*N*N + j*N + k]) / 6;   
                    }
                }
            }
        }
        for (int i = max(1,startx); i < min(endx,N-1); ++i) {
            for (int j = max(1,starty); j < min(endy,N-1) ; ++j) {
                for (int k = max(1,startz); k < min(endz,N-1); ++k) {
                    if ((i + j + k) % 2 == 1) { // Black points
                        phi_n[i*N*N + j*N + k] = (phi_o[(i + 1)*N*N + j*N + k] + phi_o[(i - 1)*N*N + j*N + k] +
                                                phi_o[i*N*N + (j + 1)*N + k] + phi_o[i*N*N + (j - 1)*N + k] +
                                                phi_o[i*N*N + j*N + k + 1] + phi_o[i*N*N + j*N + k - 1] -
                                                dx * dx * RHS[i*N*N + j*N + k]) / 6;   
                    }
                }
            }
        }

        for (int i = max(1,startx); i < min(endx,N-1); ++i) {
            for (int j = max(1,starty); j < min(endy,N-1) ; ++j) {
                for (int k = max(1,startz); k < min(endz,N-1); ++k) {
                    error += fabs(phi_o[i*N*N + j*N + k] - phi_n[i*N*N + j*N + k]);
                }
            }
        }

        global_error = 0.0;
        MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        it++;

        double *temp = phi_o;
        phi_o = phi_n;
        phi_n = temp;    
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allreduce(phi_o, phi_n, N * N * N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i=0;i<N*N*N;i++)
        // F is not zero at y boundaries
        if(((((i-i%N)/N)%N)==0)||((((i-i%N)/N)%N))==(N-1)){
            phi_o[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
            phi_n[i]=RHS[i]/(-(m*m+n*n+k*k)*Pi*Pi);
        }

    if (rankID == 0) {

        // for(int i=0; i<N;i++){
        //     cout<<"Z= "<<i*dx<< " it #"<< it<< endl;
        //     for(int j=0; j<N; j++){
        //         for (int k = 0; k < N ; k++) {
        //             cout << phi_n[k*N*N + j*N + i] << " " ;
        //         }
        //         cout << endl ; 
        //     }
        //     cout << endl ;
        // }

        long long t2 = currentTimeInMicroseconds()-t1;
        cout << "MPI with N = " << N << endl;
        cout << "Number of iterations: " << it << endl;
        cout << "Global error = " << global_error << endl;
        cout << "Time in microseconds: " << t2 << endl;
        cout << global_error << endl;
    }

    MPI_Finalize();
    return 0;
}