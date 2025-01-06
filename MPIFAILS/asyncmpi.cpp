#include <mpi.h>
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
    int rankID = 0, num_ranks = 1;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);

    //constants
    int n, m, k; // eigen solutions
    n=1;m=1;k=1;
    int N=4; // dimension
    int local_N=N/2;
    double L=1;
    double dx=L/(N-1);
    double Pi=  4*atan(1);
    //vector solutions 
    double ***phi_o =allocate_3d_array(N);
    double ***phi_n=allocate_3d_array(N);
    double ***RHS=allocate_3d_array(N);
    // double ***error=allocate_3d_array(N);
    // soultios

    double error = 20;
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

    // Define neighbor processors (north, south, east, west, in, out)
    int north = (rankID >= 2 && rankID != 4 && rankID != 5) ? rankID - 2 : MPI_PROC_NULL ; // P 2,3,6,7 
    int south = (rankID < 6 && rankID != 2 && rankID != 3) ? rankID + 2 : MPI_PROC_NULL ; // P 0,1,4,5  
    int west = (rankID % 2 == 0) ? MPI_PROC_NULL : rankID - 1; // P 1,3,5,7
    int east = (rankID % 2 == 1) ? MPI_PROC_NULL : rankID + 1; // P 0,2,4,6
    int in = (rankID >= 4) ? MPI_PROC_NULL : rankID + 4; // P 0,1,2,3
    int out = (rankID < 4) ? MPI_PROC_NULL : rankID - 4; // P 4,5,6,7

    int startx = (rankID % 2 == 0) ? 0 : local_N; // columns
    int endx = (rankID % 2 == 0) ? local_N-1 : N;
    int starty = (rankID < 6 && rankID != 2 && rankID != 3) ? 0 : local_N; // rows
    int endy = (rankID < 6 && rankID != 2 && rankID != 3) ? local_N-1 : N;
    int startz = (rankID < 4) ? 0 : local_N; // depth
    int endz = (rankID < 4) ? local_N-1 : N;

    MPI_Barrier(MPI_COMM_WORLD);

    long long t1= currentTimeInMicroseconds();

    while (global_error>0.001){
        error = 0;
        MPI_Request requests[12];
        int req_idx = 0;

        // Communication of boundary values
        double *send_hplane = new double[local_N*local_N]; // N-S
        double *recv_hplane = new double[local_N*local_N];
        double *send_vplane = new double[local_N*local_N]; // E-W
        double *recv_vplane = new double[local_N*local_N];
        double *send_rplane = new double[local_N*local_N]; // I-O
        double *recv_rplane = new double[local_N*local_N];
        // double **send_vplane = new double*[local_N]; // E-W
        // double **recv_vplane = new double*[local_N];
        // double **send_hplane = new double*[local_N]; // N-S
        // double **recv_hplane = new double*[local_N];
        // double **send_rplane = new double*[local_N]; // I-O
        // double **recv_rplane = new double*[local_N];

        // for (int i = 0; i<local_N ; i++) {
        //     send_vplane[i] = new double[local_N];
        //     recv_vplane[i] = new double[local_N];
        //     send_hplane[i] = new double[local_N];
        //     recv_hplane[i] = new double[local_N];
        //     send_rplane[i] = new double[local_N];
        //     recv_rplane[i] = new double[local_N];
        // }

        fprintf(stderr,"start\n");
        // Communication for south boundary - down
        for (int i = startx; i < endx; i++) {
            for (int k = startz; k < endz; k++) {
                send_hplane[i*local_N+k] = phi_o[local_N-1][i][k]; 
                // send_hplane[i][k] = phi_o[local_N-1][i][k]; 
            }
        }
        // MPI_Sendrecv(send_hplane, local_N*local_N, MPI_DOUBLE, south, 0,
        //                  recv_hplane, local_N*local_N, MPI_DOUBLE, north, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_hplane, local_N * local_N, MPI_DOUBLE, south, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_hplane, local_N * local_N, MPI_DOUBLE, north, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (north != MPI_PROC_NULL) {
            for (int i = startx; i < endx; i++) {
                for (int k = startz; k < endz; k++) {
                    phi_o[local_N][i][k] = recv_hplane[i*local_N+k]; 
                    // phi_o[local_N][i][k] = recv_hplane[i][k]; 
                }
            }
        }
        

        // Communication for north boundary - up
        for (int i = startx; i < endx; i++) {
            for (int k = startz; k < endz; k++) {
                send_hplane[i*local_N+k] = phi_o[local_N][i][k]; 
            }
        }
        // MPI_Sendrecv(send_hplane, local_N*local_N, MPI_DOUBLE, north, 0,
        //                  recv_hplane, local_N*local_N, MPI_DOUBLE, south, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_hplane, local_N * local_N, MPI_DOUBLE, north, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_hplane, local_N * local_N, MPI_DOUBLE, south, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (south != MPI_PROC_NULL) {
            for (int i = startx; i < endx; i++) {
                for (int k = startz; k < endz; k++) {
                    phi_o[local_N-1][i][k] = recv_hplane[i*local_N+k]; 
                }
            }
        }


        // Communication for east boundary - right
        for (int i = starty; i < endy; i++) {
            for (int k = startz; k < endz; k++) {
                send_vplane[i*local_N+k] = phi_o[i][local_N-1][k]; 
            }
        }
        // MPI_Sendrecv(send_vplane, local_N*local_N, MPI_DOUBLE, east, 0,
        //                  recv_vplane, local_N*local_N, MPI_DOUBLE, west, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_vplane, local_N * local_N, MPI_DOUBLE, east, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_vplane, local_N * local_N, MPI_DOUBLE, west, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (west != MPI_PROC_NULL) {
            for (int i = starty; i < endy; i++) {
                for (int k = startz; k < endz; k++) {
                    phi_o[i][local_N][k] = recv_vplane[i*local_N+k]; 
                }
            }
        }


        // Communication for west boundary - left
        for (int i = starty; i < endy; i++) {
            for (int k = startz; k < endz; k++) {
                send_vplane[i*local_N+k] = phi_o[i][local_N][k]; 
            }
        }
        // MPI_Sendrecv(send_vplane, local_N*local_N, MPI_DOUBLE, west, 0,
        //                  recv_vplane, local_N*local_N, MPI_DOUBLE, east, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_vplane, local_N * local_N, MPI_DOUBLE, west, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_vplane, local_N * local_N, MPI_DOUBLE, east, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (east != MPI_PROC_NULL) {
            for (int i = starty; i < endy; i++) {
                for (int k = startz; k < endz; k++) {
                    phi_o[i][local_N-1][k] = recv_vplane[i*local_N+k]; 
                }
            }
        }


        // Communication for in boundary 
        for (int i = startx; i < endx; i++) {
            for (int k = starty; k < endy; k++) {
                send_rplane[i*local_N+k] = phi_o[i][k][local_N-1]; 
            }
        }
        // MPI_Sendrecv(send_rplane, local_N*local_N, MPI_DOUBLE, in, 0,
        //                  recv_rplane, local_N*local_N, MPI_DOUBLE, out, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_rplane, local_N * local_N, MPI_DOUBLE, in, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_rplane, local_N * local_N, MPI_DOUBLE, out, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (out != MPI_PROC_NULL) {
            for (int i = startx; i < endx; i++) {
                for (int k = starty; k < endy; k++) {
                    phi_o[i][k][local_N] = recv_rplane[i*local_N+k]; 
                }
            }
        }


        // Communication for out boundary 
        for (int i = startx; i < endx; i++) {
            for (int k = starty; k < endy; k++) {
                send_rplane[i*local_N+k] = phi_o[i][k][local_N]; 
            }
        }
        // MPI_Sendrecv(send_rplane, local_N*local_N, MPI_DOUBLE, out, 0,
        //                  recv_rplane, local_N*local_N, MPI_DOUBLE, in, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Isend(send_rplane, local_N * local_N, MPI_DOUBLE, out, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        MPI_Irecv(recv_rplane, local_N * local_N, MPI_DOUBLE, in, 0, MPI_COMM_WORLD, &requests[req_idx++]);
        if (in != MPI_PROC_NULL) {
            for (int i = startx; i < endx; i++) {
                for (int k = starty; k < endy; k++) {
                    phi_o[i][k][local_N-1] = recv_rplane[i*local_N+k]; 
                }
            }
        }
        // MPI_Barrier(MPI_COMM_WORLD);
        MPI_Waitall(req_idx, requests, MPI_STATUSES_IGNORE);

        delete[] send_hplane;
        delete[] recv_hplane;
        delete[] send_vplane;
        delete[] recv_vplane;
        delete[] send_rplane;
        delete[] recv_rplane;

        fprintf(stderr,"end halo exchange\n");

        // problem with comm
        if (rankID == 0) {
            for(int i=0; i<N;i++){
                cout<<"Z= "<<i*dx<< " it #"<< it<< endl;
                for(int j=0; j<N; j++){
                    for (int k = 0; k < N ; k++) {
                        cout << phi_n[k][j][i] << " " ;
                        // cout<<RHS[i][j][k]<<" ";                
                    }
                    cout << endl ; 
                }
                cout << endl ;
            }
        }

        fprintf(stderr,"start phi computation\n");

        // Update logic for phi based on boundary values and other values

        // for (int i = max(1,startx); i < min(endx,N-1); ++i) {
        //     for (int j = max(1,starty); j < min(endy,N-1) ; ++j) {
        //         for (int k = max(1,startz); k < min(endz,N-1); ++k) {
        for (int i = 1; i < N-1; ++i) {
            for (int j = 1; j < N-1 ; ++j) {
                for (int k = 1; k < N-1; ++k) {
                    if ((i + j + k) % 2 == 0) { // Red points
                        phi_n[i][j][k] = (phi_o[i + 1][j][k] + phi_o[i - 1][j][k] +
                                            phi_o[i][j + 1][k] + phi_o[i][j - 1][k] +
                                            phi_o[i][j][k + 1] + phi_o[i][j][k - 1] -
                                            dx * dx * RHS[i][j][k]) / 6;   
                }
            }
            }
        }
        for (int i = 1; i < N-1; ++i) {
            for (int j = 1; j < N-1 ; ++j) {
                for (int k = 1; k < N-1; ++k) {
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
                    error += fabs(phi_o[i][j][k] - phi_n[i][j][k]);
                    phi_o[i][j][k] = phi_n[i][j][k];
                }
            }
        }
        fprintf(stderr,"done with comp\n");

        global_error = 0.0;
        MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (rankID == 0) cout << "Global error: " << global_error << endl;
        it++;
    }

    MPI_Barrier(MPI_COMM_WORLD);


    if (rankID == 0) {
        long long t2 = currentTimeInMicroseconds()-t1;
        cout << "Baseline with N = " << N << endl;
        cout << "Number of iterations: " << it << endl;
        cout << "Global error = " << global_error << endl;
        cout << "Time in microseconds: " << t2 << endl;
    }

    MPI_Finalize();
    return 0;
}