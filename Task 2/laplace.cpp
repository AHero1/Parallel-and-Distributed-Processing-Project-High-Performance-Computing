#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#include <omp.h>

const double EPSILON = 0.001;

// Grid sizes and fixed epochs
const int GRID_SIZES[] = {100, 500, 700, 1000, 2000, 4000};
const int MAX_EPOCHS = 5000; // Fixed for all grid sizes
const int NUM_SIZES = 4;

// Allocate 2D array
double** alloc_grid(int rows, int cols) {
    double** grid = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        grid[i] = new double[cols](); // Initialize to zero
    }
    return grid;
}

// Free 2D array
void free_grid(double** grid, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] grid[i];
    }
    delete[] grid;
}

// Initialize grid: top row = 50 (rank 0), bottom row = -50 (last rank or global grid)
void initialize_grid(double** grid, int rows, int cols, int rank = 0, int size = 1) {
    // Initialize all to zero
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            grid[i][j] = 0.0;
        }
    }

    // Set top row for rank 0 (first row of global grid)
    if (rank == 0) {
        for (int j = 0; j < cols; ++j) {
            grid[0][j] = 50.0;
        }
    }

    // Set bottom row for last rank (last row of global grid)
    if (rank == size - 1) {
        int last_row = (size == 1) ? rows - 1 : rows - 2; // Adjust for ghost row in MPI
        for (int j = 0; j < cols; ++j) {
            grid[last_row][j] = -50.0;
        }
    }
}

// Save grid to CSV (for n=100)
void save_grid(double** grid, int rows, int cols, const char* filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error opening " << filename << " for writing\n";
        return;
    }
    out << std::fixed << std::setprecision(3);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            out << grid[i][j];
            if (j < cols - 1) out << ",";
        }
        out << "\n";
    }
    out.close();
}

// Serial Laplace solver
double solve_laplace_serial(double** grid, double** new_grid, int n, int max_epochs) {
    double start_time = MPI_Wtime();
    for (int epoch = 0; epoch < max_epochs; ++epoch) {
        double max_diff = 0.0;

        // Update interior points
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                new_grid[i][j] = 0.25 * (
                    grid[i+1][j] + grid[i-1][j] +
                    grid[i][j+1] + grid[i][j-1]
                );
                double diff = std::fabs(new_grid[i][j] - grid[i][j]);
                max_diff = std::max(max_diff, diff);
            }
        }

        // Copy interior points to grid
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                grid[i][j] = new_grid[i][j];
            }
        }

        if (max_diff < EPSILON) {
            break;
        }
    }
    return MPI_Wtime() - start_time;
}

// OpenMP Laplace solver
double solve_laplace_omp(double** grid, double** new_grid, int n, int max_epochs, int num_threads) {
    double start_time = omp_get_wtime();
    omp_set_num_threads(num_threads);
    for (int epoch = 0; epoch < max_epochs; ++epoch) {
        double global_max_diff = 0.0;

        #pragma omp parallel
        {
            double local_max_diff = 0.0;
            #pragma omp for collapse(2)
            for (int i = 1; i < n - 1; ++i) {
                for (int j = 1; j < n - 1; ++j) {
                    new_grid[i][j] = 0.25 * (
                        grid[i+1][j] + grid[i-1][j] +
                        grid[i][j+1] + grid[i][j-1]
                    );
                    double diff = std::fabs(new_grid[i][j] - grid[i][j]);
                    local_max_diff = std::max(local_max_diff, diff);
                }
            }
            #pragma omp critical
            {
                global_max_diff = std::max(global_max_diff, local_max_diff);
            }
        }

        // Copy interior points to grid
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < n - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                grid[i][j] = new_grid[i][j];
            }
        }

        if (global_max_diff < EPSILON) {
            break;
        }
    }
    return omp_get_wtime() - start_time;
}

// MPI Laplace solver
double solve_laplace_mpi(double** grid, double** new_grid, int local_rows, int n, int max_epochs, int rank, int size, MPI_Comm comm) {
    double start_time = MPI_Wtime();
    for (int epoch = 0; epoch < max_epochs; ++epoch) {
        double local_max_diff = 0.0;

        // Exchange boundary rows (non-blocking)
        MPI_Request reqs[4];
        int req_count = 0;

        if (rank > 0) {
            MPI_Isend(grid[1], n, MPI_DOUBLE, rank - 1, 0, comm, &reqs[req_count++]);
            MPI_Irecv(grid[0], n, MPI_DOUBLE, rank - 1, 0, comm, &reqs[req_count++]);
        }
        if (rank < size - 1) {
            MPI_Isend(grid[local_rows - 2], n, MPI_DOUBLE, rank + 1, 0, comm, &reqs[req_count++]);
            MPI_Irecv(grid[local_rows - 1], n, MPI_DOUBLE, rank + 1, 0, comm, &reqs[req_count++]);
        }
        MPI_Waitall(req_count, reqs, MPI_STATUSES_IGNORE);

        // Update interior points (excluding ghost rows)
        for (int i = 1; i < local_rows - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                new_grid[i][j] = 0.25 * (
                    grid[i+1][j] + grid[i-1][j] +
                    grid[i][j+1] + grid[i][j-1]
                );
                double diff = std::fabs(new_grid[i][j] - grid[i][j]);
                local_max_diff = std::max(local_max_diff, diff);
            }
        }

        // Copy interior points to grid
        for (int i = 1; i < local_rows - 1; ++i) {
            for (int j = 1; j < n - 1; ++j) {
                grid[i][j] = new_grid[i][j];
            }
        }

        // Global max diff for convergence
        double global_max_diff;
        MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, comm);
        if (global_max_diff < EPSILON) {
            break;
        }
    }
    return MPI_Wtime() - start_time;
}

// Gather and save grid for MPI (n=100)
void gather_and_save_mpi(double** grid, int local_rows, int n, int rank, int size, MPI_Comm comm) {
    // Prepare send buffer (excluding ghost rows)
    double* sendbuf = new double[(local_rows - 2) * n];
    for (int i = 1; i < local_rows - 1; ++i) {
        for (int j = 0; j < n; ++j) {
            sendbuf[(i - 1) * n + j] = grid[i][j];
        }
    }

    int* recvcounts = NULL;
    int* displs = NULL;
    double* fullgrid = NULL;

    if (rank == 0) {
        recvcounts = new int[size];
        displs = new int[size];
        int offset = 0;
        for (int i = 0; i < size; ++i) {
            int base_rows = n / size;
            int extra = (i < n % size) ? 1 : 0;
            int rows = base_rows + extra;
            recvcounts[i] = rows * n;
            displs[i] = offset;
            offset += rows * n;
        }
        fullgrid = new double[n * n];
    }

    MPI_Gatherv(sendbuf, (local_rows - 2) * n, MPI_DOUBLE,
                fullgrid, recvcounts, displs, MPI_DOUBLE, 0, comm);

    if (rank == 0) {
        double** grid_2d = alloc_grid(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                grid_2d[i][j] = fullgrid[i * n + j];
            }
        }
        save_grid(grid_2d, n, n, "grid_mpi.csv");
        free_grid(grid_2d, n);
        delete[] fullgrid;
        delete[] recvcounts;
        delete[] displs;
    }

    delete[] sendbuf;
}

int main() {
    MPI_Init(NULL, NULL);
    int rank, size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    // Default mode and threads
    int mode = 0;
    int num_threads = omp_get_max_threads();

    // Prompt for mode on rank 0
    if (rank == 0) {
        std::cout << "Select execution mode (0: Serial, 1: OpenMP, 2: MPI): ";
        std::cin >> mode;
        while (mode < 0 || mode > 2) {
            std::cout << "Invalid mode. Choose 0 (Serial), 1 (OpenMP), or 2 (MPI): ";
            std::cin >> mode;
        }
        if (mode == 1) {
            std::cout << "Enter number of threads for OpenMP: ";
            std::cin >> num_threads;
            while (num_threads < 1) {
                std::cout << "Invalid number of threads. Enter a positive integer: ";
                std::cin >> num_threads;
            }
        }
    }

    // Broadcast mode and num_threads to all ranks
    MPI_Bcast(&mode, 1, MPI_INT, 0, comm);
    MPI_Bcast(&num_threads, 1, MPI_INT, 0, comm);

    // Initialize timings CSV
    if (rank == 0) {
        if (mode == 2) {
            std::ofstream out("timings_mpi.csv");
            out << "GridSize,MaxEpochs,Time\n";
            out.close();
        } else if (mode == 1) {
            std::ofstream out("timings_omp.csv");
            out << "GridSize,MaxEpochs,Threads,Time\n";
            out.close();
        } else {
            std::ofstream out("timings_serial.csv");
            out << "GridSize,MaxEpochs,Time\n";
            out.close();
        }
    }

    // Store timings for comparison
    double timings[NUM_SIZES][3] = {{0.0}};

    for (int s = 0; s < NUM_SIZES; ++s) {
        int n = GRID_SIZES[s];

        if (mode == 2) { // MPI
            // Calculate row distribution
            int base_rows = n / size;
            int extra = (rank < n % size) ? 1 : 0;
            int local_rows = base_rows + extra + 2; // Include ghost rows

            double** grid = alloc_grid(local_rows, n);
            double** new_grid = alloc_grid(local_rows, n);
            initialize_grid(grid, local_rows, n, rank, size);

            double time = solve_laplace_mpi(grid, new_grid, local_rows, n, MAX_EPOCHS, rank, size, comm);

            if (rank == 0) {
                std::ofstream out("timings_mpi.csv", std::ios::app);
                out << n << "," << MAX_EPOCHS << "," << time << "\n";
                out.close();
                timings[s][2] = time;
                std::cout << "MPI GridSize: " << n << ", Time: " << time << " seconds\n";
            }

            if (n == 100) {
                gather_and_save_mpi(grid, local_rows, n, rank, size, comm);
            }

            free_grid(grid, local_rows);
            free_grid(new_grid, local_rows);
        } else if (rank == 0) { // Serial or OpenMP
            if (mode == 0) { // Serial
                double** grid = alloc_grid(n, n);
                double** new_grid = alloc_grid(n, n);
                initialize_grid(grid, n, n);
                double time = solve_laplace_serial(grid, new_grid, n, MAX_EPOCHS);
                std::ofstream out("timings_serial.csv", std::ios::app);
                out << n << "," << MAX_EPOCHS << "," << time << "\n";
                out.close();
                timings[s][0] = time;
                std::cout << "Serial GridSize: " << n << ", Time: " << time << " seconds\n";
                if (n == 100) {
                    save_grid(grid, n, n, "grid_serial.csv");
                }
                free_grid(grid, n);
                free_grid(new_grid, n);
            } else if (mode == 1) { // OpenMP
                double** grid = alloc_grid(n, n);
                double** new_grid = alloc_grid(n, n);
                initialize_grid(grid, n, n);
                double time = solve_laplace_omp(grid, new_grid, n, MAX_EPOCHS, num_threads);
                std::ofstream out("timings_omp.csv", std::ios::app);
                out << n << "," << MAX_EPOCHS << "," << num_threads << "," << time << "\n";
                out.close();
                timings[s][1] = time;
                std::cout << "OpenMP GridSize: " << n << ", Time: " << time << " seconds\n";
                if (n == 100) {
                    save_grid(grid, n, n, "grid_omp.csv");
                }
                free_grid(grid, n);
                free_grid(new_grid, n);
            }
        }
    }

    // Broadcast timings to all ranks for output
    MPI_Bcast(timings, NUM_SIZES * 3, MPI_DOUBLE, 0, comm);
    MPI_Finalize();
    return 0;
}

