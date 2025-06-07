import matplotlib.pyplot as plt

# Grid sizes
grid_sizes = [100, 500, 700, 1000]

# Execution times from second image
serial_times = [0.730089, 31.6808, 62.3564, 128.899]
openmp_times = [0.424631, 17.6183, 34.6309, 70.5187]
mpi_times = [0.323383, 8.2036, 16.0857, 32.7119]

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(grid_sizes, serial_times, marker='o', label='Serial', color='red')
plt.plot(grid_sizes, openmp_times, marker='s', label='OpenMP (2 threads)', color='green')
plt.plot(grid_sizes, mpi_times, marker='^', label='MPI (4 processes)', color='blue')

# Labels and title
plt.xlabel('Grid Size')
plt.ylabel('Execution Time (seconds)')
plt.title('Performance Comparison: Serial vs OpenMP (2T) vs MPI (4P)')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Show plot
plt.show()
