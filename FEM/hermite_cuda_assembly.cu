

#include <iostream>
#include <iomanip>


__global__ void CSR_Assembly_Hermite(const int num_nodes, const int degree, const double *local_stiffness, const double *local_mass, int *row_offset_stiffness, int *col_indices_stiffness, double *stiffness_values, int *row_offset_mass, int *col_indices_mass, double *mass_values) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int dof_per_node = 2 * degree; // Each node has 2n degrees of freedom for Hermite polynomials of degree n

    if (idx < num_nodes) {
        int start_idx_stiffness = dof_per_node * dof_per_node * idx; // Start index for this element's stiffness entries
        int start_idx_mass = dof_per_node * dof_per_node * idx; // Start index for this element's mass entries

        int nextIdx = (idx + 1) % num_nodes;  // Handle periodic boundary conditions
        int prevIdx = (idx - 1 + num_nodes) % num_nodes;

        // Set the row offsets for the stiffness and mass matrices
        row_offset_stiffness[idx] = start_idx_stiffness;
        row_offset_mass[idx] = start_idx_mass;

        // Each element contributes to a dof_per_node x dof_per_node block in the global matrix
        for (int i = 0; i < dof_per_node; i++) {
            for (int j = 0; j < dof_per_node; j++) {
                int row_idx = start_idx_stiffness + i * dof_per_node + j;
                col_indices_stiffness[row_idx] = prevIdx * dof_per_node + j;
                stiffness_values[row_idx] = local_stiffness[i * dof_per_node + j];

                row_idx = start_idx_mass + i * dof_per_node + j;
                col_indices_mass[row_idx] = prevIdx * dof_per_node + j;
                mass_values[row_idx] = local_mass[i * dof_per_node + j];
            }
        }

        for (int i = 0; i < dof_per_node; i++) {
            for (int j = 0; j < dof_per_node; j++) {
                int row_idx = start_idx_stiffness + dof_per_node * dof_per_node + i * dof_per_node + j;
                col_indices_stiffness[row_idx] = idx * dof_per_node + j;
                stiffness_values[row_idx] = local_stiffness[(dof_per_node + i) * dof_per_node + j];

                row_idx = start_idx_mass + dof_per_node * dof_per_node + i * dof_per_node + j;
                col_indices_mass[row_idx] = idx * dof_per_node + j;
                mass_values[row_idx] = local_mass[(dof_per_node + i) * dof_per_node + j];
            }
        }

        for (int i = 0; i < dof_per_node; i++) {
            for (int j = 0; j < dof_per_node; j++) {
                int row_idx = start_idx_stiffness + 2 * dof_per_node * dof_per_node + i * dof_per_node + j;
                col_indices_stiffness[row_idx] = nextIdx * dof_per_node + j;
                stiffness_values[row_idx] = local_stiffness[(2 * dof_per_node + i) * dof_per_node + j];

                row_idx = start_idx_mass + 2 * dof_per_node * dof_per_node + i * dof_per_node + j;
                col_indices_mass[row_idx] = nextIdx * dof_per_node + j;
                mass_values[row_idx] = local_mass[(2 * dof_per_node + i) * dof_per_node + j];
            }
        }
    }

    // Handle the last element in row_offsets
    if (idx == num_nodes) {
        row_offset_stiffness[num_nodes] = dof_per_node * dof_per_node * num_nodes * 3;
        row_offset_mass[num_nodes] = dof_per_node * dof_per_node * num_nodes * 3;
    }
}


int main() {
    // Parameters
    int num_nodes = 5;  // Example number of nodes
    int degree = 2;      // Degree of Hermite polynomials
    int dof_per_node = 2 * degree;

    // Allocate memory on host
    double *h_local_stiffness = new double[dof_per_node * dof_per_node];
    double *h_local_mass = new double[dof_per_node * dof_per_node];
    int *h_row_offset_stiffness = new int[num_nodes + 1];
    int *h_col_indices_stiffness = new int[num_nodes * dof_per_node * dof_per_node * 3];
    double *h_stiffness_values = new double[num_nodes * dof_per_node * dof_per_node * 3];
    int *h_row_offset_mass = new int[num_nodes + 1];
    int *h_col_indices_mass = new int[num_nodes * dof_per_node * dof_per_node * 3];
    double *h_mass_values = new double[num_nodes * dof_per_node * dof_per_node * 3];

    // Initialize local matrices (example values, should be replaced with actual values)
    for (int i = 0; i < dof_per_node * dof_per_node; i++) {
        h_local_stiffness[i] = 1.0;  // Example stiffness values
        h_local_mass[i] = 1.0;       // Example mass values
    }

    // Allocate memory on device
    double *d_local_stiffness, *d_local_mass, *d_stiffness_values, *d_mass_values;
    int *d_row_offset_stiffness, *d_col_indices_stiffness, *d_row_offset_mass, *d_col_indices_mass;

    cudaMalloc(&d_local_stiffness, dof_per_node * dof_per_node * sizeof(double));
    cudaMalloc(&d_local_mass, dof_per_node * dof_per_node * sizeof(double));
    cudaMalloc(&d_row_offset_stiffness, (num_nodes + 1) * sizeof(int));
    cudaMalloc(&d_col_indices_stiffness, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(int));
    cudaMalloc(&d_stiffness_values, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(double));
    cudaMalloc(&d_row_offset_mass, (num_nodes + 1) * sizeof(int));
    cudaMalloc(&d_col_indices_mass, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(int));
    cudaMalloc(&d_mass_values, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(double));

    // Copy data from host to device
    cudaMemcpy(d_local_stiffness, h_local_stiffness, dof_per_node * dof_per_node * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_local_mass, h_local_mass, dof_per_node * dof_per_node * sizeof(double), cudaMemcpyHostToDevice);

    // Kernel launch
    int blockSize = 256;
    int numBlocks = (num_nodes + blockSize - 1) / blockSize;
    CSR_Assembly_Hermite<<<numBlocks, blockSize>>>(num_nodes, degree, d_local_stiffness, d_local_mass, d_row_offset_stiffness, d_col_indices_stiffness, d_stiffness_values, d_row_offset_mass, d_col_indices_mass, d_mass_values);

    // Copy results back to host
    cudaMemcpy(h_row_offset_stiffness, d_row_offset_stiffness, (num_nodes + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_col_indices_stiffness, d_col_indices_stiffness, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_stiffness_values, d_stiffness_values, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_row_offset_mass, d_row_offset_mass, (num_nodes + 1) * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_col_indices_mass, d_col_indices_mass, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_mass_values, d_mass_values, num_nodes * dof_per_node * dof_per_node * 3 * sizeof(double), cudaMemcpyDeviceToHost);

    // Cleanup
    cudaFree(d_local_stiffness);
    cudaFree(d_local_mass);
    cudaFree(d_row_offset_stiffness);
    cudaFree(d_col_indices_stiffness);
    cudaFree(d_stiffness_values);
    cudaFree(d_row_offset_mass);
    cudaFree(d_col_indices_mass);
    cudaFree(d_mass_values);

    std::cout << "Row Offsets:\n";
    for (int i = 0; i <= num_nodes; ++i) {
        std::cout << h_row_offset_stiffness[i] << " ";
    }
    std::cout << "\n\nColumn Indices:\n";
    //int total_entries = h_row_offset_stiffness[num_nodes];  // The last element of row_offsets gives the total number of entries
    for (int i = 0; i < 3*num_nodes; ++i) {
        std::cout << h_col_indices_stiffness[i] << " ";
    }
    std::cout << "\n\nValues:\n";
    for (int i = 0; i < 3*num_nodes; ++i) {
        std::cout << std::fixed << std::setprecision(2) << h_stiffness_values[i] << " ";
    }
    std::cout << std::endl;

    delete[] h_local_stiffness;
    delete[] h_local_mass;
    delete[] h_row_offset_stiffness;
    delete[] h_col_indices_stiffness;
    delete[] h_stiffness_values;
    delete[] h_row_offset_mass;
    delete[] h_col_indices_mass;
    delete[] h_mass_values;

    return 0;
}