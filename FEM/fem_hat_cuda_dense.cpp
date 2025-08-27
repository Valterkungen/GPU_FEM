#include <iostream>
#include <cmath>
#include <vector>

void setupConnectivity(int numElements, std::vector<int>& connectivity) {
    for (int i = 0; i < numElements; ++i) {
        int node1 = i;               // current node
        int node2 = (i + 1) % numElements; // next node, wraps around
        connectivity[2 * i] = node1;
        connectivity[2 * i + 1] = node2;
    }
}

void setMatrices(int matrix_size, std::vector<float>& matrix) {
    for (int i = 0; i < matrix_size; ++i) {
        matrix[i] = 0.f;
    }
}

void assembleStiffnessMatrix(int numElements, int numNodes, float h, const std::vector<float>& localMass, const std::vector<float>& localStiffness, const std::vector<int>& connectivity, std::vector<float>& globalMass, std::vector<float>& globalStiffness) {
    for (int elemIdx = 0; elemIdx < numElements; ++elemIdx) {
        int node1 = connectivity[2 * elemIdx];     // Node i
        int node2 = connectivity[2 * elemIdx + 1]; // Node j

        // Add local stiffness to global stiffness
        globalStiffness[node1 * numNodes + node1] += (1.f/h) * localStiffness[0];
        globalStiffness[node1 * numNodes + node2] += (1.f/h) * localStiffness[1];
        globalStiffness[node2 * numNodes + node1] += (1.f/h) * localStiffness[2];
        globalStiffness[node2 * numNodes + node2] += (1.f/h) * localStiffness[3];

        // Add local mass to global mass
        globalMass[node1 * numNodes + node1] += (1.f/h) * localMass[0];
        globalMass[node1 * numNodes + node2] += (1.f/h) * localMass[1];
        globalMass[node2 * numNodes + node1] += (1.f/h) * localMass[2];
        globalMass[node2 * numNodes + node2] += (1.f/h) * localMass[3];
    }
}

int main() {
    int numNodes = 5;
    int numElements = numNodes;  
    float h = 1.f / numNodes;

    std::vector<float> localMass = {1.0f/3, 1.0f/6, 1.0f/6, 1.0f/3};
    std::vector<float> localStiffness = {1, -1, -1, 1};
    std::vector<float> mass_matrix(numNodes * numNodes, 0);
    std::vector<float> stiff_matrix(numNodes * numNodes, 0);
    std::vector<int> connectivity(2 * numElements);

    setupConnectivity(numElements, connectivity);
    setMatrices(numNodes * numNodes, mass_matrix);
    setMatrices(numNodes * numNodes, stiff_matrix);

    assembleStiffnessMatrix(numElements, numNodes, h, localMass, localStiffness, connectivity, mass_matrix, stiff_matrix);

    std::cout << "Stiffness Matrix:" << std::endl;
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            std::cout << stiff_matrix[i * numNodes + j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Mass Matrix:" << std::endl;
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            std::cout << mass_matrix[i * numNodes + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
