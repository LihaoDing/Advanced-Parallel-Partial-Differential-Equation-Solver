#define _USE_MATH_DEFINES

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <mpi.h>

using namespace std;

//Define a matrix to map a 2D array to a 1D array and it will use the contiguous memory to increase the performance.
class Cmatrix
{
public:
	double* mat_1D;
	double** mat_2D;
	int n, m;

	Cmatrix() {
		mat_1D = nullptr;
		mat_2D = nullptr;
	}

	Cmatrix(int imax, int jmax) {
		n = imax;
		m = jmax;
		mat_1D = new double[n * m];
		mat_2D = new double* [n];
		for (int i = 0; i < n; i++)
			mat_2D[i] = &mat_1D[i * m];
	}

	~Cmatrix() {
		delete[] mat_1D;
		delete[] mat_2D;
	}
};


// Decomposition: Function to calculate the dimensions of the local domain for each processor
void calculate_domain_dimensions(int p, int &rows, int &cols, int id, int imax, int jmax, double x_max, double y_max,
                                 Cmatrix*& local_C1, Cmatrix*& local_C2, Cmatrix*& local_C1_old, Cmatrix*& local_C2_old, int &start_i, int &start_j, int &end_i, int &end_j, int bc) {
    // Determine the optimal row and column count to maintain as square a grid as possible
    int sqrt_p = (int)std::sqrt(p);
    while (p % sqrt_p != 0) {
        sqrt_p--;
    }
    cols = sqrt_p;  // Prefer more columns than rows to spread from left to right
    rows = p / sqrt_p;

    // Calculate local dimensions based on process ID ordered from left to right, top to bottom
    int my_row = id / cols;
    int my_col = id % cols;
    // Determine the range of rows and columns for this process
    start_i = my_row * (imax / rows) + (my_row < imax % rows ? my_row : imax % rows);
    start_j = my_col * (jmax / cols) + (my_col < jmax % cols ? my_col : jmax % cols);
    end_i = (my_row + 1) * (imax / rows) + (my_row + 1 < imax % rows ? my_row + 1 : imax % rows);
    end_j = (my_col + 1) * (jmax / cols) + (my_col + 1 < jmax % cols ? my_col + 1 : jmax % cols);

    int local_imax, local_jmax;
    
    if (bc == 0) {
        // Basic division of work
        local_imax = imax / rows + (my_row < imax % rows ? 1 : 0);
        local_jmax = jmax / cols + (my_col < jmax % cols ? 1 : 0);

        // Adjust for ghost cells, not added at the global domain's boundaries
        if (my_row > 0) {// top ghost cell
            local_imax += 1;
            start_i -= 1;
        } 
        if (my_row < rows - 1) {// bottom ghost cell
            local_imax += 1;
            end_i += 1;
        }
        if (my_col > 0) {// left ghost cell
            local_jmax += 1;
            start_j -= 1;
        }
        if (my_col < cols - 1) {// right ghost cell
            local_jmax += 1;
            end_j += 1;
        }

    } else {
        // Calculate local dimensions with ghost cells
        local_imax = imax / rows + (my_row < imax % rows ? 1 : 0) + 2; // Adding two ghost cells, one on each side
        local_jmax = jmax / cols + (my_col < jmax % cols ? 1 : 0) + 2;

        start_i -= 1;
        start_j -= 1;
        end_i += 1;
        end_j += 1;
    }

    // Allocate local matrices with conditional ghost cells
    local_C1 = new Cmatrix(local_imax, local_jmax);
    local_C2 = new Cmatrix(local_imax, local_jmax);
    local_C1_old = new Cmatrix(local_imax, local_jmax);
    local_C2_old = new Cmatrix(local_imax, local_jmax);
}


// set up simulation parameters
int imax = 301, jmax = 301;
double t_max = 30.0;
double t, t_out = 0.0, dt_out = 0.1, dt;
double y_max = 30.0, x_max = 30.0, dx, dy;

//set up simulation constants
const double f = 2.0, q = 0.002, epsilon = 0.03, D1 = 1.0, D2 = 0.6;


// Function for checking the start and end indexes for the local domain
void check_start_end(int &start, int &end, int rows, int cols, int pos, int id, int imax, int jmax, int mode) {
    int my_row = id / cols;
    int my_col = id % cols;
    if (pos == 0 || pos == 2) { //represent iterate the top or bottom row
        if (my_col == 0) {
            start = 0;
        } else {
            start = 1;
        }
        if (my_col == cols - 1) {
            end = jmax;
        } else {
            end = jmax - 1;
        }
        if (mode == 1) {
            if (my_col == 0) {
                start += 1;
            }
            if (my_col == cols - 1) {
                end -= 1;
            }
        }
    }
    if (pos == 1 || pos == 3) { //represent iterate the left or right column
        if (my_row == 0) {
            start = 0;
        } else {
            start = 1;
        }
        if (my_row == rows - 1) {
            end = imax;
        } else {
            end = imax - 1;
        }
        if (mode == 1) {
            if (my_row == 0) {
                start += 1;
            }
            if (my_row == rows - 1) {
                end -= 1;
            }
        }
    }
}


// Function to check the start and end indexes for the local domain with periodic boundaries
void check_special_start_end(int &start, int &end, int rows, int cols, int pos, int id, int imax, int jmax) {
    int my_row = id / cols;
    int my_col = id % cols;
    if (pos == 0 || pos == 2) {
        if (my_col == 0) {
            start = 1;
        } else {
            start = 0;
        }
        if (my_col == cols - 1) {
            end = jmax - 1;
        } else {
            end = jmax;
        }
    }
    if (pos == 1 || pos == 3) {
        if (my_row == 0) {
            start = 1;
        } else {
            start = 0;
        }
        if (my_row == rows - 1) {
            end = imax - 1;
        } else {
            end = imax;
        }
    }
}


// Function to initialize the local matrices for each processor
void initialize_local_matrices(int id, int rows, int cols, int imax, int jmax, double x_max, double y_max,
                               Cmatrix*& local_C1, Cmatrix*& local_C2, int start_i, int start_j, int end_i, int end_j, int bc) {
    // Calculate local indexes and dimensions
    int my_row = id / cols;
    int my_col = id % cols;

    int local_imax, local_jmax;

    if (bc == 0) {
        // Start indexes without adding ghost cells directly
        local_imax = imax / rows + (my_row < imax % rows ? 1 : 0);
        local_jmax = jmax / cols + (my_col < jmax % cols ? 1 : 0);

        // Adjust local_imax and local_jmax for ghost cells where necessary
        if (my_row > 0) local_imax++; // Top ghost cell
        if (my_row < rows - 1) local_imax++; // Bottom ghost cell
        if (my_col > 0) local_jmax++; // Left ghost cell
        if (my_col < cols - 1) local_jmax++; // Right ghost cell
    } else {
        local_imax = imax / rows + (my_row < imax % rows ? 1 : 0) + 2;  // Adding two ghost cells, one on each side
        local_jmax = jmax / cols + (my_col < jmax % cols ? 1 : 0) + 2;  // Same for columns
    }
    

    double origin_x = x_max / 2.0, origin_y = y_max / 2.0;
    double dx = x_max / (double)(imax - 1);
    double dy = y_max / (double)(jmax - 1);

    // Initialize matrices
    for (int i = 0; i < local_imax; i++) {
        for (int j = 0; j < local_jmax; j++) {
            // Calculate global indexes, accounting for ghost cells
            int global_i = start_i + i;
            int global_j = start_j + j;

            if (global_i < imax && global_j < jmax) {
                double x = global_i * dx, y = global_j * dy;
                double angle = atan2(y - origin_y, x - origin_x);

                if (angle > 0.0 && angle < 0.5)
                    local_C1->mat_2D[i][j] = 0.8;
                else
                    local_C1->mat_2D[i][j] = q * (f + 1) / (f - 1);

                local_C2->mat_2D[i][j] = q * (f + 1) / (f - 1) + angle / (8 * M_PI * f);
            } else {
                // Initialize ghost cells or external boundary cells with some default or boundary-specific values if needed
                local_C1->mat_2D[i][j] = 0; // Example default value
                local_C2->mat_2D[i][j] = 0; // Example default value
            }
        }
    }

}


// Peer to Peer Communication: Function to update the ghost cells of the local domain
void update_ghost_cells_nonblocking(int id, int rows, int cols, int local_imax, int local_jmax, 
                                    Cmatrix* local_C1, Cmatrix* local_C2, int bc) {
    MPI_Status status[16];
    MPI_Request request[16];
    int request_count = 0;
    int my_row = id / cols;
    int my_col = id % cols;
    int start, end;

    // Define communication tags to prevent conflicts
    enum { TAG_TOP_DOWN, TAG_BOTTOM_UP, TAG_LEFT_RIGHT, TAG_RIGHT_LEFT };

    // Sending and receiving updates to/from neighbors
    // Top to bottom (send bottom row, receive into bottom ghost row)
    if (my_row < rows - 1) {
        check_start_end(start, end, rows, cols, 0, id, local_imax, local_jmax, bc);
        // cout << "---processor " << id << " sending to " << id - cols << " start: " << start << " end: " << end << "Top to bottom" << " end - start: " << end - start << endl;
        MPI_Isend(&local_C1->mat_2D[local_imax - 2][start], end - start, MPI_DOUBLE, id + cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Isend(&local_C2->mat_2D[local_imax - 2][start], end - start, MPI_DOUBLE, id + cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Irecv(&local_C1->mat_2D[local_imax - 1][start], end - start, MPI_DOUBLE, id + cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Irecv(&local_C2->mat_2D[local_imax - 1][start], end - start, MPI_DOUBLE, id + cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
    }

    // Bottom to top (send top row, receive into top ghost row)
    if (my_row > 0) {
        check_start_end(start, end, rows, cols, 2, id, local_imax, local_jmax, bc);
        // cout << "---processor " << id << " sending to " << id + cols << " start: " << start << " end: " << end << "Bottom to top" << " end - start: " << end - start << endl;
        MPI_Isend(&local_C1->mat_2D[1][start], end - start, MPI_DOUBLE, id - cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Isend(&local_C2->mat_2D[1][start], end - start, MPI_DOUBLE, id - cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Irecv(&local_C1->mat_2D[0][start], end - start, MPI_DOUBLE, id - cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Irecv(&local_C2->mat_2D[0][start], end - start, MPI_DOUBLE, id - cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
    }
    if (my_col < cols - 1) {
        check_start_end(start, end, rows, cols, 1, id, local_imax, local_jmax, bc);
        
        // Create a structure type to represent a subsection of a column
        int count = end - start; // Number of blocks
        std::vector<int> blocklengths(count, 1); // All blocks are length 1
        std::vector<MPI_Aint> displacements(count);
        MPI_Aint base_address;
        MPI_Get_address(&local_C1->mat_2D[start][1], &base_address);
        for (int k = 0; k < count; k++) {
            MPI_Aint address;
            MPI_Get_address(&local_C1->mat_2D[start + k][1], &address);
            displacements[k] = address - base_address;
        }
        MPI_Datatype column_subsection_type;
        MPI_Type_create_struct(count, blocklengths.data(), displacements.data(), std::vector<MPI_Datatype>(count, MPI_DOUBLE).data(), &column_subsection_type);
        MPI_Type_commit(&column_subsection_type);

        MPI_Isend(&local_C1->mat_2D[start][local_jmax - 2], 1, column_subsection_type, id + 1, TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Isend(&local_C2->mat_2D[start][local_jmax - 2], 1, column_subsection_type, id + 1, TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Irecv(&local_C1->mat_2D[start][local_jmax - 1], 1, column_subsection_type, id + 1, TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Irecv(&local_C2->mat_2D[start][local_jmax - 1], 1, column_subsection_type, id + 1, TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Type_free(&column_subsection_type);
    }

    // Right to left (send leftmost column, receive into left ghost column)
    if (my_col > 0) {
        check_start_end(start, end, rows, cols, 3, id, local_imax, local_jmax, bc);
        
        // Create a structure type to represent a subsection of a column
        int count = end - start; // Number of blocks
        std::vector<int> blocklengths(count, 1); // All blocks are length 1
        std::vector<MPI_Aint> displacements(count);
        MPI_Aint base_address;
        MPI_Get_address(&local_C1->mat_2D[start][1], &base_address);
        for (int k = 0; k < count; k++) {
            MPI_Aint address;
            MPI_Get_address(&local_C1->mat_2D[start + k][1], &address);
            displacements[k] = address - base_address;
        }
        MPI_Datatype column_subsection_type;
        MPI_Type_create_struct(count, blocklengths.data(), displacements.data(), std::vector<MPI_Datatype>(count, MPI_DOUBLE).data(), &column_subsection_type);
        MPI_Type_commit(&column_subsection_type);

        MPI_Isend(&local_C1->mat_2D[start][1], 1, column_subsection_type, id - 1, TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Isend(&local_C2->mat_2D[start][1], 1, column_subsection_type, id - 1, TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Irecv(&local_C1->mat_2D[start][0], 1, column_subsection_type, id - 1, TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);
        MPI_Irecv(&local_C2->mat_2D[start][0], 1, column_subsection_type, id - 1, TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);

        MPI_Type_free(&column_subsection_type);
    }

    if (bc == 1) {
        if (my_row == 0) {  // Topmost row should connect to the bottommost row
            check_special_start_end(start, end, rows, cols, 2, id, local_imax, local_jmax);
            MPI_Isend(&local_C1->mat_2D[1][start], end - start, MPI_DOUBLE, id + (rows - 1) * cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Isend(&local_C2->mat_2D[1][start], end - start, MPI_DOUBLE, id + (rows - 1) * cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);

            MPI_Irecv(&local_C1->mat_2D[0][start], end - start, MPI_DOUBLE, id + (rows - 1) * cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Irecv(&local_C2->mat_2D[0][start], end - start, MPI_DOUBLE, id + (rows - 1) * cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
        }

        if (my_row == rows - 1) {  // Bottommost row should connect to the topmost row
            check_special_start_end(start, end, rows, cols, 0, id, local_imax, local_jmax);
            MPI_Isend(&local_C1->mat_2D[local_imax - 2][start], end - start, MPI_DOUBLE, id - (rows - 1) * cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Isend(&local_C2->mat_2D[local_imax - 2][start], end - start, MPI_DOUBLE, id - (rows - 1) * cols, TAG_TOP_DOWN, MPI_COMM_WORLD, &request[request_count++]);

            MPI_Irecv(&local_C1->mat_2D[local_imax - 1][start], end - start, MPI_DOUBLE, id - (rows - 1) * cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Irecv(&local_C2->mat_2D[local_imax - 1][start], end - start, MPI_DOUBLE, id - (rows - 1) * cols, TAG_BOTTOM_UP, MPI_COMM_WORLD, &request[request_count++]);
        }

        // Leftmost and rightmost columns for periodic boundaries
        if (my_col == 0) {  // Leftmost column should connect to the rightmost column
            check_special_start_end(start, end, rows, cols, 3, id, local_imax, local_jmax);
            
            // Create a structure type to represent a subsection of a column
            int count = end - start; // Number of blocks
            std::vector<int> blocklengths(count, 1); // All blocks are length 1
            std::vector<MPI_Aint> displacements(count);
            MPI_Aint base_address;
            MPI_Get_address(&local_C1->mat_2D[start][1], &base_address);
            for (int k = 0; k < count; k++) {
                MPI_Aint address;
                MPI_Get_address(&local_C1->mat_2D[start + k][1], &address);
                displacements[k] = address - base_address;
            }
            MPI_Datatype column_type;
            MPI_Type_create_struct(count, blocklengths.data(), displacements.data(), std::vector<MPI_Datatype>(count, MPI_DOUBLE).data(), &column_type);
            MPI_Type_commit(&column_type);
            
            MPI_Isend(&local_C1->mat_2D[start][1], 1, column_type, id + (cols - 1), TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Isend(&local_C2->mat_2D[start][1], 1, column_type, id + (cols - 1), TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);

            MPI_Irecv(&local_C1->mat_2D[start][0], 1, column_type, id + (cols - 1), TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Irecv(&local_C2->mat_2D[start][0], 1, column_type, id + (cols - 1), TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Type_free(&column_type);
        }

        if (my_col == cols - 1) {  // Rightmost column should connect to the leftmost column
            check_special_start_end(start, end, rows, cols, 1, id, local_imax, local_jmax);

            // Create a structure type to represent a subsection of a column
            int count = end - start; // Number of blocks
            std::vector<int> blocklengths(count, 1); // All blocks are length 1
            std::vector<MPI_Aint> displacements(count);
            MPI_Aint base_address;
            MPI_Get_address(&local_C1->mat_2D[start][1], &base_address);
            for (int k = 0; k < count; k++) {
                MPI_Aint address;
                MPI_Get_address(&local_C1->mat_2D[start + k][1], &address);
                displacements[k] = address - base_address;
            }
            MPI_Datatype column_type;
            MPI_Type_create_struct(count, blocklengths.data(), displacements.data(), std::vector<MPI_Datatype>(count, MPI_DOUBLE).data(), &column_type);
            MPI_Type_commit(&column_type);

            MPI_Isend(&local_C1->mat_2D[start][local_jmax - 2], 1, column_type, id - (cols - 1), TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Isend(&local_C2->mat_2D[start][local_jmax - 2], 1, column_type, id - (cols - 1), TAG_LEFT_RIGHT, MPI_COMM_WORLD, &request[request_count++]);

            MPI_Irecv(&local_C1->mat_2D[start][local_jmax - 1], 1, column_type, id - (cols - 1), TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Irecv(&local_C2->mat_2D[start][local_jmax - 1], 1, column_type, id - (cols - 1), TAG_RIGHT_LEFT, MPI_COMM_WORLD, &request[request_count++]);
            MPI_Type_free(&column_type);
        }
    }

    // Wait for all non-blocking operations to complete
    MPI_Waitall(request_count, request, status);
}


// Iteration function for local domain with Neumann boundary conditions
void do_local_iteration_neumann(int id, int rows, int cols, int imax, int jmax, int local_imax, int local_jmax,
                                Cmatrix*& local_C1, Cmatrix*& local_C1_old, Cmatrix*& local_C2, Cmatrix*& local_C2_old,
                                double dx, double dy, double dt, double D1, double D2, double f, double q, double epsilon, int start_i, int start_j, int end_i, int end_j) {
    // Swap pointers for local matrices
    Cmatrix* temp = local_C1;
    local_C1 = local_C1_old;
    local_C1_old = temp;
    temp = local_C2;
    local_C2 = local_C2_old;
    local_C2_old = temp;

    // Calculate new concentrations for the internal area of the local domain
    for (int i = 0; i < local_imax; i++) {
        for (int j = 0; j < local_jmax; j++) {
            int boundaryAdjust = (i == 0 || i == local_imax - 1 || j == 0 || j == local_jmax - 1) ? 0 : 1;
            if (boundaryAdjust) {
                local_C1->mat_2D[i][j] = local_C1_old->mat_2D[i][j] + dt * (
                    (local_C1_old->mat_2D[i][j] * (1.0 - local_C1_old->mat_2D[i][j]) - f * local_C2_old->mat_2D[i][j] * (local_C1_old->mat_2D[i][j] - q) / (local_C1_old->mat_2D[i][j] + q)) / epsilon
                    + D1 * ((local_C1_old->mat_2D[i + 1][j] + local_C1_old->mat_2D[i - 1][j] - 2.0 * local_C1_old->mat_2D[i][j]) / (dx * dx) + (local_C1_old->mat_2D[i][j + 1] + local_C1_old->mat_2D[i][j - 1] - 2.0 * local_C1_old->mat_2D[i][j]) / (dy * dy)));

                local_C2->mat_2D[i][j] = local_C2_old->mat_2D[i][j] + dt * (
                    local_C1_old->mat_2D[i][j] - local_C2_old->mat_2D[i][j]
                    + D2 * ((local_C2_old->mat_2D[i + 1][j] + local_C2_old->mat_2D[i - 1][j] - 2.0 * local_C2_old->mat_2D[i][j]) / (dx * dx) + (local_C2_old->mat_2D[i][j + 1] + local_C2_old->mat_2D[i][j - 1] - 2.0 * local_C2_old->mat_2D[i][j]) / (dy * dy)));
            }
        }
    }

    update_ghost_cells_nonblocking(id, rows, cols, local_imax, local_jmax, local_C1, local_C2, 0);

    // Implement Neumann boundary conditions only if the local domain is at the global domain boundary
    int my_row = id / cols;
    int my_col = id % cols;

    // Add boundary conditions for the local domain
    if (start_i == 0) {
        for (int j = 0; j < local_jmax; j++) {
            local_C1->mat_2D[0][j] = local_C1->mat_2D[1][j];
            local_C2->mat_2D[0][j] = local_C2->mat_2D[1][j];
            
        }
    }
    if (end_i == imax) {
        for (int j = 0; j < local_jmax; j++) {
            local_C1->mat_2D[local_imax - 1][j] = local_C1->mat_2D[local_imax - 2][j];
            local_C2->mat_2D[local_imax - 1][j] = local_C2->mat_2D[local_imax - 2][j];
        }
    }
    if (start_j == 0) {
        for (int i = 0; i < local_imax; i++) {
            local_C1->mat_2D[i][0] = local_C1->mat_2D[i][1];
            local_C2->mat_2D[i][0] = local_C2->mat_2D[i][1];
            
        }
    }
    if (end_j == jmax) {
        for (int i = 0; i < local_imax; i++) {
            local_C1->mat_2D[i][local_jmax - 1] = local_C1->mat_2D[i][local_jmax - 2];
            local_C2->mat_2D[i][local_jmax - 1] = local_C2->mat_2D[i][local_jmax - 2];
        }
    }

    t += dt;
}


// Iteration function for local domain with periodic boundaries
void do_local_iteration_periodic(int id, int rows, int cols, int local_imax, int local_jmax,
                                 Cmatrix*& local_C1, Cmatrix*& local_C1_old, Cmatrix*& local_C2, Cmatrix*& local_C2_old,
                                 double dx, double dy, double dt, double D1, double D2, double f, double q, double epsilon) {
    // Swap pointers for local matrices
    Cmatrix* temp = local_C1;
    local_C1 = local_C1_old;
    local_C1_old = temp;
    temp = local_C2;
    local_C2 = local_C2_old;
    local_C2_old = temp;

    // Calculate new concentrations for all the points in the local domain
    for (int i = 1; i < local_imax - 1; i++) {
        for (int j = 1; j < local_jmax - 1; j++) {
            local_C1->mat_2D[i][j] = local_C1_old->mat_2D[i][j] + dt * (
                (local_C1_old->mat_2D[i][j] * (1.0 - local_C1_old->mat_2D[i][j]) - f * local_C2_old->mat_2D[i][j] * (local_C1_old->mat_2D[i][j] - q) / (local_C1_old->mat_2D[i][j] + q)) / epsilon
                + D1 * ((local_C1_old->mat_2D[i + 1][j] + local_C1_old->mat_2D[i - 1][j] - 2.0 * local_C1_old->mat_2D[i][j]) / (dx * dx)
                + (local_C1_old->mat_2D[i][j + 1] + local_C1_old->mat_2D[i][j - 1] - 2.0 * local_C1_old->mat_2D[i][j]) / (dy * dy)));

            local_C2->mat_2D[i][j] = local_C2_old->mat_2D[i][j] + dt * (
                local_C1_old->mat_2D[i][j] - local_C2_old->mat_2D[i][j]
                + D2 * ((local_C2_old->mat_2D[i + 1][j] + local_C2_old->mat_2D[i - 1][j] - 2.0 * local_C2_old->mat_2D[i][j]) / (dx * dx)
                + (local_C2_old->mat_2D[i][j + 1] + local_C2_old->mat_2D[i][j - 1] - 2.0 * local_C2_old->mat_2D[i][j]) / (dy * dy)));
        }
    }
    update_ghost_cells_nonblocking(id, rows, cols, local_imax, local_jmax, local_C1, local_C2, 1);
    t += dt;
}


void calc_constants()
{
	dx = x_max / ((double)imax - 1);
	dy = y_max / ((double)imax - 1);

	t = 0.0;

	dt = 0.1 * pow(min(dx, dy), 2.0) / (2.0 * max(D1, D2));

    dt = min(dt, 0.001);

	// cout << "dt = " << dt << " dx = " << dx << " dy = " << dy << endl;
}


// Write the local matrix to a file
void mpi_write_filtered_matrix(Cmatrix* local_matrix, const std::string& filename, int global_imax, int global_jmax, int start_i, int start_j, int local_imax, int local_jmax, MPI_Comm comm, int bc) {
    std::vector<double> filtered_data;
    std::vector<int> positions;

    if (bc == 0) {
        for (int i = 0; i < local_imax; i++) {
            for (int j = 0; j < local_jmax; j++) {
                double value = local_matrix->mat_2D[i][j];
                int global_i = start_i + i;
                int global_j = start_j + j;

                if ((i == 0 && j == 0) || (i == 0 && j == local_jmax - 1) || (i == local_imax - 1 && j == 0) || (i == local_imax - 1 && j == local_jmax - 1)) {
                    if (!(global_i == 0) && !(global_i == global_imax - 1) && !(global_j == 0) && !(global_j == global_jmax - 1)) {
                        continue;
                    }
                }

                
                filtered_data.push_back(value);
                positions.push_back((start_i + i) * global_jmax + (start_j + j));
            }
        }
    } else {
        for (int i = 1; i < local_imax - 1; i++) {
            for (int j = 1; j < local_jmax - 1; j++) {
                double value = local_matrix->mat_2D[i][j];
                int global_i = start_i + 1 + (i - 1);
                int global_j = start_j + 1 + (j - 1);

                filtered_data.push_back(value);
                positions.push_back(global_i * global_jmax + global_j);
            }
        }
    }
    

    MPI_File fh;
    MPI_Datatype index_type;
    MPI_File_open(comm, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_Type_indexed(filtered_data.size(), std::vector<int>(filtered_data.size(), 1).data(), positions.data(), MPI_DOUBLE, &index_type);
    MPI_Type_commit(&index_type);

    MPI_File_set_view(fh, 0, MPI_DOUBLE, index_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(fh, filtered_data.data(), filtered_data.size(), MPI_DOUBLE, MPI_STATUS_IGNORE);

    MPI_Type_free(&index_type);
    MPI_File_close(&fh);
}


void write_to_file(Cmatrix* local_matrix, Cmatrix* local_matrix2, int global_imax, int global_jmax, int start_i, int start_j, int local_imax, int local_jmax, MPI_Comm comm, int bc, int out_cnt) {
    std::string fname, fname2;
    if (bc == 0) {
        fname = "./result_neu/output_C1_" + std::to_string(out_cnt) + ".dat";
        fname2 = "./result_neu/output_C2_" + std::to_string(out_cnt) + ".dat";
    } else {
        fname = "./result_per/output_C1_" + std::to_string(out_cnt) + ".dat";
        fname2 = "./result_per/output_C2_" + std::to_string(out_cnt) + ".dat";
    }
    mpi_write_filtered_matrix(local_matrix, fname, global_imax, global_jmax, start_i, start_j, local_imax, local_jmax, comm, bc);
    mpi_write_filtered_matrix(local_matrix2, fname2, global_imax, global_jmax, start_i, start_j, local_imax, local_jmax, comm, bc);
}


// ADDITIONAL FUNCTION: Save the current state of the simulation to a checkpoint file
void save_checkpoint(Cmatrix* local_C1, Cmatrix* local_C2, double current_time, int local_imax, int local_jmax, int id) {
    std::ostringstream fname;
    fname << "checkpoint_" << id << ".dat";
    std::ofstream out_file(fname.str(), ios::binary);

    // Save the current time
    out_file.write(reinterpret_cast<char*>(&current_time), sizeof(double));

    // Save the dimensions of the local matrices
    out_file.write(reinterpret_cast<char*>(&local_imax), sizeof(int));
    out_file.write(reinterpret_cast<char*>(&local_jmax), sizeof(int));

    // Save the matrix data
    out_file.write(reinterpret_cast<char*>(local_C1->mat_1D), sizeof(double) * local_imax * local_jmax);
    out_file.write(reinterpret_cast<char*>(local_C2->mat_1D), sizeof(double) * local_imax * local_jmax);

    out_file.close();
}


// ADDITIONAL FUNCTION: Load the state of the simulation from a checkpoint file
bool load_checkpoint(Cmatrix*& local_C1, Cmatrix*& local_C2, double& current_time, int local_imax, int local_jmax, int id) {
    std::ostringstream fname;
    fname << "checkpoint_" << id << ".dat";
    std::ifstream in_file(fname.str(), ios::binary);
    if (!in_file.is_open()) {
        return false;
    }

    // Load the current time
    in_file.read(reinterpret_cast<char*>(&current_time), sizeof(double));

    // Load the dimensions of the local matrices
    in_file.read(reinterpret_cast<char*>(local_C1->mat_1D), sizeof(double) * local_imax * local_jmax);
    in_file.read(reinterpret_cast<char*>(local_C2->mat_1D), sizeof(double) * local_imax * local_jmax);

    in_file.close();
    return true;
}


int main(int argc, char *argv[]) {

    // Please SET the boundary condition that you want to use! 0: Neumann, 1: Periodic
    int bc = 1;

    // Please SET the checkpoint interval
    double checkpoint_interval = 1.0;

    double last_checkpoint_time = 0;
    int out_cnt = 0, it = 0;

    MPI_Init(&argc, &argv);
    int id, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    calc_constants();

    int rows, cols;
    int start_i, start_j, end_i, end_j;
    Cmatrix *local_C1, *local_C2, *local_C1_old, *local_C2_old;
    if (bc == 0) {
        calculate_domain_dimensions(p, rows, cols, id, imax, jmax, x_max, y_max, local_C1, local_C2, local_C1_old, local_C2_old, start_i, start_j, end_i, end_j, bc);
        initialize_local_matrices(id, rows, cols, imax, jmax, x_max, y_max, local_C1, local_C2, start_i, start_j, end_i, end_j, bc);
        update_ghost_cells_nonblocking(id, rows, cols, local_C1->n, local_C1->m, local_C1, local_C2, bc);
        write_to_file(local_C1, local_C2, imax, jmax, start_i, start_j, local_C1->n, local_C1->m, MPI_COMM_WORLD, bc, out_cnt);
    } else {
        calculate_domain_dimensions(p, rows, cols, id, imax, jmax, x_max, y_max, local_C1, local_C2, local_C1_old, local_C2_old, start_i, start_j, end_i, end_j, bc);
        initialize_local_matrices(id, rows, cols, imax, jmax, x_max, y_max, local_C1, local_C2, start_i, start_j, end_i, end_j, bc);
        update_ghost_cells_nonblocking(id, rows, cols, local_C1->n, local_C1->m, local_C1, local_C2, bc);
        write_to_file(local_C1, local_C2, imax, jmax, start_i, start_j, local_C1->n, local_C1->m, MPI_COMM_WORLD, bc, out_cnt);
    }

    // ADDITIONAL OPTIMIZE: Load checkpoint if available
    if (load_checkpoint(local_C1, local_C2, t, local_C1->n, local_C1->m, id)) {
        std::cout << "Successfully loaded checkpoint." << std::endl;
    } else {
        std::cout << "No checkpoint found or failed to load." << std::endl;
    }

    out_cnt++;
    t_out += dt_out;

    while (t < t_max)
	{
        if (bc == 0) {
            do_local_iteration_neumann(id, rows, cols, imax, jmax, local_C1->n, local_C1->m, local_C1, local_C1_old, local_C2, local_C2_old, dx, dy, dt, D1, D2, f, q, epsilon, start_i, start_j, end_i, end_j);
        } else {
            do_local_iteration_periodic(id, rows, cols, local_C1->n, local_C1->m, local_C1, local_C1_old, local_C2, local_C2_old, dx, dy, dt, D1, D2, f, q, epsilon);
        }

        if (t_out <= t) {
            cout << "Writing output at t = " << t << endl;
            if (bc == 0) {
                write_to_file(local_C1, local_C2, imax, jmax, start_i, start_j, local_C1->n, local_C1->m, MPI_COMM_WORLD, bc, out_cnt);
            } else {
                write_to_file(local_C1, local_C2, imax, jmax, start_i, start_j, local_C1->n, local_C1->m, MPI_COMM_WORLD, bc, out_cnt);
            }

            out_cnt++;
            t_out += dt_out;
        }

        // ADDITIONAL OPTIMIZE: Save checkpoint if needed
        if (t - last_checkpoint_time >= checkpoint_interval) {
            save_checkpoint(local_C1, local_C2, t, local_C1->n, local_C1->m, id);
            last_checkpoint_time = t;
        }

		it++;
	}

    delete local_C1;
    delete local_C2;
    delete local_C1_old;
    delete local_C2_old;
    MPI_Finalize();
    return 0;
}
