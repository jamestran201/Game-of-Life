#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <string.h>

// Struct to hold neighbor ranks
struct Neighbor_ranks {
	int rank_top;
	int rank_bottom;
	int rank_right;
	int rank_left;
	int rank_top_right;
	int rank_top_left;
	int rank_bot_right;
	int rank_bot_left;
};

int* init_cell_grid(int r, int c);
void init_local_grid(int* grid, int* data, int r, int c);
struct Neighbor_ranks get_neighbor_ranks(MPI_Comm grid_comm, int proc_rows);
void print_grid(int* grid, int r, int c, int rank);
void print_2d(int* grid, int r, int c, int rank);
void gather_data(int* grid_with_ghost, int grid_ghost_rows, int grid_ghost_cols, int local_size[2], MPI_Comm grid_comm, struct Neighbor_ranks neighbor_ranks);
int* scatter_grid(int my_rank, int cell_rows, int cell_cols, int proc_rows, int proc_cols, int base_size[2], int remain[2], 
	int row_rank, int col_rank, int coordinates[2], int* original_data, MPI_Comm row_comm, MPI_Comm col_comm);

int main(int argc, char* argv[]) {
	int p, my_rank, proc_rows, proc_cols, cell_rows, cell_cols, max_iter;
	int dim_sizes[2];
	int wrap_around[2];
	int coordinates[2];
	int reorder = 1;
	int grid_rank;
	int* cell_grid = NULL;
	MPI_Comm grid_comm;	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	double start_time = MPI_Wtime();

	proc_rows = atoi(argv[1]); // number of rows in grid of processes
	proc_cols = atoi(argv[2]); // number of columns in grid of processes
	cell_rows = atoi(argv[3]); // number of rows in grid of cell
	cell_cols = atoi(argv[4]); // number of columns in grid of cell
	max_iter = atoi(argv[5]); // Max stage of the simulation	

	if (my_rank == 0) {
		cell_grid = init_cell_grid(cell_rows, cell_cols);
		
		// glider configuration
		cell_grid[2*cell_cols + 2] = 1;
		cell_grid[3*cell_cols + 3] = 1;
		cell_grid[4*cell_cols + 3] = 1;
		cell_grid[4*cell_cols + 2] = 1;
		cell_grid[4*cell_cols + 1] = 1;
	
		// DEBUG
		// print_grid(cell_grid, cell_rows, cell_cols, my_rank);
	}

	// Create Cartesian communicator
	dim_sizes[0] = proc_rows;
	dim_sizes[1] = proc_cols;
	wrap_around[0] = wrap_around[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes, wrap_around, reorder, &grid_comm);

	MPI_Comm_rank(grid_comm, &grid_rank);
	MPI_Cart_coords(grid_comm, grid_rank, 2, coordinates);

	// Create row and columns communicators
	MPI_Comm row_comm;
	MPI_Comm col_comm;
	int free_col_coords[2];
	int free_row_coords[2];

	free_col_coords[0] = 1;
	free_col_coords[1] = 0;
	free_row_coords[0] = 0;
	free_row_coords[1] = 1;

	MPI_Cart_sub(grid_comm, free_row_coords, &row_comm);
	MPI_Cart_sub(grid_comm, free_col_coords, &col_comm);

	// Obtain rank in row and column communicator
	int row_rank;
	int col_rank;
	MPI_Comm_rank(row_comm, &row_rank);
	MPI_Comm_rank(col_comm, &col_rank);

	// The least number of cell rows and columns that each process will have
	int base_rows = cell_rows / proc_rows;
	int base_cols = cell_cols / proc_cols;

	// The extra cell rows and columns
	int remain_rows = cell_rows % proc_rows;
	int remain_cols = cell_cols % proc_cols;
	int remain[2] = {remain_rows, remain_cols};

	// Define base size of local grid
	int base_size[2] = {base_rows, base_cols};

	// Dimensions of local cell grid
	int local_size[2]; 
	local_size[0] = (col_rank - remain[0] < 0) ? base_size[0] + 1 : base_size[0];
	local_size[1] = (row_rank - remain[1] < 0) ? base_size[1] + 1 : base_size[1];

	// Local cell grid of each process
	// Distribute the data to each process
	int* local_data = scatter_grid(my_rank, cell_rows, cell_cols, proc_rows, proc_cols, base_size, remain,row_rank, col_rank, coordinates, 
		cell_grid, row_comm, col_comm);

	// DEBUG
	//print_grid(local_data, local_size[0], local_size[1], my_rank);

	// Initialize local grid with ghost points
	int grid_ghost_rows = local_size[0] + 2;
	int grid_ghost_cols = local_size[1] + 2;
	int grid_with_ghost[grid_ghost_rows][grid_ghost_cols];

	init_local_grid(&grid_with_ghost, local_data, grid_ghost_rows, grid_ghost_cols);

	// DEBUG
	//puts("Before");
	//print_2d(&grid_with_ghost, grid_ghost_rows, grid_ghost_cols, my_rank);

	// Obtain ranks of neighbors
	struct Neighbor_ranks neighbor_ranks = get_neighbor_ranks(grid_comm, proc_rows);

	// DEBUG
	// printf("Neighbors of %d clockwise starting at top are: %d %d %d %d %d %d %d %d\n", grid_rank, neighbor_ranks.rank_top, neighbor_ranks.rank_top_right, 
	// 	neighbor_ranks.rank_right, neighbor_ranks.rank_bot_right, neighbor_ranks.rank_bottom, neighbor_ranks.rank_bot_left, neighbor_ranks.rank_left, 
	// 	neighbor_ranks.rank_top_left);

	//gather_data(&grid_with_ghost, grid_ghost_rows, grid_ghost_cols, local_size, grid_comm, neighbor_ranks);

	// Create datatype for interior grid without ghost points
	MPI_Datatype grid;
	// starting coord of subarray
	int start[2] = {1 , 1};
	//number of elements in each dimesion (rows,cols) in (full) local array
	int arrsize[2] = {grid_ghost_rows, grid_ghost_cols};

	MPI_Type_create_subarray(2, arrsize, local_size, start, MPI_ORDER_C, MPI_INT, &grid);
	MPI_Type_commit(&grid);

	// Create datatype for global grid
	MPI_Datatype view;
	//x,y starting coords of local grid in the global grid
	int startV[2];

	//Calculate the starting position of this grid in the global cell grid
	int tmp_remain_row = remain_rows;
	int tmp_remain_col = remain_cols;
	int g;
	for(g = 0; g < col_rank; g++) {
		startV[0] = startV[0] + base_rows;
		if (tmp_remain_row > 0) {
			tmp_remain_row--;
			startV[0]++;
		}
	}

	for(g = 0; g < row_rank; g++) {
		startV[1] = startV[1] + base_cols;
		if (tmp_remain_col > 0) {
			tmp_remain_col--;
			startV[1]++;
		}
	}

	//number of elements in each dimensions in full Global array
	int arrsizeV[2] = { cell_rows, cell_cols };
	//dimensions of local data to output
	MPI_Type_create_subarray(2, arrsizeV, local_size, startV, MPI_ORDER_C, MPI_INT, &view);
	MPI_Type_commit(&view);



	int iter = 0;
	int tmp_grid[local_size[0]][local_size[1]];
	while(iter < max_iter) {

		// Get and send data to neighboring processes
		gather_data(&grid_with_ghost, grid_ghost_rows, grid_ghost_cols, local_size, grid_comm, neighbor_ranks);

		// Write current state to file
		//MPI IO
		char filename[50];
		MPI_File fh;

		sprintf(filename, "stage%d", iter);
		MPI_File_open(grid_comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
		MPI_File_set_view(fh, 0, MPI_INT, view, "native", MPI_INFO_NULL);
		MPI_File_write_all(fh, &(grid_with_ghost[0][0]), 1, grid, MPI_STATUS_IGNORE);
		MPI_File_close(&fh);

		//DEBUG
		//printf("Stage %d - Process %d:\n", iter, grid_rank);
		//print_2d(grid_with_ghost, grid_ghost_rows, grid_ghost_cols, grid_rank);
	
		// Decide the new state of cells
		int i;
		for(i = 1; i < local_size[0]+1; i++) {
			int j;
			for(j = 1; j < local_size[1]+1; j++) {
				int alive_neighbors = 0;
				
				if (grid_with_ghost[i][j+1] == 1) alive_neighbors++;
				if (grid_with_ghost[i][j-1] == 1) alive_neighbors++;
				if (grid_with_ghost[i-1][j] == 1) alive_neighbors++;
				if (grid_with_ghost[i+1][j] == 1) alive_neighbors++;
				if (grid_with_ghost[i+1][j-1] == 1) alive_neighbors++;
				if (grid_with_ghost[i+1][j+1] == 1) alive_neighbors++;
				if (grid_with_ghost[i-1][j-1] == 1) alive_neighbors++;
				if (grid_with_ghost[i-1][j+1] == 1) alive_neighbors++;

				if (grid_with_ghost[i][j] == 0) {
					if (alive_neighbors == 3) tmp_grid[i-1][j-1] = 1;
					else tmp_grid[i-1][j-1] = 0;
				} else {
					if (alive_neighbors < 2 || alive_neighbors > 3) tmp_grid[i-1][j-1] = 0;
					else if (alive_neighbors == 2 || alive_neighbors == 3) tmp_grid[i-1][j-1] = 1;
				}
			}
		}	
		
		// Copy new state to local grid
		int t;
		for(t = 1; t < local_size[0] + 1; t++) {
			memcpy(&(grid_with_ghost[t][1]), tmp_grid[t-1], sizeof(int)*local_size[1]);
		}

		iter++;
	}

	// DEBUG
	//puts("After");
	//print_2d(&grid_with_ghost, grid_ghost_rows, grid_ghost_cols, my_rank);

	//Print time elapsed
	double end_time = MPI_Wtime();
	double time_taken = end_time - start_time;
	if (my_rank == 0) {
		printf("Simulation with %d processes, process grid %d x %d, cell grid %d x %d and %d stages: %.2fs\n", p, proc_rows, proc_cols, cell_rows, cell_cols, 
			max_iter, time_taken);
	}

	// Free communicators
	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&col_comm);
	MPI_Comm_free(&grid_comm);
	MPI_Type_free(&view);
	MPI_Type_free(&grid);
	//MPI_File_close(&fh);
	MPI_Finalize();

	// Deallocate memory
	free(cell_grid);
	free(local_data);

	return 0;
}

// Scatter the cell grid among the processes
int* scatter_grid(int my_rank, int cell_rows, int cell_cols, int proc_rows, int proc_cols, int base_size[2], int remain[2], 
	int row_rank, int col_rank, int coordinates[2], int* original_data, MPI_Comm row_comm, MPI_Comm col_comm) {
	
	int* row_data;
	int local_size[2]; // Dimensions of local cell grid
	local_size[0] = (col_rank - remain[0] < 0) ? base_size[0] + 1 : base_size[0];
	local_size[1] = (row_rank - remain[1] < 0) ? base_size[1] + 1 : base_size[1];


	// Segment the original grid into rows and distribute them along the first column
	if (coordinates[1] == 0) {
		// Arguments for MPI_Scatterv
		int *sendcounts = (int*) malloc(sizeof(int) * proc_rows);
		int *displs = (int*) malloc(sizeof(int) * proc_rows);
		int sum = 0;
		
		int i;
		
		// Calculate send counts and displacements
		for(i = 0; i < proc_rows; i++) {
			sendcounts[i] = base_size[0] * cell_cols;
	
			if (remain[0] > 0) {
				remain[0]--;
				sendcounts[i] = sendcounts[i] + cell_cols;
			}

			// DEBUG
			//printf("Process %d: row %d | send count %d\n", my_rank, i, sendcounts[i]);

			displs[i] = sum;

			sum = sum + sendcounts[i];
		}

		// Allocate memory for each process in the column
		row_data = init_cell_grid(local_size[0], cell_cols);

		// DEBUG
		// printf("Process %d has grid size %d x %d\n", my_rank, local_size[0], cell_cols);
		// printf("Process %d receiving %d cells\n", my_rank, sendcounts[col_rank]);
		// printf("Process %d with displacement %d\n", my_rank, displs[col_rank]);

		// Distribute the data
		MPI_Scatterv(original_data, sendcounts, displs, MPI_INT, row_data, sendcounts[col_rank], MPI_INT, 0, col_comm);

		// DEBUG
		//print_grid(row_data, local_size[0], cell_cols, col_rank);
	}

	//DEBUG
	//printf("Process %d was assigned %d x %d grid\n", my_rank, local_size[0], local_size[1]);

	MPI_Datatype original_column;
	MPI_Datatype local_column;

	// Create MPI column type for data from the first column
	MPI_Type_vector(local_size[0], 1, cell_cols, MPI_INT, &original_column);
	MPI_Type_create_resized(original_column, 0, sizeof(int), &original_column);
	MPI_Type_commit(&original_column);

	// Create MPI column type for data in each process
	MPI_Type_vector(local_size[0], 1, local_size[1], MPI_INT, &local_column);
	MPI_Type_create_resized(local_column, 0, sizeof(int), &local_column);
	MPI_Type_commit(&local_column);

	//DEBUG
	//printf("Process %d: type created successfully\n", my_rank);

	// Arguments for MPI_Scatterv
	int *sendcounts_col = (int*) malloc(sizeof(int) * proc_cols);
	int *displs_col = (int*) malloc(sizeof(int) * proc_cols);
	int sum_col = 0;
	int stride;

	//DEBUG
	//printf("Offset to the beginning of row %d is %d\n", coordinates[0], sum_col);

	// Calculate send counts and displacements
	int k;
	for(k = 0; k < proc_cols; k++) {
		sendcounts_col[k] = base_size[1];
	
		if (remain[1] > 0) {
			remain[1]--;
			sendcounts_col[k] = sendcounts_col[k] + 1;
		}

		// DEBUG
		//printf("Process %d: column %d | send count %d\n", my_rank, k, sendcounts_col[k]);

		displs_col[k] = sum_col;

		// DEBUG
		//printf("Process %d: column %d | displacements %d\n", my_rank, k, displs_col[k]);

		sum_col = sum_col + sendcounts_col[k];
	}

	// Allocate memory for local cell grid
	int* local_data = init_cell_grid(local_size[0], local_size[1]);

	// Distribute the data
	MPI_Scatterv(row_data, sendcounts_col, displs_col, original_column, local_data, sendcounts_col[row_rank], local_column, 0, row_comm);

	// DEBUG
	//print_grid(local_data, local_size[0], local_size[1], my_rank);

	if (coordinates[1] == 0) free(row_data);
	MPI_Type_free(&local_column);
	MPI_Type_free(&original_column);

	return local_data;
}

// Create a 2D cell grid
int* init_cell_grid(int r, int c) {
	int* grid = (int*) malloc(sizeof(int) * r * c);
	int i;
	for(i = 0; i < r*c; i++) {
		grid[i] = 0;
	}

	return grid;
}

// Print grid
void print_grid(int* grid, int r, int c, int rank) {
	int j;
	printf("(%d) :  ", rank);
	for(j = 0; j < r*c; j++) {
		printf("%d ", grid[j]);
	}
	puts("");
}

// Print 2D array
void print_2d(int* grid, int r, int c, int rank) {
	int i;
	for(i = 0; i < r; i++) {
		printf("Process %d: ", rank);

		int j;
		for(j = 0; j < c; j++) printf("%d ", *grid++);

		puts("");
	}
}

// Initialize local grid with ghost points
void init_local_grid(int* grid, int* data, int r, int c) {
	int x;
	for(x = 0; x < r; x++) {
		int y;
		for(y = 0; y < c; y++) {
			if (x > 0 && x < (r - 1) && y > 0 && y < (c - 1)) 
				*grid++ = *data++;
			else 
				*grid++ = 0;
		}
	}
}

// Get ranks of neighbors
struct Neighbor_ranks get_neighbor_ranks(MPI_Comm grid_comm, int proc_rows) {
	struct Neighbor_ranks neighbor_ranks;

	// Get ranks of top, right, bottom and left neighbors
	MPI_Cart_shift(grid_comm, 1, 1, &neighbor_ranks.rank_left, &neighbor_ranks.rank_right);
	MPI_Cart_shift(grid_comm, 0, 1, &neighbor_ranks.rank_top, &neighbor_ranks.rank_bottom);

	// Convert rank of right and left neighbors to coordinates to obtain ranks of remaining neighbors
	int tmp_coords_right[2], tmp_coords_left[2], tmp_coords[2];

	MPI_Cart_coords(grid_comm, neighbor_ranks.rank_right, 2, &tmp_coords_right);
	MPI_Cart_coords(grid_comm, neighbor_ranks.rank_left, 2, &tmp_coords_left);

	// Get rank of bottom right neighbor
	tmp_coords[0] = (tmp_coords_right[0] + 1) % proc_rows;
	tmp_coords[1] = tmp_coords_right[1];
	MPI_Cart_rank(grid_comm, tmp_coords, &neighbor_ranks.rank_bot_right);

	// Get rank of top right neighbor
	tmp_coords[0] = (tmp_coords_right[0] - 1 < 0) ? (proc_rows-1) : (tmp_coords_right[0] - 1);
	tmp_coords[1] = tmp_coords_right[1];
	MPI_Cart_rank(grid_comm, tmp_coords, &neighbor_ranks.rank_top_right);

	// Get rank of bottom left neighbor
	tmp_coords[0] = (tmp_coords_left[0] + 1) % proc_rows;
	tmp_coords[1] = tmp_coords_left[1];
	MPI_Cart_rank(grid_comm, tmp_coords, &neighbor_ranks.rank_bot_left);

	// Get rank of top left neighbor
	tmp_coords[0] = (tmp_coords_left[0] - 1 < 0) ? (proc_rows-1) : (tmp_coords_left[0] - 1);
	tmp_coords[1] = tmp_coords_left[1];
	MPI_Cart_rank(grid_comm, tmp_coords, &neighbor_ranks.rank_top_left);

	return neighbor_ranks;
}

// Get data from neighbors
void gather_data(int* grid_with_ghost, int grid_ghost_rows, int grid_ghost_cols, int local_size[2], MPI_Comm grid_comm, struct Neighbor_ranks neighbor_ranks) {
	// Create datatype for ghost columns
	MPI_Datatype local_col_type;
	MPI_Request requests[16];
	MPI_Status statuses[16];

	MPI_Type_vector(local_size[0], 1, grid_ghost_cols, MPI_INT, &local_col_type);
	MPI_Type_commit(&local_col_type);

	// Communicate with top neighbor
	MPI_Isend(grid_with_ghost + grid_ghost_cols + 1, local_size[1], MPI_INT, neighbor_ranks.rank_top, 0, grid_comm, &requests[0]);
	MPI_Irecv(grid_with_ghost + 1, local_size[1], MPI_INT, neighbor_ranks.rank_top, 1, grid_comm, &requests[1]);
	
	// Communicate with bottom neighbor

	MPI_Isend(grid_with_ghost + (grid_ghost_rows-2)*grid_ghost_cols + 1, local_size[1], MPI_INT, neighbor_ranks.rank_bottom, 1, grid_comm, &requests[2]);
	MPI_Irecv(grid_with_ghost + (grid_ghost_rows-1)*grid_ghost_cols + 1, local_size[1], MPI_INT, neighbor_ranks.rank_bottom, 0, grid_comm, &requests[3]);

	// Communicate with left neighbor
	MPI_Isend(grid_with_ghost + grid_ghost_cols + 1, 1, local_col_type, neighbor_ranks.rank_left, 2, grid_comm, &requests[4]);
	MPI_Irecv(grid_with_ghost + grid_ghost_cols, 1, local_col_type, neighbor_ranks.rank_left, 3, grid_comm, &requests[5]);

	// Communicate with right neighbor
	MPI_Isend(grid_with_ghost + grid_ghost_cols + (grid_ghost_cols-2), 1, local_col_type, neighbor_ranks.rank_right, 3, grid_comm, &requests[6]);
	MPI_Irecv(grid_with_ghost + grid_ghost_cols + (grid_ghost_cols-1), 1, local_col_type, neighbor_ranks.rank_right, 2, grid_comm, &requests[7]);

	// Communicate with neighbors at the corners

	//Top left
	MPI_Isend(grid_with_ghost + grid_ghost_cols + 1, 1, MPI_INT, neighbor_ranks.rank_top_left, 4, grid_comm, &requests[8]);
	MPI_Irecv(grid_with_ghost, 1, MPI_INT, neighbor_ranks.rank_top_left, 5, grid_comm, &requests[9]);

	//Bot right
	MPI_Isend(grid_with_ghost + (grid_ghost_rows-2)*grid_ghost_cols + (grid_ghost_cols-2), 1, MPI_INT, neighbor_ranks.rank_bot_right, 5, grid_comm, &requests[10]);
	MPI_Irecv(grid_with_ghost + (grid_ghost_rows-1)*grid_ghost_cols + (grid_ghost_cols-1), 1, MPI_INT, neighbor_ranks.rank_bot_right, 4, grid_comm, &requests[11]);

	//Top right
	MPI_Isend(grid_with_ghost + grid_ghost_cols + (grid_ghost_cols-2), 1, MPI_INT, neighbor_ranks.rank_top_right, 6, grid_comm, &requests[12]);
	MPI_Irecv(grid_with_ghost + (grid_ghost_cols-1), 1, MPI_INT, neighbor_ranks.rank_top_right, 7, grid_comm, &requests[13]);

	//Bot left
	MPI_Isend(grid_with_ghost + (grid_ghost_rows-2)*grid_ghost_cols + 1, 1, MPI_INT, neighbor_ranks.rank_bot_left, 7, grid_comm, &requests[14]);
	MPI_Irecv(grid_with_ghost + (grid_ghost_rows-1)*grid_ghost_cols, 1, MPI_INT, neighbor_ranks.rank_bot_left, 6, grid_comm, &requests[15]);
	
	MPI_Waitall(16, requests, statuses);
	MPI_Type_free(&local_col_type);
}
