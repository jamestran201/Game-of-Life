#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char* argv[]) {
	if (argc < 5) {
		puts("Usage: ./read_bin max_stage row col sleep_mode");
		puts("Sleep mode: 1(ON) or 0(OFF). Indicates whether the program sleeps after printing out each state");
		return -1;
	}

	int counter;
	int row = atoi(argv[3]);
	int col = atoi(argv[2]);
	int max_cell = row * col;
	int max_stage = atoi(argv[1]);
	int sleep_mode = atoi(argv[4]);
	FILE *fptr;
	int i, j;
	
	for(j = 0; j < max_stage; j++) {
		char filename[50];
		sprintf(filename, "stage%d", j);
		fptr = fopen(&filename, "rb");
	
		printf("Stage %d\n", j);
		for(counter=0; counter < max_cell; counter++) {
			fread(&i, sizeof(int), 1, fptr);
			if (i == 0) printf("-");
			else printf("O");
			if (counter % col == (col-1)) puts("");
		}
		
		if (sleep_mode == 1) sleep(1);

		puts("");
		fclose(fptr);
	}
	return 0;
}
