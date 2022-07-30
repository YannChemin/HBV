#include <stdio.h>
#include <stdlib.h>
#include "arrays.h"

#define ARRAYSIZE(x)  (sizeof(x))
/*http://cboard.cprogramming.com/c-programming/47105-how-read-csv-file.html*/
float ** readcsv(char *filename, int row, int col){
	int x,y;
	FILE *file = fopen(filename, "r");
	//float array[1910][8];//array[row][col]
	float **array	= af2d(col,row);//array[row][col] Transposed for reading
	size_t i, j;
	char buffer[BUFSIZ], *ptr;
	for ( i = 0; fgets(buffer, sizeof buffer, file); ++i ){
		for ( j = 0, ptr = buffer; j < ARRAYSIZE(*array); ++j, ++ptr ){
			array[i][j] = (float)strtof(ptr, &ptr);
		}
	}
	fclose(file);
	float **arrayout= af2d(row,col);//array[row][col]
	for (x=0;x<row;x++)
		for (y=0;y<col;y++)
			arrayout[x][y]=array[y][x];
// 	size_t k;
// 	for ( j = 0; j < i; ++j ){
// 		printf("array[%lu]: ", (long unsigned)j);
// 		for ( k = 0; k < ARRAYSIZE(*array); ++k )  printf("%.2f ", array[j][k]);
// 		putchar('\n');
// 	}
	return &arrayout[0];
}