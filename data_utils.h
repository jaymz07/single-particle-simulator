#ifndef DATA_UTILS_H
#define DATA_UTILS_H

double* linspace(double start, double end, int numSteps);
double* initializeArray(int numSteps);

void writeArray_bin(double* inputArray, int size, char* filename);

void writeFile(double* inputArray, double** dataArray, int numSteps, int numColumns, const char* filename);
void writeFile_bin(double* inputArray, double** dataArray, int numSteps, int numColumns, const char* filename);

int split_doubles(const char *str, double *arr, int max_size);

#endif
