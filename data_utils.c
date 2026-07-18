#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data_utils.h"

// ------------Convenient memory management functions------------

double* linspace(double start, double end, int numSteps) {
    double* array = (double*)malloc(numSteps * sizeof(double));
    if (array == NULL) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
    }
    double step = (end - start) / (numSteps - 1);
    for (int i = 0; i < numSteps; i++) {
        array[i] = start + i * step;
    }
    return array;
}

double* initializeArray(int numSteps) {
    double* array = (double*)malloc(numSteps * sizeof(double));
    if (array == NULL) {
        perror("Memory allocation error");
        exit(EXIT_FAILURE);
    }
    //for(int i=0; i<numSteps; i++)
    //    array[i]=0;
    return array;
}

// -----------------String Functions----------------

int split_doubles(const char *str, double *arr, int max_size) {
    char *token;
    char *endptr;
    int count = 0;

    // Make a copy of the string to avoid modifying the original
    char *str_copy = strdup(str);
    if (str_copy == NULL) {
        return -1;  // Memory allocation failed
    }

    token = strtok(str_copy, ",");
    while (token != NULL && count < max_size) {
        double val = strtof(token, &endptr);
        // Check for conversion errors
        if (token == endptr) {
            free(str_copy);
            return -1; // Conversion error
        }
        arr[count++] = val;
        token = strtok(NULL, ",");
    }

    free(str_copy);
    return count; // Number of floats stored in the array
}

// -----------------File Utilities-------------------

void writeArray_bin(double* inputArray, int size, char* filename) {
    FILE *file = fopen(filename, "wb"); // Binary mode
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    fwrite(inputArray, sizeof(double), size, file);
    fclose(file);
}

void writeFile(double* inputArray, double** dataArray, int numSteps, int numColumns, const char* filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    // Write data to the file.
    // We use exponential format to preserve precision.
    for (int i = 0; i < numSteps; i++) {
        fprintf(file, "%0.8e", inputArray[i]);
        for (int j = 0; j < numColumns; j++) {
            fprintf(file, "\t%0.8e", dataArray[j][i]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
}

void writeFile_bin(double* inputArray, double** dataArray, int numSteps, int numColumns, const char* filename) {
    FILE *file = fopen(filename, "wb");  // Open the file in binary write mode
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }
    // Write the number of steps and the number of columns first for easier data reading.
    fwrite(&numSteps, sizeof(int), 1, file);
    fwrite(&numColumns, sizeof(int), 1, file);

    // Write the input array (time) to the file
    fwrite(inputArray, sizeof(double), numSteps, file);

    // Write each data array to the file
    for (int j = 0; j < numColumns; j++) {
        fwrite(dataArray[j], sizeof(double), numSteps, file);
    }
    fclose(file);
}

