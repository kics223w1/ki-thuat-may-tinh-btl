#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Matrix structure for basic linear algebra operations
typedef struct {
    double** data;
    int rows, cols;
} Matrix;

// Matrix constructor
Matrix* matrix_create(int rows, int cols) {
    Matrix* m = (Matrix*)malloc(sizeof(Matrix));
    m->rows = rows;
    m->cols = cols;
    m->data = (double**)malloc(rows * sizeof(double*));
    for(int i = 0; i < rows; i++) {
        m->data[i] = (double*)malloc(cols * sizeof(double));
        memset(m->data[i], 0, cols * sizeof(double));
    }
    return m;
}

// Matrix destructor
void matrix_free(Matrix* m) {
    if(m) {
        for(int i = 0; i < m->rows; i++) {
            free(m->data[i]);
        }
        free(m->data);
        free(m);
    }
}

// Matrix copy
Matrix* matrix_copy(const Matrix* src) {
    Matrix* dst = matrix_create(src->rows, src->cols);
    for(int i = 0; i < src->rows; i++) {
        for(int j = 0; j < src->cols; j++) {
            dst->data[i][j] = src->data[i][j];
        }
    }
    return dst;
}

// Matrix transpose
Matrix* matrix_transpose(const Matrix* m) {
    Matrix* res = matrix_create(m->cols, m->rows);
    for(int i = 0; i < m->rows; i++)
        for(int j = 0; j < m->cols; j++)
            res->data[j][i] = m->data[i][j];
    return res;
}

// Matrix multiplication
Matrix* matrix_multiply(const Matrix* a, const Matrix* b) {
    Matrix* res = matrix_create(a->rows, b->cols);
    for(int i = 0; i < a->rows; i++)
        for(int j = 0; j < b->cols; j++)
            for(int k = 0; k < a->cols; k++)
                res->data[i][j] += a->data[i][k] * b->data[k][j];
    return res;
}

// Matrix subtraction
Matrix* matrix_subtract(const Matrix* a, const Matrix* b) {
    Matrix* res = matrix_create(a->rows, a->cols);
    for(int i = 0; i < a->rows; i++)
        for(int j = 0; j < a->cols; j++)
            res->data[i][j] = a->data[i][j] - b->data[i][j];
    return res;
}

// Matrix print
void matrix_print(const Matrix* m) {
    for(int i = 0; i < m->rows; i++) {
        for(int j = 0; j < m->cols; j++)
            printf("%.6f ", m->data[i][j]);
        printf("\n");
    }
}

// Vector operations
double vector_dot(const double* a, const double* b, int size) {
    double sum = 0.0;
    for(int i = 0; i < size; i++)
        sum += a[i] * b[i];
    return sum;
}

double* matrix_vector_multiply(const Matrix* m, const double* v, int* result_size) {
    double* res = (double*)malloc(m->rows * sizeof(double));
    memset(res, 0, m->rows * sizeof(double));
    *result_size = m->rows;

    for(int i = 0; i < m->rows; i++)
        for(int j = 0; j < m->cols; j++)
            res[i] += m->data[i][j] * v[j];
    return res;
}

double* vector_subtract(const double* a, const double* b, int size) {
    double* res = (double*)malloc(size * sizeof(double));
    for(int i = 0; i < size; i++)
        res[i] = a[i] - b[i];
    return res;
}

// Swap function for doubles
void swap_double(double* a, double* b) {
    double temp = *a;
    *a = *b;
    *b = temp;
}

// Swap function for matrix rows
void swap_matrix_rows(Matrix* A, int row1, int row2) {
    double* temp = A->data[row1];
    A->data[row1] = A->data[row2];
    A->data[row2] = temp;
}

// Gaussian elimination for solving linear system Ax = b
double* solve_linear_system(Matrix* A, const double* b, int n, int* result_size) {
    // Create a copy of A to avoid modifying the original
    Matrix* A_copy = matrix_copy(A);
    double* b_copy = (double*)malloc(n * sizeof(double));
    memcpy(b_copy, b, n * sizeof(double));

    // Forward elimination
    for(int p = 0; p < n; p++) {
        // Find pivot row
        int max = p;
        for(int i = p+1; i < n; i++)
            if(fabs(A_copy->data[i][p]) > fabs(A_copy->data[max][p]))
                max = i;

        // Swap rows
        swap_matrix_rows(A_copy, p, max);
        swap_double(&b_copy[p], &b_copy[max]);

        // Check for singular matrix
        if(fabs(A_copy->data[p][p]) < 1e-10) {
            printf("Error: Singular matrix\n");
            matrix_free(A_copy);
            free(b_copy);
            *result_size = 0;
            return NULL;
        }

        // Eliminate column
        for(int i = p+1; i < n; i++) {
            double alpha = A_copy->data[i][p] / A_copy->data[p][p];
            b_copy[i] -= alpha * b_copy[p];
            for(int j = p; j < n; j++)
                A_copy->data[i][j] -= alpha * A_copy->data[p][j];
        }
    }

    // Back substitution
    double* x = (double*)malloc(n * sizeof(double));
    for(int i = n-1; i >= 0; i--) {
        double sum = 0.0;
        for(int j = i+1; j < n; j++)
            sum += A_copy->data[i][j] * x[j];
        x[i] = (b_copy[i] - sum) / A_copy->data[i][i];
    }

    matrix_free(A_copy);
    free(b_copy);
    *result_size = n;
    return x;
}

// Read signal from file
double* read_signal(const char* filename, int* signal_size) {
    FILE* file = fopen(filename, "r");
    if(!file) {
        printf("Error: Cannot open file %s\n", filename);
        *signal_size = 0;
        return NULL;
    }

    // Count numbers first
    int count = 0;
    double temp;
    while(fscanf(file, "%lf", &temp) == 1) {
        count++;
    }

    // Reset file pointer
    fseek(file, 0, SEEK_SET);

    // Allocate memory
    double* signal = (double*)malloc(count * sizeof(double));
    *signal_size = count;

    // Read numbers
    for(int i = 0; i < count; i++) {
        fscanf(file, "%lf", &signal[i]);
    }

    fclose(file);
    return signal;
}

// Compute autocorrelation matrix R (Toeplitz matrix)
Matrix* compute_autocorrelation(const double* x, int N, int M) {
    Matrix* R = matrix_create(M, M);

    // First compute autocorrelation coefficients
    double* rxx = (double*)malloc(M * sizeof(double));
    memset(rxx, 0, M * sizeof(double));

    for(int k = 0; k < M; k++) {
        double sum = 0.0;
        for(int n = k; n < N; n++) {
            sum += x[n] * x[n - k];
        }
        rxx[k] = sum / N;
    }

    // Build Toeplitz matrix
    for(int i = 0; i < M; i++) {
        for(int j = 0; j < M; j++) {
            int lag = abs(i - j);
            R->data[i][j] = rxx[lag];
        }
    }

    free(rxx);
    return R;
}

// Compute cross-correlation vector γ_d
double* compute_cross_correlation(const double* d, const double* x, int N, int M, int* result_size) {
    double* gamma_d = (double*)malloc(M * sizeof(double));
    memset(gamma_d, 0, M * sizeof(double));
    *result_size = M;

    for(int k = 0; k < M; k++) {
        double sum = 0.0;
        for(int n = k; n < N; n++) {
            sum += d[n] * x[n - k];
        }
        gamma_d[k] = sum / N;
    }
    return gamma_d;
}

// Apply Wiener filter
double* apply_filter(const double* x, const double* h, int N, int M, int* result_size) {
    double* y = (double*)malloc(N * sizeof(double));
    memset(y, 0, N * sizeof(double));
    *result_size = N;

    for(int n = 0; n < N; n++) {
        for(int k = 0; k < M; k++) {
            int idx = n - k;
            if(idx >= 0 && idx < N) {
                y[n] += h[k] * x[idx];
            }
        }
    }
    return y;
}

// Calculate variance
double calculate_variance(const double* signal, int size) {
    double mean = 0.0;
    for(int i = 0; i < size; i++) mean += signal[i];
    mean /= size;

    double var = 0.0;
    for(int i = 0; i < size; i++) var += (signal[i] - mean) * (signal[i] - mean);
    return var / size;
}

// Calculate MMSE
double calculate_mmse(const double* desired, const double* output, int size) {
    double mse = 0.0;
    for(int i = 0; i < size; i++) {
        double error = desired[i] - output[i];
        mse += error * error;
    }
    return mse / size;
}

int main() {
    // Read input files
    int desired_size, input_size;
    double* desired = read_signal("desired.txt", &desired_size);
    double* input = read_signal("input.txt", &input_size);

    if(!desired || !input || desired_size != 10 || input_size != 10) {
        printf("Error: size not match\n");
        FILE* outfile = fopen("output.txt", "w");
        if(outfile) {
            fprintf(outfile, "Error: size not match\n");
            fclose(outfile);
        }
        free(desired);
        free(input);
        return 1;
    }

    // Filter order (number of taps)
    int M = 10;  // Use 10 taps for the Wiener filter
    int N = input_size;

    // Compute autocorrelation matrix
    Matrix* R = compute_autocorrelation(input, N, M);

    // Compute cross-correlation vector
    int gamma_size;
    double* gamma_d = compute_cross_correlation(desired, input, N, M, &gamma_size);

    // Solve for optimal filter coefficients
    int h_opt_size;
    double* h_opt = solve_linear_system(R, gamma_d, M, &h_opt_size);

    if(!h_opt) {
        printf("Error: Cannot solve linear system\n");
        matrix_free(R);
        free(desired);
        free(input);
        free(gamma_d);
        return 1;
    }

    // Apply filter
    int output_size;
    double* output = apply_filter(input, h_opt, N, M, &output_size);

    // Calculate MMSE
    double mmse = calculate_mmse(desired, output, input_size);

    // Round output values to 1 decimal place
    for(int i = 0; i < output_size; i++) {
        output[i] = round(output[i] * 10.0) / 10.0;
    }
    
    // Calculate MMSE after rounding output
    mmse = calculate_mmse(desired, output, input_size);
    mmse = round(mmse * 10.0) / 10.0;

    // Output results to console
    printf("Filtered output:");
    for(int i = 0; i < output_size; i++) {
        printf(" %9.4f", output[i]);
    }
    printf("\n");

    printf("MMSE: %.4f\n", mmse);

    // Write to output file
    FILE* outfile = fopen("output.txt", "w");
    if(outfile) {
        fprintf(outfile, "Filtered output:");
        for(int i = 0; i < output_size; i++) {
            fprintf(outfile, " %9.4f", output[i]);
        }
        fprintf(outfile, "\n");
        fprintf(outfile, "MMSE: %.4f\n", mmse);
        fclose(outfile);
    }

    // Free memory
    matrix_free(R);
    free(desired);
    free(input);
    free(gamma_d);
    free(h_opt);
    free(output);

    return 0;
}