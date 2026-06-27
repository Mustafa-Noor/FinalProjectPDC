#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_PENALTY -1

typedef enum {
    DIR_NONE = 0,
    DIR_DIAG,
    DIR_UP,
    DIR_LEFT
} Direction;

typedef struct {
    size_t rows;
    size_t cols;
    int *scores;
    Direction *directions;
} Matrix;

typedef struct {
    char *aligned_a;
    char *aligned_b;
} Alignment;

static void free_matrix(Matrix *matrix);
static void free_alignment(Alignment *alignment);
static int allocate_matrix(Matrix *matrix, size_t len_a, size_t len_b);
static void initialize_matrix(Matrix *matrix);
static int score_pair(char a, char b);
static int max3_with_direction(int diag, int up, int left, Direction *direction);
static void fill_matrix(Matrix *matrix, const char *seq_a, const char *seq_b);
static int traceback(const Matrix *matrix, const char *seq_a, const char *seq_b, Alignment *alignment);
static int calculate_score(const char *aligned_a, const char *aligned_b);
static void print_alignment(const Alignment *alignment);
static char *read_sequence(const char *prompt);
static char *generate_random_sequence(size_t length);
static int parse_size(const char *text, size_t *value);
static void print_usage(const char *program_name);

static size_t cell_index(const Matrix *matrix, size_t row, size_t col)
{
    return row * matrix->cols + col;
}

static void free_matrix(Matrix *matrix)
{
    if (matrix == NULL) {
        return;
    }

    free(matrix->scores);
    free(matrix->directions);
    matrix->scores = NULL;
    matrix->directions = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

static void free_alignment(Alignment *alignment)
{
    if (alignment == NULL) {
        return;
    }

    free(alignment->aligned_a);
    free(alignment->aligned_b);
    alignment->aligned_a = NULL;
    alignment->aligned_b = NULL;
}

static int allocate_matrix(Matrix *matrix, size_t len_a, size_t len_b)
{
    size_t total_cells;

    if (matrix == NULL) {
        return 0;
    }

    matrix->rows = len_a + 1;
    matrix->cols = len_b + 1;
    matrix->scores = NULL;
    matrix->directions = NULL;

    if (len_a == SIZE_MAX || len_b == SIZE_MAX) {
        fprintf(stderr, "Error: sequence length is too large.\n");
        return 0;
    }

    if (len_a > (size_t)INT_MAX || len_b > (size_t)INT_MAX) {
        fprintf(stderr, "Error: sequence length exceeds supported score range.\n");
        return 0;
    }

    if (matrix->rows != 0 && matrix->cols > SIZE_MAX / matrix->rows) {
        fprintf(stderr, "Error: matrix dimensions are too large.\n");
        return 0;
    }

    total_cells = matrix->rows * matrix->cols;

    if (total_cells > SIZE_MAX / sizeof(int) ||
        total_cells > SIZE_MAX / sizeof(Direction)) {
        fprintf(stderr, "Error: matrix allocation size overflow.\n");
        return 0;
    }

    matrix->scores = (int *)malloc(total_cells * sizeof(int));
    matrix->directions = (Direction *)malloc(total_cells * sizeof(Direction));

    if (matrix->scores == NULL || matrix->directions == NULL) {
        fprintf(stderr, "Error: failed to allocate DP matrix.\n");
        free_matrix(matrix);
        return 0;
    }

    return 1;
}

static void initialize_matrix(Matrix *matrix)
{
    size_t i;
    size_t j;

    /* Cell (0, 0) aligns two empty prefixes, so its score is zero. */
    matrix->scores[cell_index(matrix, 0, 0)] = 0;
    matrix->directions[cell_index(matrix, 0, 0)] = DIR_NONE;

    /* First column: align prefixes of A against an empty B using gaps. */
    for (i = 1; i < matrix->rows; i++) {
        matrix->scores[cell_index(matrix, i, 0)] = (int)i * GAP_PENALTY;
        matrix->directions[cell_index(matrix, i, 0)] = DIR_UP;
    }

    /* First row: align an empty A against prefixes of B using gaps. */
    for (j = 1; j < matrix->cols; j++) {
        matrix->scores[cell_index(matrix, 0, j)] = (int)j * GAP_PENALTY;
        matrix->directions[cell_index(matrix, 0, j)] = DIR_LEFT;
    }
}

static int score_pair(char a, char b)
{
    return a == b ? MATCH_SCORE : MISMATCH_SCORE;
}

static int max3_with_direction(int diag, int up, int left, Direction *direction)
{
    /*
     * Ties are intentionally resolved as DIAG > UP > LEFT. Keeping this
     * deterministic makes sequential and future parallel outputs comparable.
     */
    int best = diag;
    *direction = DIR_DIAG;

    if (up > best) {
        best = up;
        *direction = DIR_UP;
    }

    if (left > best) {
        best = left;
        *direction = DIR_LEFT;
    }

    return best;
}

static void fill_matrix(Matrix *matrix, const char *seq_a, const char *seq_b)
{
    size_t i;
    size_t j;

    for (i = 1; i < matrix->rows; i++) {
        for (j = 1; j < matrix->cols; j++) {
            Direction direction;
            int diag = matrix->scores[cell_index(matrix, i - 1, j - 1)] +
                       score_pair(seq_a[i - 1], seq_b[j - 1]);
            int up = matrix->scores[cell_index(matrix, i - 1, j)] + GAP_PENALTY;
            int left = matrix->scores[cell_index(matrix, i, j - 1)] + GAP_PENALTY;

            matrix->scores[cell_index(matrix, i, j)] =
                max3_with_direction(diag, up, left, &direction);
            matrix->directions[cell_index(matrix, i, j)] = direction;
        }
    }
}

static int traceback(const Matrix *matrix, const char *seq_a, const char *seq_b, Alignment *alignment)
{
    size_t i = matrix->rows - 1;
    size_t j = matrix->cols - 1;
    size_t capacity = i + j + 1;
    size_t length = 0;
    char *reverse_a = NULL;
    char *reverse_b = NULL;
    size_t k;

    alignment->aligned_a = NULL;
    alignment->aligned_b = NULL;

    reverse_a = (char *)malloc(capacity * sizeof(char));
    reverse_b = (char *)malloc(capacity * sizeof(char));
    alignment->aligned_a = (char *)malloc(capacity * sizeof(char));
    alignment->aligned_b = (char *)malloc(capacity * sizeof(char));

    if (reverse_a == NULL || reverse_b == NULL ||
        alignment->aligned_a == NULL || alignment->aligned_b == NULL) {
        fprintf(stderr, "Error: failed to allocate alignment buffers.\n");
        free(reverse_a);
        free(reverse_b);
        free_alignment(alignment);
        return 0;
    }

    while (i > 0 || j > 0) {
        Direction direction = matrix->directions[cell_index(matrix, i, j)];

        if (direction == DIR_DIAG) {
            reverse_a[length] = seq_a[i - 1];
            reverse_b[length] = seq_b[j - 1];
            i--;
            j--;
        } else if (direction == DIR_UP) {
            reverse_a[length] = seq_a[i - 1];
            reverse_b[length] = '-';
            i--;
        } else if (direction == DIR_LEFT) {
            reverse_a[length] = '-';
            reverse_b[length] = seq_b[j - 1];
            j--;
        } else {
            fprintf(stderr, "Error: invalid traceback direction at cell (%lu, %lu).\n",
                    (unsigned long)i, (unsigned long)j);
            free(reverse_a);
            free(reverse_b);
            free_alignment(alignment);
            return 0;
        }

        length++;
    }

    /* Traceback builds the alignment backwards, so reverse it for printing. */
    for (k = 0; k < length; k++) {
        alignment->aligned_a[k] = reverse_a[length - 1 - k];
        alignment->aligned_b[k] = reverse_b[length - 1 - k];
    }

    alignment->aligned_a[length] = '\0';
    alignment->aligned_b[length] = '\0';

    free(reverse_a);
    free(reverse_b);
    return 1;
}

static int calculate_score(const char *aligned_a, const char *aligned_b)
{
    int total = 0;
    size_t i;

    for (i = 0; aligned_a[i] != '\0' && aligned_b[i] != '\0'; i++) {
        if (aligned_a[i] == '-' || aligned_b[i] == '-') {
            total += GAP_PENALTY;
        } else {
            total += score_pair(aligned_a[i], aligned_b[i]);
        }
    }

    return total;
}

static void print_alignment(const Alignment *alignment)
{
    printf("Aligned sequence A: %s\n", alignment->aligned_a);
    printf("Aligned sequence B: %s\n", alignment->aligned_b);
}

static char *read_sequence(const char *prompt)
{
    size_t capacity = 64;
    size_t length = 0;
    char *buffer = (char *)malloc(capacity);
    int ch;

    if (buffer == NULL) {
        fprintf(stderr, "Error: failed to allocate input buffer.\n");
        return NULL;
    }

    printf("%s", prompt);
    fflush(stdout);

    while ((ch = getchar()) != EOF && ch != '\n') {
        if (ch == '\r') {
            continue;
        }

        if (length + 1 >= capacity) {
            char *new_buffer;

            if (capacity > SIZE_MAX / 2) {
                fprintf(stderr, "Error: input sequence is too long.\n");
                free(buffer);
                return NULL;
            }

            capacity *= 2;
            new_buffer = (char *)realloc(buffer, capacity);
            if (new_buffer == NULL) {
                fprintf(stderr, "Error: failed to grow input buffer.\n");
                free(buffer);
                return NULL;
            }
            buffer = new_buffer;
        }

        buffer[length++] = (char)toupper((unsigned char)ch);
    }

    buffer[length] = '\0';
    return buffer;
}

static char *generate_random_sequence(size_t length)
{
    static const char alphabet[] = "ACGT";
    char *sequence;
    size_t i;

    if (length == SIZE_MAX) {
        fprintf(stderr, "Error: random sequence length is too large.\n");
        return NULL;
    }

    sequence = (char *)malloc(length + 1);

    if (sequence == NULL) {
        fprintf(stderr, "Error: failed to allocate random sequence.\n");
        return NULL;
    }

    for (i = 0; i < length; i++) {
        sequence[i] = alphabet[rand() % 4];
    }
    sequence[length] = '\0';

    return sequence;
}

static int parse_size(const char *text, size_t *value)
{
    char *end = NULL;
    unsigned long parsed;

    if (text == NULL || text[0] == '\0') {
        return 0;
    }

    parsed = strtoul(text, &end, 10);
    if (*end != '\0' || parsed == ULONG_MAX) {
        return 0;
    }

    *value = (size_t)parsed;
    return 1;
}

static void print_usage(const char *program_name)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s\n", program_name);
    fprintf(stderr, "  %s --random <length_a> <length_b> [seed]\n", program_name);
}

int main(int argc, char **argv)
{
    char *seq_a = NULL;
    char *seq_b = NULL;
    Matrix matrix;
    Alignment alignment;
    clock_t start_time;
    clock_t end_time;
    double elapsed_seconds;
    int matrix_score;
    int alignment_score;

    matrix.rows = 0;
    matrix.cols = 0;
    matrix.scores = NULL;
    matrix.directions = NULL;
    alignment.aligned_a = NULL;
    alignment.aligned_b = NULL;

    if (argc == 1) {
        seq_a = read_sequence("Enter sequence A: ");
        seq_b = read_sequence("Enter sequence B: ");
    } else if ((argc == 4 || argc == 5) && strcmp(argv[1], "--random") == 0) {
        size_t len_a;
        size_t len_b;
        unsigned int seed = (unsigned int)time(NULL);

        if (!parse_size(argv[2], &len_a) || !parse_size(argv[3], &len_b)) {
            fprintf(stderr, "Error: random sequence lengths must be non-negative integers.\n");
            print_usage(argv[0]);
            return EXIT_FAILURE;
        }

        if (argc == 5) {
            size_t parsed_seed;
            if (!parse_size(argv[4], &parsed_seed) || parsed_seed > UINT_MAX) {
                fprintf(stderr, "Error: seed must be an unsigned integer.\n");
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
            seed = (unsigned int)parsed_seed;
        }

        srand(seed);
        seq_a = generate_random_sequence(len_a);
        seq_b = generate_random_sequence(len_b);
    } else {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    if (seq_a == NULL || seq_b == NULL) {
        free(seq_a);
        free(seq_b);
        return EXIT_FAILURE;
    }

    if (!allocate_matrix(&matrix, strlen(seq_a), strlen(seq_b))) {
        free(seq_a);
        free(seq_b);
        return EXIT_FAILURE;
    }

    start_time = clock();
    initialize_matrix(&matrix);
    fill_matrix(&matrix, seq_a, seq_b);

    if (!traceback(&matrix, seq_a, seq_b, &alignment)) {
        free_matrix(&matrix);
        free(seq_a);
        free(seq_b);
        return EXIT_FAILURE;
    }
    end_time = clock();

    matrix_score = matrix.scores[cell_index(&matrix, matrix.rows - 1, matrix.cols - 1)];
    alignment_score = calculate_score(alignment.aligned_a, alignment.aligned_b);
    elapsed_seconds = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    if (matrix_score != alignment_score) {
        fprintf(stderr, "Error: traceback score does not match DP matrix score.\n");
        free_alignment(&alignment);
        free_matrix(&matrix);
        free(seq_a);
        free(seq_b);
        return EXIT_FAILURE;
    }

    printf("\nSequence A: %s\n", seq_a);
    printf("Sequence B: %s\n", seq_b);
    printf("Matrix dimensions: %lu x %lu\n",
           (unsigned long)matrix.rows, (unsigned long)matrix.cols);
    printf("Alignment score: %d\n", matrix_score);
    print_alignment(&alignment);
    printf("Execution time: %.6f seconds\n", elapsed_seconds);

    free_alignment(&alignment);
    free_matrix(&matrix);
    free(seq_a);
    free(seq_b);

    return EXIT_SUCCESS;
}
