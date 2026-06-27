#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_PENALTY -1

#define MASTER_RANK 0
#define TAG_WORK 100
#define TAG_RESULT 200
#define TAG_STOP 300

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
    char *seq_a;
    char *seq_b;
} SequencePair;

typedef struct {
    int task_id;
    int ok;
    int score;
    double worker_time;
    char *aligned_a;
    char *aligned_b;
} AlignmentResult;

static size_t cell_index(const Matrix *matrix, size_t row, size_t col);
static void free_matrix(Matrix *matrix);
static void free_sequence_pairs(SequencePair *pairs, int count);
static void free_alignment_result(AlignmentResult *result);
static int allocate_matrix(Matrix *matrix, size_t len_a, size_t len_b);
static void initialize_matrix(Matrix *matrix);
static int score_pair(char a, char b);
static int max3_with_direction(int diag, int up, int left, Direction *direction);
static void fill_matrix(Matrix *matrix, const char *seq_a, const char *seq_b);
static int traceback(const Matrix *matrix, const char *seq_a, const char *seq_b,
                     char **aligned_a, char **aligned_b);
static int calculate_score(const char *aligned_a, const char *aligned_b);
static int run_needleman_wunsch(const char *seq_a, const char *seq_b, AlignmentResult *result);
static char *read_line(FILE *stream);
static void uppercase_sequence(char *sequence);
static int read_sequence_pairs(SequencePair **pairs_out, int *count_out);
static int send_work(int worker_rank, int task_id, const SequencePair *pair);
static int send_stop(int worker_rank);
static int receive_work(int *task_id, char **seq_a, char **seq_b, MPI_Status *status);
static int send_result(int master_rank, const AlignmentResult *result);
static int receive_result(AlignmentResult *result, MPI_Status *status);
static void master_process(int world_size);
static void worker_process(int rank);

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

static void free_sequence_pairs(SequencePair *pairs, int count)
{
    int i;

    if (pairs == NULL) {
        return;
    }

    for (i = 0; i < count; i++) {
        free(pairs[i].seq_a);
        free(pairs[i].seq_b);
    }
    free(pairs);
}

static void free_alignment_result(AlignmentResult *result)
{
    if (result == NULL) {
        return;
    }

    free(result->aligned_a);
    free(result->aligned_b);
    result->aligned_a = NULL;
    result->aligned_b = NULL;
}

static int allocate_matrix(Matrix *matrix, size_t len_a, size_t len_b)
{
    size_t total_cells;

    matrix->rows = len_a + 1;
    matrix->cols = len_b + 1;
    matrix->scores = NULL;
    matrix->directions = NULL;

    if (len_a == SIZE_MAX || len_b == SIZE_MAX ||
        len_a > (size_t)INT_MAX || len_b > (size_t)INT_MAX) {
        fprintf(stderr, "Error: sequence length exceeds supported range.\n");
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

    matrix->scores[cell_index(matrix, 0, 0)] = 0;
    matrix->directions[cell_index(matrix, 0, 0)] = DIR_NONE;

    for (i = 1; i < matrix->rows; i++) {
        matrix->scores[cell_index(matrix, i, 0)] = (int)i * GAP_PENALTY;
        matrix->directions[cell_index(matrix, i, 0)] = DIR_UP;
    }

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
     * Same tie policy as the sequential baseline: DIAG > UP > LEFT.
     * MPI distributes independent jobs only; it does not change scoring.
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

static int traceback(const Matrix *matrix, const char *seq_a, const char *seq_b,
                     char **aligned_a, char **aligned_b)
{
    size_t i = matrix->rows - 1;
    size_t j = matrix->cols - 1;
    size_t capacity = i + j + 1;
    size_t length = 0;
    char *reverse_a = NULL;
    char *reverse_b = NULL;
    size_t k;

    *aligned_a = NULL;
    *aligned_b = NULL;

    reverse_a = (char *)malloc(capacity);
    reverse_b = (char *)malloc(capacity);
    *aligned_a = (char *)malloc(capacity);
    *aligned_b = (char *)malloc(capacity);

    if (reverse_a == NULL || reverse_b == NULL ||
        *aligned_a == NULL || *aligned_b == NULL) {
        fprintf(stderr, "Error: failed to allocate alignment buffers.\n");
        free(reverse_a);
        free(reverse_b);
        free(*aligned_a);
        free(*aligned_b);
        *aligned_a = NULL;
        *aligned_b = NULL;
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
            fprintf(stderr, "Error: invalid traceback direction.\n");
            free(reverse_a);
            free(reverse_b);
            free(*aligned_a);
            free(*aligned_b);
            *aligned_a = NULL;
            *aligned_b = NULL;
            return 0;
        }

        length++;
    }

    for (k = 0; k < length; k++) {
        (*aligned_a)[k] = reverse_a[length - 1 - k];
        (*aligned_b)[k] = reverse_b[length - 1 - k];
    }

    (*aligned_a)[length] = '\0';
    (*aligned_b)[length] = '\0';

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

static int run_needleman_wunsch(const char *seq_a, const char *seq_b, AlignmentResult *result)
{
    Matrix matrix;
    int traceback_score;
    size_t len_a = strlen(seq_a);
    size_t len_b = strlen(seq_b);

    matrix.rows = 0;
    matrix.cols = 0;
    matrix.scores = NULL;
    matrix.directions = NULL;

    result->ok = 0;
    result->score = 0;
    result->aligned_a = NULL;
    result->aligned_b = NULL;

    if (!allocate_matrix(&matrix, len_a, len_b)) {
        return 0;
    }

    initialize_matrix(&matrix);
    fill_matrix(&matrix, seq_a, seq_b);

    if (!traceback(&matrix, seq_a, seq_b, &result->aligned_a, &result->aligned_b)) {
        free_matrix(&matrix);
        return 0;
    }

    result->score = matrix.scores[cell_index(&matrix, matrix.rows - 1, matrix.cols - 1)];
    traceback_score = calculate_score(result->aligned_a, result->aligned_b);

    if (result->score != traceback_score) {
        fprintf(stderr, "Error: traceback score does not match DP matrix score.\n");
        free_alignment_result(result);
        free_matrix(&matrix);
        return 0;
    }

    result->ok = 1;
    free_matrix(&matrix);
    return 1;
}

static char *read_line(FILE *stream)
{
    size_t capacity = 64;
    size_t length = 0;
    char *buffer = (char *)malloc(capacity);
    int ch;

    if (buffer == NULL) {
        return NULL;
    }

    while ((ch = fgetc(stream)) != EOF && ch != '\n') {
        if (ch == '\r') {
            continue;
        }

        if (length + 1 >= capacity) {
            char *new_buffer;

            if (capacity > SIZE_MAX / 2) {
                free(buffer);
                return NULL;
            }

            capacity *= 2;
            new_buffer = (char *)realloc(buffer, capacity);
            if (new_buffer == NULL) {
                free(buffer);
                return NULL;
            }
            buffer = new_buffer;
        }

        buffer[length++] = (char)ch;
    }

    if (ch == EOF && length == 0) {
        free(buffer);
        return NULL;
    }

    buffer[length] = '\0';
    return buffer;
}

static void uppercase_sequence(char *sequence)
{
    size_t i;

    for (i = 0; sequence[i] != '\0'; i++) {
        sequence[i] = (char)toupper((unsigned char)sequence[i]);
    }
}

static int read_sequence_pairs(SequencePair **pairs_out, int *count_out)
{
    char *count_line = NULL;
    char *end = NULL;
    long count_long;
    int count;
    int i;
    SequencePair *pairs = NULL;

    printf("Enter number of sequence pairs: ");
    fflush(stdout);

    count_line = read_line(stdin);
    if (count_line == NULL) {
        fprintf(stderr, "Error: failed to read number of sequence pairs.\n");
        return 0;
    }

    errno = 0;
    count_long = strtol(count_line, &end, 10);

    if (*end != '\0' || errno == ERANGE || count_long <= 0 || count_long > INT_MAX) {
        fprintf(stderr, "Error: sequence pair count must be a positive integer.\n");
        free(count_line);
        return 0;
    }
    free(count_line);

    count = (int)count_long;
    pairs = (SequencePair *)calloc((size_t)count, sizeof(SequencePair));
    if (pairs == NULL) {
        fprintf(stderr, "Error: failed to allocate sequence pair list.\n");
        return 0;
    }

    for (i = 0; i < count; i++) {
        printf("Enter sequence A for pair %d: ", i);
        fflush(stdout);
        pairs[i].seq_a = read_line(stdin);

        printf("Enter sequence B for pair %d: ", i);
        fflush(stdout);
        pairs[i].seq_b = read_line(stdin);

        if (pairs[i].seq_a == NULL || pairs[i].seq_b == NULL) {
            fprintf(stderr, "Error: failed to read sequence pair %d.\n", i);
            free_sequence_pairs(pairs, count);
            return 0;
        }

        uppercase_sequence(pairs[i].seq_a);
        uppercase_sequence(pairs[i].seq_b);
    }

    *pairs_out = pairs;
    *count_out = count;
    return 1;
}

static int send_work(int worker_rank, int task_id, const SequencePair *pair)
{
    int header[3];
    size_t len_a = strlen(pair->seq_a);
    size_t len_b = strlen(pair->seq_b);

    if (len_a > (size_t)INT_MAX - 1 || len_b > (size_t)INT_MAX - 1) {
        fprintf(stderr, "Error: sequence is too long to send with int MPI counts.\n");
        return 0;
    }

    header[0] = task_id;
    header[1] = (int)len_a + 1;
    header[2] = (int)len_b + 1;

    MPI_Send(header, 3, MPI_INT, worker_rank, TAG_WORK, MPI_COMM_WORLD);
    MPI_Send(pair->seq_a, header[1], MPI_CHAR, worker_rank, TAG_WORK, MPI_COMM_WORLD);
    MPI_Send(pair->seq_b, header[2], MPI_CHAR, worker_rank, TAG_WORK, MPI_COMM_WORLD);
    return 1;
}

static int send_stop(int worker_rank)
{
    int header[3] = {-1, 0, 0};

    MPI_Send(header, 3, MPI_INT, worker_rank, TAG_STOP, MPI_COMM_WORLD);
    return 1;
}

static int receive_work(int *task_id, char **seq_a, char **seq_b, MPI_Status *status)
{
    int header[3];

    *seq_a = NULL;
    *seq_b = NULL;

    MPI_Recv(header, 3, MPI_INT, MASTER_RANK, MPI_ANY_TAG, MPI_COMM_WORLD, status);

    if (status->MPI_TAG == TAG_STOP) {
        *task_id = -1;
        return 1;
    }

    if (status->MPI_TAG != TAG_WORK || header[1] <= 0 || header[2] <= 0) {
        fprintf(stderr, "Error: worker received invalid work message.\n");
        return 0;
    }

    *task_id = header[0];
    *seq_a = (char *)malloc((size_t)header[1]);
    *seq_b = (char *)malloc((size_t)header[2]);

    if (*seq_a == NULL || *seq_b == NULL) {
        fprintf(stderr, "Error: worker failed to allocate received sequences.\n");
        free(*seq_a);
        free(*seq_b);
        *seq_a = NULL;
        *seq_b = NULL;
        return 0;
    }

    MPI_Recv(*seq_a, header[1], MPI_CHAR, MASTER_RANK, TAG_WORK, MPI_COMM_WORLD, status);
    MPI_Recv(*seq_b, header[2], MPI_CHAR, MASTER_RANK, TAG_WORK, MPI_COMM_WORLD, status);
    return 1;
}

static int send_result(int master_rank, const AlignmentResult *result)
{
    int header[5];
    int len_a = 0;
    int len_b = 0;

    if (result->ok) {
        size_t size_a = strlen(result->aligned_a) + 1;
        size_t size_b = strlen(result->aligned_b) + 1;

        if (size_a > (size_t)INT_MAX || size_b > (size_t)INT_MAX) {
            return 0;
        }

        len_a = (int)size_a;
        len_b = (int)size_b;
    }

    header[0] = result->task_id;
    header[1] = result->ok;
    header[2] = result->score;
    header[3] = len_a;
    header[4] = len_b;

    MPI_Send(header, 5, MPI_INT, master_rank, TAG_RESULT, MPI_COMM_WORLD);
    MPI_Send((void *)&result->worker_time, 1, MPI_DOUBLE, master_rank, TAG_RESULT, MPI_COMM_WORLD);

    if (result->ok) {
        MPI_Send(result->aligned_a, len_a, MPI_CHAR, master_rank, TAG_RESULT, MPI_COMM_WORLD);
        MPI_Send(result->aligned_b, len_b, MPI_CHAR, master_rank, TAG_RESULT, MPI_COMM_WORLD);
    }

    return 1;
}

static int receive_result(AlignmentResult *result, MPI_Status *status)
{
    int header[5];

    result->task_id = -1;
    result->ok = 0;
    result->score = 0;
    result->worker_time = 0.0;
    result->aligned_a = NULL;
    result->aligned_b = NULL;

    MPI_Recv(header, 5, MPI_INT, MPI_ANY_SOURCE, TAG_RESULT, MPI_COMM_WORLD, status);
    MPI_Recv(&result->worker_time, 1, MPI_DOUBLE, status->MPI_SOURCE,
             TAG_RESULT, MPI_COMM_WORLD, status);

    result->task_id = header[0];
    result->ok = header[1];
    result->score = header[2];

    if (!result->ok) {
        return 1;
    }

    if (header[3] <= 0 || header[4] <= 0) {
        fprintf(stderr, "Error: master received invalid result lengths.\n");
        return 0;
    }

    result->aligned_a = (char *)malloc((size_t)header[3]);
    result->aligned_b = (char *)malloc((size_t)header[4]);

    if (result->aligned_a == NULL || result->aligned_b == NULL) {
        fprintf(stderr, "Error: master failed to allocate result buffers.\n");
        free_alignment_result(result);
        return 0;
    }

    MPI_Recv(result->aligned_a, header[3], MPI_CHAR, status->MPI_SOURCE,
             TAG_RESULT, MPI_COMM_WORLD, status);
    MPI_Recv(result->aligned_b, header[4], MPI_CHAR, status->MPI_SOURCE,
             TAG_RESULT, MPI_COMM_WORLD, status);
    return 1;
}

static void master_process(int world_size)
{
    SequencePair *pairs = NULL;
    AlignmentResult *results = NULL;
    int pair_count = 0;
    int next_task = 0;
    int completed = 0;
    int active_workers = 0;
    int worker;
    double start_time;
    double end_time;

    if (!read_sequence_pairs(&pairs, &pair_count)) {
        for (worker = 1; worker < world_size; worker++) {
            send_stop(worker);
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    results = (AlignmentResult *)calloc((size_t)pair_count, sizeof(AlignmentResult));
    if (results == NULL) {
        fprintf(stderr, "Error: failed to allocate result list.\n");
        free_sequence_pairs(pairs, pair_count);
        for (worker = 1; worker < world_size; worker++) {
            send_stop(worker);
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    start_time = MPI_Wtime();

    if (world_size == 1) {
        int i;

        /*
         * With only rank 0 there are no workers. This fallback keeps the
         * program useful, but real MPI speedup requires multiple processes.
         */
        for (i = 0; i < pair_count; i++) {
            results[i].task_id = i;
            results[i].worker_time = MPI_Wtime();
            run_needleman_wunsch(pairs[i].seq_a, pairs[i].seq_b, &results[i]);
            results[i].worker_time = MPI_Wtime() - results[i].worker_time;
            completed++;
        }
    } else {
        /*
         * Initial task farm seeding: give one task to each worker until we run
         * out of workers or tasks. Later tasks are sent dynamically when a
         * worker reports that it is free.
         */
        for (worker = 1; worker < world_size && next_task < pair_count; worker++) {
            send_work(worker, next_task, &pairs[next_task]);
            next_task++;
            active_workers++;
        }

        for (; worker < world_size; worker++) {
            send_stop(worker);
        }

        while (completed < pair_count) {
            MPI_Status status;
            AlignmentResult received;
            int source;

            if (!receive_result(&received, &status)) {
                fprintf(stderr, "Error: failed to receive worker result.\n");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            source = status.MPI_SOURCE;

            if (received.task_id < 0 || received.task_id >= pair_count) {
                fprintf(stderr, "Error: master received invalid task id %d.\n",
                        received.task_id);
                free_alignment_result(&received);
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }

            results[received.task_id] = received;
            completed++;

            if (next_task < pair_count) {
                send_work(source, next_task, &pairs[next_task]);
                next_task++;
            } else {
                send_stop(source);
                active_workers--;
            }
        }

        (void)active_workers;
    }

    end_time = MPI_Wtime();

    printf("\nMPI Needleman-Wunsch Results\n");
    printf("Processes: %d\n", world_size);
    printf("Sequence pairs: %d\n", pair_count);
    printf("Total execution time: %.6f seconds\n", end_time - start_time);

    for (worker = 0; worker < pair_count; worker++) {
        printf("\nPair %d\n", worker);
        printf("Sequence A: %s\n", pairs[worker].seq_a);
        printf("Sequence B: %s\n", pairs[worker].seq_b);

        if (!results[worker].ok) {
            printf("Alignment failed on worker.\n");
            continue;
        }

        printf("Alignment score: %d\n", results[worker].score);
        printf("Aligned sequence A: %s\n", results[worker].aligned_a);
        printf("Aligned sequence B: %s\n", results[worker].aligned_b);
        printf("Worker computation time: %.6f seconds\n", results[worker].worker_time);
    }

    for (worker = 0; worker < pair_count; worker++) {
        free_alignment_result(&results[worker]);
    }
    free(results);
    free_sequence_pairs(pairs, pair_count);
}

static void worker_process(int rank)
{
    (void)rank;

    while (1) {
        MPI_Status status;
        AlignmentResult result;
        char *seq_a = NULL;
        char *seq_b = NULL;
        int task_id;
        double start_time;

        if (!receive_work(&task_id, &seq_a, &seq_b, &status)) {
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        if (task_id < 0) {
            break;
        }

        result.task_id = task_id;
        result.ok = 0;
        result.score = 0;
        result.worker_time = 0.0;
        result.aligned_a = NULL;
        result.aligned_b = NULL;

        start_time = MPI_Wtime();
        run_needleman_wunsch(seq_a, seq_b, &result);
        result.task_id = task_id;
        result.worker_time = MPI_Wtime() - start_time;

        if (!send_result(MASTER_RANK, &result)) {
            fprintf(stderr, "Error: worker failed to send result.\n");
            free(seq_a);
            free(seq_b);
            free_alignment_result(&result);
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

        free(seq_a);
        free(seq_b);
        free_alignment_result(&result);
    }
}

int main(int argc, char **argv)
{
    int rank;
    int world_size;

    (void)argv;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (rank == MASTER_RANK) {
        master_process(world_size);
    } else {
        worker_process(rank);
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
