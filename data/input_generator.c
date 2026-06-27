#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const size_t DEFAULT_LENGTHS[] = {100, 500, 1000, 5000, 10000, 50000};
static const size_t DEFAULT_LENGTH_COUNT = sizeof(DEFAULT_LENGTHS) / sizeof(DEFAULT_LENGTHS[0]);
static const char DNA_ALPHABET[] = {'A', 'C', 'G', 'T'};

static int parse_positive_int(const char *text, int *value)
{
    char *end = NULL;
    long parsed;

    errno = 0;
    parsed = strtol(text, &end, 10);

    if (text == NULL || text[0] == '\0' || *end != '\0' ||
        errno == ERANGE || parsed <= 0 || parsed > INT_MAX) {
        return 0;
    }

    *value = (int)parsed;
    return 1;
}

static int parse_seed(const char *text, unsigned int *value)
{
    char *end = NULL;
    unsigned long parsed;

    errno = 0;
    parsed = strtoul(text, &end, 10);

    if (text == NULL || text[0] == '\0' || *end != '\0' ||
        errno == ERANGE || parsed > UINT_MAX) {
        return 0;
    }

    *value = (unsigned int)parsed;
    return 1;
}

static void write_random_sequence(FILE *file, size_t length)
{
    size_t i;

    for (i = 0; i < length; i++) {
        fputc(DNA_ALPHABET[rand() % 4], file);
    }
    fputc('\n', file);
}

static int generate_file(const char *output_dir, size_t length, int pair_count)
{
    char path[512];
    FILE *file;
    int pair;

    if (snprintf(path, sizeof(path), "%s/input_len_%lu_pairs_%d.txt",
                 output_dir, (unsigned long)length, pair_count) >= (int)sizeof(path)) {
        fprintf(stderr, "Error: output path is too long.\n");
        return 0;
    }

    file = fopen(path, "w");
    if (file == NULL) {
        fprintf(stderr, "Error: cannot create %s\n", path);
        return 0;
    }

    fprintf(file, "%d\n", pair_count);

    for (pair = 0; pair < pair_count; pair++) {
        write_random_sequence(file, length);
        write_random_sequence(file, length);
    }

    fclose(file);
    printf("Generated %s\n", path);
    return 1;
}

static void print_usage(const char *program_name)
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  %s <output_dir> <pairs_per_size> [seed]\n", program_name);
    fprintf(stderr, "\nGenerates lengths: 100, 500, 1000, 5000, 10000, 50000\n");
}

int main(int argc, char **argv)
{
    const char *output_dir;
    int pair_count;
    unsigned int seed = (unsigned int)time(NULL);
    size_t i;

    if (argc != 3 && argc != 4) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    output_dir = argv[1];

    if (!parse_positive_int(argv[2], &pair_count)) {
        fprintf(stderr, "Error: pairs_per_size must be a positive integer.\n");
        return EXIT_FAILURE;
    }

    if (argc == 4 && !parse_seed(argv[3], &seed)) {
        fprintf(stderr, "Error: seed must be an unsigned integer.\n");
        return EXIT_FAILURE;
    }

    srand(seed);
    printf("Seed: %u\n", seed);
    printf("Pairs per size: %d\n", pair_count);

    for (i = 0; i < DEFAULT_LENGTH_COUNT; i++) {
        if (!generate_file(output_dir, DEFAULT_LENGTHS[i], pair_count)) {
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}
