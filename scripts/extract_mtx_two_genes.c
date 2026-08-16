/* Stream Matrix Market coordinate MTX on stdin; keep columns 1057 and 11797. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    char line[4096];
    long long n_rows = 0, n_cols = 0, n_nz = 0;
    long long scanned = 0, kept = 0;
    int row, col;
    double val;

    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '%') {
            continue;
        }
        if (sscanf(line, "%lld %lld %lld", &n_rows, &n_cols, &n_nz) != 3) {
            fprintf(stderr, "bad MTX size line: %s", line);
            return 1;
        }
        break;
    }
    fprintf(stderr, "MTX header %lld x %lld nnz=%lld\n", n_rows, n_cols, n_nz);
    fputs("cell_index_1based\tgene\tcount\n", stdout);
    while (fscanf(stdin, "%d %d %lf", &row, &col, &val) == 3) {
        scanned++;
        if (col == 1057) {
            printf("%d\tTACSTD2\t%.10g\n", row, val);
            kept++;
        } else if (col == 11797) {
            printf("%d\tCLDN4\t%.10g\n", row, val);
            kept++;
        }
        if (scanned % 50000000LL == 0) {
            fprintf(stderr, "scanned %lld kept %lld\n", scanned, kept);
        }
    }
    fprintf(stderr, "done scanned=%lld kept=%lld declared=%lld\n", scanned, kept, n_nz);
    if (scanned != n_nz) {
        fprintf(stderr, "WARNING: scanned nnz != declared nnz\n");
        return 2;
    }
    return 0;
}
