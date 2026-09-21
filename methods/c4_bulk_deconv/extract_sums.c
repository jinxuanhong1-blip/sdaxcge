/* Stream a genes-by-cells gzip TSV and pool UMI sums by cell class.
   class.bin: int32 n_fields (little endian), then int16 class id per field
   (-1 = ignore). Field 0 is the gene symbol. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define NCLASS 8
#define MAX_LINE (32 * 1024 * 1024)

static int16_t *load_class(const char *path, int32_t *n_out) {
    FILE *f = fopen(path, "rb");
    int32_t n = 0;
    if (!f || fread(&n, 4, 1, f) != 1 || n < 2 || n > 5000000) {
        fprintf(stderr, "bad class file\n");
        exit(1);
    }
    int16_t *cls = malloc((size_t)n * sizeof(int16_t));
    if (!cls || fread(cls, sizeof(int16_t), (size_t)n, f) != (size_t)n) {
        fprintf(stderr, "short class file\n");
        exit(1);
    }
    fclose(f);
    *n_out = n;
    return cls;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: extract_sums matrix.txt.gz class.bin out.tsv\n");
        return 2;
    }
    int32_t n_fields = 0;
    int16_t *cls = load_class(argv[2], &n_fields);
    gzFile gz = gzopen(argv[1], "rb");
    if (!gz) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        return 1;
    }
    FILE *out = fopen(argv[3], "w");
    if (!out) return 1;
    char *line = malloc(MAX_LINE);
    if (!line) return 1;
    double lib[NCLASS];
    memset(lib, 0, sizeof(lib));
    long n_genes = 0;
    int header = 1;
    while (gzgets(gz, line, MAX_LINE)) {
        size_t len = strlen(line);
        if (len + 1 >= MAX_LINE) {
            fprintf(stderr, "line exceeds buffer\n");
            return 1;
        }
        if (line[len - 1] == '\n') line[--len] = 0;
        if (len && line[len - 1] == '\r') line[--len] = 0;
        if (header) { header = 0; continue; }
        char *p = line;
        char *tab = strchr(p, '\t');
        if (!tab) continue;
        *tab = 0;
        fputs(p, out);
        double acc[NCLASS];
        memset(acc, 0, sizeof(acc));
        int col = 1;
        p = tab + 1;
        while (*p && col < n_fields) {
            char *next = p;
            while (*next && *next != '\t') next++;
            int16_t k = cls[col];
            if (k >= 0 && k < NCLASS && next > p) {
                double v = atof(p);
                if (v > 0) {
                    acc[k] += v;
                    lib[k] += v;
                }
            }
            col++;
            if (*next == '\t') p = next + 1;
            else break;
        }
        for (int k = 0; k < NCLASS; k++) fprintf(out, "\t%.6f", acc[k]);
        fputc('\n', out);
        n_genes++;
        if (n_genes % 2000 == 0) fprintf(stderr, "genes %ld\n", n_genes);
    }
    gzclose(gz);
    fclose(out);
    fprintf(stderr, "GENES %ld\n", n_genes);
    fprintf(stderr, "LIB");
    for (int k = 0; k < NCLASS; k++) fprintf(stderr, " %.6f", lib[k]);
    fprintf(stderr, "\n");
    free(line);
    free(cls);
    return 0;
}
