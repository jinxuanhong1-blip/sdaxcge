/* Stream a genes-by-cells TSV (GEO UMI text) and emit nonzeros for selected cells.
 *
 * stdin:  TSV, row 1 = header (gene-column name, then cell ids), later rows = gene, values
 * argv:   keep_cols.txt  out_prefix
 * keep_cols.txt: one 0-based data-column index per line (column 0 is the first cell)
 *
 * Writes:
 *   out_prefix.genes.txt     one gene symbol per data row
 *   out_prefix.coo.bin       records of int32 gene_i, int32 compact_j, float32 value
 *   out_prefix.libsize.f64   float64 library size per kept cell (all genes)
 *   out_prefix.shape.txt     n_genes n_kept n_nnz
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

static int32_t *load_map(const char *path, int *n_cells_out, int *n_kept_out) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s: %s\n", path, strerror(errno));
        exit(1);
    }
    int cap = 1024, n = 0;
    int32_t *cols = malloc((size_t)cap * sizeof(int32_t));
    if (!cols) exit(1);
    int32_t maxc = -1;
    char buf[64];
    while (fgets(buf, sizeof buf, f)) {
        if (buf[0] == '\n' || buf[0] == '\0') continue;
        int32_t c = (int32_t)atoi(buf);
        if (n >= cap) {
            cap *= 2;
            cols = realloc(cols, (size_t)cap * sizeof(int32_t));
            if (!cols) exit(1);
        }
        cols[n++] = c;
        if (c > maxc) maxc = c;
    }
    fclose(f);
    if (n == 0 || maxc < 0) {
        fprintf(stderr, "empty keep list\n");
        exit(1);
    }
    int32_t *map = malloc(((size_t)maxc + 1) * sizeof(int32_t));
    if (!map) exit(1);
    for (int32_t i = 0; i <= maxc; i++) map[i] = -1;
    for (int i = 0; i < n; i++) {
        if (map[cols[i]] != -1) {
            fprintf(stderr, "duplicate keep col %d\n", cols[i]);
            exit(1);
        }
        map[cols[i]] = i;
    }
    free(cols);
    *n_cells_out = (int)maxc + 1;
    *n_kept_out = n;
    return map;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: extract_subset keep_cols.txt out_prefix < matrix.tsv\n");
        return 2;
    }
    int n_cells = 0, n_kept = 0;
    int32_t *map = load_map(argv[1], &n_cells, &n_kept);
    double *lib = calloc((size_t)n_kept, sizeof(double));
    if (!lib) return 1;

    char genes_path[4096], coo_path[4096], lib_path[4096], shape_path[4096];
    snprintf(genes_path, sizeof genes_path, "%s.genes.txt", argv[2]);
    snprintf(coo_path, sizeof coo_path, "%s.coo.bin", argv[2]);
    snprintf(lib_path, sizeof lib_path, "%s.libsize.f64", argv[2]);
    snprintf(shape_path, sizeof shape_path, "%s.shape.txt", argv[2]);
    FILE *fg = fopen(genes_path, "w");
    FILE *fc = fopen(coo_path, "wb");
    if (!fg || !fc) {
        fprintf(stderr, "cannot open outputs\n");
        return 1;
    }

    size_t cap = 1 << 20;
    char *line = malloc(cap);
    if (!line) return 1;

    /* header */
    ssize_t len = getline(&line, &cap, stdin);
    if (len < 0) {
        fprintf(stderr, "empty input\n");
        return 1;
    }
    int header_cells = 0;
    for (ssize_t i = 0; i < len; i++) if (line[i] == '\t') header_cells++;
    if (header_cells < n_cells) {
        fprintf(stderr, "header cells %d < map span %d\n", header_cells, n_cells);
        return 1;
    }

    int32_t gene_i = 0;
    int64_t nnz = 0;
    while ((len = getline(&line, &cap, stdin)) != -1) {
        if (len > 0 && line[len - 1] == '\n') line[--len] = '\0';
        if (len > 0 && line[len - 1] == '\r') line[--len] = '\0';
        char *tab = strchr(line, '\t');
        if (!tab) {
            fprintf(fg, "%s\n", line);
            gene_i++;
            continue;
        }
        *tab = '\0';
        fprintf(fg, "%s\n", line);
        char *p = tab + 1;
        int col = 0;
        while (*p) {
            /* parse a non-negative number; GEO UMIs are integers or 0.0 */
            double v = 0.0;
            if (*p == '0' && (p[1] == '\t' || p[1] == '\0')) {
                p++;
            } else {
                char *end = p;
                v = strtod(p, &end);
                if (end == p) {
                    fprintf(stderr, "parse fail gene %d col %d\n", gene_i, col);
                    return 1;
                }
                p = end;
            }
            if (*p == '\t') p++;
            if (col < n_cells) {
                int32_t j = map[col];
                if (j >= 0) {
                    if (v != 0.0) {
                        float fv = (float)v;
                        if (fwrite(&gene_i, 4, 1, fc) != 1) return 1;
                        if (fwrite(&j, 4, 1, fc) != 1) return 1;
                        if (fwrite(&fv, 4, 1, fc) != 1) return 1;
                        nnz++;
                    }
                    lib[j] += v;
                }
            }
            col++;
        }
        gene_i++;
        if (gene_i % 2000 == 0) {
            fprintf(stderr, "  streamed %d genes nnz %lld\n", gene_i, (long long)nnz);
        }
    }
    fclose(fg);
    fclose(fc);

    FILE *fl = fopen(lib_path, "wb");
    if (!fl || fwrite(lib, sizeof(double), (size_t)n_kept, fl) != (size_t)n_kept) return 1;
    fclose(fl);
    FILE *fs = fopen(shape_path, "w");
    fprintf(fs, "%d\t%d\t%lld\n", gene_i, n_kept, (long long)nnz);
    fclose(fs);
    fprintf(stderr, "done genes %d kept %d nnz %lld\n", gene_i, n_kept, (long long)nnz);
    free(line);
    free(map);
    free(lib);
    return 0;
}
