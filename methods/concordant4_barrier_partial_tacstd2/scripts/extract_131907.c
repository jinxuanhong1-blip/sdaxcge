/* Stream GSE131907 genes x cells UMI text and keep a column subset.
   Writes library size (sum of all genes) and one float32 vector per requested gene.
   Column map is int32, -1 = drop, else compact index. */
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define GENE_MAX 64
#define GNAME 64

static int lookup(char names[][GNAME], int n, const char *gene) {
    for (int i = 0; i < n; i++) {
        if (strcmp(names[i], gene) == 0) return i;
    }
    return -1;
}

static char *upper_dup(const char *s) {
    size_t n = strlen(s);
    char *o = malloc(n + 1);
    for (size_t i = 0; i < n; i++) o[i] = (char)toupper((unsigned char)s[i]);
    o[n] = 0;
    return o;
}

int main(int argc, char **argv) {
    if (argc != 6) {
        fprintf(stderr, "usage: extract_131907 <matrix.txt.gz> <genes.txt> <keep_idx.i32> <n_cells> <out_dir>\n");
        return 2;
    }
    const char *mtx = argv[1];
    const char *genes_path = argv[2];
    const char *idx_path = argv[3];
    int n_cells = atoi(argv[4]);
    const char *out_dir = argv[5];
    if (n_cells <= 0) return 2;

    char names[GENE_MAX][GNAME];
    int n_genes = 0;
    FILE *gf = fopen(genes_path, "r");
    if (!gf) { perror(genes_path); return 1; }
    char buf[GNAME];
    while (fgets(buf, sizeof buf, gf)) {
        size_t n = strlen(buf);
        while (n && (buf[n - 1] == '\n' || buf[n - 1] == '\r')) buf[--n] = 0;
        if (!n) continue;
        for (size_t i = 0; i < n; i++) buf[i] = (char)toupper((unsigned char)buf[i]);
        if (n_genes >= GENE_MAX) { fprintf(stderr, "too many genes\n"); return 1; }
        strncpy(names[n_genes], buf, GNAME - 1);
        names[n_genes][GNAME - 1] = 0;
        n_genes++;
    }
    fclose(gf);

    FILE *inf = fopen(idx_path, "rb");
    if (!inf) { perror(idx_path); return 1; }
    int32_t *map = malloc((size_t)n_cells * sizeof(int32_t));
    if (fread(map, sizeof(int32_t), (size_t)n_cells, inf) != (size_t)n_cells) {
        fprintf(stderr, "short keep map\n");
        return 1;
    }
    fclose(inf);
    int n_kept = 0;
    for (int i = 0; i < n_cells; i++) if (map[i] >= n_kept) n_kept = map[i] + 1;
    fprintf(stderr, "genes %d kept cells %d / %d\n", n_genes, n_kept, n_cells);

    double *lib = calloc((size_t)n_kept, sizeof(double));
    float **gmat = calloc((size_t)n_genes, sizeof(float *));
    for (int i = 0; i < n_genes; i++) gmat[i] = calloc((size_t)n_kept, sizeof(float));
    int *seen = calloc((size_t)n_genes, sizeof(int));

    gzFile gz = gzopen(mtx, "rb");
    if (!gz) { perror(mtx); return 1; }
    gzbuffer(gz, 1 << 20);

    char gene[GNAME];
    int gene_len = 0;
    int gene_slot = -1;
    int field = 0;
    int acc = 0;
    int in_field_nonnum = 0;
    int past_header = 0;
    long n_lines = 0;
    unsigned char *block = malloc(1 << 22);
    int nread;
    while ((nread = gzread(gz, block, 1 << 22)) > 0) {
        for (int bi = 0; bi < nread; bi++) {
            int c = block[bi];
            if (!past_header) {
                if (c == '\n') past_header = 1;
                continue;
            }
            if (c == '\t' || c == '\n') {
                if (field == 0) {
                    gene[gene_len] = 0;
                    gene_slot = lookup(names, n_genes, gene);
                } else {
                    int col = field - 1;
                    if (col >= 0 && col < n_cells && !in_field_nonnum) {
                        int slot = map[col];
                        if (slot >= 0) {
                            lib[slot] += acc;
                            if (gene_slot >= 0) {
                                gmat[gene_slot][slot] = (float)acc;
                                seen[gene_slot] = 1;
                            }
                        }
                    }
                }
                if (c == '\n') {
                    field = 0;
                    gene_len = 0;
                    gene_slot = -1;
                    n_lines++;
                    if (n_lines % 4000 == 0) fprintf(stderr, "  lines %ld\n", n_lines);
                } else {
                    field++;
                }
                acc = 0;
                in_field_nonnum = 0;
            } else if (field == 0) {
                if (gene_len < GNAME - 1) gene[gene_len++] = (char)toupper((unsigned char)c);
            } else if (c >= '0' && c <= '9' && !in_field_nonnum) {
                acc = acc * 10 + (c - '0');
            } else if (c != ' ') {
                in_field_nonnum = 1;
            }
        }
    }
    if (nread < 0) {
        fprintf(stderr, "gzread error\n");
        return 1;
    }
    gzclose(gz);
    fprintf(stderr, "parsed lines %ld\n", n_lines);

    char path[512];
    snprintf(path, sizeof path, "%s/libsize.f64", out_dir);
    FILE *out = fopen(path, "wb");
    fwrite(lib, sizeof(double), (size_t)n_kept, out);
    fclose(out);
    snprintf(path, sizeof path, "%s/genes_found.txt", out_dir);
    out = fopen(path, "w");
    for (int i = 0; i < n_genes; i++) {
        if (!seen[i]) continue;
        fprintf(out, "%s\n", names[i]);
        snprintf(path, sizeof path, "%s/%s.f32", out_dir, names[i]);
        FILE *gfout = fopen(path, "wb");
        fwrite(gmat[i], sizeof(float), (size_t)n_kept, gfout);
        fclose(gfout);
    }
    fclose(out);
    fprintf(stderr, "wrote %s\n", out_dir);
    return 0;
}
