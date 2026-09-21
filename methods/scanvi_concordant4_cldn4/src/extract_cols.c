/* Stream a genes x cells gzipped TSV and keep a barcode subset.
 *
 * Usage: extract_cols matrix.txt.gz barcodes.tsv out_prefix
 * barcodes.tsv: barcode, sub (0/1), cldn4 (0/1)
 * Writes:
 *   out_prefix.coo.bin   int32 gene, int32 cell, float32 value (nonzero subsample)
 *   out_prefix.genes.txt
 *   out_prefix.cells.txt subsample barcodes in output-column order
 *   out_prefix.cldn4.tsv barcode, UMI for cldn4-flagged cells
 *   out_prefix.tacstd2.tsv barcode, UMI for the same cells
 *   out_prefix.meta.txt  n_genes n_cells nnz
 */
#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#define MAX_BARCODE 128

typedef struct {
    char bc[MAX_BARCODE];
    int sub;
    int cldn4;
    int sub_idx;
} Keep;

static int cmp_keep(const void *a, const void *b) {
    const Keep *ka = a;
    const Keep *kb = b;
    return strcmp(ka->bc, kb->bc);
}

static char *read_line_gz(gzFile f, char **buf, size_t *cap, size_t *out_len) {
    size_t len = 0;
    int c;
    if (*buf == NULL) {
        *cap = 1 << 20;
        *buf = malloc(*cap);
        if (!*buf) return NULL;
    }
    while ((c = gzgetc(f)) != -1) {
        if (len + 2 >= *cap) {
            *cap *= 2;
            char *nb = realloc(*buf, *cap);
            if (!nb) return NULL;
            *buf = nb;
        }
        if (c == '\n') break;
        if (c != '\r') (*buf)[len++] = (char)c;
    }
    if (c == -1 && len == 0) return NULL;
    (*buf)[len] = 0;
    *out_len = len;
    return *buf;
}

/* Faster block reader that yields complete lines via a callback-style pull.
 * For the wide matrix we read in 1 MiB blocks. */
typedef struct {
    gzFile f;
    char *chunk;
    size_t chunk_cap;
    size_t begin;
    size_t end;
    int eof;
    char *line;
    size_t line_cap;
} Reader;

static void reader_init(Reader *r, gzFile f) {
    memset(r, 0, sizeof(*r));
    r->f = f;
    r->chunk_cap = 4 << 20;
    r->chunk = malloc(r->chunk_cap);
    r->line_cap = 8 << 20;
    r->line = malloc(r->line_cap);
}

static void reader_free(Reader *r) {
    free(r->chunk);
    free(r->line);
}

static char *reader_next(Reader *r, size_t *out_len) {
    size_t len = 0;
    for (;;) {
        size_t i = r->begin;
        int found = 0;
        while (i < r->end) {
            if (r->chunk[i] == '\n') {
                found = 1;
                break;
            }
            i++;
        }
        size_t take = i - r->begin;
        if (len + take + 2 >= r->line_cap) {
            while (len + take + 2 >= r->line_cap) r->line_cap *= 2;
            char *nb = realloc(r->line, r->line_cap);
            if (!nb) return NULL;
            r->line = nb;
        }
        if (take) {
            memcpy(r->line + len, r->chunk + r->begin, take);
            len += take;
        }
        if (found) {
            r->begin = i + 1;
            while (len && r->line[len - 1] == '\r') len--;
            r->line[len] = 0;
            *out_len = len;
            return r->line;
        }
        if (r->eof) {
            if (len == 0) return NULL;
            r->line[len] = 0;
            *out_len = len;
            r->begin = r->end;
            return r->line;
        }
        int n = gzread(r->f, r->chunk, (unsigned)r->chunk_cap);
        if (n < 0) return NULL;
        if (n == 0) {
            r->eof = 1;
            r->begin = 0;
            r->end = 0;
            continue;
        }
        r->begin = 0;
        r->end = (size_t)n;
    }
}

static const char *skip_field(const char *p, const char *end) {
    while (p < end && *p != '\t') p++;
    return p;
}

static const char *parse_float(const char *p, const char *end, float *out) {
    if (p >= end || *p == '\t') {
        *out = 0.f;
        return p;
    }
    char tmp[64];
    int n = 0;
    while (p < end && *p != '\t' && n < 63) tmp[n++] = *p++;
    tmp[n] = 0;
    *out = (float)atof(tmp);
    return p;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s matrix.txt.gz barcodes.tsv out_prefix\n", argv[0]);
        return 2;
    }
    FILE *bf = fopen(argv[2], "r");
    if (!bf) {
        perror(argv[2]);
        return 1;
    }
    size_t nkeep_cap = 1024, nkeep = 0;
    Keep *keeps = malloc(nkeep_cap * sizeof(Keep));
    char linebuf[512];
    while (fgets(linebuf, sizeof linebuf, bf)) {
        char bc[MAX_BARCODE];
        int sub = 0, cld = 0;
        if (sscanf(linebuf, "%127s %d %d", bc, &sub, &cld) < 3) continue;
        if (nkeep + 1 >= nkeep_cap) {
            nkeep_cap *= 2;
            keeps = realloc(keeps, nkeep_cap * sizeof(Keep));
        }
        memset(&keeps[nkeep], 0, sizeof(Keep));
        strncpy(keeps[nkeep].bc, bc, MAX_BARCODE - 1);
        keeps[nkeep].sub = sub;
        keeps[nkeep].cldn4 = cld;
        keeps[nkeep].sub_idx = -1;
        nkeep++;
    }
    fclose(bf);
    qsort(keeps, nkeep, sizeof(Keep), cmp_keep);

    int n_sub = 0;
    for (size_t i = 0; i < nkeep; i++) {
        if (keeps[i].sub) keeps[i].sub_idx = n_sub++;
    }
    fprintf(stderr, "keep barcodes=%zu subsample=%d\n", nkeep, n_sub);

    gzFile gz = gzopen(argv[1], "rb");
    if (!gz) {
        fprintf(stderr, "cannot open %s\n", argv[1]);
        return 1;
    }
    Reader rd;
    reader_init(&rd, gz);
    size_t llen = 0;
    char *header = reader_next(&rd, &llen);
    if (!header) {
        fprintf(stderr, "empty matrix\n");
        return 1;
    }
    /* Copy header because reader buffer is reused. */
    char *header_copy = malloc(llen + 1);
    memcpy(header_copy, header, llen + 1);

    size_t ncol = 0, cap = 1024;
    int *sub_of = malloc(cap * sizeof(int));
    unsigned char *want_cld = malloc(cap);
    const char *p = header_copy;
    const char *end = header_copy + llen;
    /* skip first field (Index / gene header) */
    p = skip_field(p, end);
    if (p < end && *p == '\t') p++;
    while (p < end) {
        const char *s = p;
        while (p < end && *p != '\t') p++;
        size_t bl = (size_t)(p - s);
        char bc[MAX_BARCODE];
        if (bl >= MAX_BARCODE) bl = MAX_BARCODE - 1;
        memcpy(bc, s, bl);
        bc[bl] = 0;
        if (ncol + 1 >= cap) {
            cap *= 2;
            sub_of = realloc(sub_of, cap * sizeof(int));
            want_cld = realloc(want_cld, cap);
        }
        Keep key;
        memset(&key, 0, sizeof key);
        strncpy(key.bc, bc, MAX_BARCODE - 1);
        Keep *hit = bsearch(&key, keeps, nkeep, sizeof(Keep), cmp_keep);
        if (hit && hit->sub) sub_of[ncol] = hit->sub_idx;
        else sub_of[ncol] = -1;
        want_cld[ncol] = (hit && hit->cldn4) ? 1 : 0;
        ncol++;
        if (p < end && *p == '\t') p++;
    }
    fprintf(stderr, "matrix columns=%zu\n", ncol);

    char path[1024];
    snprintf(path, sizeof path, "%s.genes.txt", argv[3]);
    FILE *gf = fopen(path, "w");
    snprintf(path, sizeof path, "%s.coo.bin", argv[3]);
    FILE *cf = fopen(path, "wb");
    snprintf(path, sizeof path, "%s.cldn4.tsv", argv[3]);
    FILE *xf = fopen(path, "w");
    snprintf(path, sizeof path, "%s.tacstd2.tsv", argv[3]);
    FILE *tf = fopen(path, "w");
    if (!gf || !cf || !xf || !tf) {
        perror("open out");
        return 1;
    }

    /* Map cldn4 columns to a compact list for the one gene. */
    size_t n_cld = 0;
    for (size_t i = 0; i < ncol; i++)
        if (want_cld[i]) n_cld++;
    size_t *cld_cols = malloc(n_cld * sizeof(size_t));
    char **cld_names = malloc(n_cld * sizeof(char *));
    n_cld = 0;
    /* Re-walk header for cldn4 names in column order. */
    p = header_copy;
    end = header_copy + strlen(header_copy);
    p = skip_field(p, end);
    if (p < end && *p == '\t') p++;
    for (size_t col = 0; col < ncol; col++) {
        const char *s = p;
        while (p < end && *p != '\t') p++;
        if (want_cld[col]) {
            size_t bl = (size_t)(p - s);
            cld_names[n_cld] = malloc(bl + 1);
            memcpy(cld_names[n_cld], s, bl);
            cld_names[n_cld][bl] = 0;
            cld_cols[n_cld] = col;
            n_cld++;
        }
        if (p < end && *p == '\t') p++;
    }

    int gene_i = 0;
    long nnz = 0;
    int cldn4_done = 0;
    int tacstd2_done = 0;
    float *cld_vals = calloc(n_cld, sizeof(float));
    while ((header = reader_next(&rd, &llen)) != NULL) {
        const char *lp = header;
        const char *le = header + llen;
        const char *gs = lp;
        while (lp < le && *lp != '\t') lp++;
        size_t gl = (size_t)(lp - gs);
        char gname[256];
        if (gl > 255) gl = 255;
        memcpy(gname, gs, gl);
        gname[gl] = 0;
        /* uppercase */
        for (char *q = gname; *q; q++) *q = (char)toupper((unsigned char)*q);
        fprintf(gf, "%s\n", gname);
        if (lp < le && *lp == '\t') lp++;
        int is_cldn4 = (strcmp(gname, "CLDN4") == 0);
        int is_tacstd2 = (strcmp(gname, "TACSTD2") == 0);
        int is_target = is_cldn4 || is_tacstd2;
        if (is_target) memset(cld_vals, 0, n_cld * sizeof(float));
        size_t col = 0;
        size_t next_cld = 0;
        while (lp < le && col < ncol) {
            int need = sub_of[col] >= 0 || (is_target && want_cld[col]);
            float v = 0.f;
            if (need) lp = parse_float(lp, le, &v);
            else lp = skip_field(lp, le);
            if (sub_of[col] >= 0 && v != 0.f) {
                int32_t rec_g = gene_i;
                int32_t rec_c = sub_of[col];
                float rec_v = v;
                fwrite(&rec_g, 4, 1, cf);
                fwrite(&rec_c, 4, 1, cf);
                fwrite(&rec_v, 4, 1, cf);
                nnz++;
            }
            if (is_target && next_cld < n_cld && cld_cols[next_cld] == col) {
                cld_vals[next_cld] = v;
                next_cld++;
            }
            col++;
            if (lp < le && *lp == '\t') lp++;
        }
        if (is_cldn4 && !cldn4_done) {
            cldn4_done = 1;
            for (size_t i = 0; i < n_cld; i++) {
                fprintf(xf, "%s\t%.8g\n", cld_names[i], cld_vals[i]);
            }
            fprintf(stderr, "CLDN4 row gene_index=%d wrote %zu cells\n", gene_i, n_cld);
        } else if (is_tacstd2 && !tacstd2_done) {
            tacstd2_done = 1;
            for (size_t i = 0; i < n_cld; i++) {
                fprintf(tf, "%s\t%.8g\n", cld_names[i], cld_vals[i]);
            }
            fprintf(stderr, "TACSTD2 row gene_index=%d wrote %zu cells\n", gene_i, n_cld);
        }
        gene_i++;
        if (gene_i % 2000 == 0) {
            fprintf(stderr, "genes %d nnz %ld\n", gene_i, nnz);
        }
    }
    fclose(gf);
    fclose(cf);
    fclose(xf);
    fclose(tf);

    snprintf(path, sizeof path, "%s.cells.txt", argv[3]);
    FILE *cells = fopen(path, "w");
    /* Emit subsample barcodes in sub_idx order. */
    char **by_idx = calloc((size_t)n_sub, sizeof(char *));
    for (size_t i = 0; i < nkeep; i++) {
        if (keeps[i].sub && keeps[i].sub_idx >= 0) by_idx[keeps[i].sub_idx] = keeps[i].bc;
    }
    for (int i = 0; i < n_sub; i++) fprintf(cells, "%s\n", by_idx[i] ? by_idx[i] : "");
    fclose(cells);

    snprintf(path, sizeof path, "%s.meta.txt", argv[3]);
    FILE *mf = fopen(path, "w");
    fprintf(mf, "n_genes\t%d\nn_cells\t%d\nnnz\t%ld\ncldn4_found\t%d\ntacstd2_found\t%d\n", gene_i, n_sub, nnz,
            cldn4_done, tacstd2_done);
    fclose(mf);
    fprintf(stderr, "done genes=%d cells=%d nnz=%ld cldn4=%d tacstd2=%d\n", gene_i, n_sub, nnz, cldn4_done,
            tacstd2_done);

    gzclose(gz);
    reader_free(&rd);
    free(header_copy);
    free(sub_of);
    free(want_cld);
    free(keeps);
    free(cld_cols);
    for (size_t i = 0; i < n_cld; i++) free(cld_names[i]);
    free(cld_names);
    free(cld_vals);
    free(by_idx);
    return (cldn4_done && tacstd2_done) ? 0 : 3;
}
