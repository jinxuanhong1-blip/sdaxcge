#!/usr/bin/env nextflow
// Draft Nextflow DSL2 — same five stages as the Snakefile.
// METHODS ONLY. Processed GEO only. No FASTQ/SRA. See playbook.md (zh+en).
//
//   nextflow run main.nf -c nextflow.config -preview
//   nextflow run main.nf -c nextflow.config -profile conda

nextflow.enable.dsl = 2

if (!params.refuse_fastq) {
    error 'params.refuse_fastq must stay true; this engine does not fetch SRA/FASTQ'
}

def script_dir = "${projectDir}/scripts"

process DOWNLOAD_GEO_PROCESSED {
    tag { cohort }
    conda "${projectDir}/envs/environment.yaml"

    input:
    tuple val(cohort), val(meta)

    output:
    tuple val(cohort), path("${cohort}/processed/matrix.bin"), path("${cohort}/processed/family.soft.gz"), path("${cohort}/download_qc.json")

    script:
    def py = "${script_dir}/download_geo.py"
    def pheno = meta.pheno_url ? "--pheno-url '${meta.pheno_url}'" : ""
    """
    set -euo pipefail
    mkdir -p ${cohort}/processed
    if [ ! -f ${py} ]; then
        echo "SKELETON: implement ${py}" >&2
        echo "  cohort=${cohort} gse=${meta.gse}" >&2
        echo "  processed_url=${meta.processed_url}" >&2
        echo "  max_bytes=${params.max_file_bytes} refuse_fastq=1" >&2
        echo "  Forbidden: prefetch fasterq-dump fastq-dump kingfisher ffq SRA RAW FASTQ" >&2
        exit 2
    fi
    python3 ${py} \\
        --cohort ${cohort} \\
        --gse ${meta.gse} \\
        --processed-url '${meta.processed_url}' \\
        --soft-url '${meta.soft_url}' \\
        ${pheno} \\
        --max-bytes ${params.max_file_bytes} \\
        --refuse-fastq \\
        --matrix-out ${cohort}/processed/matrix.bin \\
        --soft-out ${cohort}/processed/family.soft.gz \\
        --qc-out ${cohort}/download_qc.json
    """
}

process EXTRACT_GENES {
    tag { cohort }
    conda "${projectDir}/envs/environment.yaml"

    input:
    tuple val(cohort), path(matrix), path(soft), path(dl_qc), val(meta)

    output:
    tuple val(cohort), path("${cohort}/genes.tsv"), path("${cohort}/matrix_symbols.tsv.gz"), path("${cohort}/extract_qc.json"), val(meta)

    script:
    def py = "${script_dir}/extract_genes.py"
    """
    set -euo pipefail
    mkdir -p ${cohort}
    if [ ! -f ${py} ]; then
        echo "SKELETON: implement ${py}" >&2
        echo "  resolve TACSTD2=ENSG00000184292 CLDN4=ENSG00000189143" >&2
        echo "  unit=${meta.unit} panel=${meta.panel}" >&2
        echo "  if genes missing: genes_present=false; do not impute" >&2
        exit 2
    fi
    python3 ${py} \\
        --cohort ${cohort} \\
        --matrix ${matrix} \\
        --soft ${soft} \\
        --unit ${meta.unit} \\
        --panel ${meta.panel} \\
        --collapse ${params.collapse_duplicate_genes} \\
        --config ${projectDir}/config.yaml \\
        --genes-out ${cohort}/genes.tsv \\
        --symbol-matrix-out ${cohort}/matrix_symbols.tsv.gz \\
        --qc-out ${cohort}/extract_qc.json
    """
}

process DECONV {
    tag { cohort }
    conda "${projectDir}/envs/environment.yaml"

    input:
    tuple val(cohort), path(genes), path(symbol_matrix), path(extract_qc), val(meta)

    output:
    tuple val(cohort), path(genes), path("${cohort}/deconv.tsv"), path("${cohort}/deconv_qc.json")

    script:
    def rs = "${script_dir}/deconv.R"
    def methods = params.deconv_methods.join(',')
    """
    set -euo pipefail
    mkdir -p ${cohort}
    if [ "${meta.panel}" = "oncomine_immune" ] || [ "${meta.panel}" = "ncounter" ]; then
        echo "panel_too_narrow: skip immunedeconv for ${cohort}" >&2
    fi
    if [ ! -f ${rs} ]; then
        echo "SKELETON: implement ${rs}" >&2
        echo "  methods=${methods} unit=${meta.unit}" >&2
        echo "  invert log2tpm; do not treat xCell/MCP as percentages" >&2
        exit 2
    fi
    Rscript ${rs} \\
        --cohort ${cohort} \\
        --matrix ${symbol_matrix} \\
        --unit ${meta.unit} \\
        --panel ${meta.panel} \\
        --methods ${methods} \\
        --deconv-out ${cohort}/deconv.tsv \\
        --qc-out ${cohort}/deconv_qc.json
    """
}

process STATS {
    conda "${projectDir}/envs/environment.yaml"

    input:
    path gene_tables
    path deconv_tables

    output:
    path "stats/associations.tsv"

    script:
    def py = "${script_dir}/stats.py"
    """
    set -euo pipefail
    mkdir -p stats
    if [ ! -f ${py} ]; then
        echo "SKELETON: implement ${py}" >&2
        echo "  tests: MWU+rank-biserial, Spearman, Cox+median log-rank" >&2
        echo "  BH FDR within cohort family; n>=4 per arm; no fabricated p" >&2
        exit 2
    fi
    python3 ${py} \\
        --genes ${gene_tables} \\
        --deconv ${deconv_tables} \\
        --out stats/associations.tsv
    """
}

process FIGURES {
    conda "${projectDir}/envs/environment.yaml"

    input:
    path associations
    path gene_tables
    path deconv_tables

    output:
    path "figures/figure_index.tsv"

    script:
    def py = "${script_dir}/figures.py"
    """
    set -euo pipefail
    mkdir -p figures
    if [ ! -f ${py} ]; then
        echo "SKELETON: implement ${py}" >&2
        echo "  gene-by-response, deconv heatmap, TACSTD2 vs CD8, KM if survival" >&2
        exit 2
    fi
    python3 ${py} \\
        --stats ${associations} \\
        --genes ${gene_tables} \\
        --deconv ${deconv_tables} \\
        --outdir figures \\
        --index figures/figure_index.tsv
    """
}

process MANIFEST {
    input:
    val rows

    output:
    path "manifest.tsv"

    exec:
    def header = [
        "cohort", "gse", "species", "assay", "unit", "panel",
        "processed_url", "soft_url", "max_file_bytes", "refuse_fastq", "note",
    ].join("\t")
    def body = rows.collect { r ->
        [
            r.cohort, r.gse, r.species, r.assay, r.unit, r.panel,
            r.processed_url, r.soft_url, params.max_file_bytes,
            params.refuse_fastq, (r.note ?: "").replace("\t", " "),
        ].join("\t")
    }.join("\n")
    task.workDir.resolve("manifest.tsv").text = header + "\n" + body + "\n"
}

workflow {
    def cohort_ch = Channel.fromList(
        params.cohorts.collect { name, meta ->
            tuple(name, meta + [cohort: name])
        }
    )

    MANIFEST(
        Channel.value(
            params.cohorts.collect { name, meta ->
                meta + [cohort: name]
            }
        )
    )

    DOWNLOAD_GEO_PROCESSED(cohort_ch)

    extract_in = DOWNLOAD_GEO_PROCESSED.out.combine(cohort_ch, by: 0)
    EXTRACT_GENES(extract_in)

    DECONV(EXTRACT_GENES.out)

    gene_tables = DECONV.out.map { _c, genes, _d, _q -> genes }.collect()
    deconv_tables = DECONV.out.map { _c, _g, deconv, _q -> deconv }.collect()

    STATS(gene_tables, deconv_tables)
    FIGURES(STATS.out, gene_tables, deconv_tables)
}

workflow.onComplete {
    log.info "Finished: ${workflow.success}. FASTQ was never in the DAG."
}
