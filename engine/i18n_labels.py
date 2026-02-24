"""
engine/i18n_labels.py -- Pure label lookup for the engine.

Contains the TRANSLATIONS dict and stateless lookup functions.
No Streamlit dependency.

Data
----
TRANSLATIONS : dict[str, dict[str, str]]
    Nested dictionary: ``TRANSLATIONS[lang_code][key] -> str``.
    Supported languages: ``"en"`` (English), ``"es"`` (Spanish).

Functions
---------
get_label(key, lang)
    → Return the raw label string for a key in the given language.

get_label_formatted(key, lang, **kwargs)
    → Return the label with ``str.format(**kwargs)`` applied.
"""


TRANSLATIONS: dict = {

    # ──────────────────────────────────────────────────────────────────
    # ENGLISH
    # ──────────────────────────────────────────────────────────────────
    "en": {

        # ── App title & welcome ───────────────────────────────────────
        "app.title": "BeginSeq Studio — RNA-seq Differential Expression",
        "app.welcome": (
            "Welcome!\n"
            "This tool allows you to perform **RNA-seq differential expression analysis** "
            "from a raw count matrix using **DESeq2 (Python)**.\n\n"
            "**Steps:**\n"
            "1. Upload your **raw count matrix** (TSV or CSV) and **metadata** (CSV).\n"
            "2. Configure your **metadata columns** (sample ID, condition, optional filters).\n"
            "3. Click **Run differential expression** to execute the analysis.\n"
            "4. View results (Volcano, PCA, MA, Heatmap) and download the output."
        ),

        # ── Upload section ────────────────────────────────────────────
        "upload.header": "Upload your files",
        "upload.counts_label": "Upload raw count matrix (TSV, CSV, or ZIP)",
        "upload.counts_help": (
            "Tabular file with genes as rows and samples as columns. "
            "The first column should be the gene identifier. "
            "Values: raw counts (integers >= 0). "
            "ZIP files containing a CSV/TSV are also accepted."
        ),
        "upload.expander_title": "What type of count matrix does DESeq2 accept?",
        "upload.counts_info": (
            "DESeq2 performs its own internal normalization, so it **requires "
            "raw (unnormalized) integer counts**.\n\n"
            "| | Type | Examples |\n"
            "|---|---|---|\n"
            "| **Accept** | featureCounts output | Integer counts from read alignment |\n"
            "| **Accept** | HTSeq counts | `htseq-count` output |\n"
            "| **Accept** | Salmon counts (converted) | `tximport` aggregated counts |\n"
            "| **Accept** | Raw GEO matrices | Series matrix with raw counts |\n"
            "| **Reject** | FPKM | Fragments Per Kilobase per Million |\n"
            "| **Reject** | TPM | Transcripts Per Million |\n"
            "| **Reject** | CPM | Counts Per Million |\n"
            "| **Reject** | Log-normalized | `log2(CPM+1)`, `rlog`, `vst`, etc. |\n\n"
            "> **How to tell?** Raw counts are **non-negative integers** (0, 1, 5, 342, ...).\n"
            "> If your matrix has many decimal values (3.72, 0.0041, 8.156...), it is "
            "likely normalized and **should not** be used directly with DESeq2."
        ),
        "upload.metadata_label": "Upload metadata file (CSV or ZIP)",
        "upload.metadata_help": (
            "CSV file with sample information. "
            "Must have a column with sample names matching the count matrix columns, "
            "and at least one experimental condition column. "
            "You can select the columns afterwards."
        ),

        # ── File size limit ─────────────────────────────────────────────
        "upload.file_too_large": (
            "**File too large:** `{file}` is {size_mb} MB, "
            "but the limit for this tool is **{limit_mb} MB**."
        ),

        # ── File reading errors ───────────────────────────────────────
        "error.reading_files": "Error reading files: {error}",

        # ── Previews ──────────────────────────────────────────────────
        "preview.counts_title": "Preview: Count matrix",
        "preview.counts_caption": "{genes} genes x {samples} samples",
        "preview.genes_removed": (
            "{count:,} genes with zero counts across all "
            "samples were automatically removed."
        ),
        "warning.normalized_matrix": (
            "**Possible normalized matrix detected**\n\n"
            "{reason}\n\n"
            "If this matrix contains FPKM, TPM, CPM, or "
            "log-transformed values, the DESeq2 results will be "
            "**statistically invalid**. Please use raw counts.\n\n"
            '_See the "What type of count matrix does DESeq2 '
            'accept?" section above for details._'
        ),
        "preview.metadata_title": "Preview: Metadata",
        "preview.metadata_caption": "{rows} rows x {cols} columns",

        # ── Configure metadata ────────────────────────────────────────
        "config.header": "Configure metadata columns",
        "config.instructions": (
            "Select which columns to use from your metadata file:\n"
            "- **Sample column**: The column containing sample names (must match count matrix columns).\n"
            "- **Condition column**: The column defining experimental conditions (e.g., 'treatment', 'disease').\n"
            "- **Filter columns** (optional): Additional columns to filter samples before analysis."
        ),
        "config.sample_col_label": "Sample column (sample identifiers)",
        "config.sample_col_help": (
            "This column should contain the sample names that "
            "match the columns in your count matrix."
        ),
        "config.condition_col_label": "Condition column (experimental groups)",
        "config.condition_col_help": (
            "This column defines the experimental conditions "
            "for differential expression analysis."
        ),

        # ── Filters (multiple) ────────────────────────────────────────
        "filter.subheader": "Optional: Filter samples",
        "filter.checkbox_label": "Apply filters to select specific samples",
        "filter.checkbox_help": "Enable this to filter samples based on other metadata columns.",
        "filter.title": "Filter {n}",
        "filter.column_label": "Column",
        "filter.values_label": "Values to include from '{column}'",
        "filter.values_help": "Only samples with these values will be included.",
        "filter.no_values_warning": (
            "No values selected for filter {n} "
            "-- all samples will be excluded!"
        ),
        "filter.add_button": "Add another filter",
        "filter.summary": (
            "After filtering: **{count}** samples "
            "remaining ({n_filters} filter(s) applied)"
        ),
        "filter.all_removed_warning": "All samples were filtered out! Adjust your filters.",

        # ── Validation error ──────────────────────────────────────────
        "error.validation": "Validation error: {error}",

        # ── Sample mismatch ───────────────────────────────────────────
        "mismatch.warning": "**Sample mismatch detected**",
        "mismatch.metric_counts": "Samples in counts",
        "mismatch.metric_metadata": "Samples in metadata",
        "mismatch.metric_common": "Samples in common",
        "mismatch.no_common": (
            "There are NO common samples between the files. "
            "Cannot continue. Verify that sample names match."
        ),
        "mismatch.details_title": "Details: sample mismatch",
        "mismatch.only_in_counts": "**Samples in counts but NOT in metadata** ({count}):",
        "mismatch.only_in_metadata": "**Samples in metadata but NOT in counts** ({count}):",
        "mismatch.proceed_info": (
            "**{n_common}** samples are present in both "
            "files. The analysis can proceed using only these common "
            "samples. Samples without a match will be excluded. "
            "This is the standard behavior in most RNA-seq pipelines."
        ),
        "mismatch.proceed_checkbox": "I understand. Proceed with the common samples only.",

        # ── Condition levels error ────────────────────────────────────
        "error.single_condition": (
            "**Cannot run analysis:** The samples that will be "
            "analyzed only have **{n_levels}** condition "
            "level(s): **{levels}**.\n\n"
            "DESeq2 needs at least **2 different conditions** "
            "(e.g., 'tumor' and 'normal') to compute a contrast.\n\n"
            "**What to check:**\n"
            "- Does your count matrix contain samples from BOTH conditions?\n"
            "- Does your metadata cover the samples in the count matrix?\n\n"
            "Currently analyzing **{n_samples}** samples, "
            "all labeled as **'{first_level}'**."
        ),

        # ── Analysis parameters ───────────────────────────────────────
        "params.header": "Analysis parameters",
        "params.reference_label": "**Reference** condition (baseline)",
        "params.reference_help": (
            "The reference/baseline condition. "
            "Fold-changes will be calculated RELATIVE to this level. "
            "For example, if you select 'normal', a positive log2FC "
            "means the gene is MORE expressed in the test condition."
        ),
        "params.test_label": "**Test** condition (vs reference)",
        "params.test_help_single": "Only one test condition available.",
        "params.test_help_multi": (
            "The condition to compare against the reference. "
            "A positive log2FC means the gene is MORE expressed "
            "in THIS condition compared to the reference."
        ),
        "params.contrast_summary": (
            "**Contrast:** {test} vs {ref} (reference)\n\n"
            "Positive log2FC -> gene is **up-regulated** in **{test}**\n\n"
            "Negative log2FC -> gene is **down-regulated** in **{test}**"
        ),
        "params.run_button": "Run differential expression",
        "params.running": "Running DESeq2 analysis... This may take a moment.",
        "params.success": "Analysis completed! {n_genes} genes analyzed.",
        "error.unexpected": "Unexpected error during analysis: {error}",

        # ── Filter Results ────────────────────────────────────────────
        "results_filter.header": "Filter Results",
        "results_filter.stats": (
            "**Expression statistics:**\n"
            "- Minimum baseMean: **{min_bm:.2f}**\n"
            "- Median baseMean: **{med_bm:.2f}**\n"
            "- Maximum baseMean: **{max_bm:.2f}**\n\n"
            "Low-expressed genes often have unreliable statistics. "
            "It's common practice to filter them out."
        ),
        "results_filter.radio_label": "Select filtering strategy:",
        "results_filter.radio_default": "Default filters (recommended)",
        "results_filter.radio_custom": "Custom filters",
        "results_filter.radio_none": "No filtering (show all genes)",
        "results_filter.radio_help": "Choose how to filter the results before visualization.",
        "results_filter.default_info": (
            "**Default filters:**\n"
            "- baseMean >= {basemean}\n"
            "- padj < {padj}\n"
            "- |log2FoldChange| > {log2fc}"
        ),
        "results_filter.basemean_help": "Filter out genes with low expression (baseMean below this value).",
        "results_filter.padj_help": "Filter for statistically significant genes.",
        "results_filter.log2fc_help": "Filter for biologically relevant fold changes.",
        "results_filter.no_filter_warning": (
            "No filters applied. Showing all genes including those with "
            "low expression or non-significant p-values."
        ),
        "results_filter.summary": (
            "**Filtered results:** {n_filtered:,} genes "
            "({n_removed:,} removed from {n_total:,} total)"
        ),

        # ── Results Table ─────────────────────────────────────────────
        "table.subheader": "Results Table",
        "table.no_genes": "No genes passed the current filters. Try relaxing the thresholds.",
        "table.caption": "Showing top 50 of {n:,} genes (sorted by padj)",

        # ── Gene Ranking ──────────────────────────────────────────────
        "ranking.subheader": "Gene Ranking",
        "ranking.description": (
            "The ranking table contains genes sorted by their Wald statistic (`stat`). "
            "Positive values indicate up-regulation, negative values indicate down-regulation."
        ),
        "ranking.id_label": "Gene ID column name",
        "ranking.id_help": "Name for the gene identifier column in the ranking output.",
        "ranking.top_n_label": "Show top N genes",
        "ranking.top_n_help": "Number of top genes to display in the ranking.",
        "ranking.caption": "Showing top {n} of {total:,} genes (sorted by stat)",

        # ── Download ──────────────────────────────────────────────────
        "download.subheader": "Download Results",
        "download.filename_label": "File name (without extension)",
        "download.filename_help": "Enter the name for your output file.",
        "download.format_label": "Format",
        "download.format_help": "CSV uses comma separator, TSV uses tab separator.",
        "download.btn_filtered": "Filtered ({fmt})",
        "download.btn_all": "All results ({fmt})",
        "download.btn_ranking": "Ranking ({fmt})",

        # ── Visualizations ────────────────────────────────────────────
        "viz.subheader": "Visualizations",
        "viz.label_checkbox": "Label top 10 genes (by ranking stat) in Volcano Plots",
        "viz.label_help": "Show gene names for the top 10 genes based on ranking statistic.",
        "viz.tab_volcano_all": "Volcano Plot (All)",
        "viz.tab_volcano_filt": "Volcano Plot (Filtered)",
        "viz.tab_pca": "PCA",
        "viz.tab_ma": "MA Plot",
        "viz.tab_heatmap": "Heatmap",
        "viz.volcano_title": "Volcano Plot -- Differential Expression",
        "viz.no_genes_volcano": "No genes to display with current filters.",
        "viz.filtered_caption": (
            "Showing {n:,} genes after filtering "
            "(padj < {padj}, |log2FC| > {log2fc})"
        ),
        "viz.basemean_filter_label": "Min. baseMean (noise filter)",
        "viz.basemean_filter_help": (
            "Exclude genes with mean expression below this threshold. "
            "Removes low-expression noise from the volcano plot."
        ),
        "viz.legend_position_label": "Legend position",
        "viz.legend_position_help": (
            "Choose where to place the legend on volcano, PCA and MA plots."
        ),
        "viz.download_png": "Download PNG",
        "viz.download_svg": "Download SVG",

        # ── Batch correction (PCA) ──────────────────────────────────
        "pca.batch_checkbox": "Enable batch correction (ComBat) for PCA",
        "pca.batch_checkbox_help": (
            "Apply ComBat to remove batch effects from the transformed "
            "expression matrix before PCA.  Only affects visualization, "
            "NOT the DESeq2 differential expression results."
        ),
        "pca.batch_select_label": "Batch column",
        "pca.batch_select_help": (
            "Select the metadata column that identifies the batch "
            "(e.g. sequencing run, plate, center)."
        ),

        # ── Highlight ─────────────────────────────────────────────────
        "highlight.header": "Highlight Genes of Interest",
        "highlight.instructions": (
            "Enter one or more gene names/IDs to highlight them on the volcano plot. "
            "The highlighted genes will be shown as **orange diamonds** with labels, "
            "while all other genes are dimmed in the background."
        ),
        "highlight.textarea_label": "Gene names/IDs (one per line, or comma-separated)",
        "highlight.textarea_placeholder": "e.g.\nBRCA1\nTP53\nMYC",
        "highlight.textarea_help": "Enter gene identifiers exactly as they appear in your count matrix.",
        "highlight.not_found": (
            "{count} gene(s) not found in results: **{genes}**"
        ),
        "highlight.not_found_more": " and {extra} more...",
        "highlight.found": "{count} gene(s) found: **{genes}**",
        "highlight.details_title": "Details of highlighted genes",
        "highlight.download_button": "Download highlighted volcano plot (PNG)",
        "highlight.none_found": "None of the entered genes were found in the results.",

        # ── Landing page ──────────────────────────────────────────────
        "landing.title": "BeginSeq Studio",
        "landing.subtitle": "No coding required — select a tool to get started with your analysis.",
        "landing.bulk_title": "Bulk RNA-seq Classic",
        "landing.bulk_desc": "Differential expression analysis from raw count matrices using DESeq2.",
        "landing.scrna_title": "scRNA-seq Analysis",
        "landing.scrna_desc": "Single-cell RNA-seq analysis: QC, clustering, UMAP & marker genes with Scanpy.",
        "landing.dataset_title": "Dataset Creator",
        "landing.dataset_desc": "Create and prepare datasets for downstream analysis.",
        "landing.coming_soon": "Coming soon",
        "landing.open_tool": "Open tool",
        "landing.disclaimer": (
            "Designed for beginners and users without programming experience. "
            "For research and educational purposes only — not intended for clinical or diagnostic use. "
            "Provided \"as is\", without warranty of any kind. "
            "See README for full disclaimer."
        ),

        # ── Placeholder pages ─────────────────────────────────────────
        "placeholder.title_scrna": "scRNA-seq Analysis",
        "placeholder.title_dataset": "Dataset Creator",
        "placeholder.message": "This tool is under development. Check back soon!",
        "placeholder.scrna_preview": (
            "**Planned features:**\n"
            "- Quality control and filtering\n"
            "- Dimensionality reduction (PCA, UMAP, t-SNE)\n"
            "- Clustering and cell type annotation\n"
            "- Differential expression between clusters"
        ),
        "placeholder.dataset_preview": (
            "**Planned features:**\n"
            "- Merge multiple count matrices\n"
            "- Generate metadata templates\n"
            "- Format conversion (CSV, TSV, H5AD)\n"
            "- Quality control reports"
        ),

        # ── DESeq2 pipeline progress ─────────────────────────────────
        "progress.validating": "Validating data...",
        "progress.filtering_genes": "Filtering low-expression genes...",
        "progress.size_factors": "Computing size factors (normalization)...",
        "progress.genewise_disp": "Estimating gene-wise dispersions...",
        "progress.disp_trend": "Fitting dispersion trend curve...",
        "progress.map_disp": "Fitting MAP dispersions...",
        "progress.fit_lfc": "Estimating log-fold changes...",
        "progress.cooks": "Computing Cook's distances...",
        "progress.deseq2_done": "DESeq2 model fitted!",
        "progress.wald_test": "Running Wald test & p-value adjustment...",
        "progress.wald_done": "Statistical testing complete!",
        "progress.visualizations": "Generating plots (Volcano, PCA, MA, Heatmap)...",
        "progress.done": "Analysis complete! ✅",
        "progress.time_estimate": "Estimating runtime...",

        # ── Time estimation ────────────────────────────────────────────
        "time.estimate_header": "Estimated time",
        "time.estimate_detail": (
            "Dataset: **{n_samples:,}** samples × **{n_genes:,}** genes — "
            "Estimated time: **{estimate}**"
        ),
        "time.estimate_detail_filtered": (
            "Dataset: **{n_samples:,}** samples × **{n_genes_raw:,}** genes "
            "(~**{n_genes_est:,}** after filtering) — "
            "Estimated time: **{estimate}**"
        ),
        "time.elapsed": "Elapsed: {elapsed}",
        "time.remaining": "~{remaining} remaining",
        "time.step_with_eta": "{step_msg} — ~{remaining} remaining",
        "time.completed_in": "Analysis complete in **{elapsed}**! ✅",
        "time.faster_than_expected": "Completed faster than expected!",
        "time.step_breakdown": "Step timing breakdown",
        "time.step_name_size_factors": "Size factors",
        "time.step_name_genewise_disp": "Gene-wise dispersions",
        "time.step_name_disp_trend": "Dispersion trend",
        "time.step_name_map_disp": "MAP dispersions",
        "time.step_name_fit_lfc": "Log-fold changes",
        "time.step_name_cooks": "Cook's distances",
        "time.step_name_wald_test": "Wald test",
        "time.step_name_visualizations": "Visualizations",

        # ── Expression-based classification ───────────────────────────
        "classify.subheader": "Expression-based Sample Classification",
        "classify.description": (
            "Classify samples into sub-groups based on marker gene expression. "
            "This creates **compound conditions** (e.g., tumor_TNBC vs tumor_nonTNBC) "
            "for more specific differential expression analysis.\n\n"
            "**How it works:** Lightweight CPM normalization → log2(CPM+1) → "
            "classify samples by marker gene thresholds → merge into metadata."
        ),
        "classify.checkbox": "Enable expression-based classification",
        "classify.checkbox_help": (
            "Classify samples into sub-groups based on expression of marker genes. "
            "Runs BEFORE the differential expression analysis."
        ),
        "classify.gene_input_label": "Marker genes (paste list)",
        "classify.gene_input_placeholder": (
            'e.g.\nESR1, PGR, ERBB2\nor paste a comma-separated list:\n'
            '"ACTR3B","ANLN","BAG1","BCL2",...'
        ),
        "classify.gene_input_help": (
            "Enter gene names/IDs separated by commas, newlines, or semicolons. "
            "Quotes are automatically stripped. Genes not found in the count matrix "
            "will be reported but ignored."
        ),
        "classify.no_genes_warning": "Enter at least one marker gene to proceed.",
        "classify.genes_not_found_summary": (
            "**{n_found}** genes found, **{n_not_found}** not in matrix: {not_found_list}"
        ),
        "classify.genes_found_summary": "**{n}** marker genes found in the count matrix.",
        "classify.stats_expander": "Expression stats for {n} marker genes (log2 CPM+1)",
        "classify.zscore_caption": (
            "**z(thr)** = (threshold − mean) / sd — "
            "How many standard deviations the threshold is from the mean. "
            "Negative z → threshold below average; positive z → above average."
        ),
        "classify.global_params_title": "Classification threshold (applied to ALL markers)",
        "classify.global_threshold_label": "log2(CPM+1) threshold",
        "classify.global_direction_label": "Direction",
        "classify.threshold_help": (
            "log2(CPM+1) threshold applied to ALL marker genes. "
            "Expression values are in log2 scale "
            "(e.g., 1.0 ≈ 2 CPM, 3.3 ≈ 10 CPM, 6.6 ≈ 100 CPM)."
        ),
        "classify.direction_below": "Below (low expression → positive)",
        "classify.direction_above": "Above (high expression → positive)",
        "classify.direction_help": (
            "'Below': samples with ALL markers BELOW the threshold are classified as positive. "
            "'Above': samples with ALL markers ABOVE the threshold are classified as positive."
        ),
        "classify.positive_label": "Positive class label",
        "classify.positive_help": (
            "Label for samples that satisfy ALL rules "
            "(e.g., 'TNBC', 'HighExpr', 'Mutant')."
        ),
        "classify.negative_label": "Negative class label",
        "classify.negative_help": (
            "Label for samples that do NOT satisfy all rules "
            "(e.g., 'nonTNBC', 'LowExpr', 'WildType')."
        ),
        "classify.reference_keep_label": "Keep reference samples unchanged",
        "classify.reference_keep_help": (
            "If enabled, samples from the reference condition (e.g., 'control') "
            "keep their original label instead of being reclassified. "
            "Only non-reference samples get compound labels."
        ),
        "classify.preview_title": "Classification Preview",
        "classify.compound_preview": "Compound conditions after classification:",
        "classify.error": "Classification error: {error}",
        "classify.download_metadata": "Download classified metadata ({fmt})",

        # ── Pre-analysis warnings ────────────────────────────────────
        "warning.class_imbalance": (
            "**Severe class imbalance detected ({ratio}:1):** "
            "Group '**{large_group}**' has **{large_n}** samples while "
            "'**{small_group}**' only has **{small_n}**. "
            "This can reduce statistical power and bias the results."
        ),
        "warning.sparse_data": (
            "**Sparse data detected:** Every gene has at least one zero across samples. "
            "The '**poscounts**' normalization method will be used automatically "
            "instead of the default 'ratio' method (which would fail or freeze).\n\n"
            "**Implications:**\n"
            "• Changes the geometry of size factors.\n"
            "• May slightly shrink Fold Changes towards zero.\n"
            "• Required to prevent the analysis from freezing with your data."
        ),

        # ── Strict validation report ──────────────────────────────────
        "strict.report_title": "Quality Validation Report",
        "strict.errors": "Errors",
        "strict.warnings": "Warnings",
        "strict.info_items": "Info",
        "strict.details": "Details",
        "strict.download_report": "Download validation report (JSON)",

        # ── Audit log ──────────────────────────────────────────────────
        "audit.section_title": "Scientific Audit Log",
        "audit.description": (
            "A machine-readable record of all parameters, library versions, "
            "and execution details for reproducibility."
        ),
        "audit.download_json": "Download audit log (JSON)",
        "audit.download_txt": "Download audit log (TXT)",

        "strict.duplicate_samples": (
            "{count} duplicate sample name(s) found: **{samples}**. "
            "Each sample must have a unique identifier."
        ),
        "strict.duplicate_gene_ids": (
            "{count} duplicate gene ID(s) found: **{genes}**. "
            "Duplicate gene IDs will cause incorrect results."
        ),
        "strict.condition_nan": (
            "{count} sample(s) have missing or empty values in the "
            "condition column. Remove or fill these before running the analysis."
        ),
        "strict.condition_numeric": (
            "The condition column appears to contain numeric values "
            "(e.g. {sample_values}). Consider using descriptive labels "
            "like 'control' / 'treated' for clarity."
        ),
        "strict.alignment_high_loss": (
            "**{pct_lost:.0f}% of samples** were lost during alignment "
            "between counts and metadata. Check that sample names match."
        ),
        "strict.alignment_case_mismatch": (
            "Possible case mismatch detected between counts and metadata "
            "sample names: {pairs}. Names are compared case-sensitively."
        ),
        "strict.alignment_dropped": (
            "{n_dropped} sample(s) present in metadata but absent from counts: "
            "**{samples}**."
        ),
        "strict.normalized_data": (
            "Data appears to be normalized (not raw counts): {reason}. "
            "DESeq2 requires **raw integer counts** as input."
        ),
        "strict.class_imbalance": (
            "Class imbalance detected ({ratio}:1 between largest and "
            "smallest groups). This may reduce statistical power."
        ),
        "strict.sparse_data": (
            "Sparse data: a high proportion of zero counts was detected. "
            "The 'poscounts' normalization will be used automatically."
        ),
        "strict.post_filter_low": (
            "Only **{n_after}** genes remain after filtering (from {n_before}). "
            "Consider relaxing the filter thresholds."
        ),
        "strict.post_filter_critical": (
            "Only **{n_after}** genes remain after filtering (from {n_before}) — "
            "this is too few for reliable differential expression analysis."
        ),
        "strict.filter_param_aggressive": (
            "The min_samples_expressing parameter ({min_samples}) exceeds "
            "90% of total samples ({n_samples}). This may be too aggressive."
        ),
        "strict.passed": "All validation checks passed.",

        # ── Gene pre-filtering ──────────────────────────────────────
        "gene_filter.subheader": "Gene Pre-filtering",
        "gene_filter.description": (
            "Filter out lowly-expressed genes **before** running DESeq2. "
            "This is standard practice in RNA-seq analysis — it reduces memory "
            "usage, improves statistical power, and stabilises dispersion estimates."
        ),
        "gene_filter.checkbox": "Apply gene pre-filtering (recommended)",
        "gene_filter.checkbox_help": (
            "Remove genes with very low expression across all samples. "
            "Typically removes 40-60% of genes (noise) without losing "
            "biologically relevant signals."
        ),
        "gene_filter.min_total_label": "Min total count per gene",
        "gene_filter.min_total_help": (
            "A gene must have at least this many total counts "
            "(sum across ALL samples) to be kept. Default: 10."
        ),
        "gene_filter.min_samples_label": "Min samples expressing",
        "gene_filter.min_samples_help": (
            "A gene must be expressed in at least this many samples. "
            "Set to 0 for automatic: max(smallest group, 50% of samples). "
            "For large datasets (800+ samples), auto is recommended."
        ),
        "gene_filter.min_count_label": "Min count per sample",
        "gene_filter.min_count_help": (
            "Minimum raw count in a sample for the gene to be considered "
            "'expressed' in that sample. Default: 1 (any non-zero count). "
            "Higher values (5-10) give more aggressive filtering."
        ),
        "gene_filter.stats_summary": (
            "Gene filtering: {n_before:,} → {n_after:,} genes "
            "({n_removed:,} removed). Criteria: total count ≥ {min_total}, "
            "expressed (≥ {min_count}/sample) in ≥ {min_samples} samples."
        ),

        # ── Dataset Creator (GDC) ───────────────────────────────────
        "dc.title": "GDC Dataset Creator",
        "dc.subtitle": (
            "Download RNA-seq STAR-Counts datasets from the **Genomic Data Commons (GDC/TCGA)** "
            "and assemble DESeq2-ready count matrices with auto-generated metadata."
        ),
        "dc.step1_title": "Load TCGA Projects",
        "dc.step1_desc": (
            "Connect to the GDC API to retrieve the list of available TCGA projects."
        ),
        "dc.load_projects": "Load TCGA Projects",
        "dc.loading_projects": "Connecting to GDC API...",
        "dc.no_projects": "No TCGA projects found.",
        "dc.projects_loaded": "{count} TCGA projects loaded.",
        "dc.view_projects": "View all projects",
        "dc.error_api": "GDC API error: {error}",

        "dc.step2_title": "Select Project & Fetch Files",
        "dc.select_project": "Select a TCGA project",
        "dc.fetch_files": "Fetch RNA-seq Files",
        "dc.fetching_files": "Fetching RNA-seq files for {project}...",
        "dc.no_files": "No STAR-Counts RNA-seq files found for {project}.",
        "dc.files_found": "{count} RNA-seq files found for {project}.",
        "dc.metric_total": "Total files",
        "dc.metric_tumor": "Tumor",
        "dc.metric_normal": "Normal",
        "dc.metric_size": "Est. size",
        "dc.view_files": "View file details",

        "dc.step3_title": "Filter & Configure",
        "dc.filter_samples": "Filter by condition",
        "dc.select_conditions": "Select conditions to include",
        "dc.no_conditions_selected": "Select at least one condition to proceed.",
        "dc.filtered_summary": (
            "Selected: **{count}** files ({tumor} Tumor, {normal} Normal)"
        ),
        "dc.gene_id_label": "Gene identifier type",
        "dc.gene_id_help": (
            "**gene_name**: Human-readable symbol (e.g. TP53, BRCA1). "
            "**gene_id**: Ensembl ID (e.g. ENSG00000141510)."
        ),
        "dc.count_col_label": "Count column",
        "dc.count_col_help": (
            "**unstranded** (recommended): Total counts regardless of strand. "
            "Use stranded options only if your library prep was strand-specific."
        ),

        "dc.step4_title": "Download & Build Dataset",
        "dc.download_info": (
            "Ready to download **{n_files}** files (~{size} MB). "
            "Files will be downloaded, parsed, and assembled into a count matrix."
        ),
        "dc.start_download": "Download & Build Dataset",
        "dc.downloading": "Downloading files from GDC...",
        "dc.download_progress": "Downloading file {current}/{total}...",
        "dc.downloaded_count": "{count} files downloaded successfully.",
        "dc.no_files_downloaded": "No files could be downloaded. Check your connection.",
        "dc.error_download": "Download error: {error}",
        "dc.building_matrix": "Building count matrix...",
        "dc.parse_progress": "Parsing file {current}/{total}...",
        "dc.error_building": "Error building count matrix: {error}",
        "dc.done": "Done!",
        "dc.build_success": (
            "Count matrix assembled: **{genes:,}** genes x **{samples}** samples."
        ),

        "dc.step5_title": "Preview & Download Results",
        "dc.preview_counts": "Count Matrix",
        "dc.preview_counts_caption": "{genes:,} genes x {samples} samples (showing first 20 rows)",
        "dc.preview_metadata": "Metadata",
        "dc.condition_summary": "Condition Summary",
        "dc.download_results": "Download Files",
        "dc.download_counts": "Download Count Matrix (CSV)",
        "dc.download_metadata": "Download Metadata (CSV)",
        "dc.tip_bulk": (
            "**Tip:** You can use the downloaded files directly in the "
            "**Bulk RNA-seq Classic** tool for differential expression analysis!"
        ),

        # ── Sample limit controls ──────────────────────────────────
        "dc.sample_limit_title": "Limit number of samples",
        "dc.sample_limit_desc": (
            "For large projects (e.g. TCGA-BRCA with 1,000+ files), "
            "you can limit the number of samples to download. "
            "Samples are randomly selected while maintaining condition proportions."
        ),
        "dc.use_sample_limit": "Limit number of samples to download",
        "dc.limit_mode_label": "Limit mode",
        "dc.limit_mode_total": "Total samples",
        "dc.limit_mode_per_cond": "Per condition",
        "dc.limit_total_slider": "Maximum total samples",
        "dc.limit_per_cond_slider": "Maximum samples per condition",
        "dc.limited_summary": (
            "After limiting: **{count}** samples ({tumor} Tumor, {normal} Normal)"
        ),
        "dc.large_download_warning": (
            "You are about to download a large number of files. "
            "Consider using the sample limit option above to reduce download time."
        ),
        "dc.download_dir_label": "Download directory",
        "dc.download_dir_help": (
            "Directory where GDC files will be downloaded and extracted. "
            "Leave empty to use a temporary directory (files are deleted after building the matrix)."
        ),
        "dc.download_dir_persistent_info": (
            "Files will be saved to: `{path}`. "
            "They will persist after the matrix is built."
        ),

        # ── scRNA-seq Analysis ──────────────────────────────────────────
        "scrna.title": "scRNA-seq Analysis",
        "scrna.subtitle": (
            "Single-cell RNA-seq analysis pipeline powered by **Scanpy**. "
            "Upload your data and run the complete workflow: QC, normalization, "
            "dimensionality reduction, clustering, and marker gene identification."
        ),

        # ── AnnData Glossary ────────────────────────────────────────────
        "scrna.glossary_title": "AnnData Object Reference",
        "scrna.glossary_desc": (
            "All scRNA-seq data in Scanpy is stored in an **AnnData** object. "
            "This table describes its main components:"
        ),

        # ── 10x File Integrator ──────────────────────────────────────────
        "scrna.integrator_title": "10x File Integrator — Build H5AD from raw files",
        "scrna.integrator_desc": (
            "If you have raw 10x Genomics output files (**matrix.mtx**, "
            "**features.tsv** or **genes.tsv**, and **barcodes.tsv**), "
            "upload them here to generate an **.h5ad** file. "
            "Compressed files (`.gz`) are also accepted. "
            "Files with only 2 columns (gene ID + gene name) are handled automatically."
        ),
        "scrna.integrator_matrix_label": "matrix.mtx(.gz)",
        "scrna.integrator_matrix_help": "The sparse expression matrix in Matrix Market format.",
        "scrna.integrator_features_label": "features.tsv(.gz) or genes.tsv(.gz)",
        "scrna.integrator_features_help": (
            "Gene annotations file. Can have 2 columns (gene_id, gene_name) "
            "or 3 columns (gene_id, gene_name, feature_type). "
            "2-column files are fixed automatically."
        ),
        "scrna.integrator_barcodes_label": "barcodes.tsv(.gz)",
        "scrna.integrator_barcodes_help": "Cell barcode file, one barcode per line.",
        "scrna.integrator_run": "Build H5AD",
        "scrna.integrator_running": "Integrating 10x files into H5AD...",
        "scrna.integrator_success": "H5AD file generated successfully! Download it below, then upload it to run the analysis.",
        "scrna.integrator_filename_label": "File name (without extension)",
        "scrna.integrator_filename_help": "Choose a name for the output .h5ad file.",
        "scrna.integrator_download": "Download .h5ad",
        "scrna.integrator_error": "Error integrating files: {error}",

        # ── Integrator: optional metadata ──────────────────────────────
        "scrna.integrator_meta_title": "Optional: Add metadata",
        "scrna.integrator_meta_desc": (
            "Upload a CSV or TSV file with additional annotations. "
            "The system will automatically detect whether it matches "
            "cell barcodes (\u2192 adata.obs) or gene IDs (\u2192 adata.var) "
            "based on index overlap."
        ),
        "scrna.integrator_meta_label": "Metadata file (.csv, .tsv)",
        "scrna.integrator_meta_help": (
            "CSV or TSV with the first column as identifiers "
            "(cell barcodes or gene IDs). Remaining columns will be "
            "merged into the AnnData object automatically."
        ),
        "scrna.integrator_meta_auto_note": (
            "The first column will be matched against barcodes (obs) "
            "and gene names (var). The best match is selected automatically. "
            "You can override below if needed."
        ),
        "scrna.integrator_meta_target_label": "Merge metadata into:",
        "scrna.integrator_meta_target_auto": "Auto-detect (recommended)",
        "scrna.integrator_meta_target_obs": "adata.obs (cell annotations)",
        "scrna.integrator_meta_target_var": "adata.var (gene annotations)",
        "scrna.integrator_meta_merged_obs": (
            "Metadata merged into **adata.obs** (cell annotations). "
            "Matched **{overlap}** of {total} IDs ({pct}% overlap with barcodes)."
        ),
        "scrna.integrator_meta_merged_var": (
            "Metadata merged into **adata.var** (gene annotations). "
            "Matched **{overlap}** of {total} IDs ({pct}% overlap with genes)."
        ),
        "scrna.integrator_meta_no_match": (
            "Could not determine metadata target \u2014 less than 10% overlap "
            "with both barcodes and gene names. "
            "Try selecting the target manually above."
        ),

        # ── SoupX Ambient RNA Cleaner ──────────────────────────────────
        "scrna.soupx_title": "Ambient RNA Cleaner — SoupX",
        "scrna.soupx_disclaimer": (
            "**Note:** CellBender (deep-learning based) generally provides superior "
            "ambient RNA removal, but it requires a GPU and substantial computational "
            "resources that are beyond the scope of this program. "
            "SoupX is a lightweight alternative that runs entirely on CPU via R."
        ),
        "scrna.soupx_desc": (
            "Remove ambient RNA contamination from your scRNA-seq data using **SoupX** (R package via rpy2). "
            "SoupX estimates the contamination profile from the **raw (unfiltered)** droplets "
            "and subtracts it from the **filtered (cell-containing)** droplets.\n\n"
            "**You need two H5AD files:**\n"
            "- **Raw / unfiltered**: All droplets including empty ones (e.g. `raw_feature_bc_matrix`)\n"
            "- **Filtered**: Only cell-containing droplets (e.g. `filtered_feature_bc_matrix`)\n\n"
            "The cleaned output can then be uploaded below for the full analysis pipeline."
        ),
        "scrna.soupx_missing_rpy2": (
            "**rpy2** is not installed. Install it with: `pip install rpy2`"
        ),
        "scrna.soupx_missing_r": (
            "**R** is not found on this system. SoupX requires R to be installed. "
            "Download R from: https://cran.r-project.org/"
        ),
        "scrna.soupx_missing_pkg": (
            "**SoupX** R package is not installed. Open R and run: "
            "`install.packages('SoupX')`"
        ),
        "scrna.soupx_ready": "R + rpy2 + SoupX detected — ready to clean ambient RNA.",
        "scrna.soupx_raw_label": "Raw (unfiltered) H5AD",
        "scrna.soupx_raw_help": (
            "H5AD containing ALL droplets (raw_feature_bc_matrix). "
            "Includes empty droplets used to estimate the ambient RNA profile."
        ),
        "scrna.soupx_filt_label": "Filtered H5AD",
        "scrna.soupx_filt_help": (
            "H5AD containing only cell-containing droplets (filtered_feature_bc_matrix). "
            "This is the data that will be cleaned."
        ),
        "scrna.soupx_auto_label": "Automatic contamination estimation",
        "scrna.soupx_auto_help": (
            "Let SoupX automatically estimate the contamination fraction. "
            "Uncheck to set it manually (useful if auto-estimation fails)."
        ),
        "scrna.soupx_contam_label": "Contamination fraction",
        "scrna.soupx_contam_help": (
            "Manual contamination fraction override. "
            "Typical values: 0.01-0.20 (1%-20%). "
            "Higher values = more aggressive correction."
        ),
        "scrna.soupx_run": "Clean Ambient RNA (SoupX)",
        "scrna.soupx_running": "Running SoupX ambient RNA removal... This may take a few minutes.",
        "scrna.soupx_success": (
            "Ambient RNA removed successfully! Download the cleaned H5AD below, "
            "then upload it to run the analysis."
        ),
        "scrna.soupx_download": "Download soupx_cleaned.h5ad",
        "scrna.soupx_error": "SoupX error: {error}",

        # ── Upload ──────────────────────────────────────────────────────
        "scrna.upload_header": "Upload your data",
        "scrna.upload_label": "Upload scRNA-seq data (.h5ad)",
        "scrna.upload_help": (
            "Upload an **.h5ad** file (Scanpy/AnnData format).\n\n"
            "If you have raw 10x files (matrix.mtx, features.tsv, barcodes.tsv), "
            "use the **10x File Integrator** above to generate an .h5ad first."
        ),
        "scrna.data_loaded": (
            "Data loaded: **{n_cells:,}** cells x **{n_genes:,}** genes"
        ),
        "scrna.preview_title": "Data Preview",
        "scrna.preview_obs_title": "Cell Metadata (adata.obs)",
        "scrna.preview_obs_caption": "{rows:,} cells × {cols} columns",
        "scrna.preview_obs_empty": "No cell metadata found in this file.",
        "scrna.preview_var_title": "Gene Metadata (adata.var)",
        "scrna.preview_var_caption": "{rows:,} genes × {cols} columns",
        "scrna.preview_var_empty": "No gene metadata found in this file.",

        # ── QC Parameters ───────────────────────────────────────────────
        "scrna.qc_header": "Quality Control Parameters",
        "scrna.qc_description": (
            "Filter out low-quality cells and rarely detected genes. "
            "Adjust thresholds based on your QC violin plots."
        ),
        "scrna.min_genes_label": "Min genes per cell",
        "scrna.min_genes_help": "Remove cells with fewer genes than this threshold.",
        "scrna.max_genes_label": "Max genes per cell",
        "scrna.max_genes_help": (
            "Remove cells with more genes (potential doublets or multiplets)."
        ),
        "scrna.min_counts_label": "Min total counts per cell",
        "scrna.min_counts_help": "Remove cells with too few total counts.",
        "scrna.max_counts_label": "Max total counts per cell",
        "scrna.max_counts_help": "Remove cells with abnormally high counts (doublets).",
        "scrna.max_pct_mt_label": "Max % mitochondrial",
        "scrna.max_pct_mt_help": (
            "Cells with high mitochondrial gene percentage are often dying or stressed. "
            "Typical threshold: 10-20%."
        ),
        "scrna.min_cells_label": "Min cells per gene",
        "scrna.min_cells_help": "Remove genes detected in fewer cells than this threshold.",

        # ── Pipeline Parameters ──────────────────────────────────────────
        "scrna.params_header": "Analysis Parameters",
        "scrna.n_hvg_label": "Number of highly variable genes",
        "scrna.n_hvg_help": (
            "Top N variable genes used for PCA and downstream analysis. "
            "Standard: 2000. Increase for complex datasets."
        ),
        "scrna.n_pcs_label": "Number of PCs for neighbors",
        "scrna.n_pcs_help": (
            "Number of principal components used for the neighborhood graph. "
            "Check the elbow plot to choose an appropriate value."
        ),
        "scrna.n_neighbors_label": "Number of neighbors",
        "scrna.n_neighbors_help": (
            "Number of nearest neighbors for the graph. "
            "Higher = smoother clusters. Lower = finer resolution."
        ),
        "scrna.umap_min_dist_label": "UMAP min distance",
        "scrna.umap_min_dist_help": (
            "Controls how tightly UMAP packs points. "
            "Lower values = tighter clusters. Range: 0.0 - 1.0."
        ),
        "scrna.leiden_res_label": "Leiden resolution",
        "scrna.leiden_res_help": (
            "Controls the granularity of clustering. "
            "Higher = more clusters. Typical range: 0.1 - 2.0."
        ),
        "scrna.de_method_label": "DE method for markers",
        "scrna.de_method_help": "Statistical test for identifying marker genes.",
        "scrna.n_marker_genes_label": "Marker genes per cluster",
        "scrna.n_marker_genes_help": "Number of top marker genes to compute per cluster.",
        "scrna.doublet_checkbox": "Enable doublet detection (Scrublet)",
        "scrna.doublet_help": (
            "Detect and remove predicted doublets using Scrublet. "
            "Recommended for 10x Genomics data."
        ),

        # ── Batch Effect Correction ──────────────────────────────────────
        "scrna.batch_header": "Batch Effect Correction",
        "scrna.batch_checkbox": "Enable batch correction (Harmony)",
        "scrna.batch_help": (
            "Correct batch effects using Harmony, which adjusts the PCA embedding "
            "so cells from different batches integrate properly. "
            "Enable only if your data contains multiple batches."
        ),
        "scrna.batch_col_label": "Batch variable column",
        "scrna.batch_col_help": (
            "Select the column in adata.obs that identifies the batch "
            "(e.g. 'sample', 'donor', 'experiment')."
        ),
        "scrna.batch_preview": "Found **{n_batches}** batches in column `{col}`",
        "scrna.batch_no_columns": "No suitable batch columns found in the data.",

        # ── Annotation Columns ─────────────────────────────────────────────
        "scrna.annotation_header": "Annotation Columns",
        "scrna.celltype_col_label": "Cell type column (UMAP/PCA)",
        "scrna.celltype_col_help": (
            "Select a column from adata.obs containing cell type labels. "
            "This will be used as the default color for UMAP and PCA plots. "
            "Leave as '(none)' to use Leiden clusters."
        ),
        "scrna.marker_groupby_label": "Group markers by",
        "scrna.marker_groupby_help": (
            "Column used to group cells for marker gene analysis, "
            "DotPlots and Heatmaps. Default: Leiden clusters."
        ),

        # ── Visualization Settings ────────────────────────────────────────
        "scrna.viz_header": "Visualization Settings",
        "scrna.legend_position_label": "Legend position",
        "scrna.legend_position_help": "Choose where to place the legend on PCA and UMAP plots.",

        # ── Run ──────────────────────────────────────────────────────────
        "scrna.run_button": "Run scRNA-seq Analysis",
        "scrna.running": "Running scRNA-seq pipeline...",

        # ── Progress steps ──────────────────────────────────────────────
        "scrna.step.qc_annotation": "Annotating QC metrics (MT, Ribo, HB genes)...",
        "scrna.step.filtering": "Filtering cells and genes...",
        "scrna.step.doublet_detection": "Detecting doublets (Scrublet)...",
        "scrna.step.doublet_removal": "Removing predicted doublets...",
        "scrna.step.normalization": "Normalizing and log-transforming...",
        "scrna.step.hvg_selection": "Selecting highly variable genes...",
        "scrna.step.scaling": "Scaling data...",
        "scrna.step.pca": "Running PCA...",
        "scrna.step.batch_correction": "Correcting batch effects (Harmony)...",
        "scrna.step.neighbors": "Computing neighborhood graph...",
        "scrna.step.umap": "Computing UMAP embedding...",
        "scrna.step.clustering": "Running Leiden clustering...",
        "scrna.step.marker_genes": "Identifying marker genes...",
        "scrna.step.done": "Analysis complete!",

        # ── Results ──────────────────────────────────────────────────────
        "scrna.results_header": "Results",
        "scrna.stats_title": "Analysis Summary",
        "scrna.stats_cells": "Cells",
        "scrna.stats_genes": "Genes",
        "scrna.stats_clusters": "Clusters",
        "scrna.stats_hvg": "HVGs",

        "scrna.filter_stats": (
            "Filtering: **{cells_before:,}** -> **{cells_after:,}** cells "
            "({cells_removed:,} removed), "
            "**{genes_before:,}** -> **{genes_after:,}** genes "
            "({genes_removed:,} removed)"
        ),
        "scrna.doublet_stats": (
            "Doublets detected: **{n_doublets}** "
            "({doublet_rate:.1%} doublet rate)"
        ),
        "scrna.hvg_stats": (
            "Highly variable genes: **{n_hvg:,}** of **{n_total:,}** total"
        ),
        "scrna.harmony_stats": (
            "Batch correction (Harmony): **{n_batches}** batches "
            "from column `{batch_key}`, using **{n_pcs}** PCs"
        ),
        "scrna.leiden_stats": (
            "Leiden clustering: **{n_clusters}** clusters "
            "(resolution = {resolution})"
        ),

        # ── Visualization tabs ──────────────────────────────────────────
        "scrna.tab_qc": "QC Metrics",
        "scrna.tab_hvg": "Variable Genes",
        "scrna.tab_pca": "PCA",
        "scrna.tab_umap": "UMAP",
        "scrna.tab_markers": "Marker Genes",

        "scrna.umap_color_label": "Color UMAP by",
        "scrna.umap_color_help": "Select a metadata column or gene to color the UMAP plot.",
        "scrna.gene_search_label": "Search gene",
        "scrna.gene_search_help": "Type a gene name to visualize its expression on UMAP.",
        "scrna.gene_not_found": "Gene '{gene}' not found in the dataset.",
        "scrna.marker_cluster_label": "Select cluster",
        "scrna.marker_n_genes_label": "Top N genes to show",

        # ── Downloads ───────────────────────────────────────────────────
        "scrna.download_header": "Download Results",
        "scrna.download_h5ad": "Download H5AD (full analysis)",
        "scrna.download_h5ad_help": (
            "Download the complete AnnData object (.h5ad) with all "
            "analysis results, embeddings, and metadata."
        ),
        "scrna.download_markers_csv": "Download marker genes (CSV)",
        "scrna.download_obs_csv": "Download cell metadata (CSV)",
        "scrna.download_png": "Download PNG",
        "scrna.download_svg": "Download SVG",

        # ── Errors ──────────────────────────────────────────────────────
        "scrna.error_loading": "Error loading data: {error}",
        "scrna.error_pipeline": "Error during analysis: {error}",
        "scrna.error_lapack": (
            "⚠️ A numerical instability was detected during the analysis "
            "(near-singular matrix). This can happen with certain datasets."
        ),
        "scrna.error_lapack_hint": (
            "💡 **Suggestions:** Try reducing the number of PCs, disabling "
            "batch correction, changing the HVG count, or adjusting QC "
            "filters to remove low-quality cells/genes."
        ),
        "scrna.error_no_cells": (
            "No cells remain after filtering. Try relaxing the QC thresholds."
        ),

        # ── Remote access gate ──────────────────────────────────────────
        "scrna.remote_blocked": (
            "This feature is only available when running locally. "
            "If you would like to use this feature, visit: "
            "https://github.com/JACKNINES/beginseq-studio"
        ),

        # ── Page titles (browser tab) ─────────────────────────────────────
        "page_title.bulk": "Bulk RNA-seq",
        "page_title.scrna": "scRNA-seq Analysis",
        "page_title.dataset": "Dataset Creator",

        # ── Hardcoded UI strings — Bulk RNA-seq ───────────────────────────
        "bulk.step_caption": "Step {step}/{total}: {msg}",
        "bulk.params_updated_classification": "(updated with classification)",
        "bulk.more_genes_suffix": " (+{extra} more)",
        "bulk.timing_total": "Total: {elapsed}",
        "bulk.basemean_filter_label": "baseMean ≥",
        "bulk.padj_filter_label": "padj <",
        "bulk.log2fc_filter_label": "|log2FC| >",

        # ── Hardcoded UI strings — scRNA-seq page ─────────────────────────
        "scrna.cells_x_genes": "{n_cells:,} cells × {n_genes:,} genes",
        "scrna.step_caption": "Step {step}/{total}: {msg}",
        "scrna.option_none": "(none)",
        "scrna.no_markers_info": "No marker genes computed.",
        "scrna.gene_search_placeholder": "e.g. CD3D, MS4A1, NKG7",
        "scrna.dl_suffix_scatter": "(Scatter)",
        "scrna.dl_suffix_elbow": "(Elbow)",
        "scrna.dl_suffix_pca": "(PCA)",
        "scrna.dl_suffix_markers": "(Markers)",
        "scrna.dl_suffix_dotplot": "(Dotplot)",
        "scrna.dl_suffix_heatmap": "(Heatmap)",

        # ── AnnData glossary table (scRNA-seq page) ───────────────────────
        "scrna.glossary_col_element": "Element",
        "scrna.glossary_col_stands_for": "What it stands for",
        "scrna.glossary_col_contains": "What it contains",
        "scrna.glossary_x_name": "Expression matrix",
        "scrna.glossary_obs_name": "Observations (cells)",
        "scrna.glossary_var_name": "Variables (genes)",
        "scrna.glossary_obsm_name": "Observation multi-dimensional",
        "scrna.glossary_uns_name": "Unstructured annotations",
        "scrna.glossary_layers_name": "Alternative matrix layers",
        "scrna.glossary_x_desc": (
            "The gene expression matrix (cells x genes). "
            "Normalized/log-transformed values after preprocessing."
        ),
        "scrna.glossary_obs_desc": (
            "Metadata for each cell: cluster labels, QC metrics "
            "(n_genes, total_counts, pct_mt), sample IDs, cell types."
        ),
        "scrna.glossary_var_desc": (
            "Metadata for each gene: gene names, highly_variable flag, "
            "mean expression, dispersions."
        ),
        "scrna.glossary_obsm_desc": (
            "Cell-level embeddings and coordinates: PCA (X_pca), "
            "UMAP (X_umap), t-SNE (X_tsne)."
        ),
        "scrna.glossary_uns_desc": (
            "General analysis results: clustering parameters, marker "
            "gene rankings, color palettes, method settings."
        ),
        "scrna.glossary_layers_desc": (
            "Alternative representations of X: raw counts ('counts'), "
            "scaled data. All same shape as X."
        ),

        # ── Plot labels — scRNA visualization ──────────────────────────────
        "plot.qc_title": "Quality Control Metrics",
        "plot.qc_scatter_title": "QC Scatter",
        "plot.hvg_title": "Highly Variable Genes",
        "plot.hvg_label": "Highly Variable",
        "plot.hvg_other_label": "Other",
        "plot.hvg_xlabel": "Mean expression",
        "plot.hvg_ylabel_disp": "Normalized dispersion",
        "plot.hvg_ylabel_count": "Number of genes",
        "plot.pca_elbow_title": "PCA Variance Ratio (Elbow Plot)",
        "plot.pca_elbow_xlabel": "Principal Component",
        "plot.pca_elbow_ylabel": "Variance Ratio",
        "plot.pca_embed_title": "PCA Embedding",
        "plot.umap_xlabel": "UMAP1",
        "plot.umap_ylabel": "UMAP2",
        "plot.umap_title_prefix": "UMAP — {color}",
        "plot.marker_ranking_title": "Top Marker Genes per {groupby}",
        "plot.marker_score_label": "Score",
        "plot.marker_no_data": "No data",
        "plot.marker_no_genes": "No marker genes found",
        "plot.dotplot_xlabel": "Genes",
        "plot.dotplot_ylabel": "Cluster",
        "plot.dotplot_title": "Dot Plot — Marker Genes",
        "plot.dotplot_cbar_label": "Mean expression",
        "plot.dotplot_size_legend": "% expressing",
        "plot.heatmap_no_genes": "No marker genes to display",
        "plot.heatmap_cluster_title": "Cluster",
        "plot.heatmap_xlabel": "Cells (ordered by cluster)",
        "plot.heatmap_cbar_label": "Expression",
        "plot.heatmap_title": "Marker Gene Heatmap",

        # ── Plot labels — Bulk visualization ───────────────────────────────
        "plot.not_significant": "Not significant",
        "plot.up_regulated": "Up-regulated",
        "plot.down_regulated": "Down-regulated",
        "plot.volcano_summary": (
            "Total genes: {n_total:,}\n"
            "Significant: {n_sig:,} ({n_up:,} up · {n_down:,} down)"
        ),
        "plot.volcano_shrunk_suffix": "  [shrunk]",
        "plot.highlight_other_genes": "Other genes ({count:,})",
        "plot.highlight_highlighted": "Highlighted ({count:,})",
        "plot.highlight_title": "Volcano Plot — Highlighted Genes",
        "plot.highlight_info": "Highlighted: {n_found}",
        "plot.highlight_not_found": "Not found: {n_not_found}",
        "plot.pca_axis_variance": "PC{n} ({pct:.1f}% variance)",
        "plot.ma_not_significant": "Not significant",
        "plot.ma_up_regulated": "Up-regulated",
        "plot.ma_down_regulated": "Down-regulated",
        "plot.ma_summary": (
            "Total: {n_total:,} genes\n"
            "Significant: {n_sig:,} ({n_up:,} up · {n_down:,} down)"
        ),
        "plot.heatmap_zscore_label": "Z-score",
        "plot.heatmap_condition_title": "Condition",
    },

    # ──────────────────────────────────────────────────────────────────
    # SPANISH
    # ──────────────────────────────────────────────────────────────────
    "es": {

        # ── Titulo y bienvenida ───────────────────────────────────────
        "app.title": "BeginSeq Studio — Expresion Diferencial RNA-seq",
        "app.welcome": (
            "Bienvenido!\n"
            "Esta herramienta permite realizar un **analisis de expresion diferencial de RNA-seq** "
            "a partir de una matriz de conteos crudos usando **DESeq2 (Python)**.\n\n"
            "**Pasos:**\n"
            "1. Sube tu **matriz de conteos crudos** (TSV o CSV) y tus **metadatos** (CSV).\n"
            "2. Configura tus **columnas de metadatos** (ID de muestra, condicion, filtros opcionales).\n"
            "3. Haz clic en **Ejecutar expresion diferencial** para ejecutar el analisis.\n"
            "4. Visualiza resultados (Volcano, PCA, MA, Heatmap) y descarga la salida."
        ),

        # ── Seccion de carga ──────────────────────────────────────────
        "upload.header": "Sube tus archivos",
        "upload.counts_label": "Subir matriz de conteos crudos (TSV, CSV o ZIP)",
        "upload.counts_help": (
            "Archivo tabular con genes como filas y muestras como columnas. "
            "La primera columna debe ser el identificador del gen. "
            "Valores: conteos crudos (enteros >= 0). "
            "Tambien se acepta un archivo ZIP que contenga el CSV/TSV dentro."
        ),
        "upload.expander_title": "Que tipo de matriz de conteos acepta DESeq2?",
        "upload.counts_info": (
            "DESeq2 realiza su propia normalizacion interna, por lo que **requiere "
            "conteos crudos (enteros sin normalizar)**.\n\n"
            "| | Tipo | Ejemplos |\n"
            "|---|---|---|\n"
            "| **Aceptado** | Salida de featureCounts | Conteos enteros del alineamiento |\n"
            "| **Aceptado** | Conteos HTSeq | Salida de `htseq-count` |\n"
            "| **Aceptado** | Conteos Salmon (convertidos) | Conteos agregados con `tximport` |\n"
            "| **Aceptado** | Matrices GEO crudas | Serie de matrices con conteos crudos |\n"
            "| **Rechazado** | FPKM | Fragmentos por kilobase por millon |\n"
            "| **Rechazado** | TPM | Transcritos por millon |\n"
            "| **Rechazado** | CPM | Conteos por millon |\n"
            "| **Rechazado** | Log-normalizado | `log2(CPM+1)`, `rlog`, `vst`, etc. |\n\n"
            "> **Como saberlo?** Los conteos crudos son **enteros no negativos** (0, 1, 5, 342, ...).\n"
            "> Si tu matriz tiene muchos valores decimales (3.72, 0.0041, 8.156...), es "
            "probable que este normalizada y **no deberia** usarse directamente con DESeq2."
        ),
        "upload.metadata_label": "Subir archivo de metadatos (CSV o ZIP)",
        "upload.metadata_help": (
            "Archivo CSV con informacion de las muestras. "
            "Debe tener una columna con los nombres de muestras que coincidan "
            "con las columnas de la matriz de conteos, y al menos una columna de "
            "condicion experimental. Puedes seleccionar las columnas despues."
        ),

        # ── Limite de tamano ────────────────────────────────────────────
        "upload.file_too_large": (
            "**Archivo demasiado grande:** `{file}` pesa {size_mb} MB, "
            "pero el limite para esta herramienta es **{limit_mb} MB**."
        ),

        # ── Errores de lectura ────────────────────────────────────────
        "error.reading_files": "Error al leer archivos: {error}",

        # ── Previsualizacion ──────────────────────────────────────────
        "preview.counts_title": "Vista previa: Matriz de conteos",
        "preview.counts_caption": "{genes} genes x {samples} muestras",
        "preview.genes_removed": (
            "{count:,} genes con cero conteos en todas las "
            "muestras fueron eliminados automaticamente."
        ),
        "warning.normalized_matrix": (
            "**Posible matriz normalizada detectada**\n\n"
            "{reason}\n\n"
            "Si esta matriz contiene FPKM, TPM, CPM o "
            "valores log-transformados, los resultados de DESeq2 seran "
            "**estadisticamente invalidos**. Por favor usa conteos crudos.\n\n"
            '_Consulta la seccion "Que tipo de matriz de conteos acepta '
            'DESeq2?" arriba para mas detalles._'
        ),
        "preview.metadata_title": "Vista previa: Metadatos",
        "preview.metadata_caption": "{rows} filas x {cols} columnas",

        # ── Configurar metadatos ──────────────────────────────────────
        "config.header": "Configurar columnas de metadatos",
        "config.instructions": (
            "Selecciona que columnas usar de tu archivo de metadatos:\n"
            "- **Columna de muestras**: La columna con nombres de muestras (deben coincidir con las columnas de la matriz de conteos).\n"
            "- **Columna de condicion**: La columna que define las condiciones experimentales (ej. 'tratamiento', 'enfermedad').\n"
            "- **Columnas de filtro** (opcional): Columnas adicionales para filtrar muestras antes del analisis."
        ),
        "config.sample_col_label": "Columna de muestras (identificadores)",
        "config.sample_col_help": (
            "Esta columna debe contener los nombres de muestras que "
            "coincidan con las columnas de tu matriz de conteos."
        ),
        "config.condition_col_label": "Columna de condicion (grupos experimentales)",
        "config.condition_col_help": (
            "Esta columna define las condiciones experimentales "
            "para el analisis de expresion diferencial."
        ),

        # ── Filtros (multiples) ───────────────────────────────────────
        "filter.subheader": "Opcional: Filtrar muestras",
        "filter.checkbox_label": "Aplicar filtros para seleccionar muestras especificas",
        "filter.checkbox_help": "Activa esto para filtrar muestras segun columnas de metadatos.",
        "filter.title": "Filtro {n}",
        "filter.column_label": "Columna",
        "filter.values_label": "Valores a incluir de '{column}'",
        "filter.values_help": "Solo las muestras con estos valores seran incluidas.",
        "filter.no_values_warning": (
            "No se seleccionaron valores para el filtro {n} "
            "-- todas las muestras seran excluidas!"
        ),
        "filter.add_button": "Agregar otro filtro",
        "filter.summary": (
            "Despues del filtrado: **{count}** muestras "
            "restantes ({n_filters} filtro(s) aplicado(s))"
        ),
        "filter.all_removed_warning": "Todas las muestras fueron filtradas! Ajusta tus filtros.",

        # ── Error de validacion ───────────────────────────────────────
        "error.validation": "Error de validacion: {error}",

        # ── Discrepancia de muestras ──────────────────────────────────
        "mismatch.warning": "**Discrepancia de muestras detectada**",
        "mismatch.metric_counts": "Muestras en conteos",
        "mismatch.metric_metadata": "Muestras en metadatos",
        "mismatch.metric_common": "Muestras en comun",
        "mismatch.no_common": (
            "No hay NINGUNA muestra en comun entre los archivos. "
            "No se puede continuar. Verifica que los nombres de las "
            "muestras coincidan."
        ),
        "mismatch.details_title": "Detalles: discrepancia de muestras",
        "mismatch.only_in_counts": "**Muestras en conteos pero NO en metadatos** ({count}):",
        "mismatch.only_in_metadata": "**Muestras en metadatos pero NO en conteos** ({count}):",
        "mismatch.proceed_info": (
            "**{n_common}** muestras estan presentes en ambos "
            "archivos. El analisis puede continuar usando solo estas muestras "
            "comunes. Las muestras sin coincidencia seran excluidas. "
            "Este es el comportamiento estandar en la mayoria de los pipelines de RNA-seq."
        ),
        "mismatch.proceed_checkbox": "Entendido. Continuar solo con las muestras en comun.",

        # ── Error de niveles de condicion ─────────────────────────────
        "error.single_condition": (
            "**No se puede ejecutar el analisis:** Las muestras a "
            "analizar solo tienen **{n_levels}** nivel(es) de "
            "condicion: **{levels}**.\n\n"
            "DESeq2 necesita al menos **2 condiciones diferentes** "
            "(ej. 'tumor' y 'normal') para calcular un contraste.\n\n"
            "**Que verificar:**\n"
            "- Tu matriz de conteos contiene muestras de AMBAS condiciones?\n"
            "- Tus metadatos cubren las muestras de la matriz de conteos?\n\n"
            "Actualmente analizando **{n_samples}** muestras, "
            "todas etiquetadas como **'{first_level}'**."
        ),

        # ── Parametros de analisis ────────────────────────────────────
        "params.header": "Parametros del analisis",
        "params.reference_label": "Condicion de **referencia** (linea base)",
        "params.reference_help": (
            "La condicion de referencia/linea base. "
            "Los cambios de expresion se calculan RELATIVO a este nivel. "
            "Por ejemplo, si seleccionas 'normal', un log2FC positivo "
            "significa que el gen esta MAS expresado en la condicion de prueba."
        ),
        "params.test_label": "Condicion de **prueba** (vs referencia)",
        "params.test_help_single": "Solo hay una condicion de prueba disponible.",
        "params.test_help_multi": (
            "La condicion a comparar contra la referencia. "
            "Un log2FC positivo significa que el gen esta MAS expresado "
            "en ESTA condicion comparado con la referencia."
        ),
        "params.contrast_summary": (
            "**Contraste:** {test} vs {ref} (referencia)\n\n"
            "log2FC positivo -> gen **sobre-expresado** en **{test}**\n\n"
            "log2FC negativo -> gen **sub-expresado** en **{test}**"
        ),
        "params.run_button": "Ejecutar expresion diferencial",
        "params.running": "Ejecutando analisis DESeq2... Esto puede tomar un momento.",
        "params.success": "Analisis completado! {n_genes} genes analizados.",
        "error.unexpected": "Error inesperado durante el analisis: {error}",

        # ── Filtrar resultados ────────────────────────────────────────
        "results_filter.header": "Filtrar resultados",
        "results_filter.stats": (
            "**Estadisticas de expresion:**\n"
            "- baseMean minimo: **{min_bm:.2f}**\n"
            "- baseMean mediana: **{med_bm:.2f}**\n"
            "- baseMean maximo: **{max_bm:.2f}**\n\n"
            "Los genes con baja expresion suelen tener estadisticas poco confiables. "
            "Es practica comun filtrarlos."
        ),
        "results_filter.radio_label": "Selecciona estrategia de filtrado:",
        "results_filter.radio_default": "Filtros por defecto (recomendado)",
        "results_filter.radio_custom": "Filtros personalizados",
        "results_filter.radio_none": "Sin filtros (mostrar todos los genes)",
        "results_filter.radio_help": "Elige como filtrar los resultados antes de la visualizacion.",
        "results_filter.default_info": (
            "**Filtros por defecto:**\n"
            "- baseMean >= {basemean}\n"
            "- padj < {padj}\n"
            "- |log2FoldChange| > {log2fc}"
        ),
        "results_filter.basemean_help": "Filtrar genes con baja expresion (baseMean debajo de este valor).",
        "results_filter.padj_help": "Filtrar por genes estadisticamente significativos.",
        "results_filter.log2fc_help": "Filtrar por cambios de expresion biologicamente relevantes.",
        "results_filter.no_filter_warning": (
            "Sin filtros aplicados. Mostrando todos los genes incluyendo aquellos con "
            "baja expresion o p-valores no significativos."
        ),
        "results_filter.summary": (
            "**Resultados filtrados:** {n_filtered:,} genes "
            "({n_removed:,} eliminados de {n_total:,} total)"
        ),

        # ── Tabla de resultados ───────────────────────────────────────
        "table.subheader": "Tabla de resultados",
        "table.no_genes": "Ningun gen paso los filtros actuales. Intenta relajar los umbrales.",
        "table.caption": "Mostrando los 50 principales de {n:,} genes (ordenados por padj)",

        # ── Ranking de genes ──────────────────────────────────────────
        "ranking.subheader": "Ranking de genes",
        "ranking.description": (
            "La tabla de ranking contiene genes ordenados por su estadistico de Wald (`stat`). "
            "Valores positivos indican sobre-expresion, valores negativos indican sub-expresion."
        ),
        "ranking.id_label": "Nombre columna ID del gen",
        "ranking.id_help": "Nombre para la columna identificadora del gen en la salida del ranking.",
        "ranking.top_n_label": "Mostrar top N genes",
        "ranking.top_n_help": "Numero de genes top a mostrar en el ranking.",
        "ranking.caption": "Mostrando top {n} de {total:,} genes (ordenados por stat)",

        # ── Descarga ──────────────────────────────────────────────────
        "download.subheader": "Descargar resultados",
        "download.filename_label": "Nombre del archivo (sin extension)",
        "download.filename_help": "Ingresa el nombre para tu archivo de salida.",
        "download.format_label": "Formato",
        "download.format_help": "CSV usa separador de coma, TSV usa separador de tabulacion.",
        "download.btn_filtered": "Filtrados ({fmt})",
        "download.btn_all": "Todos los resultados ({fmt})",
        "download.btn_ranking": "Ranking ({fmt})",

        # ── Visualizaciones ───────────────────────────────────────────
        "viz.subheader": "Visualizaciones",
        "viz.label_checkbox": "Etiquetar top 10 genes (por ranking stat) en Volcano Plots",
        "viz.label_help": "Mostrar nombres de genes para los 10 principales segun estadistico de ranking.",
        "viz.tab_volcano_all": "Volcano Plot (Todos)",
        "viz.tab_volcano_filt": "Volcano Plot (Filtrado)",
        "viz.tab_pca": "PCA",
        "viz.tab_ma": "MA Plot",
        "viz.tab_heatmap": "Heatmap",
        "viz.volcano_title": "Volcano Plot -- Expresion Diferencial",
        "viz.no_genes_volcano": "No hay genes para mostrar con los filtros actuales.",
        "viz.filtered_caption": (
            "Mostrando {n:,} genes despues del filtrado "
            "(padj < {padj}, |log2FC| > {log2fc})"
        ),
        "viz.basemean_filter_label": "baseMean minimo (filtro de ruido)",
        "viz.basemean_filter_help": (
            "Excluir genes con expresion media debajo de este umbral. "
            "Elimina ruido de baja expresion del volcano plot."
        ),
        "viz.legend_position_label": "Posicion de leyenda",
        "viz.legend_position_help": (
            "Elige donde colocar la leyenda en los plots volcano, PCA y MA."
        ),
        "viz.download_png": "Descargar PNG",
        "viz.download_svg": "Descargar SVG",

        # ── Correccion de batch (PCA) ────────────────────────────────
        "pca.batch_checkbox": "Habilitar correccion de batch (ComBat) para PCA",
        "pca.batch_checkbox_help": (
            "Aplica ComBat para remover efectos de batch de la matriz "
            "de expresion transformada antes del PCA.  Solo afecta la "
            "visualizacion, NO los resultados de expresion diferencial."
        ),
        "pca.batch_select_label": "Columna de batch",
        "pca.batch_select_help": (
            "Selecciona la columna del metadata que identifica el batch "
            "(ej. corrida de secuenciacion, placa, centro)."
        ),

        # ── Highlight ─────────────────────────────────────────────────
        "highlight.header": "Resaltar genes de interes",
        "highlight.instructions": (
            "Ingresa uno o mas nombres/IDs de genes para resaltarlos en el volcano plot. "
            "Los genes resaltados se mostraran como **diamantes naranjas** con etiquetas, "
            "mientras que el resto de genes se atenuan en el fondo."
        ),
        "highlight.textarea_label": "Nombres/IDs de genes (uno por linea, o separados por coma)",
        "highlight.textarea_placeholder": "ej.\nBRCA1\nTP53\nMYC",
        "highlight.textarea_help": "Ingresa identificadores de genes exactamente como aparecen en tu matriz de conteos.",
        "highlight.not_found": (
            "{count} gen(es) no encontrado(s) en los resultados: **{genes}**"
        ),
        "highlight.not_found_more": " y {extra} mas...",
        "highlight.found": "{count} gen(es) encontrado(s): **{genes}**",
        "highlight.details_title": "Detalles de los genes resaltados",
        "highlight.download_button": "Descargar volcano plot resaltado (PNG)",
        "highlight.none_found": "Ninguno de los genes ingresados fue encontrado en los resultados.",

        # ── Landing page ──────────────────────────────────────────────
        "landing.title": "BeginSeq Studio",
        "landing.subtitle": "Sin necesidad de programar — selecciona una herramienta para comenzar tu analisis.",
        "landing.bulk_title": "Bulk RNA-seq Clasico",
        "landing.bulk_desc": "Analisis de expresion diferencial a partir de matrices de conteos crudos usando DESeq2.",
        "landing.scrna_title": "Analisis scRNA-seq",
        "landing.scrna_desc": "Analisis scRNA-seq: QC, clustering, UMAP y genes marcadores con Scanpy.",
        "landing.dataset_title": "Creador de Datasets",
        "landing.dataset_desc": "Crea y prepara datasets para analisis posteriores.",
        "landing.coming_soon": "Proximamente",
        "landing.open_tool": "Abrir herramienta",
        "landing.disclaimer": (
            "Disenado para principiantes y usuarios sin experiencia en programacion. "
            "Solo para fines de investigacion y educacion — no destinado a uso clinico ni diagnostico. "
            "Se proporciona \"tal cual\", sin garantia de ningun tipo. "
            "Consulta el README para el aviso legal completo."
        ),

        # ── Paginas placeholder ───────────────────────────────────────
        "placeholder.title_scrna": "Analisis scRNA-seq",
        "placeholder.title_dataset": "Creador de Datasets",
        "placeholder.message": "Esta herramienta esta en desarrollo. Vuelve pronto!",
        "placeholder.scrna_preview": (
            "**Funcionalidades planeadas:**\n"
            "- Control de calidad y filtrado\n"
            "- Reduccion de dimensionalidad (PCA, UMAP, t-SNE)\n"
            "- Clustering y anotacion de tipos celulares\n"
            "- Expresion diferencial entre clusters"
        ),
        "placeholder.dataset_preview": (
            "**Funcionalidades planeadas:**\n"
            "- Combinar multiples matrices de conteos\n"
            "- Generar plantillas de metadatos\n"
            "- Conversion de formatos (CSV, TSV, H5AD)\n"
            "- Reportes de control de calidad"
        ),

        # ── Progreso del pipeline DESeq2 ──────────────────────────────
        "progress.validating": "Validando datos...",
        "progress.filtering_genes": "Filtrando genes de baja expresion...",
        "progress.size_factors": "Calculando factores de normalizacion...",
        "progress.genewise_disp": "Estimando dispersiones gen a gen...",
        "progress.disp_trend": "Ajustando curva de tendencia de dispersiones...",
        "progress.map_disp": "Ajustando dispersiones MAP...",
        "progress.fit_lfc": "Estimando log-fold changes...",
        "progress.cooks": "Calculando distancias de Cook...",
        "progress.deseq2_done": "Modelo DESeq2 ajustado!",
        "progress.wald_test": "Ejecutando test de Wald y ajuste de p-valores...",
        "progress.wald_done": "Test estadistico completado!",
        "progress.visualizations": "Generando graficos (Volcano, PCA, MA, Heatmap)...",
        "progress.done": "Analisis completo! ✅",
        "progress.time_estimate": "Estimando tiempo de ejecucion...",

        # ── Estimacion de tiempo ───────────────────────────────────────
        "time.estimate_header": "Tiempo estimado",
        "time.estimate_detail": (
            "Dataset: **{n_samples:,}** muestras × **{n_genes:,}** genes — "
            "Tiempo estimado: **{estimate}**"
        ),
        "time.estimate_detail_filtered": (
            "Dataset: **{n_samples:,}** muestras × **{n_genes_raw:,}** genes "
            "(~**{n_genes_est:,}** despues del filtrado) — "
            "Tiempo estimado: **{estimate}**"
        ),
        "time.elapsed": "Transcurrido: {elapsed}",
        "time.remaining": "~{remaining} restante",
        "time.step_with_eta": "{step_msg} — ~{remaining} restante",
        "time.completed_in": "Analisis completado en **{elapsed}**! ✅",
        "time.faster_than_expected": "Completado mas rapido de lo esperado!",
        "time.step_breakdown": "Desglose de tiempos por paso",
        "time.step_name_size_factors": "Factores de normalizacion",
        "time.step_name_genewise_disp": "Dispersiones gen a gen",
        "time.step_name_disp_trend": "Tendencia de dispersiones",
        "time.step_name_map_disp": "Dispersiones MAP",
        "time.step_name_fit_lfc": "Log-fold changes",
        "time.step_name_cooks": "Distancias de Cook",
        "time.step_name_wald_test": "Test de Wald",
        "time.step_name_visualizations": "Visualizaciones",

        # ── Clasificación basada en expresión ─────────────────────────
        "classify.subheader": "Clasificacion de Muestras por Expresion",
        "classify.description": (
            "Clasifica muestras en sub-grupos basado en la expresion de genes marcadores. "
            "Esto crea **condiciones compuestas** (ej., tumor_TNBC vs tumor_nonTNBC) "
            "para un analisis de expresion diferencial mas especifico.\n\n"
            "**Como funciona:** Normalizacion ligera CPM → log2(CPM+1) → "
            "clasificar muestras por umbrales de genes marcadores → fusionar en metadatos."
        ),
        "classify.checkbox": "Habilitar clasificacion por expresion",
        "classify.checkbox_help": (
            "Clasifica muestras en sub-grupos basado en expresion de genes marcadores. "
            "Se ejecuta ANTES del analisis de expresion diferencial."
        ),
        "classify.gene_input_label": "Genes marcadores (pegar lista)",
        "classify.gene_input_placeholder": (
            'ej.\nESR1, PGR, ERBB2\no pega una lista separada por comas:\n'
            '"ACTR3B","ANLN","BAG1","BCL2",...'
        ),
        "classify.gene_input_help": (
            "Ingresa nombres/IDs de genes separados por comas, saltos de linea o punto y coma. "
            "Las comillas se eliminan automaticamente. Los genes no encontrados en la matriz "
            "seran reportados pero ignorados."
        ),
        "classify.no_genes_warning": "Ingresa al menos un gen marcador para continuar.",
        "classify.genes_not_found_summary": (
            "**{n_found}** genes encontrados, **{n_not_found}** no estan en la matriz: {not_found_list}"
        ),
        "classify.genes_found_summary": "**{n}** genes marcadores encontrados en la matriz de conteos.",
        "classify.stats_expander": "Estadisticas de expresion de {n} genes marcadores (log2 CPM+1)",
        "classify.zscore_caption": (
            "**z(thr)** = (umbral − media) / sd — "
            "Cuantas desviaciones estandar esta el umbral respecto a la media. "
            "z negativo → umbral por debajo del promedio; z positivo → por encima."
        ),
        "classify.global_params_title": "Umbral de clasificacion (aplicado a TODOS los marcadores)",
        "classify.global_threshold_label": "Umbral log2(CPM+1)",
        "classify.global_direction_label": "Direccion",
        "classify.threshold_help": (
            "Umbral en log2(CPM+1) aplicado a TODOS los genes marcadores. "
            "Los valores de expresion estan en escala log2 "
            "(ej., 1.0 ≈ 2 CPM, 3.3 ≈ 10 CPM, 6.6 ≈ 100 CPM)."
        ),
        "classify.direction_below": "Debajo (baja expresion → positivo)",
        "classify.direction_above": "Arriba (alta expresion → positivo)",
        "classify.direction_help": (
            "'Debajo': muestras con TODOS los marcadores POR DEBAJO del umbral se clasifican como positivo. "
            "'Arriba': muestras con TODOS los marcadores POR ENCIMA del umbral se clasifican como positivo."
        ),
        "classify.positive_label": "Etiqueta de clase positiva",
        "classify.positive_help": (
            "Etiqueta para muestras que cumplen TODAS las reglas "
            "(ej., 'TNBC', 'AltaExpr', 'Mutante')."
        ),
        "classify.negative_label": "Etiqueta de clase negativa",
        "classify.negative_help": (
            "Etiqueta para muestras que NO cumplen todas las reglas "
            "(ej., 'nonTNBC', 'BajaExpr', 'TipoSilvestre')."
        ),
        "classify.reference_keep_label": "Mantener muestras de referencia sin cambios",
        "classify.reference_keep_help": (
            "Si esta habilitado, las muestras de la condicion de referencia (ej., 'control') "
            "mantienen su etiqueta original en lugar de ser reclasificadas. "
            "Solo las muestras no-referencia obtienen etiquetas compuestas."
        ),
        "classify.preview_title": "Vista previa de clasificacion",
        "classify.compound_preview": "Condiciones compuestas despues de la clasificacion:",
        "classify.error": "Error de clasificacion: {error}",
        "classify.download_metadata": "Descargar metadata clasificado ({fmt})",

        # ── Advertencias pre-análisis ────────────────────────────────
        "warning.class_imbalance": (
            "**Desbalance severo detectado ({ratio}:1):** "
            "El grupo '**{large_group}**' tiene **{large_n}** muestras mientras que "
            "'**{small_group}**' solo tiene **{small_n}**. "
            "Esto puede reducir la potencia estadística y sesgar los resultados."
        ),
        "warning.sparse_data": (
            "**Datos dispersos detectados:** Cada gen tiene al menos un cero en las muestras. "
            "Se activará automáticamente la normalización '**poscounts**' "
            "en lugar del método 'ratio' por defecto (que fallaría o congelaría el análisis).\n\n"
            "**Implicaciones:**\n"
            "• Cambia la geometría de los size factors.\n"
            "• Puede sesgar ligeramente los Fold Changes hacia cero.\n"
            "• Es necesario para evitar que el análisis se congele con tus datos."
        ),

        # ── Validación estricta (reporte) ─────────────────────────────
        "strict.report_title": "Reporte de Validación de Calidad",
        "strict.errors": "Errores",
        "strict.warnings": "Advertencias",
        "strict.info_items": "Info",
        "strict.details": "Detalles",
        "strict.download_report": "Descargar reporte de validación (JSON)",

        # ── Registro de auditoría ─────────────────────────────────────
        "audit.section_title": "Registro de Auditoría Científica",
        "audit.description": (
            "Un registro legible por máquina de todos los parámetros, "
            "versiones de librerías y detalles de ejecución para reproducibilidad."
        ),
        "audit.download_json": "Descargar registro de auditoría (JSON)",
        "audit.download_txt": "Descargar registro de auditoría (TXT)",

        "strict.duplicate_samples": (
            "Se encontraron {count} nombre(s) de muestra duplicados: **{samples}**. "
            "Cada muestra debe tener un identificador único."
        ),
        "strict.duplicate_gene_ids": (
            "Se encontraron {count} ID(s) de gen duplicados: **{genes}**. "
            "Los IDs duplicados causarán resultados incorrectos."
        ),
        "strict.condition_nan": (
            "{count} muestra(s) tienen valores vacíos o faltantes en la "
            "columna de condición. Elimina o completa estos valores antes de ejecutar."
        ),
        "strict.condition_numeric": (
            "La columna de condición parece contener valores numéricos "
            "(ej. {sample_values}). Considera usar etiquetas descriptivas "
            "como 'control' / 'tratado' para mayor claridad."
        ),
        "strict.alignment_high_loss": (
            "**{pct_lost:.0f}% de las muestras** se perdieron durante la "
            "alineación entre counts y metadata. Verifica que los nombres coincidan."
        ),
        "strict.alignment_case_mismatch": (
            "Posible diferencia de mayúsculas/minúsculas entre nombres de "
            "muestras en counts y metadata: {pairs}. Los nombres se comparan "
            "respetando mayúsculas."
        ),
        "strict.alignment_dropped": (
            "{n_dropped} muestra(s) presentes en metadata pero ausentes en counts: "
            "**{samples}**."
        ),
        "strict.normalized_data": (
            "Los datos parecen estar normalizados (no son counts crudos): {reason}. "
            "DESeq2 requiere **counts enteros crudos** como entrada."
        ),
        "strict.class_imbalance": (
            "Desbalance de clases detectado ({ratio}:1 entre el grupo más "
            "grande y el más pequeño). Esto puede reducir la potencia estadística."
        ),
        "strict.sparse_data": (
            "Datos dispersos: se detectó una alta proporción de ceros. "
            "Se usará automáticamente la normalización 'poscounts'."
        ),
        "strict.post_filter_low": (
            "Solo quedan **{n_after}** genes después del filtrado (de {n_before}). "
            "Considera relajar los umbrales del filtro."
        ),
        "strict.post_filter_critical": (
            "Solo quedan **{n_after}** genes después del filtrado (de {n_before}) — "
            "esto es muy poco para un análisis de expresión diferencial confiable."
        ),
        "strict.filter_param_aggressive": (
            "El parámetro min_samples_expressing ({min_samples}) excede el "
            "90% del total de muestras ({n_samples}). Esto puede ser demasiado agresivo."
        ),
        "strict.passed": "Todas las validaciones pasaron correctamente.",

        # ── Pre-filtrado de genes ──────────────────────────────────
        "gene_filter.subheader": "Pre-filtrado de genes",
        "gene_filter.description": (
            "Filtra genes con baja expresion **antes** de ejecutar DESeq2. "
            "Esta es una practica estandar en analisis RNA-seq — reduce el uso "
            "de memoria, mejora la potencia estadistica y estabiliza las estimaciones de dispersion."
        ),
        "gene_filter.checkbox": "Aplicar pre-filtrado de genes (recomendado)",
        "gene_filter.checkbox_help": (
            "Elimina genes con expresion muy baja en todas las muestras. "
            "Tipicamente remueve 40-60% de los genes (ruido) sin perder "
            "señales biologicamente relevantes."
        ),
        "gene_filter.min_total_label": "Conteo total minimo por gen",
        "gene_filter.min_total_help": (
            "Un gen debe tener al menos este numero total de conteos "
            "(suma en TODAS las muestras) para ser retenido. Default: 10."
        ),
        "gene_filter.min_samples_label": "Min muestras expresando",
        "gene_filter.min_samples_help": (
            "Un gen debe estar expresado en al menos esta cantidad de muestras. "
            "Coloca 0 para automatico: max(grupo mas pequeño, 50% de muestras). "
            "Para datasets grandes (800+ muestras), se recomienda automatico."
        ),
        "gene_filter.min_count_label": "Conteo min por muestra",
        "gene_filter.min_count_help": (
            "Conteo minimo en una muestra para considerar el gen como "
            "'expresado' en esa muestra. Default: 1 (cualquier conteo >0). "
            "Valores mas altos (5-10) dan un filtrado mas agresivo."
        ),
        "gene_filter.stats_summary": (
            "Filtrado de genes: {n_before:,} → {n_after:,} genes "
            "({n_removed:,} removidos). Criterios: conteo total ≥ {min_total}, "
            "expresado (≥ {min_count}/muestra) en ≥ {min_samples} muestras."
        ),

        # ── Dataset Creator (GDC) ───────────────────────────────────
        "dc.title": "Creador de Datasets GDC",
        "dc.subtitle": (
            "Descarga datasets RNA-seq STAR-Counts del **Genomic Data Commons (GDC/TCGA)** "
            "y ensambla matrices de conteos listas para DESeq2 con metadatos auto-generados."
        ),
        "dc.step1_title": "Cargar Proyectos TCGA",
        "dc.step1_desc": (
            "Conecta con la API de GDC para obtener la lista de proyectos TCGA disponibles."
        ),
        "dc.load_projects": "Cargar Proyectos TCGA",
        "dc.loading_projects": "Conectando con la API de GDC...",
        "dc.no_projects": "No se encontraron proyectos TCGA.",
        "dc.projects_loaded": "{count} proyectos TCGA cargados.",
        "dc.view_projects": "Ver todos los proyectos",
        "dc.error_api": "Error de la API GDC: {error}",

        "dc.step2_title": "Seleccionar Proyecto y Buscar Archivos",
        "dc.select_project": "Selecciona un proyecto TCGA",
        "dc.fetch_files": "Buscar Archivos RNA-seq",
        "dc.fetching_files": "Buscando archivos RNA-seq para {project}...",
        "dc.no_files": "No se encontraron archivos STAR-Counts RNA-seq para {project}.",
        "dc.files_found": "{count} archivos RNA-seq encontrados para {project}.",
        "dc.metric_total": "Total archivos",
        "dc.metric_tumor": "Tumor",
        "dc.metric_normal": "Normal",
        "dc.metric_size": "Tamano est.",
        "dc.view_files": "Ver detalles de archivos",

        "dc.step3_title": "Filtrar y Configurar",
        "dc.filter_samples": "Filtrar por condicion",
        "dc.select_conditions": "Selecciona condiciones a incluir",
        "dc.no_conditions_selected": "Selecciona al menos una condicion para continuar.",
        "dc.filtered_summary": (
            "Seleccionados: **{count}** archivos ({tumor} Tumor, {normal} Normal)"
        ),
        "dc.gene_id_label": "Tipo de identificador de gen",
        "dc.gene_id_help": (
            "**gene_name**: Simbolo legible (ej. TP53, BRCA1). "
            "**gene_id**: ID Ensembl (ej. ENSG00000141510)."
        ),
        "dc.count_col_label": "Columna de conteos",
        "dc.count_col_help": (
            "**unstranded** (recomendado): Conteos totales sin importar la hebra. "
            "Usa las opciones stranded solo si tu libreria fue strand-specific."
        ),

        "dc.step4_title": "Descargar y Construir Dataset",
        "dc.download_info": (
            "Listo para descargar **{n_files}** archivos (~{size} MB). "
            "Los archivos seran descargados, parseados y ensamblados en una matriz de conteos."
        ),
        "dc.start_download": "Descargar y Construir Dataset",
        "dc.downloading": "Descargando archivos de GDC...",
        "dc.download_progress": "Descargando archivo {current}/{total}...",
        "dc.downloaded_count": "{count} archivos descargados exitosamente.",
        "dc.no_files_downloaded": "No se pudieron descargar archivos. Verifica tu conexion.",
        "dc.error_download": "Error de descarga: {error}",
        "dc.building_matrix": "Construyendo matriz de conteos...",
        "dc.parse_progress": "Parseando archivo {current}/{total}...",
        "dc.error_building": "Error al construir la matriz de conteos: {error}",
        "dc.done": "Listo!",
        "dc.build_success": (
            "Matriz de conteos ensamblada: **{genes:,}** genes x **{samples}** muestras."
        ),

        "dc.step5_title": "Vista Previa y Descarga de Resultados",
        "dc.preview_counts": "Matriz de Conteos",
        "dc.preview_counts_caption": "{genes:,} genes x {samples} muestras (mostrando primeras 20 filas)",
        "dc.preview_metadata": "Metadatos",
        "dc.condition_summary": "Resumen de Condiciones",
        "dc.download_results": "Descargar Archivos",
        "dc.download_counts": "Descargar Matriz de Conteos (CSV)",
        "dc.download_metadata": "Descargar Metadatos (CSV)",
        "dc.tip_bulk": (
            "**Consejo:** Puedes usar los archivos descargados directamente en la "
            "herramienta **Bulk RNA-seq Clasico** para analisis de expresion diferencial!"
        ),

        # ── Controles de limite de muestras ────────────────────────
        "dc.sample_limit_title": "Limitar numero de muestras",
        "dc.sample_limit_desc": (
            "Para proyectos grandes (ej. TCGA-BRCA con 1,000+ archivos), "
            "puedes limitar el numero de muestras a descargar. "
            "Las muestras se seleccionan aleatoriamente manteniendo las proporciones de condicion."
        ),
        "dc.use_sample_limit": "Limitar numero de muestras a descargar",
        "dc.limit_mode_label": "Modo de limite",
        "dc.limit_mode_total": "Total de muestras",
        "dc.limit_mode_per_cond": "Por condicion",
        "dc.limit_total_slider": "Maximo total de muestras",
        "dc.limit_per_cond_slider": "Maximo de muestras por condicion",
        "dc.limited_summary": (
            "Despues del limite: **{count}** muestras ({tumor} Tumor, {normal} Normal)"
        ),
        "dc.large_download_warning": (
            "Estas a punto de descargar una gran cantidad de archivos. "
            "Considera usar la opcion de limite de muestras arriba para reducir el tiempo de descarga."
        ),
        "dc.download_dir_label": "Directorio de descarga",
        "dc.download_dir_help": (
            "Directorio donde se descargaran y extraeran los archivos de GDC. "
            "Deja vacio para usar un directorio temporal (los archivos se eliminan al construir la matriz)."
        ),
        "dc.download_dir_persistent_info": (
            "Los archivos se guardaran en: `{path}`. "
            "Persistiran despues de construir la matriz."
        ),

        # ── Analisis scRNA-seq ──────────────────────────────────────────
        "scrna.title": "Analisis scRNA-seq",
        "scrna.subtitle": (
            "Pipeline de analisis de RNA-seq de celula unica con **Scanpy**. "
            "Sube tus datos y ejecuta el flujo completo: QC, normalizacion, "
            "reduccion de dimensionalidad, clustering e identificacion de genes marcadores."
        ),

        # ── Glosario AnnData ────────────────────────────────────────────
        "scrna.glossary_title": "Referencia del Objeto AnnData",
        "scrna.glossary_desc": (
            "Todos los datos de scRNA-seq en Scanpy se almacenan en un objeto **AnnData**. "
            "Esta tabla describe sus componentes principales:"
        ),

        # ── Integrador de archivos 10x ──────────────────────────────────
        "scrna.integrator_title": "Integrador de archivos 10x — Generar H5AD desde archivos crudos",
        "scrna.integrator_desc": (
            "Si tienes archivos crudos de 10x Genomics (**matrix.mtx**, "
            "**features.tsv** o **genes.tsv**, y **barcodes.tsv**), "
            "subrelos aqui para generar un archivo **.h5ad**. "
            "Se aceptan archivos comprimidos (`.gz`). "
            "Archivos con solo 2 columnas (gene ID + nombre) se manejan automaticamente."
        ),
        "scrna.integrator_matrix_label": "matrix.mtx(.gz)",
        "scrna.integrator_matrix_help": "La matriz de expresion dispersa en formato Matrix Market.",
        "scrna.integrator_features_label": "features.tsv(.gz) o genes.tsv(.gz)",
        "scrna.integrator_features_help": (
            "Archivo de anotaciones de genes. Puede tener 2 columnas (gene_id, gene_name) "
            "o 3 columnas (gene_id, gene_name, feature_type). "
            "Archivos de 2 columnas se corrigen automaticamente."
        ),
        "scrna.integrator_barcodes_label": "barcodes.tsv(.gz)",
        "scrna.integrator_barcodes_help": "Archivo de barcodes celulares, un barcode por linea.",
        "scrna.integrator_run": "Generar H5AD",
        "scrna.integrator_running": "Integrando archivos 10x en H5AD...",
        "scrna.integrator_success": "Archivo H5AD generado exitosamente! Descargalo abajo y luego subelo para ejecutar el analisis.",
        "scrna.integrator_filename_label": "Nombre del archivo (sin extensión)",
        "scrna.integrator_filename_help": "Elige un nombre para el archivo .h5ad de salida.",
        "scrna.integrator_download": "Descargar .h5ad",
        "scrna.integrator_error": "Error al integrar archivos: {error}",

        # ── Integrador: metadata opcional ──────────────────────────────
        "scrna.integrator_meta_title": "Opcional: Agregar metadata",
        "scrna.integrator_meta_desc": (
            "Sube un archivo CSV o TSV con anotaciones adicionales. "
            "El sistema detectara automaticamente si coincide con "
            "barcodes de celulas (\u2192 adata.obs) o IDs de genes (\u2192 adata.var) "
            "basandose en la superposicion de indices."
        ),
        "scrna.integrator_meta_label": "Archivo de metadata (.csv, .tsv)",
        "scrna.integrator_meta_help": (
            "CSV o TSV con la primera columna como identificadores "
            "(barcodes celulares o IDs de genes). Las columnas restantes "
            "se fusionaran en el objeto AnnData automaticamente."
        ),
        "scrna.integrator_meta_auto_note": (
            "La primera columna se comparara con los barcodes (obs) "
            "y los nombres de genes (var). Se selecciona automaticamente "
            "la mejor coincidencia. Puedes sobreescribir abajo si es necesario."
        ),
        "scrna.integrator_meta_target_label": "Fusionar metadata en:",
        "scrna.integrator_meta_target_auto": "Auto-detectar (recomendado)",
        "scrna.integrator_meta_target_obs": "adata.obs (anotaciones de celulas)",
        "scrna.integrator_meta_target_var": "adata.var (anotaciones de genes)",
        "scrna.integrator_meta_merged_obs": (
            "Metadata fusionada en **adata.obs** (anotaciones de celulas). "
            "Coincidieron **{overlap}** de {total} IDs ({pct}% superposicion con barcodes)."
        ),
        "scrna.integrator_meta_merged_var": (
            "Metadata fusionada en **adata.var** (anotaciones de genes). "
            "Coincidieron **{overlap}** de {total} IDs ({pct}% superposicion con genes)."
        ),
        "scrna.integrator_meta_no_match": (
            "No se pudo determinar el destino de la metadata \u2014 menos de 10% "
            "de superposicion con barcodes y nombres de genes. "
            "Intenta seleccionar el destino manualmente arriba."
        ),

        # ── SoupX Limpieza de RNA Ambiental ────────────────────────────
        "scrna.soupx_title": "Limpieza de RNA Ambiental — SoupX",
        "scrna.soupx_disclaimer": (
            "**Nota:** CellBender (basado en deep learning) generalmente proporciona una "
            "remocion de RNA ambiental superior, pero requiere GPU y recursos computacionales "
            "sustanciales que estan fuera del alcance de este programa. "
            "SoupX es una alternativa ligera que se ejecuta completamente en CPU via R."
        ),
        "scrna.soupx_desc": (
            "Remueve la contaminacion por RNA ambiental de tus datos scRNA-seq usando **SoupX** (paquete R via rpy2). "
            "SoupX estima el perfil de contaminacion a partir de las gotas **crudas (sin filtrar)** "
            "y lo sustrae de las gotas **filtradas (que contienen celulas)**.\n\n"
            "**Necesitas dos archivos H5AD:**\n"
            "- **Crudo / sin filtrar**: Todas las gotas incluyendo vacias (ej. `raw_feature_bc_matrix`)\n"
            "- **Filtrado**: Solo gotas con celulas (ej. `filtered_feature_bc_matrix`)\n\n"
            "La salida limpia puede subirse abajo para el pipeline completo de analisis."
        ),
        "scrna.soupx_missing_rpy2": (
            "**rpy2** no esta instalado. Instalalo con: `pip install rpy2`"
        ),
        "scrna.soupx_missing_r": (
            "**R** no se encontro en este sistema. SoupX requiere que R este instalado. "
            "Descarga R desde: https://cran.r-project.org/"
        ),
        "scrna.soupx_missing_pkg": (
            "**SoupX** paquete de R no esta instalado. Abre R y ejecuta: "
            "`install.packages('SoupX')`"
        ),
        "scrna.soupx_ready": "R + rpy2 + SoupX detectados — listo para limpiar RNA ambiental.",
        "scrna.soupx_raw_label": "H5AD crudo (sin filtrar)",
        "scrna.soupx_raw_help": (
            "H5AD con TODAS las gotas (raw_feature_bc_matrix). "
            "Incluye gotas vacias usadas para estimar el perfil de RNA ambiental."
        ),
        "scrna.soupx_filt_label": "H5AD filtrado",
        "scrna.soupx_filt_help": (
            "H5AD con solo gotas que contienen celulas (filtered_feature_bc_matrix). "
            "Estos son los datos que seran limpiados."
        ),
        "scrna.soupx_auto_label": "Estimacion automatica de contaminacion",
        "scrna.soupx_auto_help": (
            "Permite que SoupX estime automaticamente la fraccion de contaminacion. "
            "Desmarca para configurarla manualmente (util si la estimacion automatica falla)."
        ),
        "scrna.soupx_contam_label": "Fraccion de contaminacion",
        "scrna.soupx_contam_help": (
            "Override manual de la fraccion de contaminacion. "
            "Valores tipicos: 0.01-0.20 (1%-20%). "
            "Valores mas altos = correccion mas agresiva."
        ),
        "scrna.soupx_run": "Limpiar RNA Ambiental (SoupX)",
        "scrna.soupx_running": "Ejecutando remocion de RNA ambiental con SoupX... Esto puede tardar algunos minutos.",
        "scrna.soupx_success": (
            "RNA ambiental removido exitosamente! Descarga el H5AD limpio abajo, "
            "luego subelo para ejecutar el analisis."
        ),
        "scrna.soupx_download": "Descargar soupx_cleaned.h5ad",
        "scrna.soupx_error": "Error de SoupX: {error}",

        # ── Carga de datos ──────────────────────────────────────────────
        "scrna.upload_header": "Sube tus datos",
        "scrna.upload_label": "Subir datos scRNA-seq (.h5ad)",
        "scrna.upload_help": (
            "Sube un archivo **.h5ad** (formato Scanpy/AnnData).\n\n"
            "Si tienes archivos crudos de 10x (matrix.mtx, features.tsv, barcodes.tsv), "
            "usa el **Integrador de archivos 10x** de arriba para generar un .h5ad primero."
        ),
        "scrna.data_loaded": (
            "Datos cargados: **{n_cells:,}** celulas x **{n_genes:,}** genes"
        ),
        "scrna.preview_title": "Vista Previa de Datos",
        "scrna.preview_obs_title": "Metadata de Células (adata.obs)",
        "scrna.preview_obs_caption": "{rows:,} células × {cols} columnas",
        "scrna.preview_obs_empty": "No se encontró metadata de células en este archivo.",
        "scrna.preview_var_title": "Metadata de Genes (adata.var)",
        "scrna.preview_var_caption": "{rows:,} genes × {cols} columnas",
        "scrna.preview_var_empty": "No se encontró metadata de genes en este archivo.",

        # ── Parametros QC ───────────────────────────────────────────────
        "scrna.qc_header": "Parametros de Control de Calidad",
        "scrna.qc_description": (
            "Filtra celulas de baja calidad y genes poco detectados. "
            "Ajusta los umbrales segun tus violin plots de QC."
        ),
        "scrna.min_genes_label": "Min genes por celula",
        "scrna.min_genes_help": "Elimina celulas con menos genes que este umbral.",
        "scrna.max_genes_label": "Max genes por celula",
        "scrna.max_genes_help": (
            "Elimina celulas con demasiados genes (posibles dobletes o multipletes)."
        ),
        "scrna.min_counts_label": "Min conteos totales por celula",
        "scrna.min_counts_help": "Elimina celulas con muy pocos conteos totales.",
        "scrna.max_counts_label": "Max conteos totales por celula",
        "scrna.max_counts_help": "Elimina celulas con conteos anormalmente altos (dobletes).",
        "scrna.max_pct_mt_label": "Max % mitocondrial",
        "scrna.max_pct_mt_help": (
            "Celulas con alto porcentaje de genes mitocondriales suelen estar danadas o estresadas. "
            "Umbral tipico: 10-20%."
        ),
        "scrna.min_cells_label": "Min celulas por gen",
        "scrna.min_cells_help": "Elimina genes detectados en menos celulas que este umbral.",

        # ── Parametros del pipeline ──────────────────────────────────────
        "scrna.params_header": "Parametros del Analisis",
        "scrna.n_hvg_label": "Numero de genes altamente variables",
        "scrna.n_hvg_help": (
            "Top N genes variables usados para PCA y analisis posteriores. "
            "Estandar: 2000. Aumentar para datasets complejos."
        ),
        "scrna.n_pcs_label": "Numero de PCs para vecinos",
        "scrna.n_pcs_help": (
            "Numero de componentes principales usados para el grafo de vecindad. "
            "Revisa el elbow plot para elegir un valor apropiado."
        ),
        "scrna.n_neighbors_label": "Numero de vecinos",
        "scrna.n_neighbors_help": (
            "Numero de vecinos mas cercanos para el grafo. "
            "Mayor = clusters mas suaves. Menor = mayor resolucion."
        ),
        "scrna.umap_min_dist_label": "Distancia minima UMAP",
        "scrna.umap_min_dist_help": (
            "Controla que tan compactamente UMAP agrupa los puntos. "
            "Valores bajos = clusters mas compactos. Rango: 0.0 - 1.0."
        ),
        "scrna.leiden_res_label": "Resolucion de Leiden",
        "scrna.leiden_res_help": (
            "Controla la granularidad del clustering. "
            "Mayor = mas clusters. Rango tipico: 0.1 - 2.0."
        ),
        "scrna.de_method_label": "Metodo DE para marcadores",
        "scrna.de_method_help": "Test estadistico para identificar genes marcadores.",
        "scrna.n_marker_genes_label": "Genes marcadores por cluster",
        "scrna.n_marker_genes_help": "Numero de genes marcadores top a calcular por cluster.",
        "scrna.doublet_checkbox": "Habilitar deteccion de dobletes (Scrublet)",
        "scrna.doublet_help": (
            "Detecta y remueve dobletes predichos usando Scrublet. "
            "Recomendado para datos de 10x Genomics."
        ),

        # ── Corrección de Efectos de Lote ────────────────────────────────
        "scrna.batch_header": "Corrección de Efectos de Lote",
        "scrna.batch_checkbox": "Habilitar corrección de lote (Harmony)",
        "scrna.batch_help": (
            "Corrige efectos de lote usando Harmony, que ajusta el embedding PCA "
            "para que las células de diferentes lotes se integren correctamente. "
            "Habilitar solo si los datos contienen múltiples lotes."
        ),
        "scrna.batch_col_label": "Columna de variable de lote",
        "scrna.batch_col_help": (
            "Selecciona la columna en adata.obs que identifica el lote "
            "(por ejemplo, 'sample', 'donor', 'experiment')."
        ),
        "scrna.batch_preview": "Encontrados **{n_batches}** lotes en la columna `{col}`",
        "scrna.batch_no_columns": "No se encontraron columnas de lote adecuadas.",

        # ── Columnas de Anotación ──────────────────────────────────────────
        "scrna.annotation_header": "Columnas de Anotación",
        "scrna.celltype_col_label": "Columna de tipo celular (UMAP/PCA)",
        "scrna.celltype_col_help": (
            "Selecciona una columna de adata.obs con etiquetas de tipo celular. "
            "Se usará como color predeterminado en los gráficos UMAP y PCA. "
            "Deja '(none)' para usar clusters Leiden."
        ),
        "scrna.marker_groupby_label": "Agrupar marcadores por",
        "scrna.marker_groupby_help": (
            "Columna usada para agrupar células en el análisis de marcadores, "
            "DotPlots y Heatmaps. Por defecto: clusters Leiden."
        ),

        # ── Opciones de Visualización ──────────────────────────────────────
        "scrna.viz_header": "Opciones de Visualización",
        "scrna.legend_position_label": "Posición de leyenda",
        "scrna.legend_position_help": "Elige dónde colocar la leyenda en los gráficos PCA y UMAP.",

        # ── Ejecucion ───────────────────────────────────────────────────
        "scrna.run_button": "Ejecutar Analisis scRNA-seq",
        "scrna.running": "Ejecutando pipeline scRNA-seq...",

        # ── Pasos de progreso ───────────────────────────────────────────
        "scrna.step.qc_annotation": "Anotando metricas de QC (genes MT, Ribo, HB)...",
        "scrna.step.filtering": "Filtrando celulas y genes...",
        "scrna.step.doublet_detection": "Detectando dobletes (Scrublet)...",
        "scrna.step.doublet_removal": "Removiendo dobletes predichos...",
        "scrna.step.normalization": "Normalizando y log-transformando...",
        "scrna.step.hvg_selection": "Seleccionando genes altamente variables...",
        "scrna.step.scaling": "Escalando datos...",
        "scrna.step.pca": "Ejecutando PCA...",
        "scrna.step.batch_correction": "Corrigiendo efectos de lote (Harmony)...",
        "scrna.step.neighbors": "Calculando grafo de vecindad...",
        "scrna.step.umap": "Calculando embedding UMAP...",
        "scrna.step.clustering": "Ejecutando clustering Leiden...",
        "scrna.step.marker_genes": "Identificando genes marcadores...",
        "scrna.step.done": "Analisis completo!",

        # ── Resultados ──────────────────────────────────────────────────
        "scrna.results_header": "Resultados",
        "scrna.stats_title": "Resumen del Analisis",
        "scrna.stats_cells": "Celulas",
        "scrna.stats_genes": "Genes",
        "scrna.stats_clusters": "Clusters",
        "scrna.stats_hvg": "HVGs",

        "scrna.filter_stats": (
            "Filtrado: **{cells_before:,}** -> **{cells_after:,}** celulas "
            "({cells_removed:,} removidas), "
            "**{genes_before:,}** -> **{genes_after:,}** genes "
            "({genes_removed:,} removidos)"
        ),
        "scrna.doublet_stats": (
            "Dobletes detectados: **{n_doublets}** "
            "({doublet_rate:.1%} tasa de dobletes)"
        ),
        "scrna.hvg_stats": (
            "Genes altamente variables: **{n_hvg:,}** de **{n_total:,}** totales"
        ),
        "scrna.harmony_stats": (
            "Corrección de lote (Harmony): **{n_batches}** lotes "
            "de la columna `{batch_key}`, usando **{n_pcs}** PCs"
        ),
        "scrna.leiden_stats": (
            "Clustering Leiden: **{n_clusters}** clusters "
            "(resolucion = {resolution})"
        ),

        # ── Pestanas de visualizacion ───────────────────────────────────
        "scrna.tab_qc": "Metricas QC",
        "scrna.tab_hvg": "Genes Variables",
        "scrna.tab_pca": "PCA",
        "scrna.tab_umap": "UMAP",
        "scrna.tab_markers": "Genes Marcadores",

        "scrna.umap_color_label": "Colorear UMAP por",
        "scrna.umap_color_help": "Selecciona una columna de metadata o gen para colorear el UMAP.",
        "scrna.gene_search_label": "Buscar gen",
        "scrna.gene_search_help": "Escribe un nombre de gen para visualizar su expresion en UMAP.",
        "scrna.gene_not_found": "Gen '{gene}' no encontrado en el dataset.",
        "scrna.marker_cluster_label": "Seleccionar cluster",
        "scrna.marker_n_genes_label": "Top N genes a mostrar",

        # ── Descargas ───────────────────────────────────────────────────
        "scrna.download_header": "Descargar Resultados",
        "scrna.download_h5ad": "Descargar H5AD (analisis completo)",
        "scrna.download_h5ad_help": (
            "Descarga el objeto AnnData completo (.h5ad) con todos los "
            "resultados del analisis, embeddings y metadatos."
        ),
        "scrna.download_markers_csv": "Descargar genes marcadores (CSV)",
        "scrna.download_obs_csv": "Descargar metadatos celulares (CSV)",
        "scrna.download_png": "Descargar PNG",
        "scrna.download_svg": "Descargar SVG",

        # ── Errores ─────────────────────────────────────────────────────
        "scrna.error_loading": "Error al cargar datos: {error}",
        "scrna.error_pipeline": "Error durante el analisis: {error}",
        "scrna.error_lapack": (
            "⚠️ Se detectó una inestabilidad numérica durante el análisis "
            "(matriz casi singular). Esto puede ocurrir con ciertos datasets."
        ),
        "scrna.error_lapack_hint": (
            "💡 **Sugerencias:** Intenta reducir el número de PCs, desactivar "
            "la corrección de lote, cambiar la cantidad de HVGs, o ajustar "
            "los filtros de QC para eliminar células/genes de baja calidad."
        ),
        "scrna.error_no_cells": (
            "No quedan celulas despues del filtrado. Intenta relajar los umbrales de QC."
        ),

        # ── Acceso remoto ───────────────────────────────────────────────
        "scrna.remote_blocked": (
            "Esta funcion solo esta disponible cuando se ejecuta localmente. "
            "Si deseas utilizar esta funcion, visita: "
            "https://github.com/JACKNINES/beginseq-studio"
        ),

        # ── Títulos de página (pestaña del navegador) ──────────────────
        "page_title.bulk": "Bulk RNA-seq",
        "page_title.scrna": "Análisis scRNA-seq",
        "page_title.dataset": "Creador de Datasets",

        # ── Strings hardcodeados — Bulk RNA-seq ────────────────────────
        "bulk.step_caption": "Paso {step}/{total}: {msg}",
        "bulk.params_updated_classification": "(actualizado con clasificación)",
        "bulk.more_genes_suffix": " (+{extra} más)",
        "bulk.timing_total": "Total: {elapsed}",
        "bulk.basemean_filter_label": "baseMean ≥",
        "bulk.padj_filter_label": "padj <",
        "bulk.log2fc_filter_label": "|log2FC| >",

        # ── Strings hardcodeados — página scRNA-seq ────────────────────
        "scrna.cells_x_genes": "{n_cells:,} células × {n_genes:,} genes",
        "scrna.step_caption": "Paso {step}/{total}: {msg}",
        "scrna.option_none": "(ninguno)",
        "scrna.no_markers_info": "No se calcularon genes marcadores.",
        "scrna.gene_search_placeholder": "ej. CD3D, MS4A1, NKG7",
        "scrna.dl_suffix_scatter": "(Dispersión)",
        "scrna.dl_suffix_elbow": "(Codo)",
        "scrna.dl_suffix_pca": "(PCA)",
        "scrna.dl_suffix_markers": "(Marcadores)",
        "scrna.dl_suffix_dotplot": "(Dotplot)",
        "scrna.dl_suffix_heatmap": "(Mapa de calor)",

        # ── Tabla glosario AnnData (página scRNA-seq) ─────────────────
        "scrna.glossary_col_element": "Elemento",
        "scrna.glossary_col_stands_for": "Qué representa",
        "scrna.glossary_col_contains": "Qué contiene",
        "scrna.glossary_x_name": "Matriz de expresión",
        "scrna.glossary_obs_name": "Observaciones (células)",
        "scrna.glossary_var_name": "Variables (genes)",
        "scrna.glossary_obsm_name": "Multidimensional de observaciones",
        "scrna.glossary_uns_name": "Anotaciones no estructuradas",
        "scrna.glossary_layers_name": "Capas de matrices alternativas",
        "scrna.glossary_x_desc": (
            "La matriz de expresión génica (células x genes). "
            "Valores normalizados/log-transformados después del preprocesamiento."
        ),
        "scrna.glossary_obs_desc": (
            "Metadatos de cada célula: etiquetas de cluster, métricas QC "
            "(n_genes, total_counts, pct_mt), IDs de muestra, tipos celulares."
        ),
        "scrna.glossary_var_desc": (
            "Metadatos de cada gen: nombres de genes, flag highly_variable, "
            "expresión media, dispersiones."
        ),
        "scrna.glossary_obsm_desc": (
            "Embeddings y coordenadas a nivel celular: PCA (X_pca), "
            "UMAP (X_umap), t-SNE (X_tsne)."
        ),
        "scrna.glossary_uns_desc": (
            "Resultados generales del análisis: parámetros de clustering, "
            "rankings de genes marcadores, paletas de colores, configuraciones."
        ),
        "scrna.glossary_layers_desc": (
            "Representaciones alternativas de X: conteos crudos ('counts'), "
            "datos escalados. Misma dimensión que X."
        ),

        # ── Etiquetas de gráficos — visualización scRNA ───────────────
        "plot.qc_title": "Métricas de Control de Calidad",
        "plot.qc_scatter_title": "Dispersión QC",
        "plot.hvg_title": "Genes Altamente Variables",
        "plot.hvg_label": "Altamente Variable",
        "plot.hvg_other_label": "Otros",
        "plot.hvg_xlabel": "Expresión media",
        "plot.hvg_ylabel_disp": "Dispersión normalizada",
        "plot.hvg_ylabel_count": "Número de genes",
        "plot.pca_elbow_title": "Ratio de Varianza PCA (Gráfico de Codo)",
        "plot.pca_elbow_xlabel": "Componente Principal",
        "plot.pca_elbow_ylabel": "Ratio de Varianza",
        "plot.pca_embed_title": "Embedding PCA",
        "plot.umap_xlabel": "UMAP1",
        "plot.umap_ylabel": "UMAP2",
        "plot.umap_title_prefix": "UMAP — {color}",
        "plot.marker_ranking_title": "Top Genes Marcadores por {groupby}",
        "plot.marker_score_label": "Puntuación",
        "plot.marker_no_data": "Sin datos",
        "plot.marker_no_genes": "No se encontraron genes marcadores",
        "plot.dotplot_xlabel": "Genes",
        "plot.dotplot_ylabel": "Cluster",
        "plot.dotplot_title": "Dot Plot — Genes Marcadores",
        "plot.dotplot_cbar_label": "Expresión media",
        "plot.dotplot_size_legend": "% expresando",
        "plot.heatmap_no_genes": "No hay genes marcadores para mostrar",
        "plot.heatmap_cluster_title": "Cluster",
        "plot.heatmap_xlabel": "Células (ordenadas por cluster)",
        "plot.heatmap_cbar_label": "Expresión",
        "plot.heatmap_title": "Mapa de Calor de Genes Marcadores",

        # ── Etiquetas de gráficos — visualización Bulk ─────────────────
        "plot.not_significant": "No significativo",
        "plot.up_regulated": "Sobre-expresado",
        "plot.down_regulated": "Sub-expresado",
        "plot.volcano_summary": (
            "Total genes: {n_total:,}\n"
            "Significativos: {n_sig:,} ({n_up:,} sobre · {n_down:,} sub)"
        ),
        "plot.volcano_shrunk_suffix": "  [contraído]",
        "plot.highlight_other_genes": "Otros genes ({count:,})",
        "plot.highlight_highlighted": "Resaltados ({count:,})",
        "plot.highlight_title": "Volcano — Genes Resaltados",
        "plot.highlight_info": "Resaltados: {n_found}",
        "plot.highlight_not_found": "No encontrados: {n_not_found}",
        "plot.pca_axis_variance": "PC{n} ({pct:.1f}% varianza)",
        "plot.ma_not_significant": "No significativo",
        "plot.ma_up_regulated": "Sobre-expresado",
        "plot.ma_down_regulated": "Sub-expresado",
        "plot.ma_summary": (
            "Total: {n_total:,} genes\n"
            "Significativos: {n_sig:,} ({n_up:,} sobre · {n_down:,} sub)"
        ),
        "plot.heatmap_zscore_label": "Z-score",
        "plot.heatmap_condition_title": "Condición",
    },
}



# ══════════════════════════════════════════════════════════════════════
# Pure lookup functions
# ══════════════════════════════════════════════════════════════════════

def get_label(key: str, lang: str = "en") -> str:
    """Return translated string for key in the given language.

    Falls back to English, then to the key itself.
    """
    text = TRANSLATIONS.get(lang, {}).get(key)
    if text is not None:
        return text
    text = TRANSLATIONS.get("en", {}).get(key)
    if text is not None:
        return text
    return f"[{key}]"


def get_label_formatted(key: str, lang: str = "en", **kwargs) -> str:
    """Return translated string with format placeholders replaced."""
    template = get_label(key, lang)
    try:
        return template.format(**kwargs)
    except (KeyError, IndexError, ValueError):
        return template
