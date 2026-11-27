-- AREE-co Relational Export Schema (SQLite 3)

CREATE TABLE IF NOT EXISTS studies (
    study_id TEXT PRIMARY KEY,
    bioproject TEXT,
    species TEXT,
    taxonomy_id TEXT,
    stressor_class TEXT,
    stressor TEXT,
    phenotype TEXT,
    tissue TEXT,
    life_stage TEXT,
    data_status TEXT,
    biological_replication TEXT,
    study_design TEXT,
    experimental_treatment TEXT,
    control_condition TEXT,
    phenotype_direction TEXT,
    limitations_json TEXT
);

CREATE TABLE IF NOT EXISTS refseq_references (
    accession TEXT PRIMARY KEY,
    species TEXT,
    species_code TEXT,
    ncbi_taxid INTEGER,
    assembly_id TEXT,
    annotation_release TEXT,
    tx2gene_genes INTEGER,
    tx2gene_transcripts INTEGER
);

CREATE TABLE IF NOT EXISTS harmonized_evidence (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    study_id TEXT,
    feature_id_original TEXT,
    feature_id_standardized TEXT,
    feature_type TEXT,
    gene_symbol TEXT,
    description TEXT,
    sample_comparison TEXT,
    molecular_direction TEXT,
    effect_size REAL,
    effect_size_type TEXT,
    standard_error REAL,
    p_value REAL,
    adjusted_p_value REAL,
    sample_size INTEGER,
    tissue TEXT,
    life_stage TEXT,
    stressor TEXT,
    phenotype TEXT,
    resilience_classification TEXT,
    inference_status TEXT,
    ortholog_reference TEXT,
    mapping_confidence TEXT,
    quality_flags TEXT,
    workflow_version TEXT,
    genome_assembly TEXT,
    annotation_version TEXT,
    source_file TEXT,
    input_checksum TEXT
);

CREATE INDEX IF NOT EXISTS idx_ev_study ON harmonized_evidence(study_id);
CREATE INDEX IF NOT EXISTS idx_ev_feature ON harmonized_evidence(feature_id_standardized);
CREATE INDEX IF NOT EXISTS idx_ev_fdr ON harmonized_evidence(adjusted_p_value);

CREATE TABLE IF NOT EXISTS synthesis_strata (
    stratum_id TEXT PRIMARY KEY,
    title TEXT,
    phenotype TEXT,
    stressor TEXT,
    n_contrasts INTEGER,
    n_bioprojects INTEGER,
    n_independent_groups INTEGER,
    multi_group_pooled_features INTEGER,
    pooled_fdr_lt_0_05 INTEGER,
    concordant_fdr_lt_0_05 INTEGER,
    median_i2_percent REAL
);

CREATE TABLE IF NOT EXISTS synthesis_meta (
    stratum_id TEXT,
    feature_id_standardized TEXT,
    gene_symbol TEXT,
    description TEXT,
    n_studies INTEGER,
    n_independent_groups INTEGER,
    method TEXT,
    shared_control_rho REAL,
    pooled_effect REAL,
    pooled_standard_error REAL,
    pooled_ci_low REAL,
    pooled_ci_high REAL,
    p_value REAL,
    adjusted_p_value REAL,
    i2_percent REAL,
    tau2 REAL,
    direction_concordant TEXT,
    PRIMARY KEY (stratum_id, feature_id_standardized)
);

CREATE INDEX IF NOT EXISTS idx_synth_fdr ON synthesis_meta(stratum_id, adjusted_p_value);

CREATE TABLE IF NOT EXISTS ortholog_release_summary (
    species_code TEXT PRIMARY KEY,
    species TEXT,
    source_assembly TEXT,
    total_mapped_genes INTEGER,
    one_to_one INTEGER,
    many_to_one INTEGER,
    one_to_many INTEGER,
    many_to_many_inferred INTEGER,
    single_copy_8sp_orthogroups INTEGER,
    mgigas_indexed_genes INTEGER
);
