# PRJEB86646 RefSeq functional annotations

Descriptions and gene biotypes come directly from GCF_963853765.1 / RS_2024_06. This is annotation and ranking, not ontology or pathway enrichment.

## Heat Higher Fdr Lt 0.05

Genes: 4,542; RefSeq-characterized descriptions: 3,014 (66.4%).

- LOC136273757 (NCBI:GeneID:136273757): heat shock protein beta-1-like (log2FC=2.46; adjusted p=2.94e-46)
- LOC105323167 (NCBI:GeneID:105323167): small ribosomal subunit protein eS27-like (log2FC=0.58; adjusted p=2.83e-44)
- LOC136269995 (NCBI:GeneID:136269995): caspase-7-like (log2FC=1.03; adjusted p=1.07e-36)
- LOC105329404 (NCBI:GeneID:105329404): lysosomal phospholipase A and acyltransferase (log2FC=0.90; adjusted p=1.14e-35)
- LOC105327885 (NCBI:GeneID:105327885): GTP-binding nuclear protein Ran (log2FC=0.65; adjusted p=4e-35)
- LOC105335774 (NCBI:GeneID:105335774): multiple epidermal growth factor-like domains protein 10 (log2FC=3.41; adjusted p=5.87e-35)
- LOC105347892 (NCBI:GeneID:105347892): hillarin (log2FC=2.81; adjusted p=1.45e-34)
- LOC105343175 (NCBI:GeneID:105343175): F-box only protein 43 (log2FC=3.25; adjusted p=2.94e-33)
- LOC105325905 (NCBI:GeneID:105325905): developmentally-regulated GTP-binding protein 1 (log2FC=0.71; adjusted p=1.69e-28)
- LOC105326642 (NCBI:GeneID:105326642): growth hormone-inducible transmembrane protein (log2FC=0.84; adjusted p=2.25e-27)
- LOC105336253 (NCBI:GeneID:105336253): prestin (log2FC=1.35; adjusted p=4.27e-27)
- LOC105342003 (NCBI:GeneID:105342003): inactive serine/threonine-protein kinase TEX14 (log2FC=3.27; adjusted p=5.52e-27)
- LOC105345295 (NCBI:GeneID:105345295): toll-like receptor 4 (log2FC=0.81; adjusted p=1.01e-26)
- LOC105347370 (NCBI:GeneID:105347370): protein Mpv17 (log2FC=1.01; adjusted p=1.44e-26)
- LOC117680708 (NCBI:GeneID:117680708): serine/arginine-rich splicing factor 1 (log2FC=0.71; adjusted p=2.06e-26)

## Heat Lower Fdr Lt 0.05

Genes: 4,482; RefSeq-characterized descriptions: 3,027 (67.5%).

- LOC105345603 (NCBI:GeneID:105345603): putative peptidyl-prolyl cis-trans isomerase (log2FC=-1.82; adjusted p=7.56e-46)
- LOC117684197 (NCBI:GeneID:117684197): PE-PGRS family protein PE_PGRS47-like (log2FC=-3.80; adjusted p=1.33e-44)
- LOC105329390 (NCBI:GeneID:105329390): asparagine synthetase [glutamine-hydrolyzing] (log2FC=-2.10; adjusted p=8.07e-41)
- LOC105327353 (NCBI:GeneID:105327353): natterin-4-like (log2FC=-1.67; adjusted p=6.83e-31)
- LOC105341700 (NCBI:GeneID:105341700): glutamine--fructose-6-phosphate aminotransferase [isomerizing] 2 (log2FC=-0.88; adjusted p=1.16e-28)
- LOC105325018 (NCBI:GeneID:105325018): senecionine N-oxygenase (log2FC=-1.45; adjusted p=1.77e-28)
- LOC105333151 (NCBI:GeneID:105333151): 1,4-alpha-glucan-branching enzyme (log2FC=-0.96; adjusted p=8.97e-27)
- LOC105329816 (NCBI:GeneID:105329816): tyrosine-protein kinase Fyn (log2FC=-0.83; adjusted p=3.33e-26)
- LOC105323004 (NCBI:GeneID:105323004): receptor-type tyrosine-protein phosphatase epsilon (log2FC=-0.71; adjusted p=3.33e-26)
- LOC105324538 (NCBI:GeneID:105324538): aspartyl/asparaginyl beta-hydroxylase (log2FC=-0.92; adjusted p=1.85e-25)
- LOC105327045 (NCBI:GeneID:105327045): microfibril-associated glycoprotein 4-like (log2FC=-1.06; adjusted p=4.94e-25)
- LOC105326753 (NCBI:GeneID:105326753): complement C1q subcomponent subunit A (log2FC=-1.27; adjusted p=5.52e-24)
- LOC105348479 (NCBI:GeneID:105348479): low-density lipoprotein receptor (log2FC=-1.33; adjusted p=5.7e-24)
- LOC105332789 (NCBI:GeneID:105332789): cytochrome P450 26A1 (log2FC=-1.20; adjusted p=1.4e-23)
- LOC136271677 (NCBI:GeneID:136271677): uromodulin-like (log2FC=-1.23; adjusted p=1.78e-21)

## Gene-biotype over-representation

One-sided hypergeometric tests compare each significant direction with all genes having DESeq2 results; Benjamini-Hochberg adjustment covers both directions and all biotypes.

- heat higher fdr lt 0.05: protein_coding (4203 observed, 4006.5 expected; fold=1.05; adjusted p=4.38e-25)
- heat higher fdr lt 0.05: rRNA (10 observed, 2.7 expected; fold=3.73; adjusted p=0.000165)
- heat lower fdr lt 0.05: protein_coding (4157 observed, 3953.6 expected; fold=1.05; adjusted p=6.15e-27)

Many oyster loci retain automated or uncharacterized RefSeq labels. Functional claims should be based on formal ontology/pathway enrichment or orthology-aware analysis rather than names alone.
