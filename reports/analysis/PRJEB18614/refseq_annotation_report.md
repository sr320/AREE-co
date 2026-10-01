# PRJEB18614 RefSeq functional annotations

Descriptions and gene biotypes come directly from GCF_963853765.1 / RS_2024_06. This is annotation and ranking, not ontology or pathway enrichment.

## Alexandrium Higher Fdr Lt 0.05

Genes: 3,690; RefSeq-characterized descriptions: 3,132 (84.9%).

- LOC105339566 (NCBI:GeneID:105339566): dipeptidyl peptidase 1 (log2FC=1.51; adjusted p=1.99e-282)
- LOC105348980 (NCBI:GeneID:105348980): solute carrier family 15 member 4 (log2FC=1.77; adjusted p=3.63e-194)
- LOC105338718 (NCBI:GeneID:105338718): lysosome-associated membrane glycoprotein 1 (log2FC=1.00; adjusted p=4.59e-169)
- LOC105342870 (NCBI:GeneID:105342870): ATP-dependent translocase ABCB1 (log2FC=3.66; adjusted p=1.26e-153)
- LOC105345248 (NCBI:GeneID:105345248): CD63 antigen (log2FC=1.10; adjusted p=2.92e-137)
- LOC105346410 (NCBI:GeneID:105346410): arylacetamide deacetylase (log2FC=2.32; adjusted p=1.43e-131)
- LOC105337461 (NCBI:GeneID:105337461): maestro heat-like repeat-containing protein family member 1 (log2FC=1.29; adjusted p=2.09e-127)
- LOC105322438 (NCBI:GeneID:105322438): glycosylated lysosomal membrane protein (log2FC=1.31; adjusted p=6.91e-126)
- LOC105331606 (NCBI:GeneID:105331606): vacuolar protein sorting-associated protein 4B (log2FC=2.07; adjusted p=6.69e-125)
- LOC105325023 (NCBI:GeneID:105325023): nose resistant to fluoxetine protein 6 (log2FC=6.05; adjusted p=7.27e-123)
- LOC105321776 (NCBI:GeneID:105321776): syntaxin-8 (log2FC=0.92; adjusted p=1.01e-109)
- LOC105332185 (NCBI:GeneID:105332185): phospholipid-transporting ATPase ABCA3 (log2FC=1.47; adjusted p=2.81e-109)
- LOC105322704 (NCBI:GeneID:105322704): alternative oxidase, mitochondrial-like (log2FC=3.75; adjusted p=1.7e-107)
- LOC105331219 (NCBI:GeneID:105331219): vitamin D3 receptor A (log2FC=1.94; adjusted p=1.34e-104)
- LOC105329718 (NCBI:GeneID:105329718): neutral amino acid uniporter 4 (log2FC=1.82; adjusted p=1.3e-101)

## Alexandrium Lower Fdr Lt 0.05

Genes: 3,290; RefSeq-characterized descriptions: 2,032 (61.8%).

- LOC105335432 (NCBI:GeneID:105335432): ammonium transporter Rh type A (log2FC=-1.40; adjusted p=3.48e-40)
- LOC117681463 (NCBI:GeneID:117681463): toll-like receptor 4 (log2FC=-1.21; adjusted p=4.93e-38)
- LOC105343328 (NCBI:GeneID:105343328): helix-loop-helix protein delilah (log2FC=-0.67; adjusted p=8.46e-33)
- LOC105338379 (NCBI:GeneID:105338379): tumor necrosis factor receptor superfamily member 19 (log2FC=-0.81; adjusted p=1.31e-31)
- LOC105339780 (NCBI:GeneID:105339780): sodium- and chloride-dependent glycine transporter 2 (log2FC=-1.89; adjusted p=2.69e-29)
- LOC105317404 (NCBI:GeneID:105317404): sodium-coupled monocarboxylate transporter 1 (log2FC=-0.48; adjusted p=2.08e-28)
- LOC105322136 (NCBI:GeneID:105322136): caveolin-3 (log2FC=-1.01; adjusted p=6.65e-27)
- LOC105344169 (NCBI:GeneID:105344169): sodium/calcium exchanger Calx (log2FC=-0.50; adjusted p=9.75e-27)
- LOC117681556 (NCBI:GeneID:117681556): blastula protease 10-like (log2FC=-1.88; adjusted p=1.36e-26)
- LOC117680326 (NCBI:GeneID:117680326): hydroxyproline dehydrogenase-like (log2FC=-1.81; adjusted p=1.77e-25)
- LOC136271940 (NCBI:GeneID:136271940): organic cation transporter protein-like (log2FC=-2.31; adjusted p=3.64e-25)
- LOC105334801 (NCBI:GeneID:105334801): hypoxia-inducible factor 1-alpha inhibitor (log2FC=-0.42; adjusted p=3.23e-24)
- LOC105326150 (NCBI:GeneID:105326150): hillarin (log2FC=-0.46; adjusted p=1.9e-23)
- LOC105318397 (NCBI:GeneID:105318397): cilia- and flagella-associated protein 251 (log2FC=-0.96; adjusted p=5.23e-21)
- LOC105346523 (NCBI:GeneID:105346523): sodium-dependent proline transporter (log2FC=-1.59; adjusted p=6.91e-21)

## Gene-biotype over-representation

One-sided hypergeometric tests compare each significant direction with all genes having DESeq2 results; Benjamini-Hochberg adjustment covers both directions and all biotypes.

- alexandrium higher fdr lt 0.05: protein_coding (3608 observed, 3168.0 expected; fold=1.14; adjusted p=1.75e-154)
- alexandrium lower fdr lt 0.05: protein_coding (3030 observed, 2824.6 expected; fold=1.07; adjusted p=2.36e-31)

Many oyster loci retain automated or uncharacterized RefSeq labels. Functional claims should be based on formal ontology/pathway enrichment or orthology-aware analysis rather than names alone.
