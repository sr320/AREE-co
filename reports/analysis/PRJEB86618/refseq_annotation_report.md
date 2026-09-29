# PRJEB86618 RefSeq functional annotations

Descriptions and gene biotypes come directly from GCF_963853765.1 / RS_2024_06. This is annotation and ranking, not ontology or pathway enrichment.

## Starved Higher Fdr Lt 0.05

Genes: 5,941; RefSeq-characterized descriptions: 3,730 (62.8%).

- LOC105344085 (NCBI:GeneID:105344085): nucleolin (log2FC=1.33; adjusted p=1.27e-73)
- LOC105338579 (NCBI:GeneID:105338579): cartilage matrix protein (log2FC=1.61; adjusted p=2.9e-63)
- LOC105345822 (NCBI:GeneID:105345822): tyrosine aminotransferase (log2FC=1.09; adjusted p=1.04e-55)
- LOC105348704 (NCBI:GeneID:105348704): low-density lipoprotein receptor-related protein 6-like (log2FC=0.63; adjusted p=6.49e-55)
- LOC105343470 (NCBI:GeneID:105343470): secretin receptor (log2FC=1.23; adjusted p=1.97e-52)
- LOC105320770 (NCBI:GeneID:105320770): formin-J (log2FC=1.07; adjusted p=1.65e-51)
- LOC117681030 (NCBI:GeneID:117681030): tetraspanin-1 (log2FC=1.03; adjusted p=1.34e-50)
- LOC105339373 (NCBI:GeneID:105339373): monocarboxylate transporter 13-like (log2FC=2.57; adjusted p=9.12e-50)
- LOC105317849 (NCBI:GeneID:105317849): peroxisome proliferator-activated receptor alpha (log2FC=1.65; adjusted p=2.35e-49)
- LOC105342449 (NCBI:GeneID:105342449): neuronal acetylcholine receptor subunit alpha-2 (log2FC=1.41; adjusted p=4.29e-48)
- LOC105340979 (NCBI:GeneID:105340979): transforming growth factor beta receptor type 3 (log2FC=1.00; adjusted p=4.64e-45)
- LOC105323221 (NCBI:GeneID:105323221): dual specificity tyrosine-phosphorylation-regulated kinase 4 (log2FC=0.80; adjusted p=2.72e-44)
- LOC105338219 (NCBI:GeneID:105338219): CARD- and ANK-domain containing inflammasome adapter protein (log2FC=1.90; adjusted p=1.93e-43)
- LOC105318062 (NCBI:GeneID:105318062): alkaline phosphatase (log2FC=1.57; adjusted p=2.91e-43)
- LOC105331739 (NCBI:GeneID:105331739): fibropellin-1 (log2FC=0.98; adjusted p=9.84e-42)

## Starved Lower Fdr Lt 0.05

Genes: 5,989; RefSeq-characterized descriptions: 4,318 (72.1%).

- LOC105327731 (NCBI:GeneID:105327731): serine protease inhibitor Cvsi-1 (log2FC=-3.32; adjusted p=1.01e-134)
- LOC105318596 (NCBI:GeneID:105318596): NPC intracellular cholesterol transporter 2 (log2FC=-3.78; adjusted p=6.47e-127)
- LOC105317784 (NCBI:GeneID:105317784): serine protease 33 (log2FC=-2.59; adjusted p=1.09e-118)
- LOC136269432 (NCBI:GeneID:136269432): lectin BRA-3-like (log2FC=-2.70; adjusted p=8.14e-111)
- LOC105346974 (NCBI:GeneID:105346974): L-sorbose 1-dehydrogenase (log2FC=-3.07; adjusted p=2.34e-104)
- LOC105334885 (NCBI:GeneID:105334885): patched domain-containing protein 3 (log2FC=-3.36; adjusted p=1.43e-98)
- LOC105330227 (NCBI:GeneID:105330227): isocitrate dehydrogenase [NADP] cytoplasmic (log2FC=-1.52; adjusted p=4.42e-96)
- LOC105349013 (NCBI:GeneID:105349013): superoxide dismutase [Cu-Zn]-like (log2FC=-4.73; adjusted p=8.42e-81)
- LOC105329692 (NCBI:GeneID:105329692): peptidyl-prolyl cis-trans isomerase FKBP1A (log2FC=-1.03; adjusted p=8.48e-79)
- LOC105328600 (NCBI:GeneID:105328600): enoyl-CoA delta isomerase 2 (log2FC=-1.12; adjusted p=9.3e-73)
- LOC105322557 (NCBI:GeneID:105322557): putative ATP synthase subunit f, mitochondrial (log2FC=-0.81; adjusted p=3.6e-72)
- LOC105342529 (NCBI:GeneID:105342529): mitochondrial carnitine/acylcarnitine carrier protein (log2FC=-1.49; adjusted p=1.22e-71)
- LOC105328838 (NCBI:GeneID:105328838): prostaglandin reductase-3 (log2FC=-1.23; adjusted p=1.42e-71)
- LOC105330579 (NCBI:GeneID:105330579): epidermal retinol dehydrogenase 2 (log2FC=-1.28; adjusted p=2.38e-67)
- LOC105348748 (NCBI:GeneID:105348748): alpha-amylase (log2FC=-1.58; adjusted p=2.4e-64)

## Gene-biotype over-representation

One-sided hypergeometric tests compare each significant direction with all genes having DESeq2 results; Benjamini-Hochberg adjustment covers both directions and all biotypes.

- starved higher fdr lt 0.05: protein_coding (5485 observed, 5237.5 expected; fold=1.05; adjusted p=7.26e-32)
- starved lower fdr lt 0.05: protein_coding (5598 observed, 5279.8 expected; fold=1.06; adjusted p=8.97e-53)

Many oyster loci retain automated or uncharacterized RefSeq labels. Functional claims should be based on formal ontology/pathway enrichment or orthology-aware analysis rather than names alone.
