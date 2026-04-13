# Bioinformatics — Lab Portfolio

Practical assignments from the **Bioinformática** course, Biomedical Engineering degree, Universidad de Valladolid (2024–2025).

The course covers the main computational tools used in modern bioinformatics, all implemented in **R** using packages from **CRAN** and **Bioconductor**. Labs span sequence analysis, genomic data, statistical modelling and machine learning applied to biological datasets.

---

## Lab 1 · Introduction to R
**Techniques:** Vector and matrix operations, data frames, statistical hypothesis testing, data visualisation.

- Construction and manipulation of vectors, matrices and lists; indexing by position and by condition.
- Working with distributions: `rnorm()`, `rpois()`, probability, quantile and density prefixes (`r`, `p`, `q`, `d`).
- Reading and writing tabular data (`.csv`, `.txt`, `.RData`, `.xls`); working directory management.
- Data frame operations: `subset()`, `cbind()`, `rbind()`, `apply()`, `lapply()`, `sapply()`.
- Descriptive statistics: `mean()`, `sd()`, `median()`, `quantile()`, `cor()`, `cov()`, `summary()`.
- Hypothesis testing: t-Student for two independent samples, chi-squared test for contingency tables, Wilcoxon signed-rank test.
- Base R visualisation: scatter plots with regression lines (`abline()`), box plots, histograms with density overlay.

**Data:** `iris` and `sleep` built-in datasets; simulated gene expression matrices.

---

## Lab 2 · Bioconductor Annotations and PubMed Queries
**Techniques:** Gene ID conversion, GO annotation, PubMed API queries, bibliometric analysis, word clouds.

- Installation and usage of Bioconductor packages: `BiocManager`, `AnnotationDbi`, `org.Hs.eg.db`, `GO.db`.
- Conversion between gene identifier systems: Entrez ID ↔ Gene Symbol, Entrez ID ↔ Ensembl, Entrez ID ↔ GO categories, using `mget()`.
- GO annotation mapping: gene-to-GO and GO-to-gene directions; handling of `NA` and multi-mapping entries.
- PubMed querying via `RISmed` (`EUtilsSummary()`, `EUtilsGet()`): field-based search syntax, metadata extraction (PMID, title, DOI, country, publication type, date).
- Citation count retrieval with `rentrez` (`entrez_summary()`, `extract_from_esummary()`).
- Data wrangling with `dplyr`: `mutate()`, `distinct()`, `group_by()`, `summarise()`, pipe operator `%>%`.
- Visualisation with `ggplot2`: bar charts by month, polar/pie charts by publication type.
- Abstract word cloud generation with `PubMedWordcloud` and `wordcloud2`.

**Example query:** BRCA1 + Breast Cancer publications in PubMed (2024).

---

## Lab 3 · Sequence Analysis and Pairwise Alignment
**Techniques:** FASTA I/O, nucleotide composition analysis, GC content, dinucleotide statistics, pairwise alignment, statistical significance of alignments.

- Reading and writing FASTA files with `seqinr`: `read.fasta()`, `write.fasta()`, `getSequence()`, `getName()`.
- Nucleotide frequency analysis: `table()` for base counts, `seqinr::GC()` for GC content.
- Confidence interval construction for proportions (95% CI); comparison of GC content between bacterial phyla (Actinobacteria vs. Proteobacteria).
- k-mer frequency analysis with `seqinr::count()` (codon-level trinucleotides).
- Over/under-representation of dinucleotides via the ρ (rho) statistic using `seqinr::rho()`.
- Pairwise sequence alignment with `Biostrings` / `pwalign`: `nucleotideSubstitutionMatrix()`, `pairwiseAlignment()` with global (Needleman-Wunsch) and local (Smith-Waterman) modes; gap opening and extension penalties.
- Protein sequence alignment with BLOSUM50/BLOSUM62 substitution matrices.
- Statistical significance assessment of alignments via multinomial randomisation: `generateSeqsWithMultinomialModel()`, empirical p-value estimation from 1000 random alignments.

**Sequences used:** RNA polymerase of *Mycobacterium tuberculosis* and *E. coli*; toy peptide sequences PAWHEAE / HEAGAWGHEE.

---

## Lab 4 · GWAS — Genome-Wide Association Studies
**Techniques:** SNP association analysis, Hardy-Weinberg equilibrium, multiple testing, GWAS visualisation.

- Data preparation with `SNPassoc::setupSNP()`: formatting SNP genotype data, allele separator handling.
- Single-SNP association testing with `association()`: codominant, dominant, recessive, overdominant and log-additive inheritance models; OR, 95% CI, p-value and AIC reporting.
- Whole-genome association scan with `WGassociation()`: simultaneous evaluation of all SNPs against quantitative and binary phenotypes (blood pressure, protein levels, case/control status).
- Maximum statistic (`maxstat()`) for robust model-free association testing.
- Hardy-Weinberg Equilibrium (HWE) testing with `tableHWE()`: exact test per SNP in cases and controls separately; Bonferroni-corrected filtering of SNPs deviating from HWE.
- Missing genotype pattern analysis with `plotMissing()`.
- Gene–environment and gene–gene interaction modelling using interaction terms in `association()`.
- GWAS visualisation with `GWASTools`: QQ-plots (`qqPlot()`) and Manhattan plots (`manhattanPlot()`).

**Data:** `SNPs` dataset from `SNPassoc` (157 individuals, 35 SNPs); simulated 1000-SNP Manhattan plot scenario.

---

## Lab 5 · Multiple Sequence Alignment and Phylogenetics
**Techniques:** MSA with ClustalW / ClustalOmega / MUSCLE, alignment visualisation, phylogenetic tree construction.

- Reading amino acid sequences from FASTA files with `Biostrings::readAAStringSet()`.
- Multiple sequence alignment with `msa`: `msa()` function using ClustalW (default), ClustalOmega and MUSCLE algorithms; parameter tuning (gap penalties, substitution matrix, clustering method, iteration limit).
- Alignment output: `print(..., show="complete")` for full alignment display; export to FASTA with `writeXStringSet()` + `unmasked()`.
- Alignment visualisation with `ggmsa`: colour schemes (`Taylor_AA`, `Chemistry_AA`), sequence logo overlay (`geom_seqlogo()`), consensus bar chart (`geom_msaBar()`).
- Phylogenetic tree construction with `phangorn`: alignment conversion to `phyDat`, distance matrix calculation with `dist.ml()`, UPGMA tree, Neighbor-Joining (NJ) tree; `plot()` visualisation of both topologies.

**Data:** Nine PH4H (phenylalanine hydroxylase) protein sequences across species (*Homo sapiens*, *Rattus norvegicus*, *Mus musculus*, *Bos taurus*, *Chromobacterium violaceum*, *Ralstonia solanacearum*, *Caulobacter crescentus*, *Pseudomonas aeruginosa*, *Rhizobium loti*).

---

## Lab 6 · Protein Structure Analysis
**Techniques:** PDB file handling, 3D structural visualisation, UniProt annotation, multi-protein comparison.

- Reading and parsing PDB files with `bio3d::read.pdb()`: atom attributes, residue indexing, secondary structure elements.
- Atom selection with `atom.select()`: by type (`"calpha"`, `"protein"`, `"water"`), by residue ID, with `inverse=TRUE`, and combined selections via `combine.select()` (AND/OR operators).
- Multi-PDB operations: `trim.pdb()` for chain extraction, `cat.pdb()` for structural concatenation, `write.pdb()` for file export.
- Interactive 3D visualisation with `NGLVieweR`: cartoon, ball+stick, ribbon and surface representations; `setSpin()` / `setRock()` animation; region-specific colouring (`sele`, `colorScheme`); residue labelling with custom format strings.
- UniProt feature retrieval with `drawProteins::get_features()` and `feature_to_dataframe()`.
- Protein domain diagrams with `drawProteins` + `ggplot2`: `draw_canvas()`, `draw_chains()`, `draw_domains()`, `draw_regions()`, `draw_motif()`, `draw_phospho()`, `draw_folding()`.
- Multi-protein comparison: simultaneous diagram of five NF-κB family members (RelA/p65, RelB, c-Rel, p50, p52) from UniProt.

**Structures used:** KRAS (4Q21), kinesin motor domain (1BG2), NF-κB p65 (Q04206) and related family members.

---

## Lab 7 · Microarray Data Analysis
**Techniques:** GEO data retrieval, differential expression analysis with `limma`, PCA, hierarchical clustering, heatmaps, volcano plots.

- GEO dataset retrieval with `GEOquery::getGEO()`; `ExpressionSet` object structure and `pData()`, `exprs()` accessors.
- Differential expression analysis pipeline with `limma`: design matrix definition with `model.matrix()`, linear model fitting with `lmFit()`, empirical Bayes moderation with `eBayes()`, gene ranking with `topTable()` (FDR adjustment, logFC and B-statistic sorting).
- Two-condition comparison (cancer vs. normal) and multi-group design (adenoma, colorectal cancer, tumour, normal) using `makeContrasts()` and `contrasts.fit()`.
- Gene selection criteria: FDR < 0.0001 and |logFC| > 2.
- Volcano plot construction: manual (`plot()` + `abline()` + `points()`) and automated (`volcanoplot()`).
- PCA of expression matrix with `prcomp()`; 2D score plot coloured by sample subtype for pattern discovery.
- Hierarchical clustering with correlation-based distance (1 − ρ/2), complete linkage; heatmap visualisation with `gplots::heatmap.2()` and `RColorBrewer` palettes; cluster extraction with `cutree()`.

**Data:** `apColonData` from `antiProfilesData` (68 colon samples, 5339 genes, 4 subtypes: adenoma, colorectal cancer, tumour, normal); GSE24460 (Affymetrix HG-U133A 2.0, drug resistance).

---

## Lab 8 · RNA-Seq Differential Expression Analysis
**Techniques:** Count data preprocessing, TMM normalisation, exact test, GLM with batch correction, volcano plot, clustering.

- Reading raw count tables with `read.table()`; `ExpressionSet`-aware data structures.
- Exploratory visualisation of pre-filtering count distributions (log2 scale) with `ggplot2` boxplots.
- Low-count gene filtering: removal of all-zero genes with `rowSums()`.
- `edgeR` pipeline: `DGEList()` object creation, TMM normalisation with `calcNormFactors()`, CPM transformation with `cpm()`, pre/post-normalisation comparison.
- Dispersion estimation: common dispersion (`estimateCommonDisp()`), tagwise dispersion (`estimateTagwiseDisp()`).
- Exact test for two-group comparison (`exactTest()`): `topTags()` for ranked gene extraction with BH adjustment; gene selection by FDR < 0.05 and |logFC| > 1; up/downregulation classification.
- Generalised Linear Model with batch correction: design matrix incorporating replicate effect (`model.matrix()`), `estimateDisp()`, `glmFit()`, likelihood ratio test (`glmLRT()`).
- Comparison of methods (Exact Test vs. GLM) via Venn diagram (`VennDiagram`): 222 genes identified by both methods.
- Results exploration: volcano plot with `ggplot2`; hierarchical clustering and heatmap of DEGs with `mixOmics::cim()`; dendrogram cutting with `cutree()`.

**Data:** `pasilla` dataset from Bioconductor (*Drosophila melanogaster*, 7 samples: 3 treated / 4 untreated, 14,599 genes; GEO: GSE18508).

---

## Lab 9 · Statistical and Machine Learning Methods
**Techniques:** Multiple testing correction, kNN, Random Forest, SVM, PCA + k-means clustering, variable importance.

- Multiple comparisons problem: simulation of 10,000 t-tests under H₀ (5.06% false positives); p-value adjustment with `p.adjust()` using Benjamini-Hochberg (FDR) method; complete elimination of false positives post-correction.
- **kNN classification** (`class::knn()`): data standardisation with `scale()`, 80/20 stratified split, k = 10; confusion matrix and performance metrics with `caret::confusionMatrix()` (accuracy 96.67%).
- **Random Forest** (`randomForest`): formula interface, 500 trees, OOB error estimation, `predict()` on test set; variable importance with `varImpPlot()` — Petal.Length and Petal.Width identified as most informative features.
- **SVM** (`e1071::svm()`): C-classification, radial kernel (γ = 0.25); classification boundary visualisation with `plot()`; support vector identification.
- **PCA + k-means clustering**: `prcomp()` with standardisation, scree plot (`fviz_screeplot()`), optimal cluster number (`fviz_nbclust()`), k-means (`kmeans()`), cluster visualisation (`fviz_cluster()`).
- Variable importance via PCA: `fviz_pca_var()` with cos² colour coding; comparison with Random Forest importance rankings.

**Data:** `iris` (150 samples, 4 features, 3 species) for all supervised and unsupervised methods; `modencodefly_eset` (*Drosophila* expression data, 150 samples) for PCA + clustering.

---

## Skills demonstrated

| Area | Tools / Methods |
|---|---|
| R programming | Base R, `dplyr`, `ggplot2`, `apply` family, data frames, lists |
| Bioconductor | `AnnotationDbi`, `org.Hs.eg.db`, `GO.db`, `Biostrings`, `edgeR`, `limma`, `msa`, `drawProteins`, `GEOquery` |
| Sequence analysis | FASTA I/O, GC content, ρ statistic, pairwise alignment (global/local), MSA, phylogenetics |
| Genomics | GWAS, SNP association, HWE, multiple testing, microarray and RNA-Seq DEA |
| Structural bioinformatics | PDB parsing, 3D visualisation, UniProt annotation, protein domain diagrams |
| Machine learning | kNN, Random Forest, SVM, PCA, k-means clustering, feature selection |
| Statistics | Hypothesis testing, bootstrap concepts, FDR correction, ROC/Youden (theoretical) |
| Language | R (all labs) |
| Data | Real omics datasets: *Drosophila* RNA-Seq, human colon microarray, GWAS genotype data, protein structures |

---

## Context

These labs are part of the 3th-year Biomedical Engineering curriculum at Universidad de Valladolid. Several methods covered here — particularly differential expression analysis, PCA and clustering — are directly relevant to my Final Year Project on hospital occupancy prediction (Hospital Universitario Río Hortega, SACYL), where similar statistical frameworks are applied to clinical time-series data.
