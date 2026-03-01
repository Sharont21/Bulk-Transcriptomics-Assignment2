# Bulk-Transcriptomics-Assignment2
BINF 6110 Assignment 2

# Introduction

During biological wine aging, flor yeast form a biofilm (velum) at the surface of the wine that drives many of the biochemical changes affecting wine composition and sensory properties [1]. In the absence of fermentable sugars, flor yeasts undergo substantial metabolic reprogramming, shifting from fermentative to oxidative metabolism and relying on ethanol and glycerol as primary carbon sources while utilizing poor nitrogen sources under stressful conditions, including high ethanol concentrations, oxidative stress, and nutrient limitation [1,2]. Velum formation represents a key adaptive trait, enabling yeast cells to access oxygen at the wine surface while providing protection from environmental stressors [1]. This process is regulated in part by increased expression of genes such as FLO11, which enhances cell-surface hydrophobicity and promotes multicellular aggregation, allowing yeast cells to migrate to the surface and form a stable biofilm during biological aging [2,3]. Although yeast cells typically die at the end of alcoholic fermentation due to sugar and oxygen depletion, biologically aged wines are transferred to barrels with an airspace that promotes oxidative growth and biofilm formation [3]. However, comparative studies suggest that not all Saccharomyces cerevisiae strains possess the capacity to form a velum, indicating that flor yeasts represent a specialized subset adapted to biofilm-associated growth during wine aging [3].

Comparative transcriptomic analysis is a powerful approach for characterizing genome-wide changes in yeast gene expression in response to environmental, genetic, and chemical perturbations. Advances in microarray and RNA sequencing technologies have generated extensive transcriptomic datasets that have substantially improved our understanding of regulatory networks and transcriptional responses to diverse stress conditions in yeast [1]. Although many studies have focused on laboratory strains under controlled conditions, transcriptomic approaches have also been successfully applied to industrial fermentation contexts, including wine, beer, and bioethanol production, highlighting their relevance for studying gene expression dynamics in biologically and economically important environments [1].

The objective of this study was to characterize gene expression changes across three stages of velum development to identify molecular pathways associated with biofilm maturation.

Several pseudoalignment-based RNA-seq quantification tools were considered for transcript abundance estimation, including Kallisto and Salmon. Salmon was selected as the primary quantification method due to its strong balance of accuracy, robustness, and computational efficiency across diverse sequencing depths and experimental conditions [4–6]. Like Kallisto, Salmon avoids computationally intensive base-level alignments by determining read compatibility with reference transcripts; however, Salmon extends this framework through a dual-phase statistical inference procedure incorporating both online and offline expectation–maximization (EM) optimization, improving abundance estimation accuracy [4,7]. In addition, Salmon implements advanced sample-specific bias correction models that account for sequence-specific, fragment GC-content, and positional biases, reducing systematic technical artifacts that can influence transcript quantification [7]. Benchmarking studies have demonstrated that Salmon achieves higher sensitivity in downstream differential expression analyses while maintaining low false discovery rates and reducing spurious isoform switching events relative to alternative methods [4,7]. Furthermore, Salmon operates without generating large intermediate alignment files and efficiently scales across multiple CPU cores, enabling rapid quantification of large RNA-seq datasets with reduced computational demands. Based on its demonstrated accuracy, bias correction capabilities, statistical rigor, and scalability, Salmon was selected as the most suitable pseudoaligner for transcript quantification in this study [4–7].

For differential expression analysis (DEA), DESeq2 was selected over edgeR despite both being widely used and well-established Bioconductor packages for analyzing count-based sequencing data, including RNA-seq, SAGE-seq, ChIP-seq, and Hi-C datasets [8]. Both tools implement normalization strategies that have been shown to outperform alternative methods, particularly in scenarios where RNA composition varies substantially across biological conditions or when highly expressed genes are present [8]. However, benchmarking studies have demonstrated that as sample sizes increase, DESeq2 performs slightly better than many competing approaches in terms of statistical power and false discovery rate (FDR) control [9]. Specifically, for RNA-seq count data following either negative binomial or log-normal distributions, DESeq2 exhibits strong FDR control, improved power, and greater stability across varying sample sizes compared to alternative methods [9]. While other approaches such as EBSeq may be advantageous in very small sample size settings, DESeq2 is recommended for experiments with moderate or larger group sizes and remains one of the most robust and statistically reliable frameworks for RNA-seq differential expression analysis [9]. Based on its demonstrated performance, statistical rigor, and widespread validation across transcriptomic studies, DESeq2 was chosen as the most appropriate method for differential expression analysis in this study [8,9].

# Methods

## Computational Environment

All computational analyses were conducted on the Alliance Canada high-performance computing (HPC) cluster using SLURM for workload management. RNA-seq data processing was performed using command-line tools in a Linux environment. Downstream statistical analyses were performed in R (v4.5.1) using Bioconductor packages.

## Data Acquisition and Quality Control
Raw RNA-sequencing data were obtained from the NCBI Sequence Read Archive (SRA) under BioProject accession PRJNA592304. Accession numbers were retrieved and downloaded using SRA Toolkit v3.0.9 [10]. The `prefetch` utility was used to download SRA files, and `fasterq-dump` was used to convert SRA files to FASTQ format. Sequencing read quality was assessed using FastQC (v0.12.1) [11]. FastQC reports were generated for each sample to evaluate per-base sequence quality, GC content, sequence duplication levels, and adapter contamination. Reports were manually inspected to confirm that read quality was sufficient for downstream pseudoalignment-based quantification.

## Transcript Quantification
Transcript-level quantification was performed using Salmon (v1.10.2) in selective-alignment mode [12]. A reference transcriptome for _Saccharomyces cerevisiae_ (R64-1-1) was used to construct the Salmon index using `salmon index`. Quantification was performed using `salmon quant` with the `--validateMappings` option enabled to improve mapping accuracy via selective alignment. Single-end reads were supplied using the `-r` argument, and library type was explicitly specified as unstranded `-l U`. Four CPU threads were used per sample (`-p 4`). Output for each sample was written to independent directories containing transcript-level abundance estimates in `quant.sf` format. 

## Data Import and Differential Expression Analysis
All statistical analyses were conducted in R using tximport (v1.36.1) and DESeq2 (v1.48.1) [13,14]. Transcript-level abundance estimates from Salmon were imported using `tximport(type = "salmon")` and summarized to gene-level counts using a transcript-to-gene mapping derived from the reference transcriptome annotation.

Differential expression analysis was performed using DESeq2. The experimental design formula was specified as:

`design = ~ condition `

where _condition_ represented biofilm developmental stage (Stage1 (early biofilm formation), Stage2 (thin biofilm) and Stage3 (mature biofilm)). Size factor normalization and dispersion estimation were performed automatically within the `DESeq()` function. Pairwise contrasts were extracted using the Wald test for:

Stage2 vs Stage1

Stage3 vs Stage1

Stage3 vs Stage2

Log2 fold-change shrinkage was applied using `lfcShrink()` with the apeglm method (apeglm v1.30.0) to improve effect size estimation [15]. Genes with a Benjamini–Hochberg adjusted p-value (FDR) < 0.05 were considered significantly differentially expressed.
A likelihood ratio test (LRT) was additionally performed using:

`DESeq(test = "LRT", reduced = ~1)`

to identify genes exhibiting significant expression changes across all developmental stages. Genes with adjusted p-values < 0.05 were retained for clustering analysis.

## Data Visualization and Clustering
Variance stabilizing transformation (VST) was applied to normalized counts using `vst()` for visualization and clustering. Principal component analysis (PCA) was conducted using `plotPCA()` to assess global sample relationships and replicate consistency [16].
Volcano plots were generated using ggplot2 (v4.0.2), visualizing shrunken log2 fold changes versus −log10 adjusted p-values [17]. MA plots were generated using `plotMA()` to visualize mean expression versus log2 fold change [18] .
Heatmaps of the top differentially expressed genes were constructed using pheatmap (v1.0.13), with row-wise Z-score scaling applied to VST-transformed expression values to highlight stage-specific expression patterns [19].
For LRT-significant genes, Z-score scaled expression matrices were clustered using k-means clustering (k = 4). Cluster-specific expression trajectories were visualized using ggplot2 by plotting mean stage-wise expression patterns with individual gene trajectories overlaid.

## Gene Ontology Enrichment Analysis
Functional enrichment analysis was conducted using clusterProfiler (v4.16.0) and the org.Sc.sgd.db (v3.21.0) annotation database [20,21]. Over-representation analysis (ORA) was performed using `enrichGO()` for Biological Process (BP) ontology terms [22].
Significantly upregulated and downregulated genes (padj < 0.05) were analyzed separately using all expressed genes as the background universe. Multiple testing correction was performed using the Benjamini–Hochberg method, and GO terms with adjusted p-values < 0.05 were considered significantly enriched. Enrichment results were visualized using dotplots generated with the enrichplot package (v1.28.4). 

## KEGG Pathway Enrichment and Gene Set Enrichment Analysis
KEGG pathway enrichment was conducted using both over-representation analysis (`enrichKEGG()`) and gene set enrichment analysis (`gseKEGG()`) in clusterProfiler [23,24]. For KEGG ORA, significantly upregulated and downregulated genes were tested against _S. cerevisiae_ pathway annotations (organism code “sce”).For GSEA, all genes were ranked by shrunken log2 fold-change values for each pairwise contrast and analyzed using `gseKEGG()` [24]. Normalized enrichment scores (NES), adjusted p-values, and enrichment curves `(gseaplot2()`) were used to evaluate pathway directionality and coordinated transcriptional shifts across developmental stages [25]. 

# Results

## Global Transcriptomic Structure Across Biofilm Development

<p align="center">
  <img width="700" height="432" alt="image" src="https://github.com/user-attachments/assets/f6b8b96b-dd90-453f-8ff2-51ed39abcf81" />
</p>

**Figure 1. Principal component analysis (PCA) of transcriptomic profiles across biofilm developmental stages.** 
PCA was performed on variance-stabilized expression values derived from the 500 most variable genes. The first principal component (PC1) explains 72% of the total variance, while the second principal component (PC2) explains 23% of the variance. Samples cluster distinctly according to developmental stage: Stage1 (early biofilm formation), Stage2 (thin biofilm), and Stage3 (mature biofilm). Each point represents one biological replicate

Principal component analysis (PCA) was performed on variance-stabilized expression values to assess overall transcriptomic structure across biofilm development. The first principal component (PC1) accounted for 72% of the total variance, while the second principal component (PC2) explained 23% (Figure 1). Samples segregated clearly according to developmental stage along PC1, with Stage1 (early biofilm formation) clustering distinctly from Stage2 (thin biofilm) and Stage3 (mature biofilm). Stage2 samples occupied an intermediate position along PC1 and were further separated from Stage3 along PC2. Biological replicates within each stage clustered closely together, indicating consistent gene expression patterns within stages and clear differences between developmental stages.

## Differential Gene Expression Between Developmental Stages

<p align="center">
  <img src="https://github.com/user-attachments/assets/3d83f5ad-28cf-4a9c-b587-6da2a57eaa58" width="32%" />
  <img src="https://github.com/user-attachments/assets/4606e0c1-c146-41b7-a9a1-6ae284a983cb" width="32%" />
  <img src="https://github.com/user-attachments/assets/47c415f4-eb0a-4c52-860c-af7a406cfd8e" width="32%" />
</p>

**Figure 2. Volcano plots of differential gene expression across biofilm developmental stages.**
Volcano plots display shrunken log2 fold change (x-axis) versus −log10 adjusted p-value (y-axis) for each pairwise comparison: (A) Stage2 vs Stage1, (B) Stage3 vs Stage1, and (C) Stage3 vs Stage2. Genes with adjusted p-values < 0.05 are colored by direction of change (red = upregulated; blue = downregulated), while non-significant genes are shown in gray.

Differential expression analysis was performed using DESeq2 with Wald tests applied to pairwise contrasts between biofilm developmental stages. Log2 fold changes were shrunken using the apeglm method to improve effect size stability. Genes with adjusted p-values < 0.05 were considered significantly differentially expressed.
In the comparison of Stage2 vs Stage1, a total of 2,486 genes were significantly differentially expressed, including 1,209 upregulated and 1,277 downregulated genes (Figure 2A). The distribution of significant genes was relatively balanced between positive and negative fold changes.
The comparison of Stage3 vs Stage1 identified 3,154 significantly differentially expressed genes, comprising 1,522 upregulated and 1,632 downregulated genes (Figure 2B). This contrast exhibited the largest number of significant genes, indicating extensive transcriptional changes between early and mature biofilm stages.
In Stage3 vs Stage2, 2,202 genes were significantly differentially expressed, including 1,086 upregulated and 1,116 downregulated genes (Figure 2C), demonstrating continued transcriptional remodeling between the thin and mature biofilm stages.
Across all comparisons, significant genes were distributed across a broad range of effect sizes, with both moderate and large magnitude log2 fold changes observed. The volcano plots illustrate the magnitude and statistical significance of gene expression changes underlying stage-dependent biofilm development.

## Stage-Dependent Expression Patterns Identified by Likelihood Ratio Test

<p align="center">
  <img width="3000" height="2100" alt="GO_Clusters_Combined" src="https://github.com/user-attachments/assets/65f6a7ed-e245-4dfe-9d6f-09f4242faf1d" />
</p>

**Figure 3. Gene Ontology (GO) Biological Process enrichment across LRT-derived gene expression clusters.**
Genes exhibiting significant stage-dependent expression patterns identified by likelihood ratio test (LRT; adjusted p-value < 0.05) were grouped into four clusters using k-means clustering (k = 4). Over-representation analysis of GO Biological Process terms was performed for each cluster using clusterProfiler. Dot size represents gene ratio, and color indicates adjusted p-value (Benjamini–Hochberg correction)

To identify genes exhibiting significant expression changes across all developmental stages, a likelihood ratio test (LRT) was performed using a reduced model (~1) to test for stage-dependent effects. A total of 3,900 genes were identified as significantly regulated across stages (adjusted p-value < 0.05).
Variance-stabilized expression values for LRT-significant genes were Z-score scaled and subjected to k-means clustering (k = 4), resulting in four distinct expression patterns comprising 728, 1,310, 851, and 1,011 genes, respectively. These clusters represented coordinated transcriptional trajectories across biofilm development.
Functional enrichment analysis of each cluster revealed distinct biological process associations (Figure 4). Cluster 1 was enriched for translation- and ribosome-related processes, including cytoplasmic translation and ribosome biogenesis. Cluster 2 was enriched for lipid metabolic and biosynthetic processes. Cluster 3 was associated with regulatory processes, including RNA metabolic regulation and mRNA decapping. Cluster 4 showed enrichment for mitochondrial-related processes, including mitochondrial translation, mitochondrial gene expression, and mitochondrial organization.
These results indicate that genes exhibiting coordinated stage-dependent expression patterns are functionally organized into distinct biological programs during biofilm development.



# Discussion


# References

[1] Mardanov, A. V., Eldarov, M. A., Beletsky, A. V., Tanashchuk, T. N., Kishkovskaya, S. A., & Ravin, N. V. (2020). Transcriptome Profile of Yeast Strain Used for Biological Wine Aging Revealed Dynamic Changes of Gene Expression in Course of Flor Development. Frontiers in Microbiology, 11, 538. https://doi.org/10.3389/fmicb.2020.00538

[2] Carbonero-Pacheco, J., Moreno-García, J., Moreno, J., García-Martínez, T., & Mauricio, J. C. (2022). Revealing the Yeast Diversity of the Flor Biofilm Microbiota in Sherry Wines Through Internal Transcribed Spacer-Metabarcoding and Matrix-Assisted Laser Desorption/Ionization Time of Flight Mass Spectrometry. Frontiers in Microbiology, 12, 825756. https://doi.org/10.3389/fmicb.2021.825756

[3] David-Vaizant, V., & Alexandre, H. (2018). Flor Yeast Diversity and Dynamics in Biologically Aged Wines. Frontiers in Microbiology, 9, 2235. https://doi.org/10.3389/fmicb.2018.02235

[4] Zhang, C., Zhang, B., Lin, L.-L., & Zhao, S. (2017). Evaluation and comparison of computational tools for RNA-seq isoform quantification. BMC Genomics, 18(1), Article 583. https://doi.org/10.1186/s12864-017-4002-1

[5] Jin, H., Wan, Y.-W., & Liu, Z. (2017). Comprehensive evaluation of RNA-seq quantification methods for linearity. BMC Bioinformatics, 18(Suppl 4), Article 117. https://doi.org/10.1186/s12859-017-1526-y

[6] de Villarreal, J. M., Kalisz, M., Piedrafita, G., Graña-Castro, O., Chondronasiou, D., Serrano, M., & Real, F. X. (2023). Pseudoalignment tools as an efficient alternative to detect repeated transposable elements in scRNAseq data. Bioinformatics (Oxford, England), 39(1). https://doi.org/10.1093/bioinformatics/btac737

[7] Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4), 417–419. https://doi.org/10.1038/nmeth.4197

[8] Varet, H., Brillet-Guéguen, L., Coppée, J.-Y., & Dillies, M.-A. (2016). SARTools: A DESeq2- and EdgeR-Based R Pipeline for Comprehensive Differential Analysis of RNA-Seq Data. PloS One, 11(6), e0157022. https://doi.org/10.1371/journal.pone.0157022

[9] Li, D., Zand, M. S., Dye, T. D., Goniewicz, M. L., Rahman, I., & Xie, Z. (2022). An evaluation of RNA-seq differential analysis methods. PloS One, 17(9), e0264246. https://doi.org/10.1371/journal.pone.0264246

[10] Ncbi. (n.d.). NCBI/SRA-Tools: Sra tools. GitHub. https://github.com/ncbi/sra-tools 

[11} Mary Piper, R. K. (2017, September 20). Quality Control Using FASTQC. Introduction to RNA-Seq using high-performance computing - ARCHIVED. https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/qc_running_fastqc.html 

[12] COMBINE-lab. (n.d.). Combine-Lab/Salmon: 🐟 🍣 🍱 highly-accurate & wicked fast transcript-level quantification from RNA-seq reads using selective alignment. GitHub. https://github.com/COMBINE-lab/salmon 

[13] Tximport: Import transcript-level abundances and estimated counts for gene-level analysis packages. RDocumentation. (n.d.). https://www.rdocumentation.org/packages/tximport/versions/1.0.3/topics/tximport 

[14] DESEQ2 R tutorial. (n.d.). https://lashlock.github.io/compbio/R_presentation.html 

[15] azhu513. (n.d.). AZHU513/APEGLM: APEGLM provides bayesian shrinkage estimators for effect sizes for a variety of GLM models, using approximation of the posterior for individual coefficients. GitHub. https://github.com/azhu513/apeglm 

[16] PLOTPCA: Sample PCA plot for transformed data. RDocumentation. (n.d.-a). https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotPCA 

[17] Grover, J. (2024, April 21). Making volcano plots with GGPLOT2. Jeff Grover. Bioinformatics Scientist. https://groverj3.github.io/articles/2024-04-21_making-volcano-plots-with-ggplot2.html 

[18] PLOTMA: Ma-plot from base means and log fold changes. RDocumentation. (n.d.-a). https://www.rdocumentation.org/packages/DESeq2/versions/1.12.3/topics/plotMA 

[19] Pheatmap: A function to draw clustered heatmaps. RDocumentation. (n.d.-a). https://www.rdocumentation.org/packages/pheatmap/versions/1.0.13/topics/pheatmap 

[20] YuLab-SMU. (n.d.). Yulab-SMU/Clusterprofiler: :bar_chart: A universal enrichment tool for interpreting OMICS data. GitHub. https://github.com/YuLab-SMU/clusterProfiler 

[21] Org.Sc.sgd.db. Bioconductor. (n.d.). https://www.bioconductor.org/packages/release/data/annotation/html/org.Sc.sgd.db.html 

[22] EnrichGO: GO enrichment analysis of a gene set. given a vector of genes, this function will return the enrichment go categories after FDR control. RDocumentation. (n.d.-a). https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enrichGO 

[23] Enrichkegg: Kegg enrichment analysis of a gene set. given a vector of genes, this function will return the enrichment Kegg categories with FDR control. RDocumentation. (n.d.-b). https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/enrichKEGG 

[24] GSEKEGG: GSEKEGG. RDocumentation. (n.d.-c). https://www.rdocumentation.org/packages/clusterProfiler/versions/3.0.4/topics/gseKEGG 

[25] Guangchuang Yu [aut,  cre] (). (2021, January 30). GSEAPLOT2: GSEAPLOT2 in enrichplot: Visualization of functional enrichment result. gseaplot2: gseaplot2 in enrichplot: Visualization of Functional Enrichment Result. https://rdrr.io/bioc/enrichplot/man/gseaplot2.html 


