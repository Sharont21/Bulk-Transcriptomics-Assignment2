# Bulk-Transcriptomics-Assignment2
BINF 6110 Assignment 2

# Introduction

During biological wine aging, flor yeast form a biofilm (velum) at the surface of the wine that drives many of the biochemical changes affecting wine composition and sensory properties [1]. In the absence of fermentable sugars, flor yeasts undergo substantial metabolic reprogramming, shifting from fermentative to oxidative metabolism and relying on ethanol and glycerol as primary carbon sources while utilizing poor nitrogen sources under stressful conditions, including high ethanol concentrations, oxidative stress, and nutrient limitation [1,2]. Velum formation represents a key adaptive trait, enabling yeast cells to access oxygen at the wine surface while providing protection from environmental stressors [1]. This process is regulated in part by increased expression of genes such as FLO11, which enhances cell-surface hydrophobicity and promotes multicellular aggregation, allowing yeast cells to migrate to the surface and form a stable biofilm during biological aging [2,3]. Although yeast cells typically die at the end of alcoholic fermentation due to sugar and oxygen depletion, biologically aged wines are transferred to barrels with an airspace that promotes oxidative growth and biofilm formation [3]. However, comparative studies suggest that not all Saccharomyces cerevisiae strains possess the capacity to form a velum, indicating that flor yeasts represent a specialized subset adapted to biofilm-associated growth during wine aging [3].

Comparative transcriptomic analysis is a powerful approach for characterizing genome-wide changes in yeast gene expression in response to environmental, genetic, and chemical perturbations. Advances in microarray and RNA sequencing technologies have generated extensive transcriptomic datasets that have substantially improved our understanding of regulatory networks and transcriptional responses to diverse stress conditions in yeast [1]. Although many studies have focused on laboratory strains under controlled conditions, transcriptomic approaches have also been successfully applied to industrial fermentation contexts, including wine, beer, and bioethanol production, highlighting their relevance for studying gene expression dynamics in biologically and economically important environments [1].

The objective of this study was to characterize gene expression changes across three stages of velum development to identify molecular pathways associated with biofilm maturation.

Several pseudoalignment-based RNA-seq quantification tools were considered for transcript abundance estimation, including Kallisto and Salmon. Salmon was selected as the primary quantification method due to its strong balance of accuracy, robustness, and computational efficiency across diverse sequencing depths and experimental conditions [4–6]. Like Kallisto, Salmon avoids computationally intensive base-level alignments by determining read compatibility with reference transcripts; however, Salmon extends this framework through a dual-phase statistical inference procedure incorporating both online and offline expectation–maximization (EM) optimization, improving abundance estimation accuracy [4,7]. In addition, Salmon implements advanced sample-specific bias correction models that account for sequence-specific, fragment GC-content, and positional biases, reducing systematic technical artifacts that can influence transcript quantification [7]. Benchmarking studies have demonstrated that Salmon achieves higher sensitivity in downstream differential expression analyses while maintaining low false discovery rates and reducing spurious isoform switching events relative to alternative methods [4,7]. Furthermore, Salmon operates without generating large intermediate alignment files and efficiently scales across multiple CPU cores, enabling rapid quantification of large RNA-seq datasets with reduced computational demands. Based on its demonstrated accuracy, bias correction capabilities, statistical rigor, and scalability, Salmon was selected as the most suitable pseudoaligner for transcript quantification in this study [4–7].

For differential expression analysis (DEA), DESeq2 was selected over edgeR despite both being widely used and well-established Bioconductor packages for analyzing count-based sequencing data, including RNA-seq, SAGE-seq, ChIP-seq, and Hi-C datasets [8]. Both tools implement normalization strategies that have been shown to outperform alternative methods, particularly in scenarios where RNA composition varies substantially across biological conditions or when highly expressed genes are present [8]. However, benchmarking studies have demonstrated that as sample sizes increase, DESeq2 performs slightly better than many competing approaches in terms of statistical power and false discovery rate (FDR) control [9]. Specifically, for RNA-seq count data following either negative binomial or log-normal distributions, DESeq2 exhibits strong FDR control, improved power, and greater stability across varying sample sizes compared to alternative methods [9]. While other approaches such as EBSeq may be advantageous in very small sample size settings, DESeq2 is recommended for experiments with moderate or larger group sizes and remains one of the most robust and statistically reliable frameworks for RNA-seq differential expression analysis [9]. Based on its demonstrated performance, statistical rigor, and widespread validation across transcriptomic studies, DESeq2 was chosen as the most appropriate method for differential expression analysis in this study [8,9].

# Methods

# Results

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



