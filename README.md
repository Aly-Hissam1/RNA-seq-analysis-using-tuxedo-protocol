# RNA-seq-Analysis-Using-Tuxedo-Protocol

**Overview**

This project implements two protocols for RNA-seq analysis. The first utilizes the Tuxedo protocol, while the second employs state-of-the-art tools to generate a count matrix followed by differential gene expression analysis. The image provided illustrates the workflow of these protocols.

  ![image](https://github.com/Aly-Hissam1/RNA-seq-analysis-using-tuxedo-protocol/assets/117119881/e4da9851-3735-4599-a4b6-3e452c8dd330)

**Protocols Workflow**

*Quality Control*

Both protocols begin with FastP for quality control and trimming of the raw data.

*Alignment*

Tuxedo Protocol: Alignment is performed using TopHat.

State-of-the-art Tools Protocol: Alignment is achieved with Hisat2.

*Transcript Assembly*

Tuxedo Protocol: Transcript assembly is handled by Cufflinks, and the output is merged using Cuffmerge.

State-of-the-art Tools Protocol: HTSeq is used for assembly, and a custom Python script to_merge.py is utilized to merge the output into a quantification matrix.

**Differential Expression Analysis**

Tuxedo Protocol: Cuffdiff is included in the script.sh for differential expression.

State-of-the-art Tools Protocol: The quantification matrix is first TMM normalized. Subsequently, Limma Trend and five machine learning feature selection methods (Random Forest, Recursive Feature Elimination, SelectFromModel, Gradient Boosting, fdr) are applied.

**Assessment**

The intersection of differentially expressed genes by various tools were identified. 

Pathway enrichment analysis is conducted using databases such as Kaggle, Reactome, and Gene Ontology to ensure high confidence in the findings.

**Input Files**

The input files for the analysis are paired-end FastQ files, a reference genome and a GTF.

The input files are provided to the script as "Input.txt". Each line in the file corresponds to one sample and contains paths to the forward read, reverse read, the reference genome, the sample ID, and the protocol number to be used.

**Final Note**

This project is set up to streamline RNA-seq analysis using two different, yet robust protocols. By comparing the results from both, users can obtain a more comprehensive understanding of differential gene expression in their samples.

---
Feel free to contribute to the project or suggest improvements.
