# Ribo-db Approach

The <b>Ribo-db Approach</b> leverages Ribo-Seq and RNA-Seq data to create non-redundant sample-specific protein databases (Ribo-db) containing only actively translated sequences. <b>Ribo-db Approach</b> identify start codons from Ribo-seq data of cells treated with the initiation inhibitor harringtonine or lactimidomycin through a TIS-calling probabilistic algorithm that differentiates start codon positions from background considering all near-cognate start codons. From Ribo-seq elongation data collected from actively translating cells and RNA-Seq data, <b>Ribo-db Approach</b> gets sample-specific transcriptomes using StringTie software (Pertea et al., 2015). Next, to generate sample specific transcription information, high quality single-nucleotide polymorphisms (SNPs) identified from RNA-seq data are integrated to the assembled transcripts. <b>Ribo-db Approach</b> intersect genomic positions of the start codons detected to the genomic positions of the assembled transcripts to generate the set of ORFs (coupled start codon with an assembled transcript) for in silico translation. 

## Usage

<b>Ribo-db Approach</b> is a step-by-step pipeline, as follows:
```python
1_Data_Preparation              # Sample RNA-Ribo-seq preparation
2_TIS_Calling                   # TIS detection
3_Transcriptome_Assembly        # Assamblage transcriptome
4_Proteins                      # DB generation
5_Identifications               # Portein source identification
```

## Author
Maria Virginia Ruiz,
<b>IRIC</b>,
Université de Montréal
