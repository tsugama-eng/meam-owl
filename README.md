# MEAM-OWL

<p align="center">
  <img src="./MEAM-OWL.png" alt="MEAM-OWL logo" width="220">
</p>

## Overview

**MEAM-OWL (Motif Enrichment Analysis Modules on the Web, Lite)** is a
browser-based application for exploring DNA motif enrichment and associations
between motif presence or copy number and gene expression. Motifs are analyzed
as exact DNA words (k-mers). The application accepts prepared sequences or
extracts regions from a reference genome and GFF3 annotation, then performs
analysis locally in JavaScript.

The application offers two independent analysis modes:

| Sequence feature | Mode A: target vs. background enrichment analysis | Mode B: expression-based analysis |
|---|---|---|
| Motif presence or absence | Is the fraction of sequences containing a k-mer higher in the target group? | Are sequences containing a k-mer concentrated at the high or low end of a gene list ranked by the supplied expression values (e.g., TPM, FPKM, or log2FC)? |
| Motif copy number | Are per-sequence k-mer occurrence counts greater in the target group? | Are capped k-mer counts linearly associated with the supplied expression values? |

Mode A uses group membership, not expression values directly. Mode B pairs
sequences with numeric measurements by identifier. The application does not
calculate expression levels or differential expression from RNA-seq reads.

**Scope:** the current implementation is intended for exploratory analysis.
Motif enrichment or sequence compatibility alone does not establish regulatory function.

## Installation and Access

Download and extract the complete modified release package. On Windows, macOS, or Linux, keep `index.html`, `genome-stream.js`, `bh-adjust.js`, `motif_db.js`, and `pako.min.js` together, then open `index.html` in a modern desktop browser; no installation or build is needed.

## Input formats

### Sequence FASTA

```text
>seq0001
ACGTACGTGGCATGCA
>seq0002
TTGCAACGTGTCACGT
```

### Expression or other numeric table

```csv
GeneID,log2FC_condition_A,log2FC_condition_B
seq0001,1.8,-0.4
seq0002,-0.7,2.1
```

Provide plain comma-delimited (`.csv`) or tab-delimited (`.tsv`) tables. For this
example, set **ID Col Index = 1** and **Expression Col(s) = 2,3** or `2-3`.

### Genome and annotation

Use matching genome FASTA/GFF3 files with identical sequence names. GFF3 must have nine tab-separated columns and `key=value` attributes; transcript features need `ID` and `Parent`, and exon records must reference transcript IDs through `Parent`. CDS extraction also requires CDS records. Standard GTF attributes are unsupported.

## Workflow and parameters

### Mode A: Target vs. background

1. Supply **Target** and **Background** FASTA sequences.
2. Configure k-mer lengths, strand handling, and the background size limit.
3. For HOMER-Like or Mann-Whitney U, disable **Cluster** to inspect individual
   retained k-mers, or enable it to group retained candidates by a shared core.
4. Select **HOMER-Like**, **Mann-Whitney U**, or **Extended Core**, then inspect
   or export the results.

| Parameter | Default | Meaning |
|---|---|---|
| **k-mer Size** | `4-6` | Integer lengths 1–15 |
| **Forward Strand Only** | Off | When off, merges each k-mer with its reverse complement |
| **Max Background Size** | `1000` | Larger background sets are randomly subsampled without replacement |
| **Cluster** | On | Groups retained HOMER-Like or Mann-Whitney U candidates sharing a contiguous core |
| **Core length, m** | `4` | Length of the shared grouping core |

### Mode B: Expression-based analysis

1. Supply the sequence population FASTA and its corresponding measurement
   table (comma- or tab-delimited text).
2. Set **ID Col Index**, **Expression Col(s)**, k-mer lengths, strand handling, and whether to **Aggregate** results.
3. For Pearson correlation, set **Frequency Clipping Cap**.
4. Run **DRIMust Model** or **Pearson Correlation**. 

| Parameter | Default | Meaning |
|---|---|---|
| **k-mer Size** | `4-6` | Positive integer lengths up to 6; accepts a single length |
| **ID Col Index** | `1` | 1-based column index for sequence ID matching |
| **Expression Col(s)** | `2` | 1-based column indices analyzed independently, e.g., `2,3` or `2-5` |
| **Strand Specific** | Off | When off, merges reverse-complement counts |
| **Aggregate** | On | Groups retained candidates sharing a contiguous core |
| **Core length, m** | `4` | Length of the shared grouping core |
| **Frequency Clipping Cap** | `3` | Upper count bound applied for Pearson correlation only: `min(count, cap)` |

### Sequence extraction engine

1. Load matching reference genome FASTA and GFF3 files.
2. Select a region type. For **Custom**, choose an anchor and specify
   upstream/downstream distances; for **Other feature**, specify the GFF
   feature type.
3. Optionally upload an ID filter, select its ID column.
4. Optionally enable representative transcript selection (off by default).
5. Click **Extract Sequences**, review the preview, then inject the sequences into **Mode A Target**, **Mode A Background**, or **Mode B**.

## Outputs

### Result fields

| Field | Interpretation |
|---|---|
| **k-mer (Sequence/Core)** | An exact tested word, or a shared core when grouping is enabled |
| **Enrichment (×): HOMER-Like** | Target-positive fraction divided by the background rate, including its zero-count substitution |
| **Enrichment (×): U / Extended Core** | Mean target copy count divided by mean background copy count; if the background mean is zero, displays `∞x` |
| **Conc. factor (Opt. position)** | Motif-positive fraction in the selected ranked interval divided by the overall positive fraction; a minus sign denotes low-value-end enrichment |
| **Dynamic threshold detail** | High/low end, selected interval size, observed positive count, and expected positive count |
| **Coeff. (r)** | Pearson correlation between capped motif counts and the supplied measurements |
| **Clipping condition** | Count cap applied for Pearson correlation |
| **p-value / FDR** | Method-specific P value and BH-adjusted P value; retained results require FDR < 0.05 and the method-specific effect condition |
| **Similar Known Motif** | IUPAC consensus compatibility, not experimental evidence of TF binding or function |

### Downloads

- **CSV:** All stored result rows, including up to five known-motif matches per row. Matches are shown as `[database] name`; `…` indicates additional matches were omitted. The list follows database order, not confidence ranking.
- **FASTA:** Sequences produced by the extraction engine.
- **PNG:** Distribution plots for Mode B results. Counts are grouped into `0`, `1`, `2`, and `3+`.

## Methods

| Analysis                | Summary                                                      |
| ----------------------- | ------------------------------------------------------------ |
| **HOMER-Like**          | Tests whether motif presence is enriched in target sequences relative to background sequences. |
| **Mann-Whitney U**      | Compares per-sequence motif copy numbers between target and background groups. |
| **Extended Core**       | Extends enriched k-mers into longer sequences using U-test statistics. |
| **DRIMust Model**       | Searches for motif presence concentrated at either end of an expression-ranked list. |
| **Pearson Correlation** | Tests the correlation between capped motif counts and supplied measurements. |

HOMER-Like, Mann-Whitney U, and DRIMust Model use approximations in their significance calculations. The displayed methods are simplified implementations and are not direct calls to the original software packages.

Benjamini–Hochberg adjusted P values are calculated before effect-based filtering. Results must pass **FDR < 0.05** and the selected method's effect criterion to be displayed. This adjustment does not establish that the statistical assumptions are appropriate for every dataset.

## Reproducibility and validation

The repository's `test/` directory contains example inputs for trying the application. Reproduce a specific analysis using the same sequence and measurement files, reference versions, extraction settings, and analysis parameters.

## References
### Analysis Methods & Algorithms

*   **HOMER (Hypergeometric Enrichment)**
    *   Heinz S, Benner C, Spann N, et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell. 2010; 38(4): 576-589. doi:[10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)
*   **DRIMust (Rank-based Enrichment)**
    *   Leibovich L, Paz I, Yakhini Z, Mandel-Gutfreund Y. DRIMust: a web server for discovering rank imbalanced motifs using suffix trees. Nucleic Acids Res. 2013; 41(Web Server issue): W174-W179. doi:[10.1093/nar/gkt409](https://doi.org/10.1093/nar/gkt407)
*   **False Discovery Rate (Benjamini–Hochberg Adjustment)**
    *   Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B Stat Methodol*. 1995; 57(1): 289-300. doi:[10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

### Reference Transcription Factor & Cis-Regulatory Databases

*   **CIS-BP (The Catalog of Inferred Sequence Binding Preferences)**
    *   Weirauch MT, Yang A, Albu M, et al. Determination and inference of eukaryotic transcription factor sequence specificity. Cell. 2014; 158(6): 1431-1443. doi:[10.1016/j.cell.2014.08.009](https://doi.org/10.1016/j.cell.2014.08.009)
*   **JASPAR**
    *   Rauluseviciute I, Riudavets-Puig R, Blanc-Mathieu R, et al. JASPAR 2024: 20th anniversary of the open-access database of transcription factor binding profiles. Nucleic Acids Res. 2024; 52(D1): D174-D182. doi:[10.1093/nar/gkad1059](https://doi.org/10.1093/nar/gkad1059)
*   **PLACE (Plant Cis-acting Regulatory DNA Elements)**
    *   Higo K, Ugawa Y, Iwamoto M, Korenaga T. Plant cis-acting regulatory DNA elements (PLACE) database: 1999. Nucleic Acids Res. 1999; 27(1): 297-300. doi:[10.1093/nar/27.1.297](https://doi.org/10.1093/nar/27.1.297)

*Note: The statistical methods implemented inside MEAM-OWL are simplified, browser-optimized native JavaScript versions designed for fast exploratory analysis and do not make direct server-side calls to the original command-line tools or web servers.*


## License

MEAM-OWL is released under the [MIT License](LICENSE). Third-party components and database content remain subject to their respective licenses.
