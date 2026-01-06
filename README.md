# Custom Protein Domain Identification Pipeline
## Anchor-Based Analysis of Transcription Factor Architecture

A computational pipeline for systematic identification and characterization of protein domain architectures using anchor-based coordinate normalization. Designed for analyzing conserved and lineage-specific motifs in transcription factor families across multiple species.

---

## 📋 Overview

This pipeline implements an **anchor-based framework** for analyzing protein domain organization relative to a conserved reference domain. The workflow combines *de novo* motif discovery (MEME), genome-wide motif scanning (FIMO), statistical filtering, and phylogenetic context to map domain architectures across divergent protein sequences.

**Key Features:**
- Anchor-based coordinate normalization for cross-species comparison
- De novo short linear motif (SLiM) discovery
- Genome-wide motif scanning with statistical controls
- Automated low-complexity and false positive filtering
- Phylogenetic context-aware visualization

**Application Example:** Analysis of achaete-scute complex (ASC) transcription factors (bHLH family) across arthropod lineages

---

## 🔬 Background

### Transcription Factor Domain Organization

Many transcription factor families contain:
- **Conserved DNA-binding domain** (e.g., bHLH, zinc finger, homeobox)
- **Variable flanking regions** with regulatory and interaction domains
- **Lineage-specific motifs** that confer functional specialization

### Challenge

While DNA-binding domains are well-characterized, flanking regulatory domains remain poorly understood, particularly in non-model organisms. Identifying these domains is crucial for understanding functional divergence and evolutionary adaptation.

---

## 🛠️ Pipeline Architecture

### Workflow Overview

```
Input FASTA → bHLH Detection → Anchor Definition → Clade Partitioning
     ↓
De Novo Motif Discovery (MEME) → Motif Scanning (FIMO)
     ↓
FDR Correction → bHLH Masking → Low-Complexity Filtering
     ↓
Architecture Visualization → Sequence Logo Generation → Analysis
```

---

## 📊 Methodology

### 1. Anchor Definition (bHLH Domain)

#### 1.1 Domain Detection
```bash
hmmscan --domtblout domain_table.txt \
  Pfam-A.hmm \
  ASC_targets.fasta
```

**Key Parameters:**
- Model: PF00010 (bHLH domain)
- Top non-overlapping hit per protein defines anchor
- Record: domain start, end, midpoint

#### 1.2 Anchor-Centered Coordinates
```
start_rel = start - mid_bHLH
end_rel = end - mid_bHLH
```

This transformation enables cross-protein comparison of motif placement relative to the conserved bHLH region.

---

### 2. Clade Partitioning & Background Models

#### 2.1 Sequence Partitioning
Sequences partitioned by:
- **Family:** ASCa, ASCb, ASCc
- **Order:** Diptera, Lepidoptera, Coleoptera, Chelicerata
- **Subclade:** TrueSpiders A-F, ASH, ase

#### 2.2 Markov Background Computation
```bash
fasta-get-markov ASC_targets.fasta > background.meme
```

Clade-specific backgrounds mitigate amino acid composition bias during motif discovery.

---

### 3. De Novo Motif Discovery

#### 3.1 MEME Discovery
```bash
meme clade.fasta \
  -mod zoops \
  -nmotifs 10 \
  -minw 6 \
  -maxw 30 \
  -bfile background.meme \
  -oc meme_out/
```

**MEME Parameters:**
- **Mode:** ZOOPS (Zero or One Occurrence Per Sequence)
- **Motif width:** 6-30 amino acids (optimized for SLiMs)
- **Number of motifs:** Up to 10 per clade
- **Background:** Clade-specific Markov model

**Rationale:** ZOOPS mode is standard for short linear motif discovery in proteins; width range captures functional SLiMs while avoiding overly long, non-specific matches.

#### 3.2 Motif Library Assembly
```bash
cat */meme_out/meme.txt > combined_all_memes.txt
```

All MEME outputs concatenated for global scanning and cross-clade comparison.

---

### 4. Genome-Wide Motif Scanning

#### 4.1 FIMO Scanning
```bash
fimo --thresh 1e-5 \
  --max-stored-scores 20000000 \
  --bgfile clade/bg.meme \
  clade/meme_out/meme.txt \
  ASC_targets.fasta
```

**FIMO Parameters:**
- **Initial threshold:** p ≤ 1×10⁻⁵ (conservative for short motifs)
- **Background:** Clade-specific model
- **Output:** p-values and q-values for FDR control

---

### 5. Post-Processing & Quality Control

#### 5.1 FDR Correction
```python
# Benjamini-Hochberg FDR control
filtered_hits = hits[
    (hits['p_value'] <= 1e-4) & 
    (hits['q_value'] <= 0.05)
]
```

**Quality Filters:**
- **p-value:** ≤ 1×10⁻⁴
- **q-value:** ≤ 0.05 (Benjamini-Hochberg)
- **Best site per motif-sequence pair:** Retain only lowest p-value

#### 5.2 bHLH Core Masking
```python
# Remove sites overlapping PF00010 ± 8 residues
NON_BHLH = remove_overlap(hits, anchor_coords, buffer=8)
```

**Rationale:** Focus on flanking functional domains; ±8 aa buffer excludes conserved helix core while retaining adjacent flanks.

#### 5.3 Low-Complexity Filtering
```python
# NON_BHLH_SOFT: Remove generic low-complexity matches
filtered = NON_BHLH[
    ~((shannon_entropy < 1.6) & 
      (occurrence_in_families >= 2))
]
```

**Criteria for Removal:**
- Shannon entropy < 1.6 bits
- **AND** present in ≥2 ASC families

**Important:** Lineage-restricted low-entropy motifs retained (functional SLiMs often in compositionally biased regions).

---

### 6. Specificity Analysis

#### 6.1 Enrichment Calculation
```
enrichment = hit_rate_group / hit_rate_global
```

Where:
- hit_rate = (sequences with ≥1 site) / (total sequences in group)

#### 6.2 Family-Specific Candidate Shortlisting
**Criteria:**
- Present in ≥3 sequences within one ASC family
- Absent from other families
- Enrichment index > 1.5

---

### 7. Architecture Visualization

#### 7.1 Domain Mapping
Architecture figures generated directly from FIMO hits without clustering.

**Coordinate Transformation:**
```python
rel_start = start - anchor_pos
rel_end = stop - anchor_pos
```

**Visualization Elements:**
- **Gray boxes:** bHLH domain (position 0 = anchor midpoint)
- **Colored rectangles:** Individual PWM motifs (widths = true amino acid length)
- **Y-axis:** Protein sequences (color-coded by family)
- **X-axis:** Position relative to bHLH anchor (-300 to +300 aa)

#### 7.2 Legend Construction
**Representative tokens:** First 6 amino acids from matched sequence (or MEME consensus)

**Per-PWM Statistics:**
- Median position relative to bHLH
- Median motif width
- Sequence coverage (one-vote-per-sequence rule)

---

### 8. Sequence Logo Generation

#### 8.1 bHLH Core Extraction
```bash
# Extract PF00010 coordinates for proteins with motif hits
# Partition by lineage: Chelicerata vs Insecta
```

#### 8.2 Multiple Sequence Alignment
```bash
mafft --auto extracted_bHLH.fasta > aligned_bHLH.fasta
```

#### 8.3 Logo Rendering
```python
import logomaker

# Generate sequence logos for:
# - All taxa combined
# - Insects only
# - Spiders (Chelicerata) only
```

---

## 📈 Pipeline Outputs

The pipeline generates several types of outputs:

### 1. Domain Architecture Visualizations
- Anchor-normalized position plots showing motif distribution
- Family/clade-specific domain organization maps
- Per-sequence architecture diagrams with motif annotations

### 2. Sequence Conservation Analysis
- Sequence logos for conserved domains
- Lineage-specific conservation patterns
- Position-weight matrices (PWMs) for identified motifs

### 3. Statistical Summaries
- Motif presence/absence matrices
- Family-specific enrichment scores
- Positional constraint analysis
- Coverage and occurrence statistics

### 4. Quality Metrics
- FDR-controlled significance values
- Shannon entropy scores for complexity filtering
- Per-motif confidence scores

---

## 💻 Technical Requirements

### Software Dependencies

```bash
# Core bioinformatics tools
hmmer >= 3.4
meme-suite >= 5.0
mafft >= 7.0

# Python environment
python >= 3.11
pandas >= 1.3.0
numpy >= 1.21.0
biopython >= 1.79
logomaker >= 0.8
matplotlib >= 3.5.0
seaborn >= 0.11.0
```

### Python Package Installation
```bash
pip install pandas numpy biopython logomaker matplotlib seaborn
```

### Required Reference Files
```
reference_data/
├── Pfam-A.hmm              # Pfam HMM database
├── PF00010.hmm             # bHLH domain model
└── ASC_targets.fasta       # Curated ASC protein sequences
```

---

## 📂 Project Structure

```
project/
├── data/
│   ├── ASC_targets.fasta           # Input protein sequences
│   └── metadata/
│       └── clade_assignments.csv   # Phylogenetic classifications
├── hmmer_scans/
│   └── bHLH_domains.domtblout     # PF00010 domain coordinates
├── clades/
│   ├── ASCa/
│   │   ├── bg.meme                 # Background model
│   │   └── meme_out/
│   │       └── meme.txt            # Discovered motifs
│   ├── ASCb/
│   └── ASCc/
├── fimo_results/
│   ├── ALL.fimo_hits.tsv          # Raw FIMO hits
│   ├── NON_BHLH.fimo_hits.tsv     # bHLH-masked hits
│   └── NON_BHLH_SOFT.fimo_hits.tsv # Final filtered set
├── analysis/
│   ├── motif_presence_matrix.tsv
│   ├── enrichment_summary.tsv
│   └── per_sequence_hits.tsv
├── visualizations/
│   ├── domain_architecture_all.png
│   ├── ASCa_architecture.png
│   ├── ASCb_architecture.png
│   ├── ASCc_architecture.png
│   └── bHLH_logos/
│       ├── all_taxa.png
│       ├── insects.png
│       └── spiders.png
└── scripts/
    ├── 01_hmmer_scan.sh
    ├── 02_clade_partition.py
    ├── 03_meme_discovery.sh
    ├── 04_fimo_scan.sh
    ├── 05_post_processing.py
    ├── 06_architecture_plots.py
    └── 07_sequence_logos.py
```

---

## 🎯 Usage Example

### Complete Pipeline Execution

#### Step 1: bHLH Domain Detection
```bash
# Scan all proteins for PF00010 domain
hmmscan --domtblout hmmer_scans/bHLH_domains.domtblout \
  reference_data/Pfam-A.hmm \
  data/ASC_targets.fasta

# Extract anchor coordinates
python scripts/extract_anchors.py \
  --domtblout hmmer_scans/bHLH_domains.domtblout \
  --output analysis/anchor_coordinates.csv
```

#### Step 2: Clade Partitioning
```bash
# Partition sequences by clade/family
python scripts/02_clade_partition.py \
  --fasta data/ASC_targets.fasta \
  --metadata data/metadata/clade_assignments.csv \
  --outdir clades/

# Generate background models
for clade in clades/*/; do
  fasta-get-markov ${clade}sequences.fasta > ${clade}bg.meme
done
```

#### Step 3: De Novo Motif Discovery
```bash
# Run MEME for each clade
for clade in clades/*/; do
  meme ${clade}sequences.fasta \
    -mod zoops \
    -nmotifs 10 \
    -minw 6 \
    -maxw 30 \
    -bfile ${clade}bg.meme \
    -oc ${clade}meme_out/
done

# Combine all motifs
cat clades/*/meme_out/meme.txt > combined_all_memes.txt
```

#### Step 4: Motif Scanning
```bash
# Scan with FIMO
fimo --thresh 1e-5 \
  --max-stored-scores 20000000 \
  --bgfile clades/ASCa/bg.meme \
  combined_all_memes.txt \
  data/ASC_targets.fasta \
  --oc fimo_results/

# Repeat for each clade-specific background
```

#### Step 5: Post-Processing
```bash
python scripts/05_post_processing.py \
  --fimo fimo_results/fimo.tsv \
  --anchors analysis/anchor_coordinates.csv \
  --p-threshold 1e-4 \
  --q-threshold 0.05 \
  --bhlh-buffer 8 \
  --entropy-threshold 1.6 \
  --output-dir analysis/
```

#### Step 6: Visualization
```bash
# Generate architecture plots
python scripts/06_architecture_plots.py \
  --hits analysis/NON_BHLH_SOFT.fimo_hits.tsv \
  --anchors analysis/anchor_coordinates.csv \
  --metadata data/metadata/clade_assignments.csv \
  --output-dir visualizations/

# Generate sequence logos
python scripts/07_sequence_logos.py \
  --fasta data/ASC_targets.fasta \
  --anchors analysis/anchor_coordinates.csv \
  --output-dir visualizations/bHLH_logos/
```

---

## 🔧 Parameter Optimization Guide

### MEME Parameters

| Parameter | Default | Range | Rationale |
|-----------|---------|-------|-----------|
| `-minw` | 6 | 4-10 | Minimum for functional SLiMs |
| `-maxw` | 30 | 20-50 | Balance specificity vs sensitivity |
| `-nmotifs` | 10 | 5-15 | Capture diversity without over-fitting |
| `-mod` | zoops | - | Standard for SLiM discovery |

### FIMO Parameters

| Parameter | Default | Range | Rationale |
|-----------|---------|-------|-----------|
| `--thresh` | 1e-5 | 1e-4 to 1e-6 | Conservative for short motifs |
| Post-filter p | 1e-4 | 1e-3 to 1e-5 | FDR control |
| Post-filter q | 0.05 | 0.01 to 0.1 | Benjamini-Hochberg standard |

### Filtering Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| bHLH buffer | ±8 aa | Exclude helix core, retain flanks |
| Shannon entropy | <1.6 bits | Empirical boundary for low complexity |
| Multi-family occurrence | ≥2 families | Remove generic matches |

---

## 📊 Output File Formats

### 1. Anchor Coordinates File
```csv
sequence_id,anchor_start,anchor_end,anchor_mid,domain_width
Ptep_g1234.t1,145,195,170,50
Dmel_achaete,230,280,255,50
```

### 2. FIMO Hits (Filtered)
```tsv
motif_id	sequence_id	start	stop	strand	score	p_value	q_value	matched_sequence	anchor_rel_start	anchor_rel_end
KKDPFS	Ptep_g1234.t1	120	142	+	15.2	1.2e-5	0.003	KKDPFSLMQRTVWH	-50	-28
```

### 3. Motif Presence Matrix
```tsv
motif_id	ASCa	ASCb	ASCc	family_specific
KKDPFS	0	15	0	ASCb
EASSPY	12	0	0	ASCa
VNKENE	0	0	10	ASCc
```

### 4. Enrichment Summary
```tsv
motif_id	family	hit_count	hit_rate	enrichment	median_position
KKDPFS	ASCb	15	1.00	15.2	-55
EASSPY	ASCa	12	0.75	8.4	+45
```

---

## 🔬 Applications

This pipeline can be applied to:

1. **Transcription Factor Families**
   - bHLH, zinc finger, homeobox, and other TF families
   - Cross-species comparative analysis
   - Identification of family-specific regulatory domains

2. **Evolutionary Studies**
   - Motif gain/loss events across phylogeny
   - Lineage-specific domain elaboration
   - Ancestral domain reconstruction

3. **Functional Annotation**
   - Prediction of interaction surfaces
   - Identification of regulatory elements
   - Characterization of disordered regions

4. **Structural Biology**
   - Mapping functionally important regions
   - Positional constraint analysis
   - Modular domain organization

---

## 📚 Key References

### Pipeline Development
This pipeline was developed as part of a research project at Turetzek Lab (Evolutionary Developmental Biology), Ludwig-Maximilians-Universität München.

### Methodological References

**Motif Discovery:**
- Bailey et al. (2009, 2015) - MEME Suite tools  
  *Nucleic Acids Research*
- Grant et al. (2011) - FIMO: Scanning for occurrences of a given motif  
  *Bioinformatics*

**Statistical Methods:**
- Benjamini & Hochberg (1995) - Controlling the false discovery rate  
  *Journal of the Royal Statistical Society*

**Domain Databases:**
- Mistry et al. (2021) - Pfam: The protein families database  
  *Nucleic Acids Research*

**Short Linear Motifs:**
- Tompa et al. (2014) - A million peptide motifs for the molecular biologist  
  *Molecular Cell*
- Van Roey et al. (2014) - Short linear motifs: Ubiquitous and functionally diverse  
  *Chemical Reviews*

**Sequence Analysis:**
- Katoh & Standley (2013) - MAFFT multiple sequence alignment  
  *Molecular Biology and Evolution*
- Tareen & Kinney (2020) - Logomaker: Beautiful sequence logos in Python  
  *Bioinformatics*

---

## 🤝 Acknowledgments

**Research conducted at:**
- Turetzek Lab - Evolutionary Developmental Biology
- Ludwig-Maximilians-Universität München, Faculty of Biology

---

## 📧 Contact

**Author:** Merve Görkem Durmaz  
**Email:** m.gorkemdurmaz@gmail.com  
**LinkedIn:** [linkedin.com/in/merve-görkem-durmaz-50902616a](https://www.linkedin.com/in/merve-goerkem-durmaz-50902616a)  
**GitHub:** [github.com/gorkem8d](https://github.com/gorkem8d)  
**ORCID:** [0000-0003-2106-2860](https://orcid.org/0000-0003-2106-2860)

---

## 📝 License

This pipeline and methodology are available for academic and research purposes. For commercial use or collaboration inquiries, please contact the author.

---

## ⚠️ Important Notes

1. **Input Requirements:** Curated protein sequences with proper phylogenetic annotations
2. **Computational Resources:** HPC recommended for large-scale motif scanning
3. **Random Seeds:** Fixed at 42 for reproducibility
4. **Version Control:** All software versions should be pinned for reproducibility
5. **Clade-Specific Backgrounds:** Essential for reducing composition-driven artifacts
6. **Motif Interpretation:** Always consider phylogenetic context
7. **Research Status:** This pipeline is part of active research; specific biological findings are not yet published

---

## 🔮 Future Enhancements

Potential improvements to the pipeline:

1. **Integration with Additional Databases:**
   - ELM (Eukaryotic Linear Motif) resource
   - PROSITE patterns
   - InterPro domain annotations

2. **Machine Learning Integration:**
   - Deep learning for motif prediction
   - Supervised classification of domain types
   - Transfer learning from related protein families

3. **Interactive Visualization:**
   - Web-based interface for results exploration
   - Dynamic filtering and annotation
   - Exportable publication-quality figures

4. **Expanded Analysis Features:**
   - 3D structure prediction integration (AlphaFold)
   - Co-evolution analysis
   - Functional site prediction
   - Disorder region prediction

5. **Scalability Improvements:**
   - Parallelization for large datasets
   - Cloud deployment options
   - Optimized memory usage

---

**Last Updated:** January 2025  
**Pipeline Version:** 1.0  
**HMMER Version:** 3.4+  
**MEME Suite Version:** 5.x
