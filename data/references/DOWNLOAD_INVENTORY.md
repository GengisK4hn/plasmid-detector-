# Reference Data Download Inventory

**Date:** 2026-04-17
**Purpose:** SV40 and promoter element detection tool development
**Location:** `/home/dad/supracode-tool/enhanced-plasmid-finder/data/references/`

---

## ✅ Successfully Downloaded

### 1. Vaccine Plasmid Reference Sequences

#### Pfizer BNT162b2 Production Plasmid
- **Accession:** OR134577.1
- **Source:** NCBI GenBank
- **URL:** https://www.ncbi.nlm.nih.gov/nuccore/OR134577.1
- **Files:**
  - `pfizer_bnt162b2_OR134577.1.fasta` (14.8 KB)
  - `pfizer_bnt162b2_OR134577.1.gb` (14.8 KB GenBank format)
- **Length:** 7,810 bp
- **Contains:** SV40 enhancer, TATA boxes, GC boxes, Sp1 sites, ColE1 origin

#### Moderna mRNA-1273 Production Plasmid
- **Accession:** OR134578.1
- **Source:** NCBI GenBank
- **URL:** https://www.ncbi.nlm.nih.gov/nuccore/OR134578.1
- **Files:**
  - `moderna_mrna_1273_OR134578.1.fasta` (13.0 KB)
  - `moderna_mrna_1273_OR134578.1.gb` (13.0 KB GenBank format)
- **Length:** ~7,000 bp
- **Contains:** Vaccine production plasmid sequence

### 2. Research Papers and Data

#### McKernan 2023 - Bivalent Sequencing
- **Title:** Sequencing of bivalent COVID-19 vaccines
- **Date:** April 11, 2023
- **URL:** https://osf.io/preprints/osf/b9t7m_v1
- **Direct PDF:** https://cdn-ceo-ca.s3.amazonaws.com/1i4tp3q-Sequencing%20of%20bivalent_4-11-23.pdf
- **File:** `McKernan_2023_Bivalent_Sequencing.pdf` (2.8 MB)
- **Contains:** BAM files, fastq, Megahit assemblies for Pfizer & Moderna vials

#### Speicher et al. 2025 - SV40 Quantification
- **Title:** Quantification of residual plasmid DNA and SV40 promoter-enhancer sequences in Pfizer BioNTech and Moderna modRNA COVID-19 vaccines from Ontario, Canada
- **Journal:** Journal of Immunotoxicology
- **Date:** 2025
- **DOI:** 10.1080/08916934.2025.2551517
- **URL:** https://www.tandfonline.com/doi/full/10.1080/08916934.2025.2551517
- **Files:**
  - `Speicher_2025_SV40_quantification.pdf` (9.0 KB)
  - `Speicher_2025_Supplemental_Data.zip` (0 bytes - failed download)
- **Contains:** qPCR, fluorometry, ONT sequencing data, SV40 quantitation

---

## ⚠️ Failed/Requires Manual Download

### TGA FOI 25-0070
- **Title:** Australian TGA Freedom of Information response on DNA testing
- **URL:** https://www.tga.gov.au/sites/default/files/2024-12/FOI%2025-0070.pdf
- **Status:** HTTP/2 stream error - requires manual download
- **Action Required:** Manual download from TGA website

### Speicher 2025 Supplemental Data
- **URL:** Figshare dataset (direct link failed)
- **Alternative:** Access from paper page: https://tandf.figshare.com/articles/dataset/30068064
- **Status:** Requires manual download from Figshare

---

## 📋 Still To Download

### JASPAR Transcription Factor Database
- **Purpose:** Automated TATA/Sp1/GC-box scanning
- **URL:** https://jaspar.elixir.no/download
- **Recommended:** Non-redundant CORE vertebrates
- **Status:** Not yet downloaded

### McKernan 2023 Raw Data Files
- **Location:** OSF page with Mega.nz links
- **URL:** https://osf.io/preprints/osf/b9t7m_v1
- **Contains:** BAM files, fastq, assemblies
- **Status:** Requires Mega.nz download

### SV40 Promoter Reference
- **Source:** SnapGene
- **URL:** https://www.snapgene.com/plasmids/basic_cloning_vectors/SV40_promoter
- **Format:** .dna file
- **Status:** Not yet downloaded

---

## 🔍 Verification

### Sequences Downloaded
- ✅ Pfizer OR134577.1 (complete)
- ✅ Moderna OR134578.1 (complete)

### Papers Downloaded
- ✅ McKernan 2023 (complete)
- ✅ Speicher 2025 paper (complete)
- ⚠️ Speicher 2025 supplemental data (failed)
- ⚠️ TGA FOI 25-0070 (failed)

---

## 📊 Usage Notes

### For SV40 Detection
Use the downloaded FASTA files as reference sequences:
```bash
python verify_sv40_findings.py data/references/pfizer_bnt162b2_OR134577.1.fasta
python verify_sv40_findings.py data/references/moderna_mrna_1273_OR134578.1.fasta
```

### For Plasmid Detection
Use the GenBank files for complete feature annotations:
- Contains promoter locations
- Contains origin of replication positions
- Contains enhancer element coordinates

### For Cross-Reference
Use the research papers to validate findings:
- McKernan 2023: Independent sequencing methodology
- Speicher 2025: qPCR quantification methods
- TGA FOI: Regulatory perspective on DNA testing

---

## 🚀 Next Steps

1. Test SV40 detection on newly downloaded references
2. Create comprehensive promoter map from GenBank annotations
3. Integrate reference sequences into automated scanning pipeline
4. Update documentation with verified detection results

---

**Last Updated:** 2026-04-17
**Status:** Ready for analysis
**Total Files:** 8 files (5 successful, 2 failed, 1 pending)
