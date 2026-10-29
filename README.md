# Mitochondrial DNA, the gut microbiome and type 2 diabetes in AWI-Gen 2

R and PLINK code for a thesis chapter asking whether mitochondrial DNA (mtDNA) variation, measured as
haplogroups (L0–L3) and as individual imputed mtDNA SNPs, is associated with gut microbiome
composition and with type 2 diabetes (T2D) in African cohorts.

## Questions and scripts

| Script | Question it answers | Method |
|---|---|---|
| **mtDNA quality control** |  |  |
| **Haplogroup assignment** | |  |
| **Haplogroup significance to overall microbiome** | Does haplogroup explain variation in overall microbiome composition? | PERMANOVA (`vegan::adonis2`) on Bray–Curtis distances, shuffling within site |
| **Abundance between haplogroups** | For each species or genus, does abundance differ between L1, L2 or L3 and L0? | MaAsLin2 (site as a fixed effect), with and without ancestry PCs, plus an Agincourt-only run |
| **Power Analysis** | Given that we found no significant contribution of haplogroup to species abundance/genus abundance, is the reason due to no existing associations or limited power? |      |
| **PCA** | Does including PCA as a covariate in my analysis change anything |  |
| **mtSNP and Species abundance** | Is any individual mtDNA SNP associated with the abundance of a species or genus? | Linear regression per tag SNP × taxon, pooled and per-site with meta-analysis, permutation check |
| **T2D and Haplogroup association** | Is T2D associated with the microbiome, and with haplogroup? | Composite T2D definition, MaAsLin2 microbiome vs T2D (plus treatment analyses), logistic regression of T2D on haplogroup |
| **mtSNP and T2D** | Is any individual mtDNA SNP associated with T2D? | Logistic regression per tag SNP, pooled and per-site with meta-analysis, permutation check |

Suggested running order: **mtDNA quality control** first (it produces the genotypes and PCs the others use),
then the remaining scripts in any order.

## What each script does

### mtDNA quality control
- Aligns the imputed mtDNA genotypes to the revised Cambridge Reference Sequence (rCRS), so that the
  reference allele is coded 0 and the other allele 1 (mtDNA is haploid, so genotypes are 0/1, and mixed
  calls are set to missing).


### Haplogroup significance to overall microbiome
- Keeps haplogroups L0–L3, converts counts to relative abundance and measures dissimilarity between
  every pair of people (Bray–Curtis).
- Tests haplogroup with PERMANOVA after age, BMI, antibiotic use and site, testing each term after all
  others (`by = "margin"`). Haplogroup labels are shuffled within site, so site differences cannot create a
  false result.
- Reports pseudo-F, R² and the permutation p-value for species and genus.

### Abundance between haplogroups
- Runs MaAsLin2 (total-sum scaling, log transformation, linear model) for each species and genus, comparing
  L1, L2 and L3 with L0 and adjusting for age, BMI, antibiotic use and site (site as a fixed effect).
- Recomputes Benjamini–Hochberg q-values over the haplogroup rows only, because MaAsLin2's own q-value is
  pooled across all model terms.
- Repeats the analysis with ancestry PCs added, and in Agincourt alone, to check that results do not change.
- Power analysis: reports confidence intervals for the top results, the smallest effect detectable with
  80% power, and the sample size needed to detect given fold-changes.

### mtSNP and Species abundance
- Regresses log2 relative abundance on each tag SNP, adjusting for age, BMI, antibiotic use, PC1–PC3
  (and site).
- **Run 1:** Agincourt only. **Run 2:** each site separately, then a fixed-effect inverse-variance
  meta-analysis, with I² (between-site disagreement) and k (number of sites contributing).
  **Run 3:** all sites pooled with site as a fixed effect (dummy variables).
- Reports BH q-values, a permutation check (outcome shuffled within site) and the genomic inflation factor
  (lambda) as a diagnostic.

### T2D and Haplogroup association
- Defines T2D as recorded diabetes status, fasting glucose > 7 mmol/L, random glucose > 11.1 mmol/L, or
  diabetes medication.
- Tests species abundance against T2D with MaAsLin2 (site fixed), refits it with treatment as a covariate,
  and compares treated and untreated cases separately with controls.
- Tests haplogroup against T2D with logistic regression (T2D ~ haplogroup + age + BMI + site).

### mtSNP and T2D
- Fits a logistic regression of T2D on each tag SNP, adjusting for age, BMI, PC1–PC3 and site.
- **Run A:** all sites pooled. **Run B:** each site separately, then combined by inverse-variance
  meta-analysis.
- Reports odds ratios with 95% confidence intervals, BH q-values and a permutation check.

## Conventions used throughout
- **Site** is a fixed effect, coded as dummy variables (never as one numeric code).
- **Multiple testing:** Benjamini–Hochberg q-values, reported alongside raw p-values. Permutation p-values
  (smallest-p across all tests) are used as a robustness check.
- **Effect sizes** are reported with 95% confidence intervals. Power is described through detectable effect
  sizes, not observed power.
- Reporting is at **species and genus** level only.

## Requirements
- R (4.x) with `vegan`, `permute`, `Maaslin2`, `lme4`
- PLINK 1.9 and PLINK2 (for the quality-control step)
- Random seeds are set inside the scripts so shuffle-based results are reproducible.

## Data
The input files are not distributed here:
- Imputed mtDNA genotypes (PLINK format) and nuclear array data
- Species and genus count tables (mOTUs)
- Clinical metadata (age, BMI, antibiotic use, site, glucose, diabetes status and treatment)
