# Fever Cohort Analysis - Cross-Modal Autoencoder for Disease State Discovery

**Research Hypothesis**: Cross-modal embedding via autoencoders can capture "true" disease states. For patients with community-acquired fever in ER/hospital settings, different etiologies (e.g., viral vs bacterial pneumonia) will naturally separate in the learned embedding space because they generate systematically different patterns across modalities (medications, labs, vitals, diagnoses).

**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Analysis Type**: Unsupervised learning, cross-modal autoencoder, disease state discovery
**Primary Methods**: OMOP data extraction, temporal matching, autoencoder embedding, cluster analysis

---

## Project Overview

### Core Hypothesis
Different fever etiologies produce distinct multi-modal signatures that can be learned by an autoencoder without explicit labels:
- **Bacterial pneumonia**: ↑WBC, ↑CRP/procalcitonin, antibiotics, specific vital patterns
- **Viral infection**: Different lab profiles, antipyretics primarily, distinct clinical trajectory
- **Other causes**: Unique combinations across modalities

The autoencoder will learn to compress these high-dimensional multi-modal data into a lower-dimensional embedding where similar disease states cluster together.

### Cohort Definition
**Inclusion criteria:**
- Elevated temperature (>38°C) within first 24 hours of ER or hospital visit 
- During ER or inpatient hospital visit
- Complete multi-modal data available (meds, labs, vitals, diagnoses)

**No control group needed** - This is not a case-control study. All participants have fever; the goal is to discover natural clusters within the fever cohort.

### Multi-Modal Feature Extraction
For each fever episode, extract:
1. **Medications** - All
2. **Lab tests** - Commonly obtained labs (e.g., CBC/diff, CRP, procalcitonin, blood cultures, other infection markers, CMP)
3. **Vital signs** - Temperature, heart rate, blood pressure (SBP, DBP, MAP, PP), respiratory rate, etc. 
4. **Diagnoses** - Pre-existing conditions (ICD codes for pneumonia, leukemia, immunosuppression, UTI, etc.)

### Ground Truth Labels (for validation)
Extract actual clinical diagnoses AFTER admission and BEFORE discharge to validate if autoencoder clusters align with:
- Bacterial pneumonia
- Viral pneumonia
- Sepsis
- Urinary tract infection
- Other fever etiologies

---

## Project Workflow

```
Stage 1: Cohort Definition
  ↓ Identify all community-acquired fever (only fever in first 24 hours of admission or in ER) cases (elevated vitals ONLY)
  ↓ Link to ER/hospital visits
  ↓ Extract ground truth diagnoses (for validation)
  ↓ Filter for data completeness

Stage 2: Multi-Modal Feature Extraction
  ↓ Extract medications (temporal patterns around fever)
  ↓ Extract labs (infection markers, temporal trends)
  ↓ Extract vitals (temperature trajectory + other vitals)
  ↓ Extract pre-existing diagnoses (no diagnoses after admission due to data leakage)
  ↓ Temporal alignment to presentation index date (i.e. ER or inpatient)

Stage 3: Feature Engineering & Representation
  ↓ Encode medications (binary, frequency, temporal patterns)
  ↓ Encode labs (continuous, normalization, missingness)
  ↓ Encode vitals (time series features, trajectories)
  ↓ Encode diagnoses (binary indicators, counts)
  ↓ Create unified multi-modal feature matrix

Stage 4: Autoencoder Training
  ↓ Design cross-modal autoencoder architecture
  ↓ Handle missing data (modality-aware)
  ↓ Train encoder-decoder on multi-modal features
  ↓ Learn low-dimensional embedding

Stage 5: Cluster Discovery & Validation
  ↓ Visualize embedding space (UMAP, t-SNE)
  ↓ Identify natural clusters
  ↓ Characterize clusters by features
  ↓ Validate against ground truth diagnoses
  ↓ Assess clinical interpretability

Stage 6: Analysis & Interpretation
  ↓ Compare clusters to known etiologies
  ↓ Identify novel disease state patterns
  ↓ Feature importance analysis
  ↓ Clinical validation
```

---

## File Map

### Stage 1: Cohort Definition
- **01_fever_cohort.py/.ipynb** - Identify all fever cases (vitals only)
- **01b_ground_truth_labels.py/.ipynb** - Extract validation diagnoses (pneumonia, sepsis, etc.) after presentation and before discharge/death

### Stage 2: Multi-Modal Feature Extraction
- **02_extract_medications.py/.ipynb** - Medication data for first days of presentation
- **03_extract_labs.py/.ipynb** - Lab tests for first days of presentation
- **04_extract_vitals.py/.ipynb** - Temperature and other vitals for first days of presentation
- **05_extract_diagnoses.py/.ipynb** - Pre-existing conditions only

### Stage 3: Feature Engineering
- **06_feature_engineering.py/.ipynb** - Create unified multi-modal feature matrix
- **07_descriptive_stats.py/.ipynb** - Data completeness, distributions

### Stage 4: Autoencoder Training
- **08_autoencoder_model.py/.ipynb** - Model architecture and training
- **09_embedding_extraction.py/.ipynb** - Generate embeddings for all patients

### Stage 5: Cluster Discovery
- **10_clustering_analysis.py/.ipynb** - Identify clusters, UMAP visualization
- **11_cluster_validation.py/.ipynb** - Compare to ground truth diagnoses

### Stage 6: Interpretation
- **12_cluster_characterization.py/.ipynb** - Feature profiles of each cluster
- **13_clinical_interpretation.py/.ipynb** - Clinical validation and insights

---

## Key Data Sources

### OMOP Tables Used
- `condition_occurrence` and `observation` - (ICD diagnostic codes)
- `measurement` - Temperature vitals, lab values
- `visit_occurrence` - ER/hospital encounters
- `drug_exposure` - Medications
- `person` - Demographics
- `concept` - OMOP vocabulary lookups

### Reference Files
Located in `_reference/`:
- `All of Us Controlled Tier Dataset v8 CDR Data Dictionary (C2024Q3R8) - OMOP-Compatible Tables.tsv` - Data Dictionary
- `concept_domains.tsv` - concept domains
- `condition_vocabularies.tsv` - ICD-9, ICD-10, and SNOMED counts
- `measurement_concepts.tsv` - Common measurement concept IDs
- `table_row_counts.tsv` - row counts for each table
- `table_schemas.tsv` - OMOP table structures
- `visit_concepts.tsv` - Visit type concepts
- `vocabulary_structure.tsv` - Available vocabularies

---

## Fever Cohort Definition

### Inclusion Criteria

**Primary inclusion:**

#### Vital Sign-Based
**Temperature measurements:**
- Concept IDs: [TBD - extract from measurement_concepts.csv]
- Threshold: >38.0°C (100.4°F)
- Context: Must be within first 24 hours of presentation (ER or inpatient visit)

#### Combined Inclusion
Include if patient has:
- Elevated temperature (>38°C)
- **AND** encounter is ER or inpatient visit
- **AND** sufficient multi-modal data available (defined below)

### Data Completeness Requirements

For autoencoder training, require minimum data across modalities:
- **Must have**: At least 2 of 4 modalities with complete data
  - Medications (any prescription during episode)
  - Labs (at least 1 lab)
  - Vitals (temperature measurement)

Rationale: Autoencoder needs sufficient signal across modalities to learn meaningful representations

---

## Ground Truth Labels (for Validation)

Extract concurrent diagnoses for validating autoencoder clusters. Use all phecodes with appropriate specificity.

**Phecodes:**
- RE_468.1 - Viral pneumonia
- RE_468.2 - Bacterial pneumonia
- GE_981.8 - Malignant hyperthermia

These labels are **NOT used for training** the autoencoder - only for validation to see if discovered clusters align with clinical diagnoses.

---

## Temporal Windows

Define time windows relative to fever index date:

| Data Type | Prior Window | Post Window | Rationale |
|-----------|--------------|-------------|-----------|
| **Temperature vitals** | visit onset | +7 days | Capture acute fever episode |
| **Antipyretics** | visit onset | +7 days | Immediate symptom management |
| **Antibiotics** | visit onset | +14 days | Treatment for bacterial infection |
| **Labs (WBC, CRP)** | visit onset | +14 days | Infection workup and monitoring |
| **Pre-existing conditions** | All prior | - | Historical diagnoses before fever |

---

## Autoencoder Architecture Considerations

### Cross-Modal Representation Learning

**Challenge**: Different modalities have different:
- Data types (binary, continuous, categorical)
- Missingness patterns
- Scales and distributions
- Temporal characteristics

**Approach options:**

#### Option 1: Concatenated Input (Simple)
```
Input: [medications | labs | vitals | diagnoses] → Encoder → Embedding → Decoder → Output
```
- Pros: Simple architecture, single training loop
- Cons: Assumes equal importance of modalities, struggles with missing modalities

#### Option 2: Modality-Specific Encoders (Recommended)
```
Medications → Med_Encoder ─┐
Labs → Lab_Encoder ─────────┤
Vitals → Vital_Encoder ─────┼→ Fusion → Shared_Embedding → Modality_Decoders
Diagnoses → Dx_Encoder ─────┘
```
- Pros: Learns modality-specific representations, handles missing data better
- Cons: More complex architecture, requires careful fusion strategy

#### Option 3: Variational Autoencoder (VAE)
- Learn probabilistic embeddings (mean, variance)
- Better for generation and uncertainty quantification
- More complex training (KL divergence term)

### Handling Missing Data

**Reality**: Not all patients have all modalities

**Strategies:**
1. **Imputation**: Fill with modality-specific defaults (risky)
2. **Masking**: Use attention masks to ignore missing modalities (recommended)
3. **Indicator variables**: Add "missingness" as features
4. **Multi-task learning**: Train to predict missing modalities

### Feature Engineering by Modality

#### Medications
- **Binary indicators**: Presence/absence of drug class (antipyretics, antibiotics, etc.)
- **Counts**: Number of unique medications
- **Temporal**: Days from fever to first antibiotic
- **Route**: Oral vs IV (severity indicator)

#### Labs
- **Continuous normalized**: Z-scores relative to normal ranges
- **Categorical bins**: Normal/elevated/critically elevated
- **Temporal trends**: Slope of WBC over time
- **Missingness**: Binary indicator for each lab type

#### Vitals
- **Temperature trajectory**: Max, min, mean, variance, duration >38°C
- **Other vitals**: Heart rate, BP, respiratory rate at fever peak
- **Time series features**: Area under curve, rate of change
- **Normalization**: Age-adjusted Z-scores

### Embedding Dimension

**Trade-offs:**
- **Too small (2-8 dims)**: Good for visualization, may lose information
- **Too large (>128 dims)**: Captures detail, harder to interpret
- **Sweet spot (16-64 dims)**: Balance between information and interpretability

**Recommendation**: Start with 32 dimensions, evaluate reconstruction loss vs interpretability

### Loss Function Design

**Multi-modal reconstruction loss:**
```
Total_Loss = α₁ * Med_Loss + α₂ * Lab_Loss + α₃ * Vital_Loss
```

Where:
- Med_Loss: Binary cross-entropy (binary indicators)
- Lab_Loss: MSE (continuous values)
- Vital_Loss: MSE (continuous time series)
- α weights: Tune to balance modality importance

**Regularization:**
- L2 penalty on weights
- Dropout in encoder/decoder
- VAE: KL divergence term

---

## Expected Multi-Modal Data Matrix

### Dimensions
- **Rows**: N fever episodes (~50K-100K estimated)
- **Columns**: M features across modalities (~500-2000 estimated)

### Feature Breakdown Estimate
| Modality | Feature Types | Estimated Features |
|----------|---------------|-------------------|
| Medications | Drug classes, routes, timing | 50-100 |
| Labs | WBC, CRP, procalcitonin, cultures, chemistry | 20-50 |
| Vitals | Temperature trajectory, HR, BP, RR | 20-40 |
| Diagnoses | Pre-existing ICD codes | 200-500 |
| Metadata | Age, sex, visit type, time features | 10-20 |
| **Total** | | **~300-700 features** |

### Sparsity Considerations
- Medications: ~30% patients have antibiotics
- Labs: ~50% have infection workup
- Vitals: ~80% have temperature data (by definition)
- Diagnoses: ~90% have concurrent diagnoses

**Implication**: Expect 30-70% missing data per feature, requiring robust missing data handling

---

## Reusable Components

### From Sarcoid Project
```python
from sarcoid.cohort import AouQueries

# Flexible OMOP concept searching
aou = AouQueries()
aou.find_diagnosis_codes(
    search_terms=["fever", "pyrexia"],
    pattern_codes={'ICD10CM': ['R50.%']},
    exact_codes={'ICD9CM': ['780.60', '780.61']}
)
```

### From Medications Project
```python
# RxNorm ingredient extraction pattern
# Use cb_criteria_ancestor for drug hierarchy
```

---

## Implementation Notes

### Cohort Size Estimate
- All of Us CDR v8: ~400K total people with genomic data
- Estimated fever cases: 50K-100K (10-25% prevalence, rough estimate)
- With sufficient multi-modal data: ~25K-50K (50% expected after completeness filters)

### Data Volume Estimates
- Condition_occurrence: ?
- Drug_exposure: ?
- Measurement: ?
- Expected fever cohort data: ?

### Unit Conversion Needs
**Temperature:**
- Celsius ↔ Fahrenheit conversion
- Standard: Report in Celsius

**Lab values:**
- WBC: cells/μL (standard)
- CRP: mg/L or mg/dL (need conversion)

---

## Resource Usage & Machine Types

### Observed Usage

**Cohort Size**: [To be filled in during execution]

#### Phase 1: Data Extraction & Initial Processing

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| BigQuery extraction (fever diagnoses) | | | | | | | |
| BigQuery extraction (temperature vitals) | | | | | | | |
| BigQuery extraction (medications) | | | | | | | |
| BigQuery extraction (labs) | | | | | | | |
| Link to ER/hospital visits | | | | | | | |

**Notes:**
- [Add observations about what worked well or caused issues]
- [Note any out-of-memory errors or performance bottlenecks]

#### Phase 2: Data Cleaning & Transformation

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| Temperature unit conversion | | | | | | | |
| Fever threshold filtering | | | | | | | |
| Temporal matching (medications) | | | | | | | |
| Temporal matching (labs) | | | | | | | |
| Temporal matching (vitals) | | | | | | | |

**Notes:**
- [Document memory spikes during joins]
- [Note if operations were CPU or memory bound]

#### Phase 3: Analysis & Visualization

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| Frequency tables (medications) | | | | | | | |
| Frequency tables (labs) | | | | | | | |
| Frequency tables (vitals) | | | | | | | |
| Frequency tables (diagnoses) | | | | | | | |
| Table 1 generation | | | | | | | |
| Plotting | | | | | | | |

**Notes:**
- [Any visualization performance issues]

---

### Total Project Cost

| Phase | Machine Hours | Est. Cost |
|-------|---------------|-----------|
| Phase 1: Extraction | | |
| Phase 2: Processing | | |
| Phase 3: Analysis | | |
| **Total** | | |

---

### Lessons Learned

**What worked well:**
- [Machine types that handled workloads efficiently]
- [Operations that were faster than expected]

**What didn't work:**
- [Memory issues, out-of-memory errors]
- [Performance bottlenecks]
- [Over-provisioned or under-provisioned resources]

**Recommendations for next time:**
- [Specific machine types for specific operations]
- [Workflow optimizations]
- [Cost-saving opportunities]

---

### Machine Type Decision Record

Document why you chose specific machine types:

**Interactive Notebook:**
- **Chose**: [e.g., n2-standard-64, 64 vCPU, 256GB RAM]
- **Why**: [e.g., "400K cohort with 5 large datasets in memory simultaneously"]
- **Result**: [e.g., "Peak usage 180GB, good fit"]
- **Next time**: [e.g., "Could try n2-standard-32 (128GB) for smaller cohorts"]

**Dsub Jobs (if used):**
- **Chose**: [e.g., n2-standard-8]
- **Why**: [e.g., "Processing per-chromosome or parallelizable tasks"]
- **Result**: [e.g., "Completed successfully"]
- **Next time**: [e.g., "Good sizing, keep same"]

---

### Quick Reference: Cost per Machine Type

For reference (us-central1 pricing as of 2026):

| Machine Type | vCPU | RAM | Cost/hour | Good For |
|--------------|------|-----|-----------|----------|
| n2-standard-4 | 4 | 16GB | $0.19 | Light processing, small cohorts |
| n2-standard-8 | 8 | 32GB | $0.39 | Medium processing |
| n2-standard-16 | 16 | 64GB | $0.77 | Medium-large processing |
| n2-standard-32 | 32 | 128GB | $1.55 | Large cohort notebooks |
| n2-standard-64 | 64 | 256GB | $3.10 | Very large cohorts, multiple datasets |
| n2-highmem-8 | 8 | 64GB | $0.52 | Memory-heavy, low CPU |
| n2-highmem-16 | 16 | 128GB | $1.04 | Large memory workloads |
| n2-highmem-32 | 32 | 256GB | $2.08 | Very memory-intensive processing |
| n2-highcpu-8 | 8 | 16GB | $0.30 | CPU-bound, modest memory |

*Preemptible instances: 60-80% cheaper, but can be interrupted*

---

## Data Quality Checks

### Medication Data Quality
- [ ] Verify RxNorm ingredient mapping completeness
- [ ] Check for missing route information
- [ ] Validate temporal alignment with fever episodes

### Lab Data Quality
- [ ] Verify unit standardization (WBC, CRP)
- [ ] Check for missing lab values
- [ ] Validate normal ranges

### Visit Context Quality
- [ ] Verify ER vs inpatient classification
- [ ] Check for missing visit linkage
- [ ] Validate visit dates align with fever dates

---

## Future Directions

### Autoencoder Extensions
- **Temporal autoencoders**: Incorporate time-series dynamics (LSTM/GRU encoders)
- **Attention mechanisms**: Learn which modalities/features are most important per patient
- **Contrastive learning**: Pre-train on unlabeled data, fine-tune for specific tasks
- **Transfer learning**: Apply learned embeddings to other fever-related outcomes

### Clinical Applications
- **Prognostic prediction**: Use embeddings to predict outcomes (hospitalization, complications)
- **Treatment response**: Predict antibiotic effectiveness from initial presentation
- **Severity stratification**: Identify high-risk fever phenotypes
- **Differential diagnosis support**: Suggest likely etiologies based on cluster membership

### Integration Opportunities
- **Genomic data**: Test if genetic variants associate with discovered fever clusters
- **Wearables data**: Add Fitbit heart rate, activity, sleep patterns to multi-modal input
- **Longitudinal analysis**: Track patients across multiple fever episodes
- **External validation**: Test learned representations on independent datasets

### Research Questions
- Do autoencoder-discovered clusters map to known fever etiologies?
- Can we identify novel disease state subtypes not captured by ICD codes?
- Which modality contributes most to cluster separation?
- How stable are cluster assignments over time (same patient, different episodes)?

---

## Metadata
- **Author**: Bennett Waxse
- **Created**: January 2026
- **Platform**: All of Us Research Workbench
- **Data Version**: CDR v8
- **Status**: In development