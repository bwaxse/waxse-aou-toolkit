# waxse-aou-toolkit

Production-ready code for genomic and epidemiological analyses using All of Us Research Workbench data.

**Platform**: All of Us Research Workbench (OMOP CDM v8)
**Author**: Bennett Waxse
**License**: MIT

---

## 🚀 Quick Start for Claude Code

### I want to...

| **Task** | **Read First** | **Then Use** |
|----------|---------------|--------------|
| Run **GWAS** on a new phenotype | [`hpv/CLAUDE.md`](hpv/CLAUDE.md) | HPV B01-B03 files as template |
| Harmonize **BMI** or similar covariate | [`bmi/CLAUDE.md`](bmi/CLAUDE.md) | `bmi_functions.py` library |
| Understand **architectural patterns** | [`CLAUDE_ARCHITECTURE.md`](CLAUDE_ARCHITECTURE.md) | Cross-project reusable components |
| Do **PheWAS** analysis | [`hpv/CLAUDE.md`](hpv/CLAUDE.md) → B01 | PheTK integration example |
| Create **HLA visualizations** | [`hpv/CLAUDE.md`](hpv/CLAUDE.md) → B06 | Manhattan/forest plot functions |
| Match **time-stamped data** to cohorts | [`bmi/CLAUDE.md`](bmi/CLAUDE.md) | `hierarchical_temporal_matching()` |

---

## 📁 Documentation Structure

```
├── README.md (you are here)           ← Start here - quick decision guide
├── CLAUDE_ARCHITECTURE.md             ← Cross-project patterns & reusable functions
├── hpv/
│   ├── CLAUDE.md                      ← HPV GWAS: Complete 7-stage pipeline
│   ├── B01.a HPV Cohort v0.py         ← Cohort definition template
│   ├── B02.a SAIGE GWAS v0.py         ← SAIGE workflow template
│   ├── B03 METAL Meta-analysis.py     ← Meta-analysis template
│   └── B06.1 Plot HLA Results.py      ← HLA visualization template
└── bmi/
    ├── CLAUDE.md                      ← BMI harmonization: Complete function library
    ├── bmi_functions.py               ← Importable functions (USE THESE!)
    └── bmi_harmonization.ipynb        ← Usage examples & demonstrations
```

---

## 🎯 Available Projects

### 1. **HPV** - Multi-Ancestry GWAS Pipeline
**Type**: Genome-wide association study with meta-analysis
**Status**: Production
**Documentation**: [`hpv/CLAUDE.md`](hpv/CLAUDE.md)

**What it does**:
- Extracts HPV case/control cohorts from OMOP data
- Runs SAIGE GWAS across 5 ancestries (EUR, AFR, AMR, EAS, SAS)
- Meta-analyzes results with METAL
- Creates HLA region visualizations
- Distributed computing via Google Cloud Batch (dsub)

**Reusable for**: Any binary or quantitative trait GWAS

---

### 2. **BMI** - Covariate Harmonization Library
**Type**: Data harmonization with temporal matching
**Status**: Production
**Documentation**: [`bmi/CLAUDE.md`](bmi/CLAUDE.md)

**What it does**:
- Extracts BMI/weight/height from OMOP measurement table
- Sophisticated cleaning (unit conversion, outlier detection)
- Temporal matching with hierarchical preference for prior measurements
- Quality assessment and visualization

**Reusable for**: Any time-stamped measurement (labs, vitals, anthropometrics)

---

## 🛠️ Adding a New Project

### Step 1: Choose Template
- **Genetic association** → Use **HPV** as template
- **Covariate harmonization** → Use **BMI** as template

### Step 2: Copy & Modify
```bash
# For GWAS:
cp hpv/B01.a\ HPV\ Cohort\ v0.py my_project/B01_my_phenotype_cohort.py
# Then swap concept IDs, keep infrastructure code

# For harmonization:
cp bmi/bmi_functions.py my_project/my_measurement_functions.py
# Then swap measurement concepts, keep cleaning/matching logic
```

### Step 3: Read Relevant CLAUDE.md
Each project has detailed documentation explaining:
- What each file does
- Which functions to reuse
- How to adapt for your use case
- Common pitfalls to avoid

See [`CLAUDE_ARCHITECTURE.md`](CLAUDE_ARCHITECTURE.md) → "Adding a New Project" for detailed workflow.

---

## 📚 Key Reusable Functions

| Function | Location | Use For |
|----------|----------|---------|
| `polars_gbq()` | `bmi/bmi_functions.py:36` | BigQuery → polars DataFrame |
| `hierarchical_temporal_matching()` | `bmi/bmi_functions.py:307` | Match measurements to time points |
| `dsub_script()` | `hpv/B02...py:141` | Distributed computing (Cloud Batch) |
| `run_saige_null()` + `run_saige_test()` | `hpv/B02...py` | Complete SAIGE GWAS workflow |
| `run_metal()` | `hpv/B03...py:525` | METAL meta-analysis |
| `update_covariates()` | `hpv/B01...py:762` | Demographic covariate engineering |

Full catalog: [`CLAUDE_ARCHITECTURE.md`](CLAUDE_ARCHITECTURE.md) → "Reusable Function Library"

---

## 🔬 Requirements

- Access to All of Us Researcher Workbench
- Python 3.8+ with: pandas, polars, numpy, google-cloud-bigquery
- For GWAS: Docker containers (SAIGE, METAL) via dsub

Code works with other OMOP implementations with adjustments to concept IDs.

---

## ⚠️ Important Notes

### All of Us Data Policies
- **Suppress small counts**: Display "<20" for counts 1-19 in all outputs
- **Genomic data**: Filter to `has_whole_genome_variant = 1` before GWAS
- **Workspace variables**: Always use `os.environ['WORKSPACE_CDR']` and `os.environ['WORKSPACE_BUCKET']`

### Best Practices
- Always validate results against your research questions
- Review data use agreements before sharing derived datasets
- Use `polars_gbq()` for BigQuery (faster than pandas)
- Use `dsub` for parallelizable tasks (GWAS per chromosome)
- Store outputs in GCS, not local filesystem

---

## 🤝 Contributing

Found a bug? Have a useful snippet to share? Please open an issue or submit a pull request.

When adding new projects:
1. Create `{project}/CLAUDE.md` following the HPV/BMI template
2. Update `CLAUDE_ARCHITECTURE.md` with new reusable functions
3. Update this README's "Available Projects" section

---

## 📖 Support

- **All of Us Help**: [Researcher Help Desk](https://support.researchallofus.org/hc/en-us)
- **Documentation**: Start with project-specific `CLAUDE.md` files
- **Architecture**: See [`CLAUDE_ARCHITECTURE.md`](CLAUDE_ARCHITECTURE.md) for patterns

---

**Last updated**: December 2025
