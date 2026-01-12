---
name: check-ancestry-consistency
enabled: true
event: file
action: warn
conditions:
  - field: new_text
    operator: regex_match
    pattern: (EUR|AFR|AMR|EAS|SAS)
  - field: file_path
    operator: regex_match
    pattern: (GWAS|SAIGE|B02|B03).*\.(py|ipynb)
---

🧬 **Ancestry code detected in GWAS file**

You're working with ancestry-specific code in a GWAS pipeline file.

**Verify consistency across your pipeline**:

1. **Cohort Definition** (B01 files):
   ```python
   # Check SQL WHERE clause for ancestry filtering
   WHERE ancestry_pred_other IN ('EUR', 'AFR', 'AMR', ...)
   ```

2. **SAIGE GWAS** (B02 files):
   ```bash
   # Verify --ancestry parameter or input file filtering
   --ancestry=EUR
   # OR ancestry-specific phenotype file
   ```

3. **METAL Meta-Analysis** (B03 files):
   ```python
   # Ensure all ancestries are included in input files
   ancestries = ['EUR', 'AFR', 'AMR']
   for anc in ancestries:
       assert results_exist(anc), f"Missing results for {anc}"
   ```

4. **Output Paths**:
   ```python
   # Confirm ancestry in filenames for traceability
   output = f"gs://{bucket}/gwas/{phenotype}_{ancestry}_results.txt"
   ```

**All of Us ancestry distributions** (approximate WGS counts):
- `EUR` (European): ~190,000 participants
- `AFR` (African): ~40,000 participants
- `AMR` (Admixed American): ~10,000 participants
- `EAS` (East Asian): ~4,000 participants
- `SAS` (South Asian): ~1,000 participants

**Common mistakes**:
- Filtering cohort to EUR but running SAIGE on all ancestries
- Hardcoding ancestry list in meta-analysis but missing one GWAS run
- Inconsistent ancestry labels across files (e.g., 'eur' vs 'EUR')

**Best practice**: Define ancestry list once at top of pipeline:
```python
ANCESTRIES = ['EUR', 'AFR', 'AMR']  # Define once, use everywhere
```
