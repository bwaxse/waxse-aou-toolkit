# Fever Cohort todos for BENNETT
# Don't edit, Claude

### Action Items
- [ ] Decide on medication concepts and update SELECTED_MED_CONCEPTS
- [ ] Decide on lab concepts and update SELECTED_LAB_CONCEPTS
- [ ] Decide on vital concepts and update SELECTED_VITAL_CONCEPTS
- [ ] Update total cohort filter to require data from these specific concepts

- [ ] Review temperature concept output, but adjust for F to C conversion
- [ ] Unifier update to recent package
- [ ] Compare Core and Peripheral temperature distributions -- different threshold for fever?

## Ground Truth Labels (for Validation)

- [ ] Extract ground truth ICD diagnosis codes for validation
  - **Temporal rule**: Only use diagnoses WITHIN the macrovisit (between `macrovisit_start_date` and `macrovisit_end_date`)
  - **No pre-existing diagnoses**: Exclude ICD codes from before admission
  - **Source tables**: `condition_occurrence` joined to macrovisits
  - **Vocabularies**: ICD9CM, ICD10CM, SNOMED