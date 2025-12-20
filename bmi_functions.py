"""
Enhanced BMI Harmonization Functions for All of Us Research
===========================================================

Standalone functions for BMI data extraction, harmonization, and temporal matching.
These functions provide sophisticated quality control and temporal logic for 
All of Us Research Program data.

Author: Bennett Waxse
Created: June 2025
Compatible: All of Us CDR v8
Dependencies: pandas, polars, numpy, matplotlib, seaborn, google-cloud-bigquery

Usage:
    from bmi_functions import harmonize_bmi_complete, hierarchical_temporal_matching
    
    # Complete pipeline
    bmi_data, matched_cohorts, quality = harmonize_bmi_complete(
        cohort_dfs, version, matching_strategy='standard'
    )
"""

import pandas as pd
import polars as pl
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from google.cloud import bigquery
from typing import List, Dict, Tuple, Optional, Union

# BMI-related concept IDs for All of Us
WEIGHT_CONCEPTS = [3010220, 3013762, 3025315, 3027492, 3023166]
HEIGHT_CONCEPTS = [3036798, 3023540, 3019171, 3036277, 40655804]
BMI_CONCEPTS = [3038553, 4245997]

def polars_gbq(query: str) -> pl.DataFrame:
    """
    Execute BigQuery SQL and return result as polars dataframe
    
    Args:
        query: BigQuery SQL query string
    
    Returns:
        pl.DataFrame: Query results
    """
    client = bigquery.Client()
    query_job = client.query(query)
    rows = query_job.result()
    df = pl.from_arrow(rows.to_arrow())
    return df

def extract_bmi_data(person_ids: Union[List, pd.Series], version: str) -> pl.DataFrame:
    """
    Extract BMI-related measurements for specified persons
    
    Args:
        person_ids: List or pandas Series of person_ids, or list of DataFrames
        version: BigQuery dataset name
    
    Returns:
        pl.DataFrame: Raw BMI measurements with units
    """
    # Handle multiple DataFrame inputs
    if isinstance(person_ids, list) and all(isinstance(x, pd.DataFrame) for x in person_ids):
        person_id_series = [df['person_id'] for df in person_ids]
        combined_person_ids = pd.concat(person_id_series)
        unique_person_ids = combined_person_ids.drop_duplicates()
    else:
        unique_person_ids = pd.Series(person_ids).drop_duplicates()
    
    person_ids_str = ', '.join(map(str, unique_person_ids))
    
    query = f"""
    WITH combined AS (
        SELECT 
            person_id, 
            measurement_date AS date, 
            MAX(IF(measurement_concept_id IN ({','.join(map(str, WEIGHT_CONCEPTS))}), 
                value_as_number, NULL)) AS wt,
            MAX(IF(measurement_concept_id IN ({','.join(map(str, WEIGHT_CONCEPTS))}), 
                c2.concept_name, NULL)) AS wt_units,
            MAX(IF(measurement_concept_id IN ({','.join(map(str, HEIGHT_CONCEPTS))}), 
                value_as_number, NULL)) AS ht,
            MAX(IF(measurement_concept_id IN ({','.join(map(str, HEIGHT_CONCEPTS))}), 
                c2.concept_name, NULL)) AS ht_units,
            MAX(IF(measurement_concept_id IN ({','.join(map(str, BMI_CONCEPTS))}), 
                value_as_number, NULL)) AS bmi,
            MAX(IF(measurement_concept_id IN ({','.join(map(str, BMI_CONCEPTS))}), 
                c2.concept_name, NULL)) AS bmi_units
        FROM 
            {version}.measurement m
        INNER JOIN 
            {version}.concept c2 ON unit_concept_id = c2.concept_id
        WHERE 
            measurement_concept_id IN ({','.join(map(str, WEIGHT_CONCEPTS + HEIGHT_CONCEPTS + BMI_CONCEPTS))})
            AND person_id IN ({person_ids_str})
        GROUP BY 
            person_id, measurement_date
    )
    SELECT * FROM combined ORDER BY person_id, date
    """
    
    return polars_gbq(query)

def clean_units(df: pl.DataFrame) -> pl.DataFrame:
    """
    Clean unit inconsistencies and invalid measurements
    
    Args:
        df: polars.DataFrame with measurement data
    
    Returns:
        pl.DataFrame: Cleaned data
    """
    # Clean weight units
    df = df.with_columns(
        pl.when(pl.col('wt_units') == 'No matching concept')
          .then(pl.lit(None))
          .otherwise(pl.col('wt_units'))
          .alias('wt_units')
    )
    
    # Clean height units - remove ambiguous measurements
    df = df.with_columns([
        pl.when(pl.col('ht_units').is_in(["percent", "No matching concept"]))
          .then(pl.lit(None))
          .otherwise(pl.col('ht_units'))
          .alias('ht_units'),
        pl.when(pl.col('ht_units').is_in(["percent", "No matching concept"]))
          .then(pl.lit(None))
          .otherwise(pl.col('ht'))
          .alias('ht')
    ])
    
    # Clean BMI units - remove ambiguous measurements
    df = df.with_columns([
        pl.when(pl.col('bmi_units').is_in(["ratio", "no value", "No matching concept"]))
          .then(pl.lit(None))
          .otherwise(pl.col('bmi_units'))
          .alias('bmi_units'),
        pl.when(pl.col('bmi_units').is_in(["ratio", "no value", "No matching concept"]))
          .then(pl.lit(None))
          .otherwise(pl.col('bmi'))
          .alias('bmi')
    ])
    
    return df

def apply_clean_outliers(df: pl.DataFrame, column: str, unit: Optional[str] = None, 
                        sigma_threshold: float = 4) -> pl.DataFrame:
    """
    Remove outliers using sigma-based thresholds by unit type
    
    Args:
        df: polars.DataFrame
        column: Column name to clean
        unit: Unit column name (optional)
        sigma_threshold: Number of standard deviations for outlier threshold
    
    Returns:
        pl.DataFrame: Data with outliers removed
    """
    if unit:
        # Calculate stats by unit type
        stats = df.filter(pl.col(unit).is_not_null()).group_by(unit).agg([
            pl.col(column).mean().alias(f'{column}_mean'),
            pl.col(column).std().alias(f'{column}_std')
        ])
        df = df.join(stats, on=unit, how='left')
        
        # Calculate bounds and remove outliers
        lower_bound = pl.col(f'{column}_mean') - sigma_threshold * pl.col(f'{column}_std')
        upper_bound = pl.col(f'{column}_mean') + sigma_threshold * pl.col(f'{column}_std')
        
        df = df.with_columns(
            pl.when((pl.col(column) < lower_bound) | (pl.col(column) > upper_bound))
              .then(None)
              .otherwise(pl.col(column))
              .alias(column)
        ).drop([f'{column}_mean', f'{column}_std'])
    else:
        # Global outlier removal
        mean = df[column].mean()
        std = df[column].std()
        lower_bound = mean - sigma_threshold * std
        upper_bound = mean + sigma_threshold * std
        
        df = df.with_columns(
            pl.when((pl.col(column) < lower_bound) | (pl.col(column) > upper_bound))
              .then(None)
              .otherwise(pl.col(column))
              .alias(column)
        )
    
    return df

def convert_units(df: pl.DataFrame) -> pl.DataFrame:
    """
    Convert measurements to standard units (kg, cm)
    
    Args:
        df: polars.DataFrame with measurements
    
    Returns:
        pl.DataFrame: Standardized units
    """
    # Convert height: inches to centimeters
    df = df.with_columns(
        pl.when(pl.col("ht_units") == "inch (US)")
          .then(pl.col("ht") * 2.54)
          .otherwise(pl.col("ht"))
          .alias("ht")
    )
    
    # Convert weight: pounds to kilograms
    df = df.with_columns(
        pl.when(pl.col("wt_units") == "pound (US)")
          .then(pl.col("wt") * 0.45359237)
          .otherwise(pl.col("wt"))
          .alias("wt")
    )
    
    return df

def enhanced_height_carryforward(df: pl.DataFrame, max_carryforward_days: int = 1095) -> pl.DataFrame:
    """
    Carry forward height measurements with configurable time limits
    
    Args:
        df: polars DataFrame with height measurements
        max_carryforward_days: Maximum days to carry forward height (default: 3 years)
    
    Returns:
        polars DataFrame with carried forward heights and source flags
    """
    df = df.sort("person_id", "date")
    
    # Create forward-filled height
    df = df.with_columns([
        pl.col("date").cast(pl.Date),
        pl.col("ht_cm").forward_fill().over("person_id").alias("ht_cm_filled")
    ])
    
    # Calculate days between measurements
    df = df.with_columns([
        pl.col("date").diff().over("person_id").dt.total_days().alias("days_since_last")
    ])
    
    # Use carried forward height only within time limit
    df = df.with_columns([
        pl.when(
            (pl.col("ht_cm").is_null()) & 
            (pl.col("days_since_last") <= max_carryforward_days)
        )
        .then(pl.col("ht_cm_filled"))
        .otherwise(pl.col("ht_cm"))
        .alias("ht_cm_final"),
        
        # Add source flag
        pl.when(pl.col("ht_cm").is_not_null())
        .then(pl.lit("measured"))
        .when(
            (pl.col("ht_cm").is_null()) & 
            (pl.col("days_since_last") <= max_carryforward_days)
        )
        .then(pl.lit("carried_forward"))
        .otherwise(pl.lit("missing"))
        .alias("height_source")
    ])
    
    # Clean up temporary columns
    df = df.drop(["ht_cm_filled", "days_since_last"])
    df = df.rename({"ht_cm_final": "ht_cm"})
    
    return df

def calculate_bmi_with_validation(df: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate BMI from height/weight and validate against recorded BMI
    
    Args:
        df: polars.DataFrame with height and weight
    
    Returns:
        pl.DataFrame: With calculated BMI and quality metrics
    """
    # Calculate BMI from height/weight
    df = df.with_columns(
        (pl.col("wt_kg") / (pl.col("ht_cm") / 100) ** 2).alias("bmi_calc")
    )
    
    # Calculate difference between recorded and calculated BMI
    df = df.with_columns(
        (pl.col("bmi") - pl.col("bmi_calc")).alias("bmi_diff")
    )
    
    # Use calculated BMI where recorded BMI is missing
    df = df.with_columns(
        pl.when(pl.col("bmi").is_null())
          .then(pl.col("bmi_calc"))
          .otherwise(pl.col("bmi"))
          .alias("bmi")
    )
    
    return df

def hierarchical_temporal_matching(cohort_df: pd.DataFrame, bmi_df: pd.DataFrame, 
                                 time_col: str = 'time_zero', include_post: bool = True, 
                                 max_prior_days: int = 365, max_post_days: int = 90) -> pd.DataFrame:
    """
    Hierarchical temporal matching with preference for prior measurements
    
    Args:
        cohort_df: DataFrame with cohort and time reference column
        bmi_df: DataFrame with BMI measurements and date column
        time_col: Column name for reference time in cohort_df
        include_post: Whether to include post-time_zero measurements
        max_prior_days: Maximum days before time_zero to consider
        max_post_days: Maximum days after time_zero to consider
    
    Returns:
        DataFrame: Cohort with matched BMI and timing flags
    """
    # Ensure datetime format
    cohort_df = cohort_df.copy()
    bmi_df = bmi_df.copy()
    
    cohort_df[time_col] = pd.to_datetime(cohort_df[time_col]).dt.tz_localize(None)
    bmi_df['date'] = pd.to_datetime(bmi_df['date']).dt.tz_localize(None)
    
    # Merge all possible matches
    merged = pd.merge(cohort_df, bmi_df, on='person_id', how='left')
    
    # Calculate temporal relationships
    merged['days_diff'] = (merged['date'] - merged[time_col]).dt.days
    merged['abs_days_diff'] = merged['days_diff'].abs()
    
    # Apply time window filters
    if include_post:
        valid_measurements = (
            (merged['days_diff'] <= 0) & (merged['abs_days_diff'] <= max_prior_days) |
            (merged['days_diff'] > 0) & (merged['days_diff'] <= max_post_days)
        )
    else:
        valid_measurements = (
            (merged['days_diff'] <= 0) & (merged['abs_days_diff'] <= max_prior_days)
        )
    
    merged = merged[valid_measurements]
    
    # Hierarchical priority scoring (lower = better)
    merged['priority'] = np.where(
        merged['days_diff'] <= 0,  # Prior measurements
        merged['abs_days_diff'],   # Prefer closer to time_zero
        1000 + merged['abs_days_diff']  # Heavily penalize post-time_zero
    )
    
    # Keep best match per person
    merged = merged.sort_values(['person_id', 'priority'])
    result = merged.drop_duplicates('person_id', keep='first')
    
    # Add interpretive flags
    result['bmi_timing'] = np.where(
        result['days_diff'] <= 0, 'prior', 'post'
    )
    
    result['bmi_quality'] = np.where(
        result['abs_days_diff'] <= 30, 'excellent',
        np.where(result['abs_days_diff'] <= 90, 'good',
                np.where(result['abs_days_diff'] <= 180, 'acceptable', 'poor'))
    )
    
    # Select final columns
    original_cols = cohort_df.columns.tolist()
    bmi_cols = ['date', 'bmi', 'wt_kg', 'ht_cm', 'days_diff', 'bmi_timing', 'bmi_quality']
    available_bmi_cols = [col for col in bmi_cols if col in result.columns]
    result = result[original_cols + available_bmi_cols]
    
    return result

def analyze_temporal_matching_quality(merged_cohorts: List[pd.DataFrame]) -> Dict:
    """
    Analyze the quality of temporal matching across cohorts
    
    Args:
        merged_cohorts: List of DataFrames from hierarchical_temporal_matching
    
    Returns:
        dict: Summary statistics of matching quality
    """
    all_matches = pd.concat(merged_cohorts, ignore_index=True)
    
    summary = {
        'total_participants': len(all_matches),
        'participants_with_bmi': len(all_matches.dropna(subset=['bmi'])),
        'matching_rate': len(all_matches.dropna(subset=['bmi'])) / len(all_matches),
        
        # Timing distribution
        'timing_distribution': all_matches['bmi_timing'].value_counts().to_dict() if 'bmi_timing' in all_matches.columns else {},
        'timing_percentages': all_matches['bmi_timing'].value_counts(normalize=True).to_dict() if 'bmi_timing' in all_matches.columns else {},
        
        # Quality distribution
        'quality_distribution': all_matches['bmi_quality'].value_counts().to_dict() if 'bmi_quality' in all_matches.columns else {},
        'quality_percentages': all_matches['bmi_quality'].value_counts(normalize=True).to_dict() if 'bmi_quality' in all_matches.columns else {},
        
        # Temporal statistics
        'days_diff_stats': {
            'mean': all_matches['days_diff'].mean() if 'days_diff' in all_matches.columns else None,
            'median': all_matches['days_diff'].median() if 'days_diff' in all_matches.columns else None,
            'std': all_matches['days_diff'].std() if 'days_diff' in all_matches.columns else None,
            'min': all_matches['days_diff'].min() if 'days_diff' in all_matches.columns else None,
            'max': all_matches['days_diff'].max() if 'days_diff' in all_matches.columns else None
        }
    }
    
    return summary

def plot_temporal_matching_quality(merged_cohorts: List[pd.DataFrame], figsize: Tuple[int, int] = (15, 5)):
    """
    Create visualizations of temporal matching quality
    
    Args:
        merged_cohorts: List of DataFrames from hierarchical_temporal_matching
        figsize: Figure size tuple
    
    Returns:
        matplotlib.figure.Figure: The created figure
    """
    all_matches = pd.concat(merged_cohorts, ignore_index=True)
    all_matches = all_matches.dropna(subset=['bmi'])
    
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    
    # Days difference distribution
    if 'days_diff' in all_matches.columns:
        axes[0].hist(all_matches['days_diff'], bins=50, alpha=0.7, edgecolor='black')
        axes[0].axvline(0, color='red', linestyle='--', label='Time Zero')
        axes[0].set_xlabel('Days from Time Zero')
        axes[0].set_ylabel('Count')
        axes[0].set_title('BMI Measurement Timing Distribution')
        axes[0].legend()
    
    # Timing categories
    if 'bmi_timing' in all_matches.columns:
        timing_counts = all_matches['bmi_timing'].value_counts()
        axes[1].pie(timing_counts.values, labels=timing_counts.index, autopct='%1.1f%%')
        axes[1].set_title('Prior vs Post Time Zero')
    
    # Quality categories
    if 'bmi_quality' in all_matches.columns:
        quality_counts = all_matches['bmi_quality'].value_counts()
        axes[2].bar(quality_counts.index, quality_counts.values)
        axes[2].set_xlabel('Temporal Quality')
        axes[2].set_ylabel('Count')
        axes[2].set_title('BMI Measurement Quality')
        axes[2].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.show()
    
    return fig

def harmonize_bmi_complete(cohort_dfs: List[pd.DataFrame], version: str, 
                          matching_strategy: str = 'standard', time_col: str = 'time_zero',
                          height_carryforward_days: int = 1095, **matching_kwargs) -> Tuple[pl.DataFrame, List[pd.DataFrame], Dict]:
    """
    Complete BMI harmonization pipeline with enhanced temporal matching
    
    Args:
        cohort_dfs: List of DataFrames with person_id columns
        version: All of Us CDR version (workspace variable)
        matching_strategy: 'standard', 'strict_prior', 'flexible', or 'drug_study'
        time_col: Column name for reference time in cohort DataFrames
        height_carryforward_days: Days to carry forward height measurements
        **matching_kwargs: Additional arguments for hierarchical_temporal_matching
    
    Returns:
        tuple: (harmonized_bmi_df, merged_cohort_dfs, quality_summary)
    """
    print("=== Complete BMI Harmonization Pipeline ===")
    
    # 1. Extract BMI data
    print("Extracting BMI measurements...")
    bmi_df = extract_bmi_data(cohort_dfs, version)
    print(f"Raw measurements extracted: {len(bmi_df):,}")
    
    # 2. Clean units and remove invalid measurements
    print("Cleaning units and invalid measurements...")
    bmi_df = clean_units(bmi_df)
    
    # 3. Remove outliers by unit type
    print("Removing outliers by unit type...")
    bmi_df = apply_clean_outliers(bmi_df, "wt", "wt_units")
    bmi_df = apply_clean_outliers(bmi_df, "ht", "ht_units")
    bmi_df = apply_clean_outliers(bmi_df, "bmi")
    
    # 4. Convert to standard units
    print("Converting to standard units (kg, cm)...")
    bmi_df = convert_units(bmi_df)
    
    # 5. Filter rows with measurements and rename columns
    bmi_df = bmi_df.filter(
        pl.col("wt").is_not_null() | pl.col("ht").is_not_null() | pl.col("bmi").is_not_null()
    )
    bmi_df = bmi_df.rename({"wt": "wt_kg", "ht": "ht_cm"})
    bmi_df = bmi_df.drop(["wt_units", "ht_units", "bmi_units"])
    
    # 6. Enhanced height carry-forward
    print(f"Applying height carry-forward ({height_carryforward_days} days)...")
    bmi_df = enhanced_height_carryforward(bmi_df, max_carryforward_days=height_carryforward_days)
    
    # 7. Calculate BMI with validation
    print("Calculating BMI with validation...")
    bmi_df = calculate_bmi_with_validation(bmi_df)
    
    # 8. Final cleanup
    bmi_df = bmi_df.select(['person_id', 'date', 'wt_kg', 'ht_cm', 'bmi', 'height_source'])
    bmi_df = bmi_df.filter(pl.col('bmi').is_not_null())
    bmi_df = bmi_df.sort("person_id", "date")
    
    print(f"Final harmonized measurements: {len(bmi_df):,}")
    
    # 9. Define matching parameters based on strategy
    matching_params = {
        'standard': {
            'include_post': True,
            'max_prior_days': 365,
            'max_post_days': 90
        },
        'strict_prior': {
            'include_post': False,
            'max_prior_days': 365,
            'max_post_days': 0
        },
        'flexible': {
            'include_post': True,
            'max_prior_days': 730,
            'max_post_days': 180
        },
        'drug_study': {
            'include_post': False,
            'max_prior_days': 90,
            'max_post_days': 0
        }
    }
    
    # Override with user-provided parameters
    params = matching_params.get(matching_strategy, matching_params['standard'])
    params.update(matching_kwargs)
    
    print(f"Using {matching_strategy} matching strategy:")
    print(f"  Include post: {params['include_post']}")
    print(f"  Max prior days: {params['max_prior_days']}")
    print(f"  Max post days: {params['max_post_days']}")
    
    # 10. Temporal matching with cohorts
    bmi_pandas = bmi_df.to_pandas()
    merged_cohorts = [
        hierarchical_temporal_matching(
            df, bmi_pandas,
            time_col=time_col,
            **params
        )
        for df in cohort_dfs
    ]
    
    # 11. Quality assessment
    quality_summary = analyze_temporal_matching_quality(merged_cohorts)
    
    print(f"=== Harmonization Complete ===")
    print(f"Overall matching rate: {quality_summary['matching_rate']:.1%}")
    if quality_summary['timing_percentages']:
        print(f"Prior/Post distribution: {quality_summary['timing_percentages']}")
    
    return bmi_df, merged_cohorts, quality_summary

# Convenience functions for different research scenarios
def bmi_for_gwas(cohort_dfs: List[pd.DataFrame], version: str, time_col: str = 'enrollment_date') -> Tuple[pl.DataFrame, List[pd.DataFrame], Dict]:
    """GWAS-specific BMI harmonization (strict prior only)"""
    return harmonize_bmi_complete(
        cohort_dfs, version, 
        matching_strategy='strict_prior',
        time_col=time_col,
        max_prior_days=365
    )

def bmi_for_drug_study(cohort_dfs: List[pd.DataFrame], version: str, time_col: str = 'drug_start_date') -> Tuple[pl.DataFrame, List[pd.DataFrame], Dict]:
    """Drug study BMI harmonization (recent baseline only)"""
    return harmonize_bmi_complete(
        cohort_dfs, version,
        matching_strategy='drug_study', 
        time_col=time_col,
        max_prior_days=90
    )

def bmi_flexible(cohort_dfs: List[pd.DataFrame], version: str, time_col: str = 'time_zero') -> Tuple[pl.DataFrame, List[pd.DataFrame], Dict]:
    """Flexible BMI harmonization for longitudinal studies"""
    return harmonize_bmi_complete(
        cohort_dfs, version,
        matching_strategy='flexible',
        time_col=time_col,
        max_prior_days=730,
        max_post_days=180
    )

# Legacy function for backward compatibility
def merge_closest_bmi(cohort_df: pd.DataFrame, bmi_df: pd.DataFrame) -> pd.DataFrame:
    """Legacy function - use hierarchical_temporal_matching instead"""
    return hierarchical_temporal_matching(
        cohort_df, bmi_df,
        time_col='time_zero',
        include_post=True,
        max_prior_days=365,
        max_post_days=365
    )
