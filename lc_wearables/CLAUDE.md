# LC Wearables Analysis - Claude Code Documentation

## Project Overview
This directory contains analysis pipelines for Lincoln County (LC) wearables data from the All of Us Research Program. The project focuses on processing, analyzing, and visualizing wearable device data (Fitbit) to understand health patterns and relationships between physical activity, sleep, and health outcomes.

## Directory Structure
```
lc_wearables/
├── notebooks/
│   └── wearables_data_pull_export_20241227_cleaned.ipynb  # Main analysis notebook
├── python/
│   └── wearables_data_pull_export_20241227_cleaned.py     # Exported Python script
└── CLAUDE.md                                               # This documentation
```

## Data Sources
- **Fitbit Wearables Data**: Device-level data from All of Us participants
  - Heart rate measurements
  - Steps and activity levels
  - Sleep patterns and quality
  - Distance and calories

## Main Workflow

### 1. Environment Setup (`wearables_data_pull_export_20241227_cleaned.ipynb:1-2`)
```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime, timedelta
```

Standard scientific Python stack with visualization libraries for exploratory data analysis.

### 2. Data Loading and Exploration
The notebook performs comprehensive data exploration including:
- Loading wearables device data
- Examining data structure and completeness
- Identifying key variables and measurement types
- Temporal analysis of data collection patterns

### 3. Analysis Components

#### Activity Analysis
- **Daily step counts**: Aggregation and trending
- **Activity levels**: Sedentary, lightly active, moderately active, very active
- **Distance tracking**: Total daily distance covered
- **Calorie expenditure**: Activity-based calorie calculations

#### Heart Rate Analysis
- **Resting heart rate**: Baseline cardiovascular metrics
- **Heart rate variability**: Health and recovery indicators
- **Heart rate zones**: Time spent in different intensity zones

#### Sleep Analysis
- **Sleep duration**: Total sleep time per night
- **Sleep stages**: Deep, light, REM, and awake time
- **Sleep efficiency**: Quality metrics
- **Sleep timing**: Bedtime and wake time patterns

### 4. Visualization and Reporting
- Time series plots for temporal trends
- Distribution analyses for population-level patterns
- Correlation analyses between different health metrics
- Summary statistics and data quality assessments

## Key Functions and Reusable Components

### Data Processing Functions
The notebook likely contains reusable functions for:
- **Date/time handling**: Converting timestamps and aggregating by time periods
- **Data cleaning**: Handling missing values and outliers in wearables data
- **Aggregation**: Daily, weekly, and monthly summaries
- **Feature engineering**: Derived metrics from raw sensor data

### Visualization Functions
- **Time series plotting**: Standardized plots for temporal data
- **Distribution plots**: Histograms and density plots for health metrics
- **Correlation heatmaps**: Relationships between variables
- **Multi-panel dashboards**: Comprehensive health summaries

## Integration with Toolkit

### Shared Resources
This project uses common toolkit components:
- **AoU Utilities** (`aou_utils/`): All of Us workspace interaction functions
- **Analysis Utils** (`analysis_utils/`): Statistical and data processing tools
- **Visualization Utils** (`viz_utils/`): Standardized plotting functions

### Reusable Patterns
The wearables analysis follows established toolkit patterns:
1. **Data Loading**: Using BigQuery and workspace storage
2. **Data Validation**: Quality checks and completeness assessment
3. **Exploratory Analysis**: Systematic examination of data characteristics
4. **Visualization**: Consistent styling and formatting
5. **Export**: Saving cleaned data and results

## Usage Examples

### Loading Wearables Data
```python
# Example pattern (specific code in notebook)
# Load Fitbit device data
wearables_df = pd.read_csv('fitbit_data.csv')

# Parse dates
wearables_df['date'] = pd.to_datetime(wearables_df['date'])

# Set date as index for time series operations
wearables_df.set_index('date', inplace=True)
```

### Daily Activity Summary
```python
# Aggregate daily metrics
daily_summary = wearables_df.groupby('date').agg({
    'steps': 'sum',
    'distance': 'sum',
    'calories': 'sum',
    'heart_rate': 'mean'
})
```

### Sleep Analysis
```python
# Sleep duration and quality metrics
sleep_metrics = wearables_df.groupby('date').agg({
    'sleep_minutes': 'sum',
    'deep_sleep_minutes': 'sum',
    'light_sleep_minutes': 'sum',
    'rem_sleep_minutes': 'sum',
    'sleep_efficiency': 'mean'
})
```

## Data Quality Considerations

### Missing Data
- **Device compliance**: Not all participants wear devices consistently
- **Measurement gaps**: Technical issues or device removal
- **Sync delays**: Data may not upload immediately

### Validation Checks
- **Physiologically plausible ranges**: Filter unrealistic values
- **Temporal consistency**: Check for date/time errors
- **Completeness assessment**: Track missingness patterns

## Analysis Best Practices

1. **Temporal Alignment**: Ensure all metrics are properly aligned by date/time
2. **Individual Variation**: Account for person-specific baselines
3. **Seasonal Effects**: Consider time-of-year influences on activity
4. **Weekend vs Weekday**: Separate analysis for different day types
5. **Data Quality Filtering**: Exclude incomplete or unreliable measurements

## Next Steps and Extensions

### Potential Enhancements
1. **Longitudinal Analysis**: Track individual changes over time
2. **Health Correlations**: Link wearables data to EHR outcomes
3. **Predictive Modeling**: Use activity patterns to predict health events
4. **Clustering**: Identify participant phenotypes based on activity patterns
5. **Comparative Analysis**: Compare Lincoln County patterns to national data

### Integration Opportunities
- Link to genomics data for gene-environment interactions
- Combine with survey data for behavioral insights
- Merge with clinical data for health outcome prediction
- Connect to environmental data (weather, air quality)

## Dependencies
```
pandas>=1.3.0
numpy>=1.20.0
matplotlib>=3.4.0
seaborn>=0.11.0
scipy>=1.7.0
```

## Notes for Claude Code

### When Working with This Project
1. The notebook contains the primary analysis - reference it for specific implementations
2. Look for reusable data processing functions that could be extracted to `analysis_utils/`
3. Visualization code may be useful for other projects - consider adding to `viz_utils/`
4. Data loading patterns specific to wearables may benefit other device data analyses

### Common Tasks
- **Add new metrics**: Follow existing aggregation patterns in the notebook
- **Modify visualizations**: Use consistent color schemes and styling from toolkit
- **Export results**: Save processed data in standard formats (CSV, parquet)
- **Document findings**: Create markdown summaries for stakeholders

### Troubleshooting
- **Memory issues**: Wearables data can be large - use chunking for big datasets
- **Date parsing**: Ensure timezone awareness when working with timestamps
- **Missing values**: Use forward-fill or interpolation cautiously with time series
- **Outliers**: Apply domain knowledge for physiologically valid ranges

## References
- All of Us Research Program: https://www.researchallofus.org/
- Fitbit Data Documentation: Available in AoU Researcher Workbench
- Related Projects: See `gwas_sarcoid/` for similar analysis patterns
