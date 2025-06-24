# Utilities

Helper functions for common All of Us research tasks.

## Available Modules

### `dsub_helpers.py` ⭐
**Purpose**: Simplify dsub job management in All of Us Researcher Workbench

**Key Functions**:
- `dsub_script()` - Submit jobs with comprehensive parameter handling
- `check_dsub_status()` - Monitor job status with filtering options
- `job_details()` - Get detailed job information including failures
- `cancel_job()` - Cancel specific jobs safely
- `get_job_logs()` - Find job log locations

**Quick Start**:
```python
from utilities.dsub_helpers import dsub_script, check_dsub_status

# Submit a job
dsub_script(
    label='bmi-analysis',
    machine_type='n1-standard-4',
    script='my_analysis.sh',
    envs={'COHORT_SIZE': '10000'},
    in_params={'INPUT_DATA': 'gs://my-bucket/data.csv'},
    out_params={'RESULTS': 'gs://my-bucket/results/'}
)

# Check status
check_dsub_status(age='1d', full=True)
```

## Features

### 🚀 Job Submission
- **Environment integration** - Uses workspace variables automatically
- **Comprehensive logging** - Structured log paths with timestamps
- **Parameter validation** - Input validation with helpful error messages

### 📊 Job Monitoring
- **Flexible status checking** - Filter by age, user, job ID
- **Detailed job information** - Full job details including failure reasons
- **Log path discovery** - Easy access to job logs in workspace bucket

### 🛠 Job Management
- **Safe job cancellation** - Cancel specific jobs or all running jobs
- **User-scoped operations** - Defaults to current user, supports multi-user environments
- **Comprehensive error handling** - Clear error messages and validation

### 🎯 Convenience Functions
- `submit_standard_job()` - n1-standard-4 with sensible defaults
- `submit_highmem_job()` - n1-highmem-8 for memory-intensive tasks  
- `submit_compute_job()` - c4-standard-8 for CPU-intensive work

## Usage Examples

### Basic Job Submission
```python
from utilities.dsub_helpers import dsub_script

# Simple analysis job
dsub_script(
    label='gwas-analysis',
    machine_type='n1-standard-8',
    script='run_gwas.sh',
    envs={
        'PHENOTYPE': 'diabetes',
        'ANCESTRY': 'EUR'
    },
    in_params={
        'GENOTYPE_DATA': 'gs://bucket/genotypes.pgen',
        'PHENOTYPE_DATA': 'gs://bucket/phenotypes.csv'
    },
    out_params={
        'RESULTS': 'gs://bucket/gwas_results/',
        'LOG_FILE': 'gs://bucket/logs/gwas.log'
    }
)
```

### Job Monitoring
```python
from utilities.dsub_helpers import check_dsub_status, job_details

# Check recent jobs
check_dsub_status(age='3d')  # Last 3 days

# Get detailed information
job_details(job='my-job-id')

# Check full details for all recent jobs
check_dsub_status(age='1d', full=True)
```

### Job Management
```python
from utilities.dsub_helpers import cancel_job, get_job_logs

# Cancel a specific job
cancel_job('problematic-job-12345')

# Find job logs
get_job_logs('my-analysis-job')
```

### Convenience Functions
```python
from utilities.dsub_helpers import submit_standard_job, submit_highmem_job

# Quick standard job
submit_standard_job(
    label='quick-analysis',
    script='analysis.sh',
    envs={'SAMPLE_SIZE': '1000'}
)

# High-memory job for large datasets
submit_highmem_job(
    label='large-cohort-analysis',
    script='process_large_cohort.sh',
    in_params={'BIG_DATA': 'gs://bucket/huge_dataset.csv'}
)
```

## Machine Type Recommendations

### Compute Jobs (CPU-intensive)
```python
machine_type='c4-standard-8'  # 8 vCPUs, optimized for computation
```

### Memory Jobs (Large datasets)
```python
machine_type='n1-highmem-8'   # 8 vCPUs, 52 GB RAM
```

### Standard Jobs (Most analyses)
```python
machine_type='n1-standard-4'  # 4 vCPUs, 15 GB RAM
```

### GPU Jobs (Deep learning)
```python
machine_type='n1-standard-4'
# Add: --accelerator-type nvidia-tesla-t4 --accelerator-count 1
```

## Common Patterns

### GWAS Pipeline
```python
# Submit multiple related jobs
for ancestry in ['EUR', 'AFR', 'AMR', 'EAS']:
    dsub_script(
        label=f'gwas-{ancestry}',
        machine_type='n1-standard-8',
        script='gwas_pipeline.sh',
        envs={'ANCESTRY': ancestry},
        in_params={'GENO': f'gs://bucket/{ancestry}_genotypes.pgen'},
        out_params={'RESULTS': f'gs://bucket/results/{ancestry}/'}
    )
```

### Monitoring Pipeline
```python
# Check pipeline status
import time

def monitor_pipeline(job_prefix, check_interval=300):
    while True:
        print(f"Checking jobs with prefix: {job_prefix}")
        check_dsub_status(age='1d')
        print(f"Waiting {check_interval} seconds...")
        time.sleep(check_interval)
```

## Best Practices

1. **Use descriptive labels** - Makes job identification easier
2. **Set appropriate machine types** - Match compute needs to avoid waste
3. **Use preemptible instances** - Default is True for cost savings
4. **Monitor job logs** - Use `get_job_logs()` for debugging
5. **Clean up failed jobs** - Regular monitoring and cleanup of old jobs
6. **Test with small datasets** - Validate scripts before large runs

## Troubleshooting

### Common Issues:
- **Job fails immediately** - Check script permissions and file paths
- **Out of memory** - Use highmem machine type or optimize script
- **Slow performance** - Consider compute-optimized machines
- **Permission errors** - Verify service account permissions

### Debug Steps:
1. Check job status with `job_details()`
2. Examine logs using `get_job_logs()`
3. Test script locally with small data
4. Verify input/output paths are accessible

---

**Note**: All functions use current workspace environment variables and user credentials automatically.
