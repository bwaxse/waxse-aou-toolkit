# Resource Monitoring Template

Copy this section into your project's CLAUDE.md file to track machine usage and costs.

---

## Resource Usage & Machine Types

### Monitoring Code

Add this to notebooks to track actual resource usage:

```python
import psutil
import os
from datetime import datetime

def print_resource_usage(label=""):
    """Monitor memory and CPU usage with recommendations"""
    memory = psutil.virtual_memory()
    cpu_percent = psutil.cpu_percent(interval=1)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    print(f"\n{'='*60}")
    print(f"Resource Usage: {label}")
    print(f"Timestamp: {timestamp}")
    print(f"{'='*60}")
    print(f"Memory: {memory.used / 1e9:.1f}GB / {memory.total / 1e9:.1f}GB "
          f"({memory.percent:.1f}% used)")
    print(f"CPU: {cpu_percent:.1f}% used across {os.cpu_count()} cores")

    # Recommendations
    if memory.percent > 85:
        print("⚠️  HIGH MEMORY - Consider upgrading to highmem or optimizing")
    elif memory.percent < 40 and cpu_percent < 40:
        print("💰 LOW USAGE - Could downsize to save costs")
    else:
        print("✅ Good utilization")
    print(f"{'='*60}\n")

# Usage: Call after major operations
df = polars_gbq(query)
print_resource_usage("After loading diagnosis data")

merged = df1.join(df2, on='person_id')
print_resource_usage("After join operation")
```

---

### Observed Usage

**Cohort Size**: [e.g., 400K people, 50K cases, 350K controls]

#### Phase 1: Data Extraction & Initial Processing

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| BigQuery extraction (diagnoses) | | | | | | | |
| BigQuery extraction (medications) | | | | | | | |
| BigQuery extraction (labs) | | | | | | | |
| BigQuery extraction (vitals) | | | | | | | |

**Notes:**
- [Add observations about what worked well or caused issues]
- [Note any out-of-memory errors or performance bottlenecks]

#### Phase 2: Data Cleaning & Transformation

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| Unit conversion & outlier removal | | | | | | | |
| Temporal matching | | | | | | | |
| Multi-way joins | | | | | | | |

**Notes:**
- [Document memory spikes during joins]
- [Note if operations were CPU or memory bound]

#### Phase 3: Analysis & Visualization

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| Statistical summaries | | | | | | | |
| Table 1 generation | | | | | | | |
| Plotting | | | | | | | |

**Notes:**
- [Any visualization performance issues]

#### Phase 4: GWAS (if applicable)

| Operation | Machine Type | vCPU | RAM | Peak Memory | Peak CPU | Duration | Cost |
|-----------|--------------|------|-----|-------------|----------|----------|------|
| SAIGE null model (EUR) | | | | | | | |
| SAIGE null model (AFR) | | | | | | | |
| SAIGE chr tests (22 jobs) | | | | | | | |
| METAL meta-analysis | | | | | | | |

**Notes:**
- [Job failures, retries, preemptible interruptions]

---

### Total Project Cost

| Phase | Machine Hours | Est. Cost |
|-------|---------------|-----------|
| Phase 1: Extraction | | |
| Phase 2: Processing | | |
| Phase 3: Analysis | | |
| Phase 4: GWAS | | |
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

**SAIGE Jobs:**
- **Chose**: [e.g., n2-highmem-16 for null model]
- **Why**: [e.g., "Large sample size (400K) needs more RAM than standard"]
- **Result**: [e.g., "Completed successfully in 2.5 hours"]
- **Next time**: [e.g., "Good sizing, keep same"]

**Other Workloads:**
- [Document any other specialized jobs]

---

### Quick Reference: Cost per Machine Type

For reference (us-central1 pricing as of 2026):

| Machine Type | vCPU | RAM | Cost/hour | Good For |
|--------------|------|-----|-----------|----------|
| n2-standard-4 | 4 | 16GB | $0.19 | Light processing, small cohorts |
| n2-standard-8 | 8 | 32GB | $0.39 | SAIGE chr tests, meta-analysis |
| n2-standard-16 | 16 | 64GB | $0.77 | Medium processing, plotting |
| n2-standard-32 | 32 | 128GB | $1.55 | Large cohort notebooks |
| n2-standard-64 | 64 | 256GB | $3.10 | Very large cohorts, multiple datasets |
| n2-highmem-8 | 8 | 64GB | $0.52 | Memory-heavy, low CPU |
| n2-highmem-16 | 16 | 128GB | $1.04 | SAIGE null models (large cohorts) |
| n2-highmem-32 | 32 | 256GB | $2.08 | Very memory-intensive processing |
| n2-highcpu-4 | 4 | 8GB | $0.24 | CPU-bound, modest memory |
| n2-highcpu-8 | 8 | 16GB | $0.30 | Parallel computation |

*Preemptible instances: 60-80% cheaper, but can be interrupted*
