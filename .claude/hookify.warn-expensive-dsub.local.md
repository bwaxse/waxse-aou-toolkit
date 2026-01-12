---
name: warn-expensive-dsub
enabled: true
event: bash
pattern: dsub.*--tasks\s+\d{3,}
action: warn
---

💰 **Large-scale dsub job detected**

You're about to submit 100+ tasks via dsub.

**Cost Check**:
- How many tasks? (Check --tasks parameter)
- Instance type? (n1-standard-1 vs n1-highmem-8)
- Expected runtime per task?
- Approximate cost = tasks × hours × instance_cost

**Best Practice**: Test with --tasks 1 first, then scale up.

**Example costs**:
- 100 tasks × 2 hours × n1-standard-1 = ~$1
- 1000 tasks × 4 hours × n1-highmem-8 = ~$160

If this is intentional, you can proceed. This is just a reminder!
