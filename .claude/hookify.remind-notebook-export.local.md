---
name: remind-notebook-export
enabled: true
event: stop
action: warn
---

📓 **Notebook export reminder**

You may have worked with .ipynb files this session.

**Before finishing, ensure**:
- [ ] Notebooks are saved in Jupyter
- [ ] Exported to .py files: `jupyter nbconvert --to python notebook.ipynb`
- [ ] Both .ipynb and .py versions are current

**Why this matters**:
Your project uses a dual-file strategy (see README.md):
- `.ipynb` for interactive development and visualization
- `.py` for version control, code review, and imports

**Quick export command**:
```bash
# Single file
jupyter nbconvert --to python "path/to/notebook.ipynb"

# All notebooks in a directory
jupyter nbconvert --to python path/to/*.ipynb
```

**Verification**:
```bash
# Check which .py files are older than their .ipynb counterparts
find . -name "*.ipynb" -newer "*.py"
```

This is just a reminder - proceed when ready!
