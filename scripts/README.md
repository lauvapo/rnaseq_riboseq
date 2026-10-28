# Utility Scripts

Standalone scripts for specialized analyses not integrated into the main pipeline.

## Scripts

### `homoplasy_counter.py`
Phylogenetic homoplasy analysis script.

**Purpose**: Count the number of times each mutation appears on different branches (homoplasy index) and the number of leaves affected.

**Input files** (must be in the same directory):
- `branch_mutations.txt` — Tab-delimited file with columns: node, state1, pos, state2
- `nodes_child.tsv` — Tree structure with columns: parent, child

**Output**: `homoplasies_count.tsv`

**Usage**:
```bash
python scripts/homoplasy_counter.py
```

**Output format**:
```
mutation	total_appearances	homoplasy_index
G24965T	2	3
A45123C	1	1
...
```

- `total_appearances`: Number of branches where this mutation occurs
- `homoplasy_index`: Number of leaves descended from those branches

**Dependencies**: Python 3.7+, only built-in libraries

## Adding New Scripts

Place standalone analysis scripts here. Use clear naming (e.g., `descriptive_name.py`, `analysis_tool.R`) and include:
1. A header comment explaining purpose and usage
2. Input/output file specifications
3. Any non-standard dependencies

Keep scripts modular and self-contained so they can be run independently of the pipeline.
