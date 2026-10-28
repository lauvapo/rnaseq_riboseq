"""
Script for counting homoplasies in a phylogenetic tree based on branch_mutations.txt file.

The script reads branch_mutations.txt and counts:
- total_appearances: number of branches where each mutation appears
- homoplasy_index: number of leaves that have that mutation

Output format:
mutation	total_appearances	homoplasy_index
"""

from collections import defaultdict

branch_mutations_file = "branch_mutations.txt"
nodes_file = "nodes_child.tsv"
output_file = "homoplasies_count.tsv"

# Read branch mutations
mutations_by_branch = defaultdict(list)
all_mutations = set()

print("Reading branch mutations...")
with open(branch_mutations_file) as f:
    next(f)  # skip header
    for line in f:
        if line.strip():
            parts = line.rstrip().split('\t')
            node = parts[0]
            state1, pos, state2 = parts[1], parts[2], parts[3]
            mutation = f"{state1}{pos}{state2}"  # e.g., "G24965T"
            mutations_by_branch[node].append(mutation)
            all_mutations.add(mutation)

print(f"Found {len(all_mutations)} unique mutations in {len(mutations_by_branch)} branches")

# Build tree structure: node -> list of children
print("Building tree structure...")
children = defaultdict(list)
with open(nodes_file) as f:
    for line in f:
        if line.strip():
            parts = line.rstrip().split('\t')
            parent, child = parts[0], parts[1]
            children[parent].append(child)

# Identify all nodes and leaves
all_nodes = set(children.keys())
for child_list in children.values():
    all_nodes.update(child_list)

leaves = {node for node in all_nodes if node not in children}
print(f"Tree has {len(all_nodes)} nodes total, {len(leaves)} leaves")

# Get all descendants of a node (used to trace mutations to leaves)
def get_descendants(node, children_dict):
    """Recursively get all descendant nodes including the node itself."""
    descendants = {node}
    if node in children_dict:
        for child in children_dict[node]:
            descendants.update(get_descendants(child, children_dict))
    return descendants

# Count homoplasies for each mutation
print("Counting homoplasies...")
mutation_data = defaultdict(lambda: {"appearances": 0, "leaves": set()})

for branch, muts in mutations_by_branch.items():
    for mutation in muts:
        mutation_data[mutation]["appearances"] += 1
        # Find all leaves descended from this branch
        descendants = get_descendants(branch, children)
        leaf_descendants = descendants & leaves
        mutation_data[mutation]["leaves"].update(leaf_descendants)

# Write output
print(f"Writing results to {output_file}...")
with open(output_file, 'w') as f:
    f.write("mutation\ttotal_appearances\thomoplasy_index\n")
    for mutation in sorted(all_mutations):
        data = mutation_data[mutation]
        total_appearances = data["appearances"]
        homoplasy_index = len(data["leaves"])
        f.write(f"{mutation}\t{total_appearances}\t{homoplasy_index}\n")

print(f"Complete! Results written to {output_file}")
print(f"\nSummary statistics:")
print(f"  Total unique mutations: {len(all_mutations)}")
print(f"  Mutations appearing in multiple branches: {sum(1 for d in mutation_data.values() if d['appearances'] > 1)}")
print(f"  Max appearances of a single mutation: {max(d['appearances'] for d in mutation_data.values())}")
print(f"  Max leaves affected by a single mutation: {max(len(d['leaves']) for d in mutation_data.values())}")
