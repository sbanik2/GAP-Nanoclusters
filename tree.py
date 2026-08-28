from pathlib import Path
import sys

def print_tree(path, depth, prefix=""):
    if depth < 0:
        return

    dirs = sorted(p for p in Path(path).iterdir() if p.is_dir())

    for i, d in enumerate(dirs):
        last = i == len(dirs) - 1
        print(prefix + ("└── " if last else "├── ") + d.name)

        if depth > 0:
            print_tree(d, depth - 1, prefix + ("    " if last else "│   "))

depth = int(sys.argv[1]) if len(sys.argv) > 1 else 2
print_tree(".", depth)
