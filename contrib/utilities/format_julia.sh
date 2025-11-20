#!/bin/bash

# Check that this is run from the top-level directory
if [ ! -f "scripts/simplecell.jl" ]; then
    echo "Please run this script from the top-level directory of the repository."
    exit 1
fi

# Find all .jl files and format them
find . -name "*.jl" -type f | while read -r file; do
    julia -e "using Runic; Runic.format_file(\"$file\", inplace = true)"
done

echo "OK"
exit 0
