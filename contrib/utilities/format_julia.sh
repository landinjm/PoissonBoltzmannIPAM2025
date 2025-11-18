#!/bin/bash

# Check that this is run from the top-level directory
if [ ! -f "scripts/simplecell.jl" ]; then
    echo "Please run this script from the top-level directory of the repository."
    exit 1
fi

# Check that JuliaFormatter is installed
if ! julia -e "using Pkg; Pkg.status(\"JuliaFormatter\")" &> /dev/null; then
    echo "JuliaFormatter is not installed. Installing..."
    julia -e "using Pkg; Pkg.add(\"JuliaFormatter\")"
fi

# Find all .jl files and format them
find . -name "*.jl" -type f | while read -r file; do
    julia -e "using JuliaFormatter; format(\"$file\")"
done

echo "OK"
exit 0
