# this file needs to be sourced

# Pick package manager: prefer mamba, fallback to conda
if command -v mamba >/dev/null 2>&1; then
   package_manager="mamba"
elif command -v conda >/dev/null 2>&1; then
   package_manager="conda"
else
   echo "Error: neither mamba nor conda found in PATH."
   return 1
fi
echo "Using package manager: $package_manager"

# Initialize the package manager to use the current shell
eval "$($package_manager shell hook --shell bash)"

# Create local environment if it doesn't exist
if ! test -d condaenv ; then
   $package_manager create -y -p ./condaenv pip
fi

# Activate conda environment   
$package_manager activate condaenv/

# Install juliacall & numpy
pip install -r requirements.txt

# Force juliacall to use standard julia executable
export PYTHON_JULIACALL_EXE=`which julia`

# Force juliacall to use environment in this directory
export PYTHON_JULIACALL_PROJECT=`pwd`

# Set up julia environment
julia --project --startup-file=no -e "using Pkg; Pkg.resolve(); using JuliaMPBSolver"

# Successfully return
echo "OK"
return 0
