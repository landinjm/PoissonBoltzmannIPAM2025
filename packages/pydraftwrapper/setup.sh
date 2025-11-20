# this file needs to be sourced

# Create local conda environment
if ! test -d condaenv ; then
   conda create -p condaenv pip
fi

# Activate conda environment   
conda activate condaenv/

# install juliacall & numpy
pip install -r requirements.txt

# force juliacall to use standard julia executable
export PYTHON_JULIACALL_EXE=`which julia`

# force juliacall to use environment in this directory
export PYTHON_JULIACALL_PROJECT=`pwd`

# Set up julia environment
julia --project --startup-file=no -e "using Pkg; Pkg.resolve(); using JuliaDraftSolver"
