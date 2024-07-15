#!/bin/bash
set -e

# Initialize variables
PROJECT_NAME=""
CUTOFFS=""

# Display help message
function show_help() {
    echo "Usage: $0 -p PROJECT_NAME -c CUTOFFS"
    echo "  -p    Set the project name."
    echo "  -c    Set the cutoffs, separated by commas."
    echo "  -h    Display this help message."
    echo ""
    echo "Example:"
    echo "  bash workflow/scripts/run_bigscape_notebooks.sh -p Lactobacillus_delbrueckii -c \"0.40,0.50\""
}

# Parse command line options
while getopts ":hp:c:" opt; do
    case ${opt} in
        h )
            show_help
            exit 0
            ;;
        p )
            PROJECT_NAME=$OPTARG
            ;;
        c )
            CUTOFFS=$(echo $OPTARG | tr ',' ' ') # Replace commas with spaces
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            show_help
            exit 1
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            show_help
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Check if required options are provided
if [ -z "${PROJECT_NAME}" ] || [ -z "${CUTOFFS}" ]; then
    echo "Both project name and cutoffs are required."
    show_help
    exit 1
fi

# Your script continues here
DOCS_DIR="data/processed/$PROJECT_NAME/docs"

# Create environment file
ENV_NAME="bgcflow_notes"

if ! conda info --envs | grep -q "^${ENV_NAME}"; then
    echo "Creating Conda environment: ${ENV_NAME}"
    mamba env create -f workflow/envs/bgcflow_notes.yaml
else
    echo "Conda environment '${ENV_NAME}' already exists."
fi

for cutoff in $CUTOFFS
do
    NEW_NOTEBOOK="bigscape.$cutoff.ipynb"
    (cd $DOCS_DIR && cp bigscape.ipynb $NEW_NOTEBOOK)
    # Use double quotes for variable expansion and escape inner double quotes
    (cd $DOCS_DIR && sed -i 's/cutoff = \\"0.30\\"/cutoff = \\"'"$cutoff"'\\"/' $NEW_NOTEBOOK)
    (cd $DOCS_DIR && conda run -n bgcflow_notes jupyter nbconvert --execute $NEW_NOTEBOOK --to markdown)
done

echo "You can find the resulting BiG-SCAPE graphml network file in $DOCS_DIR/assets/data"
