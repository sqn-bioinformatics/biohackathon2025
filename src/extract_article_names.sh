#!/bin/bash
# Extract article names from downloaded papers
# Portable script that works on any computer without hardcoded paths

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Set default paths relative to the script location
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
SOURCE_DIR="${SOURCE_DIR:-$PROJECT_ROOT/data/europepmc_articles}"
OUTPUT_FILE="${OUTPUT_FILE:-$PROJECT_ROOT/data/article_names.txt}"

# Check if source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory not found: $SOURCE_DIR"
    echo "Usage: SOURCE_DIR=/path/to/articles OUTPUT_FILE=/path/to/output.txt $0"
    exit 1
fi

# Count JSON files
json_count=$(find "$SOURCE_DIR" -maxdepth 1 -name "*_metadata.json" 2>/dev/null | wc -l)

if [ "$json_count" -eq 0 ]; then
    echo "Error: No metadata JSON files found in $SOURCE_DIR"
    exit 1
fi

echo "Extracting article names from $json_count articles..."
echo "Source: $SOURCE_DIR"
echo "Output: $OUTPUT_FILE"
echo ""

# Create output directory if needed
mkdir -p "$(dirname "$OUTPUT_FILE")"

# Clear output file and start fresh
> "$OUTPUT_FILE"

# Extract titles from JSON metadata files using jq if available, otherwise Python
if command -v jq &> /dev/null; then
    # Use jq for faster processing
    for json_file in "$SOURCE_DIR"/*_metadata.json; do
        if [ -f "$json_file" ]; then
            pmcid=$(basename "$json_file" _metadata.json)
            title=$(jq -r '.title // "No title"' "$json_file" 2>/dev/null || echo "No title")
            echo " $pmcid: $title" >> "$OUTPUT_FILE"
        fi
    done
else
    # Fallback to Python
    for json_file in "$SOURCE_DIR"/*_metadata.json; do
        if [ -f "$json_file" ]; then
            pmcid=$(basename "$json_file" _metadata.json)
            title=$(python3 -c "
import json
try:
    with open('$json_file', 'r') as f:
        data = json.load(f)
        print(data.get('title', 'No title'))
except Exception:
    print('No title')
" 2>/dev/null)
            echo " $pmcid: $title" >> "$OUTPUT_FILE"
        fi
    done
fi

# Count extracted entries
entry_count=$(wc -l < "$OUTPUT_FILE")

echo ""
echo "============================================"
echo "Article names extracted successfully!"
echo "============================================"
echo "Output file: $OUTPUT_FILE"
echo "Total entries: $entry_count"
echo "============================================"
