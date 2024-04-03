#!/bin/bash

# Loop through files matching the pattern
for file in *.labels; do
    # Check if the file exists and is a regular file
    if [ -f "$file" ]; then
        # Run the mcxload command for each file
        mcxload -etc-ai "$file" -strict-tabr abc_format_graph.tab -o "out.${file%.labels}.mci"
    fi
done

