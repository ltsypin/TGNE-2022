#!/bin/bash

check_line_counts() {
    local pattern="$1"

    # Use find to expand the wildcard pattern and get the line counts for each file
    local line_counts=$(find $pattern -type f -exec wc -l {} + | sed '$d' | awk '{print $1}')
    local sorted_line_counts=$(echo "$line_counts" | sort -n)
    local frequencies=$(echo "$sorted_line_counts" | uniq -c | sort -n -r)
    local unique_counts=$(echo "$frequencies" | awk '{print $2}')
    local sorted_frequencies=$(echo "$frequencies" | awk '{print $1}')

    local num_files=$(echo "$line_counts" | wc -l | tr -d '[:space:]')

    if [ "$(echo "$unique_counts" | wc -l)" -eq 1 ]; then
        echo "All $num_files files have the same line count: $unique_counts"
    else
        echo "Different line counts and their frequencies for $num_files files:"
        (echo "Line_Count Frequency"; paste <(echo "$unique_counts") <(echo "$sorted_frequencies")) | column -t
    fi
}

# Example usage: Check line counts for the specified wildcard pattern
check_line_counts "$1"

