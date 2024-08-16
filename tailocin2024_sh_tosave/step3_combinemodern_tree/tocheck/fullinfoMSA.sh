#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

INPUT_FILE=$1
OUTPUT_FILE=$2

# Read the input file and process it to remove columns with 'N' or '-'
awk '
BEGIN {
    # Sequence ID pattern
    seq_id = "^>"
    seq_count = 0
}
{
    if ($0 ~ seq_id) {
        # If not the first sequence, print the filtered sequence
        if (seq_count > 0) {
            sequences[seq_count] = current_sequence
            current_sequence = ""
        }
        seq_count++
        sequence_ids[seq_count] = $0
    } else {
        # Concatenate sequence lines
        current_sequence = current_sequence $0
    }
}
END {
    sequences[seq_count] = current_sequence
    seq_length = length(sequences[1])

    # Identify columns to keep
    for (i = 1; i <= seq_length; i++) {
        keep_column = 1
        for (j = 1; j <= seq_count; j++) {
            char = substr(sequences[j], i, 1)
            if (char == "N" || char == "-") {
                keep_column = 0
                break
            }
        }
        if (keep_column) {
            filtered_columns[i] = 1
        }
    }

    # Print the filtered sequences
    for (j = 1; j <= seq_count; j++) {
        print sequence_ids[j]
        for (i = 1; i <= seq_length; i++) {
            if (i in filtered_columns) {
                printf "%s", substr(sequences[j], i, 1)
            }
        }
        printf "\n"
    }
}
' "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Filtered alignment saved to $OUTPUT_FILE"
