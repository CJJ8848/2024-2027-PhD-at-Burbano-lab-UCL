import re

# Input and output file paths
input_file = "lengthdiff.txt"
output_file = "length_calculations.txt"

# Open the output file for writing
with open(output_file, "w") as out_file:
    # Read each pair of lines in the input file
    with open(input_file, "r") as in_file:
        lines = in_file.readlines()
        
        # Iterate through lines in pairs
        for i in range(0, len(lines), 2):
            # Extract sample name and coordinates from pairs of lines
            sample_name = lines[i].strip()  # Sample name on the first line
            coords = lines[i+1].strip().lstrip(">")  # Coordinates on the next line, removing '>'
            
            # Extract start and end positions from the coordinates
            match = re.match(r"(\d+):(\d+)-(\d+)", coords)
            if match:
                start = int(match.group(2))
                end = int(match.group(3))
                
                # Calculate the absolute length difference
                length = abs(end - start)+1
                
                # Write the result to the output file
                out_file.write(f"{sample_name}\t{length}\n")
            else:
                print(f"Coordinate format issue: {coords}")

print(f"Calculation complete. Results saved in {output_file}.")
