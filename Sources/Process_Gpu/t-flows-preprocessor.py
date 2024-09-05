import sys

# Expecting two arguments: input .fpp file and output .f90 file
if len(sys.argv) != 3:
    print("Usage: python3 t-flows-pre.py input.fpp output.f90")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Basic example: just copy the content of .fpp to .f90 (you can modify this)
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    content = infile.read()
    outfile.write(content)

print(f"Processed {input_file} -> {output_file}")

