import sys
import re

#-------------------------------------------------------------------#
#                                                                   #
#   Expecting two arguments: input .fpp file and output .f90 file   #
#                                                                   #
#-------------------------------------------------------------------#
if len(sys.argv) != 3:
    print("Usage: python3 t-flows-preprocessor.py input.fpp output.f90")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

#----------------------------------------------------------------#
#                                                                #
#   Function to perform the transformation inside !$tfp blocks   #
#                                                                #
#----------------------------------------------------------------#
def process_tfp_block(block):
    # Add pointer declarations at the beginning of the block

    # Initialize pointers
    pointer_setup = ""

    # Initialize OpenACC clauses
    openacc_setup = (
        "      !$acc parallel loop          &\n"
        "      !$acc present(               &\n")

    #--------------------------------------------------#
    #   Variables for browsing through grid entities   #
    #--------------------------------------------------#

    if "Faces_In_Region" in block:
        pointer_setup = (
            "      grid_region_f_face => Grid % region % f_face\n"
            "      grid_region_l_face => Grid % region % l_face\n")
        openacc_setup += (
            "      !$acc   grid_region_f_face,  &\n"
            "      !$acc   grid_region_l_face,  &\n")
        block = re.sub(r'Grid % region % f_face', 'grid_region_f_face', block)
        block = re.sub(r'Grid % region % l_face', 'grid_region_l_face', block)

    #-------------------#
    #   Grid entities   #
    #-------------------#

    # Grid % faces_c
    if re.search(r'Grid % faces_c\b', block):
        pointer_setup += (
            "      grid_faces_c       => Grid % faces_c\n")
        openacc_setup += (
            "      !$acc   grid_faces_c,        &\n")
        block = re.sub(r'Grid % faces_c', 'grid_faces_c', block)

    # Grid % sx, Grid % sy, Grid % sz and Grid % s
    if re.search(r'Grid % sx\b', block):
        pointer_setup += (
            "      grid_sx            => Grid % sx\n")
        openacc_setup += (
            "      !$acc   grid_sx,             &\n")
    if re.search(r'Grid % sy\b', block):
        pointer_setup += (
            "      grid_sy            => Grid % sy\n")
        openacc_setup += (
            "      !$acc   grid_sy,             &\n")
    if re.search(r'Grid % sz\b', block):
        pointer_setup += (
            "      grid_sz            => Grid % sz\n")
        openacc_setup += (
            "      !$acc   grid_sz,             &\n")
    if re.search(r'Grid % s\b', block):
        pointer_setup += (
            "      grid_s             => Grid % s\n")
        openacc_setup += (
            "      !$acc   grid_s,              &\n")

    #-------------------------#
    #   Dependent variables   #
    #-------------------------#
    if re.search(r'Flow % v_flux % n\b', block):
        pointer_setup += (
            "      flow_v_flux_n      => Flow % v_flux % n\n")
        openacc_setup += (
            "      !$acc   flow_v_flux_n,       &\n")
        block = re.sub(r'Flow % v_flux % n', 'flow_v_flux_n', block)

    if re.search(r'b\b', block):
        pointer_setup += (
            "      b                  => Flow % Nat % b\n")
        openacc_setup += (
            "      !$acc   b,                   &\n")
        block = re.sub(r'b', 'b', block)  # no change to 'b' but kept for completeness

    openacc_setup += (
        "      !$acc )\n")

    # Replace the last comma with a space in the present clause
    last_comma_index = openacc_setup.rfind(",")
    if last_comma_index != -1:
        openacc_setup = openacc_setup[:last_comma_index] + " " + openacc_setup[last_comma_index + 1:]

    # Replace the 'do' loop with the OpenACC-parallelized version
    block = re.sub(
        r'do s = Faces_In_Region\(reg\)',
        'do s = grid_region_f_face(reg), grid_region_l_face(reg)',
        block
    )

    # Add the !$acc end parallel at the end
    block = block.replace("end do", "end do\n      !$acc end parallel")

    # Return the modified block with the pointer setup at the beginning
    return pointer_setup + openacc_setup + block

#-----------------------------------------------#
#                                               #
#   Looks as if the script really starts here   #
#    This is like the main function, I guess    #
#                                               #
#-----------------------------------------------#

#-------------------------#
#   Read the input file   #
#-------------------------#
with open(input_file, 'r') as infile:
    lines = infile.readlines()

# Process each line and look for !$tfp blocks
inside_tfp_block = False
output_lines = []

for line in lines:
    if line.strip().startswith("!$tfp begin"):
        inside_tfp_block = True
        tfp_block = []
    elif line.strip().startswith("!$tfp end"):
        inside_tfp_block = False
        # Process the whole tfp block and append it to output
        processed_block = process_tfp_block("".join(tfp_block))
        output_lines.append(processed_block)
    elif inside_tfp_block:
        tfp_block.append(line)
    else:
        output_lines.append(line)

# Write the transformed content to the output file
with open(output_file, 'w') as outfile:
    outfile.writelines(output_lines)

print(f"Processed {input_file} -> {output_file}")

