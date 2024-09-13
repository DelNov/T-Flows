#==============================================================================#
#   This script is designed to automatically insert OpenACC directives into
#   T-Flows' Process_Gpu code for the blocks marked with !$tf-acc comments.
#
#   It processes the input T-Flows code, detects specific array structures,
#   handles them according to predefined rules, and inserts OpenACC directives
#   for parallel processing.
#
#   It relies on the assumption that there is only a finite number of T-Flows'
#   loop constructes, and that these can be automated with a script like this.
#   This is not a procedure which would work for any Fortran code, but only
#   to T-Flows, and even for T-Flows it relies on its strict coding standards.
#
#   Breakdown of the Script's Workflow
#
#   1. Argument Handling: The script expects two command-line arguments:
#      - The input .fpp file containing Fortran code blocks.
#      - The output .f90 file where the transformed code will be written.
#
#      If incorrect arguments are provided, the script will exit with error.
#
#   2. Initialization: The script initializes some helper variables:
#
#      - indent: Sets the indentation level to three spaces for aligning
#        inserted OpenACC directives.
#      - excluded_macros: Lists macros (sych as Faces_In_Region,
#        Cells_In_Domain and many others) that will not be processed as arrays
#        (since they are treated differently).
#
#   3. Main Logic - Process_Tfp_Block(block) Function: This function processes
#      each !$tf-acc code block in the input file. Below is a step-by-step
#      explanation of its functionality:
#
#      a. Initialize OpenACC Pointers and Clauses:
#         - The script initializes the pointer_setup and openacc_setup
#           strings, which are used to store OpenACC directives.
#         - For each !$tf-acc block, it builds a list of variables that need
#           to be "present" in the OpenACC data region. This list will be
#           inserted into the OpenACC directives.
#
#      b. Handle Special Macros:
#         - Faces_In_Region: If the block contains Faces_In_Region, the script
#           assigns pointers to the corresponding grid structures (such as
#           grid_region_f_face and grid_region_l_face) and adds them to the
#           OpenACC present clause.
#         - Cells_In_Domain: Similarly, if Cells_In_Domain is detected,
#           pointers are created and the variables are added to the present
#           clause.
#
#      c. Array Handling:
#         - The script uses a regular expression to find arrays in the block
#           (variables followed by parentheses) and creates pointer names for
#           each array.
#         - If the array is part of a structure (Grid % cells_c, for example),
#           a pointer setup line is created to alias the array, and the array
#           is added to the OpenACC present clause.
#         - The script replaces occurrences of the full array name with its
#           pointer name within the block.
#
#      d. Detect Reduction Variables:
#         - The script checks for reduction operations (variables that
#           accumulate values in a loop, such as var = var + something) using
#           a regular expression.
#         - If reduction variables are detected, it adds a reduction clause to
#           the OpenACC directives.
#
#      e. Finalizing OpenACC Directives:
#         - After processing the arrays and reduction variables, the
#           openacc_setup string is completed by closing the present clause.
#
#      f. Insert OpenACC Parallel Directives:
#         - The script looks for "do" loops in the block:
#           > It finds the second "do" loop and inserts an OpenACC sequential
#             loop directive (!$acc loop seq) before it.
#           > The script also inserts "!$acc end parallel" after the last
#             "end do" statement and "!$acc end loop" after the second-to-last
#             "end do".
#
#   4. File Processing:
#
#      - The script reads the input file line by line, looking for "!$tf-acc"
#        blocks.  When it finds the start of such a block ("!$tf-acc loop
#        begin"), it begins accumulating lines into a buffer (tfp_block).
#      - When the end of the block ("!$tf-acc loop end") is found, the script
#        processes the entire block using the Process_Tfp_Block function,
#        transforming it into a new block with OpenACC directives.
#      - The transformed block is written to the output file.
#
#   5. Writing the Output:
#
#      - Once all "!$tf-acc" blocks are processed, the script writes the
#        modified code into the output .f90 file.
#------------------------------------------------------------------------------#
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

indent  = "   "
RED     = "\033[31m"
GREEN   = "\033[32m"
YELLOW  = "\033[33m"
BLUE    = "\033[34m"
MAGENTA = "\033[35m"
CYAN    = "\033[36m"
RESET   = "\033[0m"

copy_to_device_file = "./Gpu_Mod/Grid/Copy_To_Device.f90"

import re

#==============================================================================#
#                                                                              #
#   Check if a string is in the file                                           #
#                                                                              #
#------------------------------------------------------------------------------#
def Command_In_File(file_path, target_string):

  # Remove spaces from the target string for comparison
  cleaned_target = re.sub(r'\s+', '', target_string)

  # Open the file and read its content
  with open(file_path, 'r') as file:
    file_content = file.read()

  # Remove spaces from the file content
  cleaned_file_content = re.sub(r'\s+', '', file_content)

  # Check if the cleaned target string is in the cleaned file content
  return cleaned_target in cleaned_file_content

#==============================================================================#
#                                                                              #
#   Function to find all arrays in a block                                     #
#                                                                              #
#------------------------------------------------------------------------------#
def Find_Arrays_In_Block(block):

  #------------------------------------------------#
  #   Remove standard Fortran control structures   #
  #------------------------------------------------#

  # Remove all comments
  cleaned_block = re.sub(r'!.*', '', block)
  print(cleaned_block)

  # Replace every occurrence of spaces followed by "%" with "@%"
  # and every occurence of "%" followed by spacess with "%@"
  cleaned_block = re.sub(r'\s+%', '@%', cleaned_block)
  cleaned_block = re.sub(r'%\s+', '%@', cleaned_block)

  # Save blocks for searching, it will be used later on
  block_for_searching = re.sub(r'\s+%', ' %', block)
  block_for_searching = re.sub(r'%\s+', '% ', block_for_searching)

  # Remove "do"
  # (This will remove "do" also from "end do" statements)
  cleaned_block = re.sub(r'\bdo\b', '', cleaned_block)

  # Remove "enddo" without spaces
  # (Any "end do" is by now only "end", and "end" is treaed in the end)
  cleaned_block = re.sub(r'\benddo\b', '', cleaned_block)

  # Remove "if" followed by opening parentheses with or without spaces
  # (Unlike "do" above, this will leave "end if" combinations around)
  cleaned_block = re.sub(r'\bif\s*\(', '', cleaned_block)

  # Remove "endif" with or without spaces in between
  # (Since the command above left "end if" combinations around)
  cleaned_block = re.sub(r'\bend\s*if\b', '', cleaned_block)

  # Remove "then" preceded by a closing parenthesis with or without spaces
  cleaned_block = re.sub(r'\)\s*then\b', '', cleaned_block)

  # Remove "elseif"
  # (If there was any "else if", the "if" part was already removed above)
  cleaned_block = re.sub(r'\belseif\b', '', cleaned_block)

  # Remove "else"
  cleaned_block = re.sub(r'\belse\b', '', cleaned_block)

  # Remove "end"
  cleaned_block = re.sub(r'\bend\b', '', cleaned_block)

  # Remove "abs"
  cleaned_block = re.sub(r'\babs\b', '', cleaned_block)

  # Remove "max"
  cleaned_block = re.sub(r'\bmax\b', '', cleaned_block)

  # Remove "min"
  cleaned_block = re.sub(r'\bmin\b', '', cleaned_block)

  # Remove "maxval"
  cleaned_block = re.sub(r'\bmaxval\b', '', cleaned_block)

  # Remove "minval"
  cleaned_block = re.sub(r'\bminval\b', '', cleaned_block)

  # Remove "merge"
  cleaned_block = re.sub(r'\bmerge\b', '', cleaned_block)

  # Remove empty parentheses
  cleaned_block = re.sub(r'\(\s*\)', '', cleaned_block)

  # Remove "\n"
  cleaned_block = re.sub(r"\n", "", cleaned_block)

  # Remove "&"
  cleaned_block = re.sub(r"&", "", cleaned_block)

  # Put spaces around mathematical operators
  cleaned_block = re.sub(r"(?<! )\+(?! )", " + ", cleaned_block)
  cleaned_block = re.sub(r"(?<! )-(?! )", " - ", cleaned_block)
  cleaned_block = re.sub(r"(?<! )\*(?! )", " * ", cleaned_block)
  cleaned_block = re.sub(r"(?<! )/(?! )", " / ", cleaned_block)

  result = set()
  seen_arrays = set()  # Track arrays to avoid duplicates

  #----------------------------------------------------------#
  #   Browse through the reduced block and look for arrays   #
  #----------------------------------------------------------#
  for i, char in enumerate(cleaned_block):
    if char == '(':
      # The array name is just before the opening brace
      array_name_end = i
      # Go backwards to find where the array name starts
      j = i - 1
      while j >= 0 and (cleaned_block[j] not in [' '] 
                        and
                        cleaned_block[j] not in ['(']):
        j -= 1
      array_name_start = j + 1
      array_name = cleaned_block[array_name_start:array_name_end]

      # Only consider valid array names (non-empty) and avoid duplicates
      if array_name and array_name not in seen_arrays:
        seen_arrays.add(array_name)  # mark array as seen
        array_name = re.sub(r'@', ' ', array_name)
        result.add(array_name)

  #---------------------------------------------------------#
  #   At this point, we have found all arrays, now search   #
  #   through the original block to find their positions    #
  #---------------------------------------------------------#

  # Convert the result set to a list for further processing
  array_list = list(result)

  # Find the position of each array in the original block
  array_positions = []
  for array in array_list:
    position = block_for_searching.find(array)
    array_positions.append((array, position))

  # Sort the arrays based on their positions
  sorted_array_positions = sorted(array_positions, key=lambda x: x[1])

  # Return only the array names in the order of their appearance
  return [array[0] for array in sorted_array_positions]

#==============================================================================#
#                                                                              #
#   Function to perform the transformation inside !$tf-acc blocks              #
#                                                                              #
#------------------------------------------------------------------------------#
def Process_Tfp_Block(block):

  print("")
  print(f"{YELLOW}#======================{RESET}")
  print(f"{YELLOW}#                      {RESET}")
  print(f"{YELLOW}# Preprocessing block: {RESET}")
  print(f"{YELLOW}#                      {RESET}")
  print(f"{YELLOW}#----------------------{RESET}")
  print(block, end="")

  # Initialize pointers and OpenACC clauses
  pointer_setup = ""
  openacc_setup = ""
  if "Cells_In_Domain_And_Buffers" in block:
    openacc_setup += (indent + "!$acc parallel loop independent  &\n")
  else:
    if "Cells_In_Domain" in block:
      openacc_setup += (indent + "!$acc parallel loop  &\n")
  if "Faces_In_Region" in block:
    openacc_setup += (indent + "!$acc parallel loop  &\n")
  if "Faces_In_Domain_And_At_Buffers" in block:
    openacc_setup += (indent + "!$acc parallel loop  &\n")

  openacc_setup += (indent + "!$acc present(  &\n")

  # Set to track processed variables to avoid duplicates
  processed_vars = set()

  # Set to track reduction variables
  reduction_vars = set()

  # List of macros to exclude from processing as arrays
  excluded_macros = {"Faces_In_Region",
                     "Faces_In_Domain_And_At_Buffers",
                     "Face_Value",
                     "Cells_In_Domain",
                     "Cells_In_Domain_And_Buffers"}

  #----------------------------------------------------#
  #   Special handling for Faces_In_Region variables   #
  #----------------------------------------------------#
  if "Faces_In_Region" in block:

    print("")
    print(f"{CYAN}  #-------------------------------{RESET}")
    print(f"{CYAN}  # Block of type Faces_In_Region {RESET}")
    print(f"{CYAN}  #-------------------------------{RESET}")

    print("")
    print(f"{RED}  # Pointers used in the block{RESET}")

    commands = ("grid_region_f_face => Grid % region % f_face",
                "grid_region_l_face => Grid % region % l_face")

    for command in commands:
      if not Command_In_File(copy_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", copy_to_device_file, sep="")

    openacc_setup += (
      indent + "!$acc   grid_region_f_face,  &\n" +
      indent + "!$acc   grid_region_l_face,  &\n")
    block = re.sub(r'Grid % region % f_face', 'grid_region_f_face', block)
    block = re.sub(r'Grid % region % l_face', 'grid_region_l_face', block)

    # Add these to the processed set to avoid duplicates
    processed_vars.add('grid_region_f_face')
    processed_vars.add('grid_region_l_face')

  #-------------------------------------------------------------------#
  #   Special handling for Faces_In_Domain_And_At_Buffers variables   #
  #-------------------------------------------------------------------#
  if "Faces_In_Domain_And_At_Buffers" in block:

    print("")
    print(f"{RED}  # Pointers used in the block{RESET}")

    commands = ("grid_region_f_face => Grid % region % f_face",
                "grid_region_l_face => Grid % region % l_face")

    for command in commands:
      if not Command_In_File(copy_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", copy_to_device_file, sep="")

    openacc_setup += (
      indent + "!$acc   grid_region_f_cell,  &\n" +
      indent + "!$acc   grid_region_l_cell,  &\n")
    block = re.sub(r'Grid % n_regions',       'grid_n_regions',     block)
    block = re.sub(r'Grid % region % f_cell', 'grid_region_f_cell', block)
    block = re.sub(r'Grid % region % l_cell', 'grid_region_l_cell', block)

  #--------------------------------------------------#
  #   Special handling for Cells_In_Domain variables #
  #--------------------------------------------------#
  if "Cells_In_Domain" in block and not "Cells_In_Domain_And_Buffers" in block:

    print("")
    print(f"{CYAN}  #-------------------------------{RESET}")
    print(f"{CYAN}  # Block of type Cells_In_Domain {RESET}")
    print(f"{CYAN}  #-------------------------------{RESET}")

    print("")
    print(f"{RED}  # Pointers used in the block{RESET}")

    commands = ("grid_region_f_cell => Grid % region % f_cell",
                "grid_region_l_cell => Grid % region % l_cell")

    for command in commands:
      if not Command_In_File(copy_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", copy_to_device_file, sep="")

    openacc_setup += (
      indent + "!$acc   grid_region_f_cell,  &\n" +
      indent + "!$acc   grid_region_l_cell,  &\n")
    block = re.sub(r'Grid % region % f_cell', 'grid_region_f_cell', block)
    block = re.sub(r'Grid % region % l_cell', 'grid_region_l_cell', block)

    # Add these to the processed set to avoid duplicates
    processed_vars.add('grid_region_f_cell')
    processed_vars.add('grid_region_l_cell')

  #----------------------------------------------------------------#
  #   Special handling for Cells_In_Domain_And_Buffers variables   #
  #----------------------------------------------------------------#
  if "Cells_In_Domain_And_Buffers" in block:

    print("")
    print(f"{CYAN}  #-------------------------------------------{RESET}")
    print(f"{CYAN}  # Block of type Cells_In_Domain_And_Buffers {RESET}")
    print(f"{CYAN}  #-------------------------------------------{RESET}")

    print("")
    print(f"{RED}  # Pointers used in the block{RESET}")

    commands = ("grid_region_f_cell => Grid % region % f_cell",
                "grid_region_l_cell => Grid % region % l_cell")

    for command in commands:
      if not Command_In_File(copy_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", copy_to_device_file, sep="")

    openacc_setup += (
      indent + "!$acc   grid_region_f_cell,  &\n" +
      indent + "!$acc   grid_region_l_cell,  &\n")
    block = re.sub(r'Grid % n_regions',       'grid_n_regions',     block)
    block = re.sub(r'Grid % region % f_cell', 'grid_region_f_cell', block)
    block = re.sub(r'Grid % region % l_cell', 'grid_region_l_cell', block)

    # Add these to the processed set to avoid duplicates
    processed_vars.add('grid_region_f_cell')
    processed_vars.add('grid_region_l_cell')

  #-----------------------------------------------------#
  #   General handling for arrays in the remaining code #
  #-----------------------------------------------------#
  # Regular expression to detect arrays (anything followed by parentheses)
  array_pattern = re.compile(r'(\w+(?: % \w+)*)\s*\(([^)]*)\)')

  #----------------------------------#
  #   Find all arrays in the block   #
  #----------------------------------#
  arrays = Find_Arrays_In_Block(block)

  print("")
  print("  # Arrays found in the block ")

  #------------------------#
  #   Process each array   #
  #------------------------#
  for array in arrays:

    full_name = array

    # Create pointer variable name
    variable_name = full_name.replace(" % ", "_").lower()

    # Skip excluded macros (like Faces_In_Region)
    if any(macro in full_name for macro in excluded_macros):
      continue

    # Check if variable has already been processed
    if variable_name not in processed_vars:
      if "%" in full_name:  # If it's part of a structure
        # Add to the pointer setup
        pointer_setup += f"{indent}{variable_name} => {full_name}\n"
        print("  ", variable_name, "  (", full_name, ")", sep="")

      else:
        print("  ", variable_name, sep="")

      # Add to the OpenACC setup
      openacc_setup += f"{indent}!$acc   {variable_name},  &\n"

      # Replace occurrences of the original variable
      # name with the pointer name in the block
      block = re.sub(re.escape(full_name), variable_name, block)

      # Mark variable as processed
      processed_vars.add(variable_name)

  #----------------------------------------------------#
  #   Detect reduction variables within the do-loop    #
  #----------------------------------------------------#

  # Pattern to detect reduction operations: var = var + something
  # or var = var - something.  The regex ensures that the variable
  # being reduced is a scalar (no parentheses for indices)
  reduction_pattern = re.compile(r'(\b\w+\b)\s*=\s*\1\s*[\+\-]')

  # Find reduction variables in the block (only scalars, no arrays)
  reductions = reduction_pattern.findall(block)
  for var in reductions:
    reduction_vars.add(var)

  # Add reduction clause if any reductions are found
  if reduction_vars:
    reduction_clause = "reduction(+:"
    reduction_clause += ",".join(reduction_vars)
    reduction_clause += ")"
    openacc_setup = openacc_setup.replace(
      "parallel loop",
      f"parallel loop {reduction_clause}"
    )

  # Append the closing parenthesis for OpenACC present clause
  openacc_setup += indent + "!$acc )\n"

  # Replace the last comma with a space in the present clause
  last_comma_index = openacc_setup.rfind(",")
  if last_comma_index != -1:
    openacc_setup = (
        openacc_setup[:last_comma_index]
      + " "
      + openacc_setup[last_comma_index + 1:]
    )

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Faces_In_Region\(reg\)',
    'grid_region_f_face(reg), grid_region_l_face(reg)',
    block
  )

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Faces_In_Domain_And_At_Buffers\(\)',
    'grid_region_f_face(grid_n_regions), grid_region_l_face(grid_n_regions)',
    block
  )

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Cells_In_Domain\(\)',
    'grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)',
    block
  )

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Cells_In_Domain_And_Buffers\(\)',
    'grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions+1)',
    block
  )

  # Insert OpenACC directives before the second "do" loop
  first_do_index = block.find("do ")
  if first_do_index != -1:
    second_do_index = block.find("do ", first_do_index + 1)
    if second_do_index != -1:

      # Find the last end of line before the second do statement
      last_eol_index = block.rfind("\n",  0, second_do_index)

      block = (
          block[:last_eol_index + 1]
        + indent
        + "!$acc loop seq\n"
        + block[last_eol_index + 1:]
      )

  # Add the !$acc end parallel at the end
  last_end_do_index = block.rfind("end do")
  if last_end_do_index != -1:
    block = (
        block[:last_end_do_index]
      + "end do\n"
      + indent
      + "!$acc end parallel"
      + block[last_end_do_index + len("end do"):]
    )

  # Now, search for the second-to-last end do before the first one
  second_last_end_do_index = block[:last_end_do_index].rfind("end do")
  if second_last_end_do_index != -1:
    block = (
        block[:second_last_end_do_index]
      + "end do\n" + indent + "!$acc end loop"
      + block[second_last_end_do_index + len("end do"):]
    )

  # Return the modified block with the pointer setup at the beginning
  return pointer_setup + openacc_setup + block

#==============================================================================#
#                                                                              #
#   The script really starts here                                              #
#                                                                              #
#------------------------------------------------------------------------------#

#-------------------------------#
#                               #
#   Read the whole input file   #
#                               #
#-------------------------------#
with open(input_file, 'r') as infile:
  lines = infile.readlines()

# Process each line and look for !$tf-acc blocks
inside_tfp_block = False
output_lines = []
line_count = 0

#----------------------------------------#
#                                        #
#   Browse through all the input lines   #
#                                        #
#----------------------------------------#
for line in lines:

  line_count = line_count + 1

  #----------------------------------------#
  #  If you find the !$tf-acc loop begin   #
  #----------------------------------------#
  if line.strip().startswith("!$tf-acc loop begin"):

    # Adjust indentation
    i = line.find("!")
    indent = line[:i]

    inside_tfp_block = True
    tfp_block = []

  #--------------------------------------#
  #  If you find the !$tf-acc loop end   #
  #--------------------------------------#
  elif line.strip().startswith("!$tf-acc loop end"):

    # Adjust indentation
    i = line.find("!")
    indent = line[:i]

    inside_tfp_block = False

    # Process the whole tfp block and append it to output
    processed_block = Process_Tfp_Block("".join(tfp_block))
    output_lines.append(processed_block)

  #------------------------------------------#
  #  If you find an invalid !$tf-acc entry   #
  #------------------------------------------#
  elif line.strip().startswith("!$tf-acc"):
    print(f"{RED}#============================#{RESET}")
    print(f"{RED}#                            #{RESET}")
    print(f"{RED}#   ERROR IN PREPROCESSING   #{RESET}")
    print(f"{RED}#                            #{RESET}")
    print(f"{RED}#============================#{RESET}")
    print("Invalid !$tf-acc directive on line", line_count, ":", line.strip())
    print("Valid entries are: ")
    print(f"{GREEN}   !$tf-acc loop begin{RESET}")
    print(f"{GREEN}   !$tf-acc loop end{RESET}")
    sys.exit(1)

  elif inside_tfp_block:
    tfp_block.append(line)

  else:
    output_lines.append(line)

#------------------------------------------------------#
#                                                      #
#   Write the transformed content to the output file   #
#                                                      #
#------------------------------------------------------#
with open(output_file, 'w') as outfile:
  outfile.writelines(output_lines)

print(f"Processed {input_file} -> {output_file}")

