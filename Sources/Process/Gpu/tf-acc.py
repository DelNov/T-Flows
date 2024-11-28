import sys
import re

#-------------------------------------------------------------------#
#                                                                   #
#   Expecting two arguments: input .fpp file and output .f90 file   #
#                                                                   #
#-------------------------------------------------------------------#
if len(sys.argv) != 4:
  print("Usage: python3 tf-acc.py [openacc|openmp] input.fpp output.f90")
  sys.exit(1)

backend     = sys.argv[1]
input_file  = sys.argv[2]
output_file = sys.argv[3]

if backend not in ['openmp', 'openacc']:
  print("Error: The first argument must be 'openmp' or 'openacc'.")
  sys.exit(1)

indent  = "   "

# Regular Colors
BLACK   = "\033[30m"
RED     = "\033[31m"
GREEN   = "\033[32m"
YELLOW  = "\033[33m"
BLUE    = "\033[34m"
MAGENTA = "\033[35m"
CYAN    = "\033[36m"
WHITE   = "\033[37m"

# Bright Colors
BRIGHT_RED     = "\033[91m"
BRIGHT_GREEN   = "\033[92m"
BRIGHT_YELLOW  = "\033[93m"
BRIGHT_BLUE    = "\033[94m"
BRIGHT_MAGENTA = "\033[95m"
BRIGHT_CYAN    = "\033[96m"
BRIGHT_WHITE   = "\033[97m"

# Text Styles
BOLD      = "\033[1m"
UNDERLINE = "\033[4m"
REVERSED  = "\033[7m"

RESET = "\033[0m"

grid_to_device_file = "./Gpu_Mod/Grid/Copy_To_Device.f90"
flow_to_device_file = "./Gpu_Mod/Field/Copy_To_Device.f90"
turb_to_device_file = "./Gpu_Mod/Turb/Copy_To_Device.f90"
gpu_pointers_file   = "./Gpu_Pointers_Mod.f90"

begin_parallel_loop_independent = "!$acc parallel loop independent  &"
begin_parallel_loop             = "!$acc parallel loop  &"
end_parallel_loop               = "!$acc end parallel"
begin_sequential_loop           = "!$acc loop seq"
end_sequential_loop             = "!$acc end loop"
begin_present                   = "!$acc present(  &"
end_present                     = "!$acc )"
begin_private                   = ""
end_private                     = ""

if backend == "openmp":
  begin_parallel_loop_independent = "!$omp parallel do  &"
  begin_parallel_loop             = "!$omp parallel do  &"
  end_parallel_loop               = "!$omp end parallel do"
  begin_sequential_loop           = ""
  end_sequential_loop             = ""
  begin_present                   = ""
  end_present                     = ""
  begin_private                   = "!$omp private(  &"
  end_private                     = "!$omp )"

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

  # Remove comments from the file content
  cleaned_file_content = re.sub(r'!.*', '', file_content)

  # Remove spaces from the file content
  cleaned_file_content = re.sub(r'\s+', '', cleaned_file_content)

  # Check if the cleaned target string is in the cleaned file content
  return cleaned_target in cleaned_file_content

#==============================================================================#
#                                                                              #
#   Check if a block has inner loop                                            #
#                                                                              #
#------------------------------------------------------------------------------#
def Has_Inner_Loops(block):

  # Find all positions of standalone "do" statements
  do_positions = [match.start()
                  for match in re.finditer(r'^\s*do\b', block, re.MULTILINE)]

  # Find all positions of "end do" (allowing
  # for variations like "enddo" and spaces)
  end_do_positions = [match.end()
                      for match in re.finditer(r'\bend\s*do\b', block)]

  # If there are at least 2 "do" statements and
  # at least 1 "end do," we have an inner loop
  return len(do_positions) >= 2 and len(end_do_positions) > 0

#==============================================================================#
#                                                                              #
#   Check if a block has inner loop                                            #
#                                                                              #
#------------------------------------------------------------------------------#
def Remove_Inner_Loops(block):

  # Find all positions of standalone "do" statements
  do_positions = [match.start()
                  for match in re.finditer(r'^\s*do\b', block, re.MULTILINE)]

  # Find all positions of "end do" (allowing
  # for variations like "enddo" and spaces)
  end_do_positions = [match.end()
                      for match in re.finditer(r'\bend\s*do\b', block)]

  # If we have fewer than 2 "do" statements
  # or no "end do", there's no inner loop
  if len(do_positions) < 2 or len(end_do_positions) == 0:
    return block

  # The second "do" statement starts the inner loop
  inner_do_start = do_positions[1]

  # The first "end do" position marks the end of the inner loop
  inner_end_do_end = end_do_positions[0]

  # Remove everything from the second "do" to the first "end do"
  cleaned_block = block[:inner_do_start] + block[inner_end_do_end:]

  return cleaned_block

#==============================================================================#
#                                                                              #
#   Function to find temporary variables in a block                            #
#                                                                              #
#------------------------------------------------------------------------------#
def Find_Temporary_Variables(block):

  # Remove comments for cleaner parsing
  cleaned_block = re.sub(r'!.*', '', block)

  # Regex to detect variable assignments (excluding arrays)
  assignment_pattern = re.compile(r'(\b\w+\b)\s*=\s*')

  # Find all matches
  matches = assignment_pattern.findall(cleaned_block)

  # Filter out arrays and known global variables (arrays have parentheses)
  temporary_vars = set(var for var in matches if '(' not in var)

  return temporary_vars

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

  # Remove "exp"
  cleaned_block = re.sub(r'\bexp\b', '', cleaned_block)

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

  # Remove "sqrt"
  cleaned_block = re.sub(r'\bsqrt\b', '', cleaned_block)

  # Remove empty parentheses
  cleaned_block = re.sub(r'\(\s*\)', '', cleaned_block)

  # Remove "\n"
  cleaned_block = re.sub(r"\n", "", cleaned_block)

  # Remove "&"
  cleaned_block = re.sub(r"&", "", cleaned_block)

  # Put spaces around mathematical operators
  cleaned_block = re.sub(r'(?<!\s)\+', ' +', cleaned_block)
  cleaned_block = re.sub(r'\+(?!\s)',  '+ ', cleaned_block)
  cleaned_block = re.sub(r'\-(?!\s)',  '- ', cleaned_block)
  cleaned_block = re.sub(r'(?<!\s)\-', ' -', cleaned_block)
  cleaned_block = re.sub(r'(?<!\s)\*', ' *', cleaned_block)
  cleaned_block = re.sub(r'\*(?!\s)',  '* ', cleaned_block)
  cleaned_block = re.sub(r'(?<!\s)\/', ' /', cleaned_block)
  cleaned_block = re.sub(r'\/(?!\s)',  '/ ', cleaned_block)

  result = {}
  seen_arrays = {}  # track arrays to avoid duplicates

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
        seen_arrays[array_name] = None  # mark array as seen
        array_name = re.sub(r'@', ' ', array_name)
        result[array_name] = None

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
  print(f"{BRIGHT_YELLOW}#======================={RESET}")
  print(f"{BRIGHT_YELLOW}#                       {RESET}")
  print(f"{BRIGHT_YELLOW}# Preprocessing a block {RESET}")
  print(f"{BRIGHT_YELLOW}#                       {RESET}")
  print(f"{BRIGHT_YELLOW}#-----------------------{RESET}")
  print(f"{BRIGHT_CYAN}",  end="")
  print(indent + "# Block before preprocessing:")
  print(block,      end="")
  print(f"{RESET}", end="")

  # Initialize pointers and OpenACC clauses
  pointer_setup = ""
  present_setup = ""
  if ("Cells_At_Boundaries_In_Domain_And_Buffers" in block or
      "Cells_At_Boundaries"                       in block or
      "Cells_In_Domain"                           in block or
      "Cells_In_Domain_And_Buffers"               in block):
    present_setup += (indent + begin_parallel_loop_independent + "\n")
  elif ("Faces_In_Region"                in block or
        "Faces_In_Domain_And_At_Buffers" in block):
    present_setup += (indent + begin_parallel_loop + "\n")
  else:  # this covers loops through non-zeroes and all Linalg_Mod
    present_setup += (indent + begin_parallel_loop_independent + "\n")

  if begin_present:
    present_setup += (indent + begin_present + "\n")

  # Set to track processed variables to avoid duplicates
  processed_vars = {}

  # Set to track reduction variables
  sum_reduction_vars = {}
  max_reduction_vars = {}
  min_reduction_vars = {}

  # Set to track reduction variables
  temporary_vars = {}

  # List of macros to exclude from processing as arrays
  excluded_macros = {"Cells_At_Boundaries_In_Domain_And_Buffers",
                     "Cells_At_Boundaries",
                     "Cells_In_Domain",
                     "Cells_In_Domain_And_Buffers",
                     "Faces_In_Region",
                     "Faces_In_Domain_And_At_Buffers",
                     "Face_Value"}

  print("")
  print(f"{RED}  # Pointers used in the block{RESET}")

  #-------------------------------------------------#
  #   Special handling for Faces_In_... variables   #
  #-------------------------------------------------#
  if ("Faces_In_Region"                in block or
      "Faces_In_Domain_And_At_Buffers" in block):

    commands = ("grid_region_f_face => Grid % region % f_face",
                "grid_region_l_face => Grid % region % l_face")

    for command in commands:
      if not Command_In_File(grid_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", grid_to_device_file, sep="")

    if begin_present:
      present_setup += (
        indent + "!$acc   grid_region_f_face,  &\n" +
        indent + "!$acc   grid_region_l_face,  &\n")

    # Add these to the processed set to avoid duplicates
    processed_vars['grid_region_f_face'] = None
    processed_vars['grid_region_l_face'] = None

  #------------------------------------------------------------------#
  #   Special handling for Cells_In_... and Cells_At_... variables   #
  #------------------------------------------------------------------#
  if ("Cells_At_Boundaries_In_Domain_And_Buffers" in block or
      "Cells_At_Boundaries"                       in block or
      "Cells_In_Domain"                           in block or
      "Cells_In_Domain_And_Buffers"               in block):

    commands = ("grid_region_f_cell => Grid % region % f_cell",
                "grid_region_l_cell => Grid % region % l_cell")

    for command in commands:
      if not Command_In_File(grid_to_device_file, command):
        pointer_setup += (indent + command + "\n")
      else:
        print("  ", command, " already in ", grid_to_device_file, sep="")

    if begin_present:
      present_setup += (
        indent + "!$acc   grid_region_f_cell,  &\n" +
        indent + "!$acc   grid_region_l_cell,  &\n")

    # Add these to the processed set to avoid duplicates
    processed_vars['grid_region_f_cell'] = None
    processed_vars['grid_region_l_cell'] = None

  #----------------------------------------------------#
  #   These can, and should, be replaced in any case   #
  #----------------------------------------------------#
  block = re.sub(r'Grid % n_bnd_regions',   'grid_n_bnd_regions', block)
  block = re.sub(r'Grid % n_regions',       'grid_n_regions',     block)
  block = re.sub(r'Grid % region % f_cell', 'grid_region_f_cell', block)
  block = re.sub(r'Grid % region % l_cell', 'grid_region_l_cell', block)
  block = re.sub(r'Grid % region % f_face', 'grid_region_f_face', block)
  block = re.sub(r'Grid % region % l_face', 'grid_region_l_face', block)

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
  print(f"{RED}  # Arrays found in the block{RESET}")

  #------------------------#
  #   Process each array   #
  #------------------------#
  for array in arrays:

    full_name = array

    # Create pointer variable name
    pointer_name = full_name.replace(" % ", "_").lower()

    # Check if this pointer is declared
    if "%" in full_name:  # If it's part of a structure
      if not Command_In_File(gpu_pointers_file, pointer_name):
        print(f"{RED}#============================#{RESET}")
        print(f"{RED}#                            #{RESET}")
        print(f"{RED}#   ERROR IN PREPROCESSING   #{RESET}")
        print(f"{RED}#                            #{RESET}")
        print(f"{RED}#============================#{RESET}")
        print("  Variable", pointer_name,
              "is not defined in", gpu_pointers_file)

    # Skip excluded macros (like Faces_In_Region)
    if any(macro in full_name for macro in excluded_macros):
      continue

    # Check if variable has already been processed
    if pointer_name not in processed_vars:
      if "%" in full_name:  # If it's part of a structure
        # Add to the pointer setup
        command = f"{pointer_name} => {full_name}"

        if Command_In_File(grid_to_device_file, command):
          print("  ", command, " already in ", grid_to_device_file, sep="")
        elif Command_In_File(flow_to_device_file, command):
          print("  ", command, " already in ", flow_to_device_file, sep="")
        elif Command_In_File(turb_to_device_file, command):
          print("  ", command, " already in ", turb_to_device_file, sep="")
        else:
          print(f"{BRIGHT_RED}", end="")
          print(" ", command, " is not present in any of the three files:")
          print(" ", grid_to_device_file, "\n ",
                     flow_to_device_file, "\n ",
                     turb_to_device_file)
          print("  Adding it to this block!")
          pointer_setup += (indent + command + "\n")
          print(f"{RESET}",      end="")

      # Add to the OpenACC setup
      if begin_present:
        present_setup += f"{indent}!$acc   {pointer_name},  &\n"

      # Replace occurrences of the original variable
      # name with the pointer name in the block
      block = re.sub(re.escape(full_name), pointer_name, block)

      # Mark variable as processed
      processed_vars[pointer_name] = None

  #----------------------------------------------------#
  #   Detect reduction variables within the do-loop    #
  #----------------------------------------------------#

  # Pattern to detect scalar reduction operations: var = var + something
  # or var = var - something. The regex ensures that the variable being
  # reduced is a scalar (no parentheses for indices).
  reduction_pattern = re.compile(r'(\b\w+\b)\s*=\s*\1\s*[\+\-]\s*.+')

  # Pattern to detect max reduction operations: var = max(var, something)
  # The regex ensures that the variable being reduced is scalar (no parentheses).
  max_pattern = re.compile(r'(\b\w+\b)\s*=\s*max\s*\(\s*\1\s*,\s*.+\)')

  # Pattern to detect min reduction operations: var = min(var, something)
  # The regex ensures that the variable being reduced is scalar (no parentheses).
  min_pattern = re.compile(r'(\b\w+\b)\s*=\s*min\s*\(\s*\1\s*,\s*.+\)')

  # Block without inner loops
  cleaned_block = Remove_Inner_Loops(block)

  # Find reduction variables in the cleaned_block (only scalars, no arrays)
  reductions = reduction_pattern.findall(cleaned_block)
  max_reductions = max_pattern.findall(cleaned_block)
  min_reductions = min_pattern.findall(cleaned_block)

  # Dictionaries to track reduction variables (both +, -, max, and min)
  sum_reduction_vars = {var: None for var in reductions}
  max_reduction_vars = {var: None for var in max_reductions}
  min_reduction_vars = {var: None for var in min_reductions}

  # Add reduction clause if any reductions are found
  if sum_reduction_vars or max_reduction_vars or min_reduction_vars:
    reduction_clause = []

    if sum_reduction_vars:
      reduction_clause.append(f"reduction(+: {','.join(sum_reduction_vars)})")
    if max_reduction_vars:
      reduction_clause.append(f"reduction(max: {','.join(max_reduction_vars)})")
    if min_reduction_vars:
      reduction_clause.append(f"reduction(min: {','.join(min_reduction_vars)})")

    # Join both reduction clauses if present
    final_reduction_clause = " ".join(reduction_clause)

    # Add the appropriate reduction clause(s) at the end of first line
    # Find the index of the first newline character
    newline_index = present_setup.find("\n")

    # Initialize a variable to track the position
    # of the last letter in the first line
    last_letter_index = -1

    # If a newline character is found, search backward for the last letter
    if newline_index != -1:
      for i in range(newline_index - 1, -1, -1):
        if present_setup[i].isalpha():
          last_letter_index = i
          break

    # Insert the reduction clause after the last letter in the first line
    if last_letter_index != -1:
      present_setup = (  present_setup[:last_letter_index + 1]
                       + f" {final_reduction_clause}"
                       + present_setup[last_letter_index + 1:])

  #-----------------------------------------------#
  #   Find all temporary variables in the block   #
  #    and write them all out as a directive      #
  #        This is needed for OpenMP only!        #
  #-----------------------------------------------#
  temporary_vars = Find_Temporary_Variables(block)

  if begin_private:
    present_setup += (indent + begin_private + "\n")
    for var in temporary_vars:
      if (var not in sum_reduction_vars and
          var not in max_reduction_vars and
          var not in min_reduction_vars):
        present_setup += (indent + "!$omp   " + var + ",  &\n")
    present_setup += (indent + "!$omp )\n")

  #---------------------------------------------------------------#
  #   Append the closing parenthesis for OpenACC present clause   #
  #---------------------------------------------------------------#
  if begin_present:
    present_setup += indent + end_present + "\n"

  #---------------------------------------------------------------#
  #   Replace the last comma with a space in the present clause   #
  #---------------------------------------------------------------#
  last_comma_index = present_setup.rfind(",")
  if last_comma_index != -1:
    present_setup = (
        present_setup[:last_comma_index]
      + " "
      + present_setup[last_comma_index + 1:]
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

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Cells_At_Boundaries_In_Domain_And_Buffers\(\)',
    'grid_region_f_cell(1), grid_region_l_cell(grid_n_regions+1)',
    block
  )

  # Replace the 'do' loop with the OpenACC-parallelized version (if present)
  block = re.sub(
    r'Cells_At_Boundaries\(\)',
    'grid_region_f_cell(1), grid_region_l_cell(grid_n_bnd_regions)',
    block
  )


  # Insert accelerator directives before the second "do" loop
  # (This part is relevant only for OpenACC, since OpenMP does
  #  not accept the inner sequential loop.)
  if begin_sequential_loop:
    first_do_index = block.find("do ")
    if first_do_index != -1:
      second_do_index = block.find("do ", first_do_index + 1)
      if second_do_index != -1:

        # Find the last end of line before the second do statement
        last_eol_index = block.rfind("\n",  0, second_do_index)

        block = (
            block[:last_eol_index + 1]
          + indent
          + begin_sequential_loop + "\n"
          + block[last_eol_index + 1:]
        )

  # Add the end_parallel_loop at the end
  last_end_do_index = block.rfind("end do")
  if last_end_do_index != -1:
    block = (
        block[:last_end_do_index]
      + "end do\n"
      + indent
      + end_parallel_loop
      + block[last_end_do_index + len("end do"):]
    )

  # Now, search for the second-to-last end do before the first one
  if end_sequential_loop:
    second_last_end_do_index = block[:last_end_do_index].rfind("end do")
    if second_last_end_do_index != -1:
      block = (
          block[:second_last_end_do_index]
        + "end do\n"
        + indent
        + end_sequential_loop
        + block[second_last_end_do_index + len("end do"):]
      )

  print(f"{BRIGHT_GREEN}",  end="")
  print(indent + "# Block after preprocessing:")
  print(pointer_setup + present_setup + block, end="")
  print(f"{RESET}", end="")

  # Return the modified block with the pointer setup at the beginning
  return pointer_setup + present_setup + block

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

output_lines.append("!==============================================================================!")
output_lines.append("\n")
output_lines.append("!   DO NOT EDIT THIS FILE! IT IS AUTOMATICALLY GENERATED FROM ITS FPP ORIGIN   !")
output_lines.append("\n")
output_lines.append("!==============================================================================!")
output_lines.append("\n")
output_lines.append("\n")

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

