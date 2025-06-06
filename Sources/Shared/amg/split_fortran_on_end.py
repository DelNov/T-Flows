def split_fortran_on_end(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    block = []
    block_num = 1

    for line in lines:
        block.append(line)
        if line.strip().lower() == 'end':
            out_filename = f"{block_num}.f"
            with open(out_filename, 'w') as out:
                out.writelines(block)
            print(f"Wrote {out_filename}")
            block = []
            block_num += 1

    if block:
        print(f"Warning: final block {block_num} has no 'end' — not written.")

# Example usage:
split_fortran_on_end('amg1r5.f')

