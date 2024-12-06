import os
import subprocess
import gzip
import shutil


def decompress_gzip(file_path, output_path):
    """Decompress a gzipped file."""
    with gzip.open(file_path, 'rt') as f_in, open(output_path, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)


def compress_gzip(file_path, output_path):
    """Compress a file into gzip format."""
    with open(file_path, 'rb') as f_in, gzip.open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def convert_to_bash_path(path):
    """Convert a Windows-style path to a Unix-style path for Git Bash."""
    return path.replace("\\", "/").replace("C:", "/c")


def run_bbmerge(bbmap_dir, in1, in2, out, outu1, outu2):
    """Run the BBMerge tool on paired FASTQ files."""
    # Convert paths to Unix-style for Git Bash
    bbmap_dir_bash = convert_to_bash_path(bbmap_dir)
    in1_bash = convert_to_bash_path(in1)
    in2_bash = convert_to_bash_path(in2)
    out_bash = convert_to_bash_path(out)
    outu1_bash = convert_to_bash_path(outu1)
    outu2_bash = convert_to_bash_path(outu2)

    # Prepare command
    command = [
        "bash",
        "./bbmerge.sh",
        f"in1={in1_bash}",
        f"in2={in2_bash}",
        f"out={out_bash}",
        f"outu1={outu1_bash}",
        f"outu2={outu2_bash}",
    ]
    print(f"outu1={outu1_bash}")
    print(f"outu2={outu2_bash}")

    try:
        # Change directory to BBMap directory
        current_dir = os.getcwd()
        os.chdir(bbmap_dir)  # Change to the bbmap directory for proper execution

        # Execute the command
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print output and errors
        print(f"Processing {in1} and {in2}...")
        print("Output:", result.stdout)
        if result.stderr:
            print("Error:", result.stderr)

        # Restore the original directory
        os.chdir(current_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
        os.chdir(current_dir)  # Ensure we return to the original directory


def process_all_pairs(bbmap_dir, fastq_dir):
    """Process all R1 and R2 pairs in the directory."""
    # List all gzipped FASTQ files
    fastq_files = [f for f in os.listdir(fastq_dir) if f.endswith(".fastq.gz")]

    # Group R1 and R2 pairs by matching names
    pairs = {}
    for file in fastq_files:
        if "_R1_" in file:
            base_name = file.replace("_R1_", "_")
            pairs.setdefault(base_name, [None, None])[0] = file
        elif "_R2_" in file:
            base_name = file.replace("_R2_", "_")
            pairs.setdefault(base_name, [None, None])[1] = file

    # Process each pair
    for base_name, (r1_file, r2_file) in pairs.items():
        if not r1_file or not r2_file:
            print(f"Skipping incomplete pair: {base_name}")
            continue

        # Paths for decompressed files
        r1_decompressed = os.path.join(fastq_dir, r1_file.replace(".gz", ""))
        r2_decompressed = os.path.join(fastq_dir, r2_file.replace(".gz", ""))

        # Decompress input files
        decompress_gzip(os.path.join(fastq_dir, r1_file), r1_decompressed)
        decompress_gzip(os.path.join(fastq_dir, r2_file), r2_decompressed)

        # Paths for output files
        merged_output = os.path.join(fastq_dir, f"{base_name}_merged.fastq")
        unassembled_forward = os.path.join(fastq_dir, f"{base_name}_unassembled_forward.fastq")
        unassembled_reverse = os.path.join(fastq_dir, f"{base_name}_unassembled_reverse.fastq")

        # Run BBMerge
        run_bbmerge(
            bbmap_dir=bbmap_dir,
            in1=r1_decompressed,
            in2=r2_decompressed,
            out=merged_output,
            outu1=unassembled_forward,
            outu2=unassembled_reverse
        )

        # # Compress output files
        # compress_gzip(merged_output, f"{merged_output}.gz")
        # compress_gzip(unassembled_forward, f"{unassembled_forward}.gz")
        # compress_gzip(unassembled_reverse, f"{unassembled_reverse}.gz")

        # # Cleanup decompressed files
        # os.remove(r1_decompressed)
        # os.remove(r2_decompressed)
        # os.remove(merged_output)
        # os.remove(unassembled_forward)
        # os.remove(unassembled_reverse)

        print(f"Completed processing for pair: {base_name}")


if __name__ == '__main__':

    # Paths and directories
    bbmap_dir = r"C:/Users/User/Downloads/BBMap_39.10/bbmap/"
    fastq_dir = r"C:\Users\User\PycharmProjects\dna_storage_comb_syn_analysis\composite_analysis\ratio_2_4\data"

    # Process all R1 and R2 pairs in the folder
    process_all_pairs(bbmap_dir, fastq_dir)
