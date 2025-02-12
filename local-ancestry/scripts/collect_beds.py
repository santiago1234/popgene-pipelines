import os
import sys
import pandas as pd

def main():
    # Ensure there are enough arguments
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_dir1> <input_dir2> ... <output_dir>")
        sys.exit(1)
    
    # Extract input directories and output directory from arguments
    bed_dirs = sys.argv[1:-1]  # All arguments except the last one are input directories
    outdir = sys.argv[-1]      # The last argument is the output directory

    # Create the output directory if it doesn't exist
    os.makedirs(outdir, exist_ok=True)

    # Initialize an empty dictionary to store the mappings
    id_hap_to_files = {}

    # Loop through each directory in bed_dirs
    for bed_dir in bed_dirs:
        # Loop through each file in the directory
        for filename in os.listdir(bed_dir):
            if filename.endswith(".bed"):
                # Extract ID_HAP by removing the '.bed' extension
                id_hap = filename.split('.')[0]
                
                # If the key does not exist in the dictionary, initialize it with an empty list
                if id_hap not in id_hap_to_files:
                    id_hap_to_files[id_hap] = []
                
                # Append the file path to the list for the corresponding ID_HAP
                id_hap_to_files[id_hap].append(os.path.join(bed_dir, filename))

    # Initialize a dictionary to store the concatenated data frames
    id_hap_to_df = {}

    # Loop through each ID_HAP in the dictionary
    for id_hap, file_list in id_hap_to_files.items():
        # Read each file into a data frame and concatenate them
        df_list = [pd.read_csv(file, sep='\t') for file in file_list]
        concatenated_df = pd.concat(df_list, ignore_index=True)
        
        # Store the concatenated data frame in the dictionary
        id_hap_to_df[id_hap] = concatenated_df
        
        # Define the output file path
        output_file = os.path.join(outdir, f"{id_hap}.bed")
        
        # Save the concatenated data frame to a .bed file
        concatenated_df.to_csv(output_file, sep='\t', index=False)

    # Display the dictionary with file paths for verification
    print(f"Files saved to {outdir}:")
    for id_hap in id_hap_to_df.keys():
        print(f"{id_hap}.bed")

if __name__ == "__main__":
    main()
