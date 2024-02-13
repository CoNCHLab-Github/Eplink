#!/bin/bash

# Make sure that graham is mounted at ~/Graham (using sshfs)
# Adjust the path to fMRIprep output
fMRIprep_path="/home/ali/Graham/EpLink/eplink-p3-fMRIprep-freesurfer"
# Adjust the destination to which you wannt to copy the reports
destination_path="./fMRIprep-reports-p3"
mkdir -p "$destination_path"

# Find all subject folders
pattern="sub-"
directories=$(find "$fMRIprep_path" -maxdepth 1 -type d -name "*$pattern*")

# Copy figure folders from subject folder (source) to destination
for source_folder in $directories; do
    # Extract the sub-folder name (e.g., sub-01, sub-02)
    sub_folder_name=$(basename "${source_folder}")
    echo "Fetching reports for ${sub_folder_name}"
    # Build the destination path
    destination_folder="${destination_path}/${sub_folder_name}"
    
    # Create destination folder if it doesn't exist
    mkdir -p "${destination_folder}"

    # Copy html file
    cp "${fMRIprep_path}/${sub_folder_name}.html" "${destination_path}"
    # Copy figure folders
    cp -r "${source_folder}/figures" "${destination_folder}"
done

echo "All reports fetched successfully."

