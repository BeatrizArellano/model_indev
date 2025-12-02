#!/bin/bash

# Read JSON configuration
config_file="config_launch.json"

# Obtaining configuration values
path_program=$(jq -r '.path_program' "$config_file")
threads=$(jq -r '.threads // empty' "$config_file")
percentage=$(jq -r '.percentage // empty' "$config_file")
environment=$(jq -r '.environment' "$config_file")
dynamic=$(jq -r '.enable_dynamic // empty' "$config_file")

# Check if the program path is valid
if [[ ! -x "$path_program" ]]; then
    echo -e "Error: The program at '$path_program' is not executable or does not exist. \n\tVerify the path or that the program was compiled successfully."

    exit 1
fi

# Defining the number of threads to be used
if [[ -n "$percentage" ]]; then
	if [[ "$percentage" =~ ^[0-9]+$ ]] && [ "$percentage" -ge 0 ] && [ "$percentage" -lt 100 ]; then
		# If a valid percentage is provided, calculate the number of threads
		total_threads=$(nproc --all)  # Get the total number of available threads
		num_threads=$(( (total_threads - 1) * percentage / 100 ))  # Leave 1 thread free
	    else
		echo "Error: Invalid percentage value. It must be an integer between 0 and 99."
		exit 1
	 fi
elif [[ -n "$threads" && "$threads" -gt 0 ]]; then
    # Use the specified number of threads directly if percentage is not provided
    num_threads=$threads
else
    echo "Error: No valid thread configuration found in $config_file."
    exit 1
fi

# Print the number of threads being used
echo "Running model using $num_threads processor threads"

# Setting the number of threads
export "$environment=$num_threads"
# Set OMP_DYNAMIC based on configuration
if [[ "$dynamic" == "true" ]]; then
    export OMP_DYNAMIC=1  # Enable dynamic adjustment
fi

# Running model
"$path_program"
