#!/bin/bash

# List of files to keep
files_to_keep=(
"SF10022" "SF10484" "SF10099" "SF10099v2" "SF10441" "SF10127" "SF10565" "SF10565v2" 
"SF10432" "SF11248" "SF10592" "SF11857" "SF11082" "SF11488" "SF11344" "SF12460" 
"SF1199" "SF1343" "SF11780" "SF12243" "SF12243v2" "SF2501" "SF2628" "SF2777" 
"SF2979" "SF2990" "SF3073" "SF3076" "SF3243" "SF3391" "SF3448" "SF3996" "SF4449" 
"SF4449v2" "SF4297" "SF6621" "SF4400" "SF6186" "SF4810" "SF4810v2" "SF6098" 
"SF4849" "SF6118" "SF6118v2" "SF5581" "SF9715" "SF9715v2" "SF6809" "SF7025" 
"SF6996" "SF7062" "SF7307" "SF7307v2" "SF7388" "SF8963" "SF12165" "SF12165v2" 
"SF9259R/S" "SF9510" "SF9358" "SF9962" "SF9798" "SF9494" "SF9871" "SF10108" 
"SF9372" "SF10433" "SF10433v2" "SF4209" "SF4209v2" "SF4324" "SF11916" "SF12382" 
"SF11815" "SF12408" "SF10857" "SF12008" "SF12115" "SF12751" "SF11587" "SF11981" 
"SF12407" "SF12754" "SF11720" "SF11720v2" "SF11873" "SF11331"
)

# Create extra directory if it doesn't exist

extra='~/project_GBM/gbm_DATA/gbm_DATA_single_atlas_extra'
# Iterate over all files in the current directory
for file in *; do
    # Skip if it's a directory
    [ -d "$file" ] && continue

    # Check if file should be kept
    keep=false
    for prefix in "${files_to_keep[@]}"; do
        if [[ $file == ${prefix}* ]]; then
            keep=true
            break
        fi
    done

    # Move file to extra directory if it shouldn't be kept
    if [ "$keep" = false ]; then
        mv "$file" extra/
    fi
done
