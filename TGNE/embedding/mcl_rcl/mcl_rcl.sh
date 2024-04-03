#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <undirected, weighted graph/network file in abc-format>"
    echo
    echo "abc-format: Each line encodes an edge in terms of two labels (the 'A' and the 'B') and a numerical value (the 'C'), all separated by white space."
    exit 1  # exit the script with an error code
fi

function display_error {
  echo "Script failed. Printing additional information:"
  echo

  echo "# MCL/RCL install"
  echo "mkdir installmcl"
  echo "cd installmcl"
  echo "wget https://raw.githubusercontent.com/micans/mcl/main/install-this-mcl.sh -o install-this-mcl"
  echo "chmod u+x install-this-mcl.sh"
  echo "./install-this-mcl.sh"
  echo "mcl --version # test install"

  echo "# MOVE THE NEW INSTALL OF MCL/RCL TO BIN ON PATH"
  echo "sudo mv ~/local/bin/* /opt/local/bin/"

  echo "# ENSURE RCL USES BASH 5+ TO HAVE ACCESS TO 'mapfile' command"
  echo "cd ~/../../opt/local/bin/"
  echo "brew install bash"
  echo "which -a bash"
  echo "OUT: <*/homebrew/bin/bash> <– <PATH TO BASH 5+>"
  echo "vim rcl"
  echo ":%s@!/bin/bash@<PATH TO BASH 5+>@g"
}

function handle_error {
  if [ $? -ne 0 ]; then
    display_error
    exit 1
  fi
}

# Get the current working directory
current_directory=$(pwd)

# Check if the last part of the path is not 'mcl_rcl'
if [[ "$current_directory" != *"/mcl_rcl" ]]; then
  new_directory="TGNE/embedding/mcl_rcl"
  cd "$new_directory"
fi

# echo "$current_directory"
# echo "$(pwd)"

# cleanup working directories
rm pfx.*
rm TAG/rcl.*
rm TAG/*.mci
# find TAG -type f -not -name "out.*" -exec rm {} \;

network_file="$1"
network_file_mci="${network_file%.txt}.mci"
network_file_tab="${network_file%.txt}.tab"

# Convert network format from abc to mci and create a tab file for indexing
mcxload --stream-mirror -abc $network_file -o $network_file_mci -write-tab $network_file_tab
handle_error

# Run mcl with specified inflation values
# rcl mcl TAG -p 10 -n $network_file_mci -I "1.3 1.35 1.4 1.45 1.5 1.55 1.6 1.65 1.7 1.8 1.9 2" # FIXME
# handle_error

cd TAG

# Convert .labels clusterings to .mci format
for file in *.labels; do
    # Check if the file exists and is a regular file
    if [ -f "$file" ]; then
        # Run the mcxload command for each file
        mcxload -etc-ai "$file" -strict-tabr "../${network_file_tab}" -o "out.${file%.labels}.mci"
    fi
done

cd ..

# Run rcl setup
rcl setup TAG -F -n $network_file_mci -t $network_file_tab TAG/out.*
handle_error

# Generate quality control plots
rcl-qc qc TAG
#handle_error

rcl-qc qc2 TAG -p 10 -F
#handle_error

# Run rcl tree
rcl tree TAG -p 10 -F
handle_error

# Run rcl select with specified resolution parameters
rcl select TAG -r "100 200 400 800 1600 3200"
handle_error
