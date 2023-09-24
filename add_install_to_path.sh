#!/bin/bash

# Get the directory containing the script
script_dir=$(dirname "$(realpath "$0")")

# Check if the script is inside a Git repository
git_root=""
if git -C "$script_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  # Get the root directory of the Git repository
  git_root=$(git -C "$script_dir" rev-parse --show-toplevel)
  echo "Git root directory: $git_root"
else
  echo "Not inside a Git repository."
fi

# Get the path to INSTALL folder
ROCKABLE_INSTALL_PATH="$git_root/INSTALL"

# Check if the INSTALL directory exists
if [ -d "$ROCKABLE_INSTALL_PATH" ]; then
  # Add the INSTALL directory to the PATH
  export PATH="$ROCKABLE_INSTALL_PATH:$PATH"
  echo "Added $ROCKABLE_INSTALL_PATH to the PATH."
  export LD_LIBRARY_PATH="$ROCKABLE_INSTALL_PATH:$LD_LIBRARY_PATH"
  echo "Added $ROCKABLE_INSTALL_PATH to the LD_LIBRARY_PATH."
else
  echo "The INSTALL directory does not exist in the current directory."
fi
