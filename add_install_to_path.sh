#!/bin/bash

# Get the directory containing the script
script_dir=$(dirname "$(realpath "$0")")
echo "Script dir: ${script_dir}"

# Check if the script is inside a Git repository
git_root=""
if git -C "${script_dir}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  # Get the root directory of the Git repository
  git_root=$(git -C "${script_dir}" rev-parse --show-toplevel)
  echo "Git root directory: ${git_root}"
else
  echo "Not inside a Git repository."
fi

# Get the path to INSTALL folder
ROCKABLE_INSTALL_PATH="${git_root}/INSTALL"

# Check if the INSTALL directory exists
if [ -d "$ROCKABLE_INSTALL_PATH" ]; then
  # Add the INSTALL directory to the PATH
  export PATH="${ROCKABLE_INSTALL_PATH}/bin:${PATH}"
  echo "Added ${ROCKABLE_INSTALL_PATH}/bin to the PATH."
  export LD_LIBRARY_PATH="${ROCKABLE_INSTALL_PATH}/lib:${LD_LIBRARY_PATH}"
  echo "Added ${ROCKABLE_INSTALL_PATH}/lib to the LD_LIBRARY_PATH."
  
  # Check if the -always argument is provided
  if [[ "$*" == *"-always"* ]]; then
    # Define the lines to add
    path_line="export PATH=\"${ROCKABLE_INSTALL_PATH}/bin:\${PATH}\""
    ld_line="export LD_LIBRARY_PATH=\"${ROCKABLE_INSTALL_PATH}/lib:\${LD_LIBRARY_PATH}\""

    # Check if the lines are already in ~/.bashrc
    path_exists=$(grep -Fxq "${path_line}" ~/.bashrc; echo $?)
    ld_exists=$(grep -Fxq "${ld_line}" ~/.bashrc; echo $?)

    if [ "${path_exists}" -eq 0 ]; then
      echo "The PATH export line is already present in ~/.bashrc. Skipping."
    else
      echo "${path_line}" >> ~/.bashrc
      echo "Added ${ROCKABLE_INSTALL_PATH} to PATH in ~/.bashrc."
    fi

    if [ "${ld_exists}" -eq 0 ]; then
      echo "The LD_LIBRARY_PATH export line is already present in ~/.bashrc. Skipping."
    else
      echo "${ld_line}" >> ~/.bashrc
      echo "Added ${ROCKABLE_INSTALL_PATH} to LD_LIBRARY_PATH in ~/.bashrc."
    fi
  fi
else
  echo "The INSTALL directory does not exist in the current directory."
  exit 1
fi
