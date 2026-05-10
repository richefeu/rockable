#!/bin/sh

INPUT_LIST="$1"
COMMAND="$2"
NPROCS="$3"

if [ $# -ne 3 ]; then
    echo "Usage: $0 input_list command nprocs"
    exit 1
fi

xargs -I{} -P "$NPROCS" sh -c '
    input="$1"
    command="$2"

    dir=$(dirname "$input")
    file=$(basename "$input")

    echo "Running $file in folder $dir"

    cd "$dir" || {
        echo "ERROR: cannot cd into $dir" >&2
        exit 1
    }

    eval "$command \"$file\"" > run.log 2>&1
' _ {} "$COMMAND" < "$INPUT_LIST"