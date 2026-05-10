# `batchable.sh`

Portable shell script to launch simulation jobs in parallel from a list of input files.

Works on:

* macOS
* Linux
* HPC clusters
* any POSIX-compatible shell environment

The script uses `xargs` to run multiple jobs concurrently while limiting the maximum number of parallel processes.

---

# Features

* Portable POSIX shell (`/bin/sh`)
* Parallel execution with configurable number of workers
* Executes each job inside the directory containing the input file
* Redirects output and errors to `run.log`
* Handles directory change failures safely
* Works with arbitrary commands

---

# File Structure Example

```text
project/
├── batchable.sh
├── all_inputs.txt
├── sim1/
│   └── input.txt
├── sim2/
│   └── input.txt
└── BUILD/
    └── rockable
```

Example `all_inputs.txt`:

```text
sim1/input.txt
sim2/input.txt
```

---

# Usage

```bash
./batchable.sh input_list command nprocs
```

Arguments:

| Argument     | Description                                       |
| ------------ | ------------------------------------------------- |
| `input_list` | Text file containing one input file path per line |
| `command`    | Command to execute                                |
| `nprocs`     | Maximum number of parallel processes              |

---

# Example

```bash
./batchable.sh all_inputs.txt ../../BUILD/rockable 4
```

For each line in `all_inputs.txt`, the script will:

1. Determine the directory containing the input file
2. Change into that directory
3. Execute:

```bash
../../BUILD/rockable input.txt
```

At most 4 jobs will run simultaneously.

---

# Logs

For each simulation directory, the script creates:

```text
run.log
```

containing:

* standard output
* standard error

from the executed command.

---

# Example Workflow

Suppose:

```text
all_inputs.txt
```

contains:

```text
cases/A/input.txt
cases/B/input.txt
cases/C/input.txt
```

Running:

```bash
./batchable.sh all_inputs.txt ../../BUILD/rockable 2
```

will execute:

```bash
cd cases/A && ../../BUILD/rockable input.txt
cd cases/B && ../../BUILD/rockable input.txt
cd cases/C && ../../BUILD/rockable input.txt
```

with at most 2 concurrent jobs.

---

# Error Handling

If a directory cannot be entered, the script prints:

```text
ERROR: cannot cd into <directory>
```

and skips that job safely.

Other running jobs continue normally.

---

# Advanced Usage

Because the command is evaluated dynamically, additional arguments can be passed.

Example:

```bash
./batchable.sh all_inputs.txt "python3 run.py --mode fast" 8
```

This becomes:

```bash
python3 run.py --mode fast input.txt
```

inside each simulation directory.

---

# Permissions

Make the script executable:

```bash
chmod +x batchable.sh
```

---

# Dependencies

Required tools:

* POSIX shell (`/bin/sh`)
* `xargs`
* `dirname`
* `basename`

These are available by default on most Unix-like systems.



