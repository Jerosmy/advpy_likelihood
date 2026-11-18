# timelog

A small Python decorator that times a function and appends one row per call to a CSV file. The file is in long format (one row = one measurement), so it can be loaded directly in R or pandas to compare versions of your code while optimizing it.

## Features

- One decorator, no configuration needed: `@timelog`
- The return value and the exceptions of the decorated function are left untouched
- Failed calls are logged too (`status = error`)
- The git commit is recorded automatically, with a `-dirty` suffix if there are uncommitted changes
- An optional free label (e.g. `baseline`, `with_cache`) set through an environment variable
- Uses only the standard library

## Requirements

- Python 3.7 or later
- git (optional): without it, the `commit` column is `unknown`

## Installation

Copy `timelog.py` next to your script, or anywhere on your Python path. Nothing to install.

Add the log folder to your `.gitignore`, otherwise the CSV file itself counts as an uncommitted change and every commit will end with `-dirty`:

```
benchmark/
```

## Usage

### Decorate a function or a method

```python
from timelog import timelog

@timelog
def slowSum(n):
    return sum(i * i for i in range(n))

class Pipeline:
    @timelog
    def process(self, data):
        ...

slowSum(1_000_000)    # one call = one row in benchmark/benchmarks.csv
```

To get several measurements of the same function, call it several times (for example in a `for` loop). One measurement is not reliable: compare medians and spread, not a single run.

If you combine `@timelog` with `@staticmethod` or `@classmethod`, put `@timelog` below it (closest to the function):

```python
class Pipeline:
    @staticmethod
    @timelog
    def helper(x):
        ...
```

### Label the version being tested

Set the `BENCHMARK_LABEL` environment variable before running your script. The code does not need to change.

```bash
# Linux / macOS
BENCHMARK_LABEL=baseline python my_script.py

# Windows (cmd)
set BENCHMARK_LABEL=baseline
python my_script.py
```

Change the label after each optimization (`with_cache`, `numpy`, ...) and run again.

### Configuration

Constants at the top of `timelog.py`:

| Constant | Default | Meaning |
|---|---|---|
| `LOG_DIR` | `benchmark` | Folder of the log file (created if missing) |
| `LOG_FILE` | `benchmarks.csv` | Name of the log file |
| `COLUMNS` | see below | Column names and order |

## Output format

The file is created on the first logged call, with its header. Rows are always appended: existing rows are never modified.

| Column | Description |
|---|---|
| `timestamp` | Start of the call, UTC, format `YYYY-MM-DDTHH:MM:SS` |
| `run_id` | Random id generated once per program execution; groups all rows of one run |
| `commit` | Short git hash, `-dirty` suffix if uncommitted changes, `unknown` if git is unavailable |
| `label` | Value of `BENCHMARK_LABEL` (empty if not set) |
| `os` | Operating system (`platform.system()`) |
| `cpu_count` | Number of logical CPUs |
| `function_name` | Qualified name, e.g. `Pipeline.process` |
| `duration_s` | Wall-clock duration in seconds (float) |
| `status` | `ok` or `error` |

Example:

```
timestamp,run_id,commit,label,os,cpu_count,function_name,duration_s,status
2026-11-19T00:12:03,a3f9c1e2,4b2d7a1,baseline,Linux,8,slowSum,0.0831,ok
2026-11-19T00:12:03,a3f9c1e2,4b2d7a1,baseline,Linux,8,slowSum,0.0824,ok
```

## Analysis in R

```r
library(tidyverse)

read_csv("benchmark/benchmarks.csv") |>
  filter(status == "ok") |>
  group_by(label, function_name) |>
  summarise(median_s = median(duration_s), n = n(), .groups = "drop")
```

The column is called `function_name` rather than `function` because `function` is a reserved word in R.

## How it works

- **Once, at import:** the run id, the git commit, the label, the OS and the CPU count are computed.
- **At each call:** the wrapper takes the timestamp, starts a timer (`time.perf_counter`), runs the function, stops the timer, then writes the row. The row is written in a `finally` block, so it is logged even if the function raises. The duration is computed before writing, so the cost of logging is not part of the measure.

## Limitations

**Measurement**

- **Nested calls overlap.** If a decorated function calls another decorated function, the duration of the first includes the second. Do not sum durations across functions to get a total time. Recursive functions log one row per recursion level.
- **Overhead on fast functions.** The log file is opened and closed at every call. This is fine for functions that take more than a few milliseconds. For a very fast function called thousands of times, the logging cost can distort what you measure; accumulate results in memory and write them at the end instead.
- **Wall-clock time only.** It includes time spent waiting (I/O, network, other processes). CPU time is not recorded. Run benchmarks under comparable conditions (same data, same machine, no other heavy load).
- **Context is limited.** Input size and arguments are not recorded. A duration is hard to interpret without knowing the size of the data, so add a column if you need it (see below).

**Types of functions**

- **`async` functions** are not measured correctly: only the creation of the coroutine is timed, not its execution.
- **Generators** have the same problem: only the creation of the generator object is timed, not the iteration.

**Versioning**

- **`-dirty` detection counts untracked files too.** A new, uncommitted file in the repository makes the commit `-dirty`. Commit (or ignore) your files before benchmarking.
- **Git commands run in the current working directory**, not in the folder of the script. Launch your script from inside the repository, otherwise the commit will be `unknown`.
- **The commit and the label are fixed at import time.** Changing `BENCHMARK_LABEL` while the program is running has no effect.
- **The commit only describes the code.** Machine, Python version, library versions and data also change the timings; only the OS and CPU count are recorded.

**File handling**

- **The log path is relative** to the working directory, so running the script from another folder creates another `benchmark/` folder there.
- **Changing `COLUMNS` on an existing file** produces rows that do not match the header. Start a new file (or rename the old one) when the columns change.
- **Concurrent writes are not safe.** Several processes writing to the same file (multiprocessing, parallel runs) can interleave their rows. Use one file per process and merge them in R.
- **A failure in the logging itself** (permission denied, disk full) raises an exception from the `finally` block, which can hide the result or the original error of the decorated function.
- **Interrupting with Ctrl+C** (`KeyboardInterrupt`) logs the call with `status = ok`, because only `Exception` is caught. Discard the last row of an interrupted run.

## Adding a column

1. Add its name to `COLUMNS`.
2. Fill it in the dictionary passed to `writeLog` in the wrapper (or compute it once at the top of the module if it is run-level information).
3. Start a new log file so that all rows share the same columns.