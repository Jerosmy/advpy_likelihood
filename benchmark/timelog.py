#!/opt/anaconda3/bin/python

# Benchmark logger: a decorator that times a function and appends
# one row per call to a CSV file (long format: one row = one measurement)

import csv
import os
import platform
import subprocess
import uuid
from datetime import datetime, timezone
from functools import wraps
from time import perf_counter


def getCommit():
    """
    Gets the current git commit, to know which version of the code was measured

    Returns:
    str: short commit hash, with a '-dirty' suffix if there are uncommitted
    changes, or 'unknown' if git is not available
    """
    try:
        commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
        changes = subprocess.run(["git", "status", "--porcelain"],
                                 capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return f"{commit}-dirty" if changes else commit


# Configuration
LOG_DIR = "../benchmark"
LOG_FILE = "benchmarks.csv"
LOG_PATH = os.path.join(LOG_DIR, LOG_FILE)
COLUMNS = ["timestamp", "run_id", "commit", "label", "os", "cpu_count",
           "function_name", "duration_s", "status"]

# Run-level information: computed ONCE, when the module is imported
RUN_ID = uuid.uuid4().hex[:8]                    # groups all rows of one execution
COMMIT = getCommit()
LABEL = os.environ.get("BENCHMARK_LABEL", "")    # free label, e.g. 'baseline'
OS_NAME = platform.system()
CPU_COUNT = os.cpu_count()


def writeLog(row, columns=COLUMNS, path=LOG_PATH):
    """
    Appends one row to the log file. Creates the folder and the header if needed

    Args:
        row (dict): values to write, keys must be in columns
        columns (list of str): names of all columns, in order
        path (str): path of the csv file
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    needsHeader = not os.path.exists(path) or os.path.getsize(path) == 0

    # Append mode: existing rows are never touched
    with open(path, "a", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        if needsHeader:
            writer.writeheader()
        writer.writerow(row)


def timelog(function):
    """
    Decorator: times each call of a function and logs it in the csv file.
    The return value and the exceptions of the function are left untouched

    Args:
        function (callable): the function to measure
    """
    @wraps(function)    # keeps the name and docstring of the original function
    def wrapper(*args, **kwargs):
        # Taken before the timer starts, so it is not part of the measure
        timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
        status = "ok"

        start = perf_counter()
        try:
            return function(*args, **kwargs)
        except Exception:
            status = "error"
            raise                       # the error still reaches the caller
        finally:
            # Runs on success AND on error. The duration is computed
            # before writing, so the cost of logging is not included
            duration = perf_counter() - start
            writeLog({"timestamp": timestamp,
                      "run_id": RUN_ID,
                      "commit": COMMIT,
                      "label": LABEL,
                      "os": OS_NAME,
                      "cpu_count": CPU_COUNT,
                      "function_name": function.__qualname__,
                      "duration_s": duration,
                      "status": status})
    return wrapper