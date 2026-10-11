"""Optional launcher at the project root: forwards every argument to part1/run_all.py.

    python run_all.py                  # required tasks + the full A-H table
    python run_all.py --required-only
    python run_all.py --data path/to/other_assignment_data.py
Part 2 is validated separately:  python part2/test_hashes.py"""
import os
import runpy
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "part1"))
runpy.run_path(os.path.join(HERE, "part1", "run_all.py"), run_name="__main__")
