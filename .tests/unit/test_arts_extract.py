import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_arts_extract():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/arts_extract/data")
        expected_path = PurePosixPath(".tests/unit/arts_extract/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("data/interim/arts/antismash-7.0.0/GCA_000056065.1_arts_all_hits.json data/interim/arts/antismash-7.0.0/GCA_000056065.1_arts_bgctable_summary.json data/interim/arts/antismash-7.0.0/GCA_000056065.1_arts_coretable_summary.json data/interim/arts/antismash-7.0.0/GCA_000056065.1_arts_knownhits.json", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "data/interim/arts/antismash-7.0.0/GCA_000056065.1_arts_all_hits.json",
            "-f", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
