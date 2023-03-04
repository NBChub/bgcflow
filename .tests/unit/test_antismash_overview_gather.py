import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

import common

sys.path.insert(0, os.path.dirname(__file__))


def test_antismash_overview_gather():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/antismash_overview_gather/data")
        expected_path = PurePosixPath(".tests/unit/antismash_overview_gather/expected")

        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "data/processed/Lactobacillus_delbrueckii/tables/df_antismash_6.1.1_bgc.csv",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "data/processed/Lactobacillus_delbrueckii/tables/df_regions_antismash_6.1.1.csv",
                "-f",
                "-j1",
                "--keep-target-files",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
