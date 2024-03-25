import pytest
import subprocess

def test_create_screen():
    cmd = "bean-create-screen tests/data/var_mini_guides.csv tests/data/var_mini_samples.csv tests/data/var_mini_counts.csv"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
