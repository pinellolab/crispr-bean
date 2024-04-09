import pytest
import subprocess


def test_create_screen():
    cmd = "bean profile tests/data/var_mini_screen.h5ad "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
