import pytest
import subprocess


def test_profile_screen():
    cmd = "bean profile tests/data/var_mini_screen.h5ad "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc

def test_profile_screen_dualeditor():
    cmd = "bean profile tests/data/var_mini_screen_dual.h5ad "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc