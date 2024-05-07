import pytest
import subprocess


@pytest.mark.order(500)
def test_create_screen():
    cmd = "bean create-screen tests/data/var_mini_guides.csv tests/data/var_mini_samples.csv tests/data/var_mini_counts.csv -o tests/test_res/bean_count_var_mini_screen"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(501)
def test_downstream_create_screen():
    cmd = "bean qc tests/test_res/bean_count_var_mini_screen.h5ad -o tests/test_res/bean_count_var_mini_screen_masked.h5ad -r tests/test_res/qc_report_bean_count_var_mini_screen"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
    cmd2 = "bean run sorting variant tests/test_res/bean_count_var_mini_screen_masked.h5ad -o tests/test_res/ --uniform-edit --ignore-bcmatch"
    try:
        subprocess.check_output(
            cmd2,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
