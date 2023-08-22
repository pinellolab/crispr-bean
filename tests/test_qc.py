import pytest
import subprocess


@pytest.mark.order(3)
def test_qc():
    cmd = "bean-qc tests/data/var_mini_screen.h5ad -o tests/data/var_mini_screen_masked.h5ad -r tests/test_res/qc_report_var_mini_screen --count-correlation-thres 0.6"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(4)
def test_qc_tiling():
    cmd = "bean-qc tests/data/tiling_mini_screen.h5ad -o tests/data/tiling_mini_screen_masked.h5ad -r tests/test_res/qc_report_tiling_mini_screen --count-correlation-thres 0.6"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
