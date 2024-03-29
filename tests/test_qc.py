import pytest
import subprocess


@pytest.mark.order(7)
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


@pytest.mark.order(8)
def test_qc_tiling():
    cmd = "bean-qc tests/data/tiling_mini_screen.h5ad -o tests/data/tiling_mini_screen_masked.h5ad -r tests/test_res/qc_report_tiling_mini_screen --count-correlation-thres 0.6  --posctrl-col ''"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(9)
def test_dummy_insertion_varscreen():
    cmd = "bean-qc tests/data/var_mini_screen_missing.h5ad -o tests/data/var_mini_screen_missing_masked.h5ad -r tests/test_res/qc_report_var_mini_screen_missing --count-correlation-thres 0.6"
    try:
        subprocess.check_output(
            cmd, shell=True, universal_newlines=True, stderr=subprocess.STDOUT
        )
        raise ValueError("Filtering should fail with too small number of replicates.")
    except subprocess.CalledProcessError as exc:
        if "Too small number of replicate left after QC" not in exc.output:
            raise exc


@pytest.mark.order(10)
def test_dummy_insertion_tilingscreen():
    cmd = "bean-qc tests/data/tiling_mini_screen_missing.h5ad -o tests/data/tiling_mini_screen_missing_masked.h5ad -r tests/test_res/qc_report_tiling_mini_screen_missing --count-correlation-thres 0.6  --posctrl-col ''"
    try:
        subprocess.check_output(
            cmd, shell=True, universal_newlines=True, stderr=subprocess.STDOUT
        )
        raise ValueError("Filtering should fail with too small number of replicates.")
    except subprocess.CalledProcessError as exc:
        if "Too small number of replicate left after QC" not in exc.output:
            raise exc
