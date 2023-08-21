import subprocess


def test_qc():
    cmd = "bean-qc tests/data/bean_count_test_screen.h5ad -o bean_count_test_screen.h5ad -r tests/test_res/qc_report_bean_count_test_screen"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
