import subprocess


def test_filter_varscreen():
    cmd = "bean-filter tests/data/var_mini_screen.h5ad -o tmp -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


def test_filter_tiling_screen():
    cmd = "bean-filter tests/data/tiling_mini_screen.h5ad -o tmp -s 0 -e 19 -w -b -t -ap 0.1 -sp 0.3"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
