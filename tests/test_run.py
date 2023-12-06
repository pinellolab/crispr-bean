import pytest
import subprocess


@pytest.mark.order(13)
def test_run_variant_wacc():
    cmd = "bean-run sorting variant tests/data/var_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal_chr6.bw -o tests/test_res/var/ --repguide-mask None"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(14)
def test_run_variant_noacc():
    cmd = "bean-run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(15)
def test_run_variant_wo_negctrl_uniform():
    cmd = "bean-run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ --uniform-edit "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(16)
def test_run_tiling_wo_negctrl():
    cmd = "bean-run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal.bw -o tests/test_res/tiling/ --allele-df-key allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 --control-guide-tag None  --repguide-mask None"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(17)
def test_run_tiling_with_wo_negctrl_noacc():
    cmd = "bean-run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --allele-df-key allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 --control-guide-tag None  --repguide-mask None"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(18)
def test_run_tiling_with_wo_negctrl_uniform():
    cmd = "bean-run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --uniform-edit --allele-df-key allele_counts_spacer_0_19_A.G_translated_prop0.1_0.3 --control-guide-tag None --repguide-mask None"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
