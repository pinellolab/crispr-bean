import pytest
import subprocess


@pytest.mark.order(16)
def test_run_variant_wacc():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal_chr6.bw -o tests/test_res/var/ --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(17)
def test_run_variant_noacc():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(18)
def test_run_variant_wo_negctrl_uniform():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ --uniform-edit --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(19)
def test_run_variant_wacc_negctrl():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal_chr6.bw -o tests/test_res/var/ --repguide-mask None --n-iter 10 --fit-negctrl "
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(20)
def test_run_variant_noacc_negctrl():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ --fit-negctrl --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(21)
def test_run_variant_uniform_negctrl():
    cmd = "bean run sorting variant tests/data/var_mini_screen_annotated.h5ad -o tests/test_res/var/ --uniform-edit --fit-negctrl --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(22)
def test_run_tiling_wo_negctrl():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal.bw -o tests/test_res/tiling/ --control-guide-tag None  --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(23)
def test_run_tiling_with_wo_negctrl_noacc():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --control-guide-tag None  --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(23)
def test_run_tiling_with_wo_negctrl_uniform():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --uniform-edit --allele-df-key allele_counts_spacer_0_19_noindels_A.G_translated_prop0.1_0.3 --control-guide-tag None --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(24)
def test_run_tiling_negctrl_allelekey():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad --scale-by-acc --acc-bw-path tests/data/accessibility_signal.bw -o tests/test_res/tiling/ --fit-negctrl --negctrl-col strand --negctrl-col-value neg --control-guide-tag neg --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(25)
def test_run_tiling_with_negctrl_noacc():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --allele-df-key allele_counts_spacer_0_19_noindels_A.G_translated_prop0.1_0.3 --fit-negctrl --negctrl-col strand --negctrl-col-value neg --control-guide-tag neg --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(26)
def test_run_tiling_with_negctrl_uniform():
    cmd = "bean run sorting tiling tests/data/tiling_mini_screen_annotated.h5ad -o tests/test_res/tiling/ --uniform-edit --allele-df-key allele_counts_spacer_0_19_noindels_A.G_translated_prop0.1_0.3 --fit-negctrl --negctrl-col strand --negctrl-col-value neg --control-guide-tag neg --repguide-mask None --n-iter 10"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


# Add fit_negctrl examples


@pytest.mark.order(17)
def test_survival_run_variant_noacc():
    cmd = "bean run survival variant tests/data/survival_var_mini_screen_masked.h5ad -o tests/test_res/var/ --n-iter 10 --control-condition=D7"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(18)
def test_survival_run_variant_wo_negctrl_uniform():
    cmd = "bean run survival variant tests/data/survival_var_mini_screen_masked.h5ad -o tests/test_res/var/ --uniform-edit --n-iter 10 --control-condition=D7"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(20)
def test_survival_run_variant_noacc_negctrl():
    cmd = "bean run survival variant tests/data/survival_var_mini_screen_masked.h5ad -o tests/test_res/var/ --fit-negctrl --n-iter 10 --control-condition=D7"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc


@pytest.mark.order(21)
def test_survival_run_variant_uniform_negctrl():
    cmd = "bean run survival variant tests/data/survival_var_mini_screen_masked.h5ad -o tests/test_res/var/ --uniform-edit --fit-negctrl --n-iter 10 --control-condition=D7"
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as exc:
        raise exc
