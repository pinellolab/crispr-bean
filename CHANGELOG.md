# Changelog
## 1.2.9
* Add sgRNA information as output
* Add pseudocounts for LFC calculation (`perturb-tools>=0.3.5`)
## 1.2.8
* Change .pyx files to be compatible with more recent numpy versions
## 1.2.7
* **CRITICAL** Fix sample ordering & masking issue for survival screens
## 1.2.6
* Fix overflow in `bean run survival` and autograde error related to inplace assignment for `bean run survival tiling`.
## 1.2.5
* Allow `bean run .. tiling` for untranslated `--allele-df-key`.