.. bean documentation master file, created by
   sphinx-quickstart on Fri Mar 29 19:10:46 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
====================================
Welcome to `bean` documentation!
====================================

.. image:: https://img.shields.io/pypi/pyversions/crispr-bean
   :target: https://pypi.org/project/crispr-bean/
   :alt: PyPI pyversions

.. image:: https://img.shields.io/pypi/v/crispr-bean
   :target: https://pypi.org/project/crispr-bean/
   :alt: PyPI version

.. image:: https://github.com/pinellolab/crispr-bean/actions/workflows/CI.yml/badge.svg
   :target: https://github.com/pinellolab/crispr-bean/actions/workflows/CI.yml
   :alt: Test

.. image:: https://github.com/pinellolab/crispr-bean/actions/workflows/documentation.yml/badge.svg
   :target: https://github.com/pinellolab/crispr-bean/actions/workflows/documentation.yml
   :alt: Documentation

.. image:: https://img.shields.io/badge/License-AGPL%20v3-blue.svg
   :target: https://www.gnu.org/licenses/agpl-3.0
   :alt: License: AGPL v3


`bean` improves CRISPR pooled screen analysis by 1) unconfounding variable per-guide editing outcome by considering genotypic outcome from *reporter* sequence and 2) through accurate modeling of screen procedure.

.. image:: assets/summary.png
  :width: 700
  :alt: BEAN schematic


Workflows
--------------------------
.. toctree::
    :maxdepth: 2

    tutorials


Model description
--------------------------
.. toctree::
    :maxdepth: 2

    model


API references
--------------------------
.. toctree::
    :maxdepth: 3

    input
    subcommands


Screen data structure
--------------------------
.. toctree::

    reporterscreen


Indices and tables
--------------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
