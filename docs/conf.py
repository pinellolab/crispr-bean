# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "bean"
copyright = "2024, Pinello lab"
author = "Jayoung Ryu, Pinello lab"
release = "1.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinxarg.ext",
    "m2r",
    "sphinx.ext.extlinks",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
extlinks = {
    "git_tag": ("https://github.com/sphinx-doc/alabaster/tree/%s", "%s"),
    "bug": ("https://github.com/sphinx-doc/alabaster/issues/%s", "#%s"),
    "feature": ("https://github.com/sphinx-doc/alabaster/issues/%s", "#%s"),
    "issue": ("https://github.com/sphinx-doc/alabaster/issues/%s", "#%s"),
}


root_doc = "index"
numpydoc_show_class_members = False
source_suffix = [".rst", ".md"]
# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
html_theme = "alabaster"
html_static_path = ["_static"]
html_logo = "assets/beans.svg"
html_theme_options = {
    "description": "Activity-normalized variant effect size estimation from pooled CRISPR screens",
    "github_user": "pinellolab",
    "github_repo": "crispr-bean",
    "github_button": "true",
    "github_count": "false",
}
html_sidebars = {
    "**": [
        "about.html",
        "globaltoc.html",
        "searchbox.html",
    ]
}
html_favicon = "assets/beans_32.ico"
