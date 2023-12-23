import importlib.metadata

project = "AbStochKin"
copyright = "2024, Alex Plakantonakis"
author = "Alex Plakantonakis"
version = release = importlib.metadata.version("AbStochKin")

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "autodoc2"]
autodoc2_packages = ["../abstochkin"]
autodoc2_render_plugin = "myst"

source_suffix = [".rst", ".md"]
exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    "Thumbs.db",
    ".DS_Store",
    ".env",
    ".venv",
]

html_theme = "furo"

myst_enable_extensions = [
    "colon_fence",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

nitpick_ignore = [
    ("py:class", "_io.StringIO"),
    ("py:class", "_io.BytesIO"),
]

always_document_param_types = True
