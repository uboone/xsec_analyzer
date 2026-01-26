from datetime import datetime
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'xsec_analyzer'
copyright = f'{datetime.now().year}, Steven Gardiner, Burke Irwin, Liang Lu, et. al.'
author = 'Steven Gardiner, Burke Irwin, Liang Lu, et. al.'

show_authors = True

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.extlinks',
]

extlinks = {
    'github': ('https://github.com/%s',
               None),
    'docdb': (
        'https://microboone-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=%s',
        None
    )
}

templates_path = ['_templates']
exclude_patterns = []

# hawkmoth_root = os.path.abspath('../../include/')
# hawkmoth_clang = [
#     '-I/home/niam/phd/code/xsec_analyzer/include',
#     '-I/home/niam/Downloads/root/include'
# ]
# hawkmoth_source_uri = "https://github.com/uboone/xsec_analyzer/tree/main/include/{source}#L{line}"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
html_theme_options = {
    "source_repository": "https://github.com/uboone/xsec_analyzer/",
    "source_branch": "main",
    "source_directory": "docs/source/"
}


# -- For theme = 'alabaster' -----------
# html_theme = 'alabaster'
# html_static_path = ['_static']
# html_sidebars = {
#     '**': [
#         'about.html',
#         'navigation.html',
#         'relations.html',
#         'searchfield.html',
#     ]
# }
# html_theme_options = {
#     'description': 'A framework for MicroBooNE cross-sections',
#     'github_user': 'uboone',
#     'github_repo': 'xsec_analyzer',
#     'page_width': '60rem',
#     'sidebar_width': '15rem',
#     'extra_nav_links': {
#         "GitHub Repository": "https://github.com/uboone/xsec_analyzer"
#     }
# }
