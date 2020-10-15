import scanpy as sc

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

#sc.settings.set_figure_params(dpi=160, frameon=False, figsize=(12, 12))  # low dpi (dots per inch) yields small inline figures
sc.settings.autoshow=False
sc.settings.autosave=True

sc.settings.figdir='./figures/'
sc.settings.cachedir='./cache'

sc.settings.cache_compression="gzip"
sc.settings.n_jobs=8
sc.settings.file_format_figs="pdf"

