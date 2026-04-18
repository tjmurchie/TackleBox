"""
MetaMerge — ensemble aDNA classification by merging MEGAN and Holi/metaDMG.

MetaMerge combines a conservative BLASTn→MEGAN eukaryote count matrix with
Holi/metaDMG damage-support outputs to assign reproducible, DNA-only support
categories to each detected taxon.

Typical usage via the command line::

    metamerge check --megan-counts counts.tsv --holi metadmg.csv --linker linker.csv
    metamerge run   --megan-counts counts.tsv --holi metadmg.csv --linker linker.csv --outdir results/

See README.md for a full description of inputs, outputs, and the classification
logic used by the package.
"""

__version__ = "1.1.0"
__author__ = "Tyler Murchie"
