#!/usr/bin/env python3
"""TackleBox: Spinner — thin entrypoint.

All logic lives in the ``spinner/`` package.  This file exists so that the
bash wrapper (``./Spinner``) can call ``python3 Spinner.py "$@"`` unchanged.
"""
from spinner.cli import main

if __name__ == "__main__":
    main()
