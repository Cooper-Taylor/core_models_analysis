"""Entry modules: the six files an author actually edits.

Imported explicitly and in a fixed order by :func:`cma.registry.load` -- never
by directory scan, because ten modules under ``scripts/`` call ``sys.exit`` at
import time.
"""
