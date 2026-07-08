from __future__ import annotations

from pandora.schemas.annotation import AnnotationLayer

__all__ = ["AnnotationLayer"]

# TODO: Add shared annotation plugin protocols here. This module should stay
# small: common annotation types, plugin call signatures, and validation helpers
# belong here; concrete annotation algorithms belong in entry.py, pairwise.py,
# or future domain-specific modules.
