"""Direction-map registry and the ``DirectionMap`` value object.

A direction map is always ``{bare MSDB rxn id -> operator}``, but it lived in
five mutually incompatible on-disk shapes with no schema check anywhere. This
module gives it one type, one validator and one content hash, so two maps can
be compared and a malformed one fails loudly.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field

from .kinds import VALID_DIRECTIONS, DirectionMapSpec
from .registry import Registry, ensure_loaded

_reg: Registry[DirectionMapSpec] = Registry("direction map")
register = _reg.register


def get(key: str) -> DirectionMapSpec:
    ensure_loaded()
    return _reg.get(key)


def all() -> list[DirectionMapSpec]:
    ensure_loaded()
    return _reg.all()


def keys() -> list[str]:
    ensure_loaded()
    return _reg.keys()


#: Word spellings.
NORMALIZE = {
    "forward": ">", "fwd": ">", "f": ">",
    "reverse": "<", "rev": "<", "r": "<",
    "reversible": "=", "both": "=", "bidirectional": "=",
    "unknown": "?", "": "?", "none": "?", "null": "?", "na": "?", "nan": "?",
}

#: Arrow spellings, matched WHOLE. These must never be decided by the first
#: character: "<=>" means reversible but starts with "<", and "=>" means
#: forward but starts with "=". Reading the first character turned every
#: reversible reaction in an arrow-notation table into reverse-only.
ARROWS = {
    "<=>": "=", "<->": "=", "<-->": "=", "<==>": "=",
    "=>": ">", "->": ">", "-->": ">", "==>": ">", "→": ">",
    "<=": "<", "<-": "<", "<--": "<", "<==": "<", "←": "<",
    "↔": "=", "⇌": "=", "⇔": "=",
}


class UnknownOperator(ValueError):
    """A direction value that could not be interpreted."""


def normalize_operator(value, strict: bool = False) -> str:
    """Coerce a direction value to one of ``> < = ?``.

    ``strict=True`` raises :class:`UnknownOperator` instead of falling back to
    ``?``. The fallback matters because ``?`` is applied to a model as
    ``(-1000, 1000)`` -- fully reversible -- so a typo in a direction table
    silently *removes* a constraint rather than failing.
    """
    if value is None:
        return "?"
    v = str(value).strip()
    if v in VALID_DIRECTIONS:
        return v
    if v in ARROWS:
        return ARROWS[v]
    low = v.lower()
    if low in NORMALIZE:
        return NORMALIZE[low]
    if low in ARROWS:
        return ARROWS[low]
    if strict:
        raise UnknownOperator(
            f"cannot interpret direction {value!r}. Use one of > < = ?, a word "
            f"({', '.join(sorted(NORMALIZE))}), or an arrow "
            f"({', '.join(sorted(ARROWS))})."
        )
    return "?"


@dataclass
class DirectionMap:
    """A loaded, validated direction map."""

    key: str
    ops: dict = field(default_factory=dict)
    coverage: str = "partial"
    source: str | None = None
    data_ref: str = "msdb:live"
    masked: bool = False

    @classmethod
    def load(cls, key: str) -> DirectionMap:
        spec = get(key)
        raw = spec.loader()
        ops = {k: normalize_operator(v) for k, v in raw.items()}
        m = cls(key=key, ops=ops, coverage=spec.coverage, source=spec.source,
                data_ref=spec.data_ref, masked=spec.masked)
        m.validate()
        return m

    def validate(self) -> DirectionMap:
        bad_ops = {v for v in self.ops.values()} - VALID_DIRECTIONS
        if bad_ops:
            raise ValueError(f"{self.key}: operators outside {sorted(VALID_DIRECTIONS)}: {sorted(bad_ops)}")
        bad_keys = [k for k in self.ops if not (k.startswith("rxn") and k[3:].isdigit())]
        if bad_keys:
            raise ValueError(
                f"{self.key}: {len(bad_keys)} keys are not bare MSDB reaction ids, "
                f"e.g. {bad_keys[:5]}"
            )
        return self

    @property
    def content_hash(self) -> str:
        blob = json.dumps(self.ops, sort_keys=True, separators=(",", ":")).encode()
        return hashlib.sha256(blob).hexdigest()[:16]

    def counts(self) -> dict:
        out = {op: 0 for op in sorted(VALID_DIRECTIONS)}
        for v in self.ops.values():
            out[v] += 1
        return out

    def overlay_on(self, base: DirectionMap) -> DirectionMap:
        merged = dict(base.ops)
        merged.update(self.ops)
        return DirectionMap(key=f"{self.key}@{base.key}", ops=merged,
                            coverage=base.coverage, source=self.source)

    def diff(self, other: DirectionMap) -> dict:
        shared = set(self.ops) & set(other.ops)
        return {k: (other.ops[k], self.ops[k]) for k in shared if other.ops[k] != self.ops[k]}

    def __len__(self) -> int:
        return len(self.ops)
