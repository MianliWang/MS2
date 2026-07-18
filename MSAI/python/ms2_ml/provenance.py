"""Formatting-independent content fingerprints for shadow provenance."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from typing import Any


def canonical_json_sha256(value: Any) -> str:
    """Hash JSON content independently of whitespace and key order."""

    encoded = json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def canonical_table_sha256(fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> str:
    """Hash ordered tabular content independently of CSV byte formatting."""

    field_order = [str(field) for field in fields]
    payload = {
        "fields": field_order,
        "rows": [
            ["" if row.get(field) is None else str(row.get(field)) for field in field_order]
            for row in rows
        ],
    }
    return canonical_json_sha256(payload)


__all__ = ["canonical_json_sha256", "canonical_table_sha256"]
