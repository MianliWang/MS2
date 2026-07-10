"""Dependency-free XML binary decoding shared by MS1 and MS2 readers."""

from __future__ import annotations

import base64
import re
import struct
import xml.etree.ElementTree as ET
import zlib


def _local_name(tag: str) -> str:
    """Return an XML tag name without its optional namespace."""

    return tag.rsplit("}", 1)[-1]


def _parse_mzxml_duration(value: str) -> float | None:
    """Parse an mzXML ISO-8601 retention-time duration into seconds."""

    match = re.fullmatch(
        r"PT(?:(?P<hours>\d+(?:\.\d+)?)H)?(?:(?P<minutes>\d+(?:\.\d+)?)M)?(?:(?P<seconds>\d+(?:\.\d+)?)S)?",
        value,
    )
    if match is None:
        return None
    parts = {key: float(raw or 0) for key, raw in match.groupdict().items()}
    return parts["hours"] * 3600 + parts["minutes"] * 60 + parts["seconds"]


def _decode_mzxml_peaks(element: ET.Element) -> tuple[list[float], list[float]]:
    """Decode one mzXML ``<peaks>`` element into m/z and intensity arrays."""

    encoded = "".join((element.text or "").split())
    if not encoded:
        return [], []
    try:
        payload = base64.b64decode(encoded, validate=True)
    except ValueError as exc:
        raise ValueError("Invalid base64 peak array in mzXML") from exc
    compression = element.attrib.get("compressionType", "none").lower()
    if compression == "zlib":
        payload = zlib.decompress(payload)
    elif compression not in {"none", ""}:
        raise ValueError(f"Unsupported mzXML peak compression: {compression}")

    precision = element.attrib.get("precision", "32")
    if precision == "32":
        code, width = "f", 4
    elif precision == "64":
        code, width = "d", 8
    else:
        raise ValueError(f"Unsupported mzXML peak precision: {precision}")
    if len(payload) % (2 * width):
        raise ValueError("Malformed mzXML peak array: odd m/z-intensity payload")

    byte_order = element.attrib.get("byteOrder", "network").lower()
    endian = ">" if byte_order in {"network", "big", "big-endian"} else "<"
    values = struct.unpack(f"{endian}{len(payload) // width}{code}", payload)
    return [float(value) for value in values[0::2]], [float(value) for value in values[1::2]]


__all__ = ["_decode_mzxml_peaks", "_local_name", "_parse_mzxml_duration"]
