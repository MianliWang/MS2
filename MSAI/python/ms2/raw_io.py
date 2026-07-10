"""Streaming mzML/mzXML readers for centroided MS2 spectra."""

from __future__ import annotations

import base64
import importlib
import struct
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path
from typing import Any, cast

from .models import Spectrum
from .utils import number
from .xml_codec import _decode_mzxml_peaks, _local_name, _parse_mzxml_duration


def load_ms2_spectra(path: Path) -> list[Spectrum]:
    """Read centroided MS2 scans from mzML or mzXML.

    ``pyopenms`` is used when installed.  Conventional ProteoWizard mzML and
    mzXML files also have streaming standard-library fallbacks, so the core
    workflow has no mandatory third-party raw-reader dependency.
    """

    path = Path(path)
    suffix = path.suffix.lower()
    try:
        oms = cast(Any, importlib.import_module("pyopenms"))
    except ImportError as exc:
        if suffix == ".mzxml":
            return _load_mzxml_stdlib(path)
        if suffix == ".mzml":
            return _load_mzml_stdlib(path)
        raise ValueError(f"Unsupported raw file type: {path.name}") from exc

    experiment = oms.MSExperiment()
    if suffix == ".mzml":
        oms.MzMLFile().load(str(path), experiment)
    elif suffix == ".mzxml":
        oms.MzXMLFile().load(str(path), experiment)
    else:
        raise ValueError(f"Unsupported raw file type: {path.name}")

    spectra: list[Spectrum] = []
    for scan in experiment:
        if scan.getMSLevel() != 2:
            continue
        precursors = scan.getPrecursors()
        if not precursors:
            continue
        mz, intensity = scan.get_peaks()
        precursor = precursors[0]
        lower_offset = float(precursor.getIsolationWindowLowerOffset())
        upper_offset = float(precursor.getIsolationWindowUpperOffset())
        precursor_mz = float(precursor.getMZ())
        spectra.append(
            Spectrum(
                precursor_mz=precursor_mz,
                rt=float(scan.getRT()),
                mz=[float(value) for value in mz],
                intensity=[float(value) for value in intensity],
                isolation_lower_mz=precursor_mz - lower_offset if lower_offset > 0 else None,
                isolation_upper_mz=precursor_mz + upper_offset if upper_offset > 0 else None,
            )
        )
    return spectra


def _load_mzxml_stdlib(path: Path) -> list[Spectrum]:
    """Stream a centroided mzXML file without third-party dependencies."""

    spectra: list[Spectrum] = []
    for _event, element in ET.iterparse(path, events=("end",)):
        if _local_name(element.tag) != "scan":
            continue
        try:
            ms_level = int(element.attrib.get("msLevel", "0"))
        except ValueError:
            ms_level = 0
        if ms_level == 2:
            precursor_element = next(
                (child for child in element if _local_name(child.tag) == "precursorMz"),
                None,
            )
            peaks_element = next(
                (child for child in element if _local_name(child.tag) == "peaks"),
                None,
            )
            if precursor_element is not None and peaks_element is not None:
                precursor_mz = number(precursor_element.text)
                rt = _parse_mzxml_duration(element.attrib.get("retentionTime", ""))
                if precursor_mz is not None and rt is not None:
                    mz, intensity = _decode_mzxml_peaks(peaks_element)
                    width = number(precursor_element.attrib.get("windowWideness"))
                    spectra.append(
                        Spectrum(
                            precursor_mz=precursor_mz,
                            rt=rt,
                            mz=mz,
                            intensity=intensity,
                            isolation_lower_mz=(precursor_mz - width / 2) if width else None,
                            isolation_upper_mz=(precursor_mz + width / 2) if width else None,
                        )
                    )
        element.clear()
    return spectra


def _load_mzml_stdlib(path: Path) -> list[Spectrum]:
    """Stream centroided MS2 spectra from conventional ProteoWizard mzML."""

    spectra: list[Spectrum] = []
    for _event, element in ET.iterparse(path, events=("end",)):
        if _local_name(element.tag) != "spectrum":
            continue
        if _mzml_cv_value(element, "MS:1000511") == "2":
            precursor_mz = number(
                _mzml_cv_value(element, "MS:1000827")
                or _mzml_cv_value(element, "MS:1000744")
            )
            rt = _mzml_scan_time(element)
            arrays = _decode_mzml_arrays(element)
            if precursor_mz is not None and rt is not None and "mz" in arrays and "intensity" in arrays:
                lower_offset = number(_mzml_cv_value(element, "MS:1000828"))
                upper_offset = number(_mzml_cv_value(element, "MS:1000829"))
                spectra.append(
                    Spectrum(
                        precursor_mz=precursor_mz,
                        rt=rt,
                        mz=arrays["mz"],
                        intensity=arrays["intensity"],
                        isolation_lower_mz=(precursor_mz - lower_offset) if lower_offset else None,
                        isolation_upper_mz=(precursor_mz + upper_offset) if upper_offset else None,
                    )
                )
        element.clear()
    return spectra


def _mzml_cv_value(element: ET.Element, accession: str) -> str | None:
    for child in element.iter():
        if _local_name(child.tag) == "cvParam" and child.attrib.get("accession") == accession:
            return child.attrib.get("value", "")
    return None


def _mzml_scan_time(element: ET.Element) -> float | None:
    for child in element.iter():
        if _local_name(child.tag) != "cvParam" or child.attrib.get("accession") != "MS:1000016":
            continue
        value = number(child.attrib.get("value"))
        if value is None:
            return None
        return value * 60 if child.attrib.get("unitAccession") == "UO:0000031" else value
    return None


def _decode_mzml_arrays(element: ET.Element) -> dict[str, list[float]]:
    decoded: dict[str, list[float]] = {}
    for binary_array in element.iter():
        if _local_name(binary_array.tag) != "binaryDataArray":
            continue
        accessions = {
            child.attrib.get("accession")
            for child in binary_array
            if _local_name(child.tag) == "cvParam"
        }
        if "MS:1000514" in accessions:
            kind = "mz"
        elif "MS:1000515" in accessions:
            kind = "intensity"
        else:
            continue
        if accessions & {"MS:1002312", "MS:1002313", "MS:1002314"}:
            raise ValueError("MS-Numpress-compressed mzML arrays are not supported.")
        code, width = ("d", 8) if "MS:1000523" in accessions else ("f", 4)
        binary = next(
            (child for child in binary_array if _local_name(child.tag) == "binary"),
            None,
        )
        payload = base64.b64decode("" if binary is None else (binary.text or ""))
        if "MS:1000574" in accessions:
            payload = zlib.decompress(payload)
        if len(payload) % width:
            raise ValueError("Malformed mzML binary array length.")
        values = struct.unpack(f"<{len(payload) // width}{code}", payload)
        decoded[kind] = [float(value) for value in values]
    return decoded


__all__ = ["load_ms2_spectra"]
