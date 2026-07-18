"""Dependency-free HTML renderer for shadow disagreement review queues."""

from __future__ import annotations

import html
from typing import Any


def render_shadow_disagreement_html(
    rows: list[dict[str, Any]], manifest: list[dict[str, str]]
) -> str:
    """Render the v2-versus-shadow disagreement review as deterministic HTML."""

    if not rows:
        cards = (
            '<section class="empty"><h2>No v2/ML diagnostic disagreements</h2>'
            "<p>The CSV is still present with its complete header.</p></section>"
        )
    else:
        blocks = []
        for row, asset in zip(rows, manifest, strict=True):
            compound = str(row.get("Compound_ID") or row.get("SGC ID for Component") or "unknown")
            v2 = str(row.get("comparison_baseline_status") or "")
            ml = str(row.get("ml_diagnostic_status") or "")
            probability = str(row.get("ml_same_compound_probability") or "—")
            reason = str(row.get("ml_abstention_reason") or "—")
            flags = str(row.get("ml_review_flags") or "—")
            model = str(row.get("ml_model_id") or "UNTRAINED")
            png = str(asset.get("png_asset_path") or "")
            svg = str(asset.get("svg_asset_path") or "")
            preview = png or svg
            primary = svg or png
            links = " · ".join(
                f'<a href="{html.escape(path)}">{label}</a>'
                for label, path in (("SVG", svg), ("PNG", png))
                if path
            )
            image = (
                f'<a href="{html.escape(primary)}"><img loading="lazy" '
                f'src="{html.escape(preview)}" alt="paired spectrum for {html.escape(compound)}"></a>'
                if preview and primary
                else ""
            )
            blocks.append(
                '<article class="card">'
                f"<h2>{html.escape(compound)}</h2>"
                f'<p class="transition">{html.escape(v2)} &rarr; {html.escape(ml)}</p>'
                f"<p><code>p(same)={html.escape(probability)}</code> · "
                f"model <code>{html.escape(model)}</code></p>"
                f"<p>abstention: <code>{html.escape(reason)}</code><br>"
                f"flags: <code>{html.escape(flags)}</code></p>"
                f"{image}<p>{links}</p></article>"
            )
        cards = "".join(blocks)
    return f"""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>v2 / shadow ML disagreements</title><style>body{{font:15px Segoe UI,Arial;background:#eef1f5;color:#172033;margin:24px}}header,.empty,.card{{max-width:1240px;margin:0 auto 18px;background:#fff;border:1px solid #dfe4eb;border-radius:12px;padding:18px}}.transition{{font-size:20px;font-weight:700;color:#6941c6}}img{{display:block;width:100%;height:auto;border:1px solid #dfe4eb}}code{{overflow-wrap:anywhere}}</style></head><body><header><h1>v2 / shadow ML disagreement review</h1><p>Shadow predictions are review-only and do not replace the v2 result.</p></header>{cards}</body></html>"""


__all__ = ["render_shadow_disagreement_html"]
