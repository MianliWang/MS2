import hashlib
import unittest

from MSAI.python.ms2_ml.review_html import render_shadow_disagreement_html


class ShadowReviewHtmlTests(unittest.TestCase):
    def test_empty_report_is_deterministic_and_complete(self):
        rendered = render_shadow_disagreement_html([], [])

        self.assertEqual(rendered, render_shadow_disagreement_html([], []))
        self.assertIn("No v2/ML diagnostic disagreements", rendered)
        self.assertEqual(rendered.count("<!doctype html>"), 1)
        self.assertTrue(rendered.endswith("</body></html>"))
        self.assertEqual(
            hashlib.sha256(rendered.encode()).hexdigest(),
            "2586497aa11a36e4aeba7279149eca886d8a3fba891ffa8f1bac8b8859146452",
        )

    def test_report_preserves_existing_layout_and_escapes_untrusted_values(self):
        rows = [
            {
                "Compound_ID": 'CMP<&"1',
                "comparison_baseline_status": "conflicting_spectra",
                "ml_diagnostic_status": "supported_same_compound",
                "ml_same_compound_probability": "0.91",
                "ml_abstention_reason": '<none & "safe">',
                "ml_review_flags": "flag<a>",
                "ml_model_id": 'model"1',
            }
        ]
        manifest = [
            {
                "svg_asset_path": 'assets/svg/CMP<&"1.svg',
                "png_asset_path": 'assets/png/CMP<&"1.png',
            }
        ]

        rendered = render_shadow_disagreement_html(rows, manifest)

        expected_card = (
            '<article class="card"><h2>CMP&lt;&amp;&quot;1</h2>'
            '<p class="transition">conflicting_spectra &rarr; supported_same_compound</p>'
            "<p><code>p(same)=0.91</code> · model <code>model&quot;1</code></p>"
            "<p>abstention: <code>&lt;none &amp; &quot;safe&quot;&gt;</code><br>"
            "flags: <code>flag&lt;a&gt;</code></p>"
            '<a href="assets/svg/CMP&lt;&amp;&quot;1.svg"><img loading="lazy" '
            'src="assets/png/CMP&lt;&amp;&quot;1.png" '
            'alt="paired spectrum for CMP&lt;&amp;&quot;1"></a><p>'
            '<a href="assets/svg/CMP&lt;&amp;&quot;1.svg">SVG</a> · '
            '<a href="assets/png/CMP&lt;&amp;&quot;1.png">PNG</a></p></article>'
        )
        self.assertIn(expected_card, rendered)
        self.assertNotIn('CMP<&"1', rendered)
        self.assertEqual(
            hashlib.sha256(rendered.encode()).hexdigest(),
            "ba1ec5909b29fb6f818a7bff21cb964ba3ba068354f53808df4a5a87861c5c4e",
        )

    def test_manifest_length_mismatch_fails_closed(self):
        with self.assertRaises(ValueError):
            render_shadow_disagreement_html([{"Compound_ID": "CMP-1"}], [])


if __name__ == "__main__":
    unittest.main()
