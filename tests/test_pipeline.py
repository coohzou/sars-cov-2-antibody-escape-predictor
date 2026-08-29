"""Basic regression tests for pipeline and predictor readiness."""

import os
import sys
import unittest

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from utils.neutralization_predictor import NeutralizationPredictor
from utils.paths import PREDICTION_DATA_DIR
from utils.sequence_analyzer import SequenceAnalyzer


class PipelineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.analyzer = SequenceAnalyzer()

    def test_ba4_detects_omicron_signature_mutations(self):
        path = os.path.join(PREDICTION_DATA_DIR, "ba4_complete.fasta")
        result = self.analyzer.analyze_sequence_file(path)
        self.assertTrue(result.get("success"), result.get("error"))
        self.assertEqual(result.get("variant"), "Omicron")
        detected = set(result.get("detected_mutations", {}).keys())
        for mutation in {"L452R", "F486V", "K417N", "N501Y", "D614G"}:
            self.assertIn(mutation, detected, f"Missing {mutation} in {sorted(detected)}")

    def test_predictor_models_available(self):
        self.assertTrue(NeutralizationPredictor.models_available())

    def test_gamma_not_alpha_like(self):
        path = os.path.join(PREDICTION_DATA_DIR, "gamma_complete.fasta")
        result = self.analyzer.analyze_sequence_file(path)
        detected = set(result.get("detected_mutations", {}).keys())
        self.assertIn("E484K", detected)
        self.assertGreaterEqual(len(detected), 4)


if __name__ == "__main__":
    unittest.main()
