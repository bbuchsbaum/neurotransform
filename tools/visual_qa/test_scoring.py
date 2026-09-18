"""Acceptance-policy checks; run with python -m unittest discover -s tools/visual_qa."""
import unittest
import numpy as np
from report import score
from run import GATES


class ScoringTests(unittest.TestCase):
    def setUp(self):
        self.native = np.linspace(0, 1, 2000).reshape(10, 10, 20)
        self.coords = np.zeros(self.native.shape + (3,))
        self.mask = np.ones(self.native.shape, dtype=bool)

    def evaluate(self, image=None, coords=None):
        return score(self.native, self.native if image is None else image,
                     self.coords, self.coords if coords is None else coords,
                     self.mask, GATES)[0]

    def test_exact_agreement_passes(self):
        self.assertTrue(self.evaluate()['pass'])

    def test_one_invalid_sample_cannot_disappear_from_denominator(self):
        image = self.native.copy()
        image[0, 0, 0] = np.nan
        result = self.evaluate(image=image)
        self.assertFalse(result['pass'])
        self.assertEqual(result['samples'], 2000)
        self.assertEqual(result['invalid_samples'], 1)

    def test_coordinate_error_fails_despite_identical_images(self):
        coords = self.coords.copy()
        coords[0, 0, 0, 0] = .021
        self.assertFalse(self.evaluate(coords=coords)['pass'])

    def test_image_error_fails_despite_identical_coordinates(self):
        self.assertFalse(self.evaluate(image=self.native + .001)['pass'])

    def test_small_mask_cannot_receive_pass(self):
        self.mask.ravel()[999:] = False
        self.assertFalse(self.evaluate()['pass'])

    def test_outside_padding_does_not_change_fixed_interior_score(self):
        self.mask[0] = False
        image = self.native.copy()
        image[0] = np.nan
        self.assertTrue(self.evaluate(image=image)['pass'])


if __name__ == '__main__':
    unittest.main()
