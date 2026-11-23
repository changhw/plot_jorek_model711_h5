"""Check SI defaults and component separation in the runnable grid examples."""
from pathlib import Path
import unittest
import numpy as np

import plot_jorek_model711_h5 as plot
from example.plot_jorek_model711_grid import interpolate_prescribed_grid

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless((ROOT / 'jorek060000_3D.h5').exists(), 'local restart fixture absent')
class PrescribedGridTests(unittest.TestCase):
    def test_si_default_and_both_components(self):
        restart = plot.load_restart(ROOT / 'jorek060000_3D.h5')
        R = np.array([[1.8, 10.]])
        Z = np.array([[.1, 10.]])
        total, nonaxis = interpolate_prescribed_grid(
            restart, R, Z, include_nonaxisymmetric=True,
        )
        self.assertTrue(restart['si_units'])
        # Direct point evaluation returns normalized fields even when SI labels
        # are enabled; compare against independent physical conversion factors.
        point = plot.evaluate_at_position(restart, 1.8, .1, 30.)['values']
        np.testing.assert_allclose(total['rho'][0, 0], point['rho'] * 1.01e20)
        np.testing.assert_allclose(total['uR'][0, 0], point['uR'] / restart['t_norm'])
        np.testing.assert_allclose(total['T_i'][0, 0],
                                   point['T_i'] / (plot.ELEMENTARY_CHARGE * plot.MU_ZERO * 1.01e20))
        np.testing.assert_allclose(total['BR'][0, 0], point['BR'])
        self.assertEqual(nonaxis['Fprof_eq'][0, 0], 0.)
        self.assertEqual(nonaxis['psi_eq'][0, 0], 0.)
        for component in (total, nonaxis):
            for values in component.values():
                self.assertEqual(values.shape, R.shape)
                self.assertTrue(np.isnan(values[0, 1]))
        # Opting out also resets SI display metadata on a previously used restart.
        raw, raw_nonaxis = interpolate_prescribed_grid(
            restart, R, Z, include_nonaxisymmetric=True, si_units=False,
        )
        self.assertFalse(restart['si_units'])
        np.testing.assert_allclose(raw['rho'][0, 0], point['rho'])
        np.testing.assert_allclose(nonaxis['rho'][0, 0], raw_nonaxis['rho'][0, 0]*1.01e20)

    def test_axisymmetric_nonaxis_fields_are_zero(self):
        path = ROOT / 'jorek000000.h5'
        if not path.exists():
            self.skipTest('axisymmetric fixture absent')
        restart = plot.load_restart(path)
        _, nonaxis = interpolate_prescribed_grid(
            restart, np.array([1.8]), np.array([.1]), include_nonaxisymmetric=True,
        )
        for values in nonaxis.values():
            np.testing.assert_array_equal(values, [0.])


if __name__ == '__main__':
    unittest.main()
