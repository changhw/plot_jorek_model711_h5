"""Numerical regression checks; run with python -m unittest discover -s tests -v."""
import contextlib
import io
import math
from pathlib import Path
import unittest
from unittest.mock import patch
import warnings

import h5py
import numpy as np

import plot_jorek_model711_h5 as plot

ROOT = Path(__file__).resolve().parents[1]


class NormalizationTests(unittest.TestCase):
    def restart(self, **extra):
        return dict(central_density=1.01, central_mass=1.0, **extra)

    def test_stored_normalization_without_mass(self):
        for name in ('t_norm', 'sqrt_mu0_rho0'):
            r = self.restart(**{name: 4.6074909698208224e-7})
            r.pop('central_mass')
            scales = plot.enable_si_units(r)
            self.assertEqual(r['time_scale'], 4.6074909698208224e-7)
            self.assertEqual(r['time_scale_source'], name)
            self.assertEqual(scales['u_phi'], 1 / r['time_scale'])
            self.assertEqual(scales['rho'], 1.01e20)

    def test_amu_fallback_and_legacy_density(self):
        r = self.restart()
        scales = plot.enable_si_units(r)
        expected = math.sqrt(plot.MU_ZERO * 1.01e20 * 1.660539040e-27)
        self.assertEqual(r['time_scale'], expected)
        legacy = dict(central_density=1.01e20, central_mass=1.0)
        self.assertEqual(scales, plot.enable_si_units(legacy))
        self.assertEqual(legacy['time_scale'], expected)

    def test_invalid_and_conflicting_metadata(self):
        with self.assertRaisesRegex(ValueError, 'inconsistent'):
            plot.enable_si_units(self.restart(t_norm=1e-6, sqrt_mu0_rho0=2e-6))
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter('always')
            r = self.restart(t_norm=float('nan'), sqrt_mu0_rho0=2e-6)
            plot.enable_si_units(r)
        self.assertEqual(len(recorded), 1)
        self.assertEqual(r['time_scale'], 2e-6)
        r = self.restart()
        r['central_mass'] = 0
        with self.assertRaisesRegex(ValueError, 'central_mass'):
            plot.enable_si_units(r)


class NewtonTests(unittest.TestCase):
    @staticmethod
    def affine_geometry(restart, element, s, t):
        return 2+s, t, 1., 0., 0., 1.

    def test_final_update_is_evaluated(self):
        with patch.object(plot, 'geometry_in_element', self.affine_geometry), \
             patch.object(plot, 'NEWTON_STARTS', ((.5, .5),)), \
             patch.object(plot, 'NEWTON_MAX_ITERATIONS', 1):
            result = plot.newton_local_coordinates({}, 0, 2.7, .7, 1e-13)
            self.assertIsNotNone(result)
            np.testing.assert_allclose(result, (.7, .7), rtol=0, atol=1e-15)

    def test_stagnation_and_outside_point(self):
        with patch.object(plot, 'geometry_in_element', side_effect=self.affine_geometry) as geometry:
            # Deliberately loose candidate box includes a target outside the element.
            bounds = np.array([[1., 4., -1., 2.]])
            self.assertIsNone(plot.find_element({}, 3.1, .5, bounds))
            self.assertLess(geometry.call_count, 30)
        with patch.object(plot, 'geometry_in_element', self.affine_geometry):
            self.assertIsNone(plot.find_element({}, 3. + 1e-11, .5, bounds))
            self.assertIsNotNone(plot.find_element({}, 3., .5, bounds))

    def test_failure_message_is_not_definitive(self):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            plot.print_point_result({}, 2., 0., None)
        self.assertIn('or point location did not converge', output.getvalue())


@unittest.skipUnless((ROOT / 'jorek060000_3D.h5').exists(), 'local restart fixture absent')
class RestartTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.restart = plot.load_restart(ROOT / 'jorek060000_3D.h5')
        cls.bounds = plot.make_element_search_bounds(cls.restart)

    def test_stored_normalization_in_supplied_files(self):
        for name in ('jorek000000.h5', 'jorek050000_2D.h5', 'jorek060000_3D.h5'):
            path = ROOT / name
            if not path.exists():
                continue
            r = plot.load_restart(path)
            self.assertEqual(np.count_nonzero(r['axis_node']), 80)
            plot.enable_si_units(r)
            with h5py.File(path, "r") as h:
                self.assertEqual(r['time_scale'], float(h['t_norm'][0]))
            self.assertEqual(r['time_scale_source'], 't_norm')
            amu_scale = math.sqrt(plot.MU_ZERO * r['central_density'] * 1e20
                                  * r['central_mass'] * plot.ATOMIC_MASS_UNIT)
            self.assertGreater(abs(r['time_scale'] / amu_scale - 1), 0.003)

    def test_near_axis_points(self):
        r = self.restart
        center = r['coordinates'][:2, 0, 0, :][:, r['axis_node']].mean(axis=1)
        for radius in (1e-9, 1e-7, 1e-5):
            for angle in np.arange(16) * 2*np.pi/16:
                with self.subTest(radius=radius, angle=angle):
                    R, Z = center + radius*np.array([np.cos(angle), np.sin(angle)])
                    located = plot.find_element(r, float(R), float(Z), self.bounds)
                    self.assertIsNotNone(located)
                    e, s, t = located
                    actual_R, actual_Z, *_ = plot.geometry_in_element(r, e, s, t)
                    self.assertLessEqual(abs(actual_R-R)+abs(actual_Z-Z),
                                         1e-13*max(1, abs(R), abs(Z)))
                    self.assertTrue(0 <= s <= 1 and 0 <= t <= 1)

    def test_point_batch_and_zero_radius(self):
        r = self.restart
        basis = plot.basis_functions(r['order'], np.array([.37]), np.array([.61]))
        fourier = plot.toroidal_basis(r['n_tor'], r['n_period'], .31)
        geometry = plot.calculate_geometry(r, basis, fourier)
        for include_eq in (True, False):
            fb = fourier if include_eq else plot.nonaxisymmetric_toroidal_basis(fourier)
            fields = plot.calculate_fields(r, basis, fb, geometry, include_eq)
            for e in (150, 3000):
                point = plot.calculate_values_in_element(r, e, .37, .61, fb,
                                                        include_equilibrium=include_eq)
                for name, value in point.items():
                    np.testing.assert_allclose(fields[name][e], value, rtol=1e-9, atol=1e-10)
        geometry['R'] = np.zeros_like(geometry['R'])
        fields = plot.calculate_fields(r, basis, fourier, geometry)
        for name in ('BR', 'BZ', 'B_phi'):
            self.assertTrue(np.isnan(fields[name]).all())


if __name__ == '__main__':
    unittest.main()
