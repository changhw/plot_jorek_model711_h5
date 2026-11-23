"""Independent analytic checks requiring neither local restarts nor Fortran."""
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import warnings

import h5py
import numpy as np

import plot_jorek_model711_h5 as plot


def synthetic_restart(order=3):
    """Encode known polynomials via Bernstein control points, not plotter bases.

    The geometry is R=2+s+0.1*s*t, Z=t+0.1*s**2. Conversion from
    Bernstein controls uses explicit cubic/quintic endpoint transforms.
    Random signed element sizes check orientation and scaling independently.
    """
    q = np.linspace(0, 1, order+1)
    B = np.array([[math.comb(order, k)*x**k*(1-x)**(order-k)
                   for k in range(order+1)] for x in q])
    half = (order+1)//2
    E = np.zeros((order+1, order+1))
    block = np.array([[1, 0], [1, 1]]) if order == 3 else np.array(
        [[1, 0, 0], [1, 1, 0], [1, 2, 1]])
    E[:half, :half] = block
    E[half:, half:] = block[::-1]
    pairs = [(0, 0), (1, 0), (0, 1), (1, 1)]
    if order == 5:
        pairs += [(2, 0), (0, 2), (2, 1), (1, 2), (2, 2)]
    sides = [(0, 0), (1, 0), (1, 1), (0, 1)]
    rng = np.random.default_rng(711)
    sizes = rng.uniform(.3, 1.7, (len(pairs), 4, 1)) * rng.choice(
        [-1., 1.], (len(pairs), 4, 1))
    s, t = np.meshgrid(q, q, indexing='ij')
    R, Z = 2+s+.1*s*t, t+.1*s*s

    def encode(samples):
        C = np.linalg.solve(B, np.broadcast_to(samples, s.shape))
        C = np.linalg.solve(B, C.T).T
        D = np.linalg.solve(E, C)
        D = np.linalg.solve(E, D.T).T
        return np.array([[D[a*half+ds, b*half+dt] for a, b in sides]
                         for ds, dt in pairs]) / sizes[:, :, 0]

    coordinates = np.stack([encode(R), encode(Z)])[:, :, None, :]
    values = np.zeros((9, len(pairs), 5, 4))
    # A3=R*Z+0.3*R*cos(4phi), AR=Z*(1+cos(2phi)), AZ=R*(.2+sin(2phi)).
    values[0, :, 0] = encode(R*Z)
    values[0, :, 3] = encode(.3*R)
    values[1, :, 0] = encode(Z)
    values[1, :, 1] = encode(Z)
    values[2, :, 0] = encode(.2*R)
    values[2, :, 2] = encode(R)
    return dict(order=order, model=711, n_var=9, n_tor=5, n_period=2,
                coordinates=coordinates, values=values, size=sizes,
                vertex=np.arange(4)[:, None], psi_eq=encode(R*Z),
                fprof_eq=encode(3+.2*Z))


class AnalyticFieldTests(unittest.TestCase):
    def test_curved_cubic_and_quintic_fields(self):
        s, t = np.array([.13, .37, .81]), np.array([.72, .61, .19])
        R, Z = 2+s+.1*s*t, t+.1*s*s
        for order in (3, 5):
            r = synthetic_restart(order)
            basis = plot.basis_functions(order, s, t)
            for phi in (0., .37, 1.21):
                fourier = plot.toroidal_basis(5, 2, phi)
                geom = plot.calculate_geometry(r, basis, fourier)
                np.testing.assert_allclose(geom['R'][0], R, atol=2e-13, rtol=0)
                np.testing.assert_allclose(geom['Z'][0], Z, atol=2e-13, rtol=0)
                for total in (True, False):
                    with self.subTest(order=order, phi=phi, total=total):
                        fb = fourier if total else plot.nonaxisymmetric_toroidal_basis(fourier)
                        expected = {
                            'BR': -2*np.cos(2*phi) + np.zeros_like(R),
                            'BZ': (-2*Z*np.sin(2*phi)-.3*np.cos(4*phi))/R,
                            'B_phi': (np.sin(2*phi)-np.cos(2*phi)) + np.zeros_like(R),
                        }
                        if total:
                            expected['BR'] += 1
                            expected['BZ'] -= Z/R
                            expected['B_phi'] += -.8+(3+.2*Z)/R
                        result = plot.calculate_fields(r, basis, fb, geom, total)
                        for k, wanted in expected.items():
                            np.testing.assert_allclose(result[k], wanted, atol=3e-12, rtol=0)
                        for i in range(len(s)):
                            point = plot.calculate_values_in_element(
                                r, 0, s[i], t[i], fb, include_equilibrium=total)
                            for k, wanted in expected.items():
                                self.assertAlmostEqual(point[k], wanted[i], delta=3e-12)

    def test_singular_and_zero_radius_masking(self):
        for order in (3, 5):
            for failure in ('singular', 'zero_radius'):
                r = synthetic_restart(order)
                if failure == 'singular':
                    r['coordinates'][1] = 0  # Z is constant; J is zero.
                    expected_geometry = (2., 0., 1., 0., 0., 0.)
                else:
                    expected_geometry = (0., 0., 1., 0., 0., 1.)
                basis = plot.basis_functions(order, np.array([.37]), np.array([.61]))
                fb = plot.toroidal_basis(5, 2, .37)
                geom = plot.calculate_geometry(r, basis, fb)
                if failure == 'zero_radius':
                    geom['R'][:] = 0
                for total in (True, False):
                    result = plot.calculate_fields(r, basis, fb, geom, total)
                    with patch.object(plot, 'geometry_in_element', return_value=expected_geometry):
                        point = plot.calculate_values_in_element(r, 0, .37, .61, fb,
                                                                include_equilibrium=total)
                    for name in ('BR', 'BZ', 'B_phi'):
                        self.assertTrue(np.isnan(result[name]).all())
                        self.assertTrue(math.isnan(point[name]))


class AxisMetadataTests(unittest.TestCase):
    def load_flags(self, flags, dtype=None):
        r = synthetic_restart()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'synthetic.h5'
            with h5py.File(path, 'w') as h:
                for key, value in dict(values=r['values'], x=r['coordinates'],
                                       size=r['size'], vertex=r['vertex']+1,
                                       psi_eq=r['psi_eq'], Fprof_eq=r['fprof_eq'],
                                       jorek_model=711, n_var=9, n_order=3,
                                       n_tor=5, n_period=2).items():
                    h[key] = value
                if flags is not None:
                    h.create_dataset('axis_node', data=flags, dtype=dtype)
            return plot.load_restart(path)

    def test_optional_flags_and_encodings(self):
        self.assertNotIn('axis_node', self.load_flags(None))
        for flags, dtype in (
            (np.array([b'T', b'F', b'F', b'T']), None),
            (np.array([b' t ', b' f ', b' F ', b' T ']), None),
            (['t', 'f', 'F', 'T'], h5py.string_dtype('utf-8')),
            ([b'T', b'F', b'F', b'T'], h5py.string_dtype('ascii')),
        ):
            np.testing.assert_array_equal(self.load_flags(flags, dtype)['axis_node'],
                                          [True, False, False, True])

    def test_invalid_flags_and_shape(self):
        for flags in ([b'T', b'X', b'F', b'T'], [[b'T'], [b'F'], [b'F'], [b'T']]):
            with self.assertRaisesRegex(ValueError, 'axis_node'):
                self.load_flags(flags)


class PortableNewtonTests(unittest.TestCase):
    def test_collapsed_edge_damping(self):
        # A polar sector spanning 180 degrees. Starting on the opposite side
        # of the sector clips the first radial Newton step onto s=0, where J=0.
        def geometry(r, e, s, t):
            a = math.pi*t
            c, v = math.cos(a), math.sin(a)
            return 2+s*c, s*v, c, -math.pi*s*v, v, math.pi*s*c
        restart = dict(vertex=np.arange(4)[:, None],
                       axis_node=np.array([True, False, False, True]))
        R, Z, *_ = geometry(None, 0, .01, .1)
        with patch.object(plot, 'geometry_in_element', geometry), \
             patch.object(plot, 'NEWTON_STARTS', ((.5, .9),)):
            self.assertIsNone(plot.newton_local_coordinates({}, 0, R, Z, 1e-13))
            point = plot.newton_local_coordinates(restart, 0, R, Z, 1e-13)
            self.assertIsNotNone(point)
            np.testing.assert_allclose(point, [.01, .1], atol=2e-11, rtol=0)

    def test_extreme_scaling(self):
        for scale in (1e-160, 1e160):
            def geometry(r, e, s, t):
                return scale*s, scale*t, scale, 0., 0., scale
            with patch.object(plot, 'geometry_in_element', geometry):
                point = plot.newton_local_coordinates({}, 0, .6*scale, .4*scale, 1e-13*scale)
            np.testing.assert_allclose(point, [.6, .4], atol=1e-14, rtol=0)

    def test_invalid_queries_and_jacobians(self):
        bounds = np.array([[0., 1., 0., 1.]])
        for bad in (math.nan, math.inf, -math.inf):
            self.assertIsNone(plot.find_element({}, bad, .5, bounds))
            self.assertIsNone(plot.find_element({}, .5, bad, bounds))
        for bad in (0., -1., math.nan, math.inf):
            with self.assertRaises(ValueError):
                plot.find_element({}, .5, .5, bounds, bad)
        for derivatives in ((1., 1., 1., 1.), (1., math.nan, 0., 1.)):
            with patch.object(plot, 'geometry_in_element', return_value=(0., 0., *derivatives)):
                self.assertIsNone(plot.newton_local_coordinates({}, 0, .5, .5, 1e-13))


class NormalizationEdgeTests(unittest.TestCase):
    def test_precedence_and_invalid_stored_values(self):
        r = dict(central_density=1., t_norm=2e-6, sqrt_mu0_rho0=2e-6*(1+5e-11))
        plot.enable_si_units(r)
        self.assertEqual(r['time_scale'], 2e-6)
        for bad in (0., -1., math.inf, math.nan):
            with self.subTest(bad=bad), warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                r = dict(central_density=1., t_norm=2e-6, sqrt_mu0_rho0=bad)
                plot.enable_si_units(r)
                self.assertEqual(r['time_scale'], 2e-6)
                self.assertEqual(len(caught), 1)
                r = dict(central_density=1., central_mass=2., t_norm=bad, sqrt_mu0_rho0=bad)
                plot.enable_si_units(r)
                self.assertEqual(r['time_scale_source'], 'atomic-mass fallback')
                self.assertAlmostEqual(r['time_scale'], math.sqrt(plot.MU_ZERO*1e20*2*1.660539040e-27))
                self.assertEqual(len(caught), 3)

    def test_missing_and_invalid_density(self):
        for density in (None, 0., -1., math.nan, math.inf):
            r = dict(t_norm=2e-6)
            if density is not None:
                r['central_density'] = density
            with self.assertRaisesRegex(ValueError, 'central_density'):
                plot.enable_si_units(r)

    def test_all_si_scales_and_time(self):
        r = dict(central_density=2., t_norm=4e-7, time=5.)
        scales = plot.enable_si_units(r)
        converted = plot.convert_fields_to_si(dict.fromkeys(plot.PLOT_VARIABLES, 2.), scales)
        for name in ('uR', 'uZ', 'u_phi'):
            self.assertEqual(converted[name], 5e6)
        for name in ('T_i', 'T_e'):
            self.assertAlmostEqual(converted[name], 2/(1.602176565e-19*4e-7*math.pi*2e20))
        self.assertEqual(converted['rho'], 4e20)
        for name in ('A3/psi', 'AR', 'AZ', 'psi_eq', 'Fprof_eq', 'BR', 'BZ', 'B_phi'):
            self.assertEqual(converted[name], 2.)
        self.assertAlmostEqual(plot.displayed_time(r), .002)


if __name__ == '__main__':
    unittest.main()
