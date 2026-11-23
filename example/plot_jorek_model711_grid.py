#!/usr/bin/env python3
"""Plot total model-711 fields on a prescribed ASDEX Upgrade grid (SI by default).

Run without arguments to read jorek000000.h5 and open a figure. Use --grid
21 31 for a quick preview, or --output total.png to save instead of displaying.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

# Allow direct execution from example/ without installing the repository.
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from plot_jorek_model711_h5 import (
    PLOT_VARIABLES,
    calculate_values_in_element,
    convert_fields_to_si,
    enable_si_units,
    evaluate_at_position,
    labeled_output_path,
    load_restart,
    make_element_search_bounds,
    nonaxisymmetric_toroidal_basis,
    plot_rectangular_fields,
    prepare_point_fields,
    select_variables,
    toroidal_basis,
)


def interpolate_prescribed_grid(restart, R_grid, Z_grid, phi_degrees=30.0,
                                include_nonaxisymmetric=False, si_units=True):
    """Return (total, nonaxisymmetric) dictionaries on equal-shaped R/Z grids.

    Nonaxisymmetric is None unless requested. Each dictionary contains all 14
    fields, in SI-derived units by default (temperatures in eV). Failed point
    searches remain NaN. A point is located only once for both components.
    This function also accepts nonuniform and curvilinear coordinate arrays.
    """
    R_grid = np.asarray(R_grid, dtype=float)
    Z_grid = np.asarray(Z_grid, dtype=float)
    if R_grid.shape != Z_grid.shape:
        raise ValueError("R_grid and Z_grid must have identical shapes")
    if not math.isfinite(phi_degrees):
        raise ValueError("phi must be finite")
    scales = enable_si_units(restart) if si_units else {}
    if not si_units:
        restart['si_units'] = False
        restart['time_scale'] = 1.0
    bounds = make_element_search_bounds(restart)
    fourier = toroidal_basis(restart['n_tor'], restart['n_period'],
                            math.radians(phi_degrees))
    prepared = prepare_point_fields(restart, fourier)
    total = {name: np.full(R_grid.shape, np.nan) for name in PLOT_VARIABLES}
    nonaxis = None
    if include_nonaxisymmetric:
        nonaxis = {name: np.full(R_grid.shape, np.nan) for name in PLOT_VARIABLES}
        nonaxis_basis = nonaxisymmetric_toroidal_basis(fourier)
        nonaxis_prepared = prepare_point_fields(restart, nonaxis_basis)

    for index in tqdm(np.ndindex(R_grid.shape), total=R_grid.size,
                      desc='Interpolating prescribed grid', unit='point'):
        result = evaluate_at_position(
            restart, float(R_grid[index]), float(Z_grid[index]), phi_degrees,
            bounds, fourier, prepared,
        )
        if result is None:
            continue
        for name, value in result['values'].items():
            total[name][index] = value
        if nonaxis is not None:
            values = calculate_values_in_element(
                restart, result['element_index'], result['s'], result['t'],
                nonaxis_basis, nonaxis_prepared, include_equilibrium=False,
            )
            for name, value in values.items():
                nonaxis[name][index] = value
    if scales:
        total = convert_fields_to_si(total, scales)
        if nonaxis is not None:
            nonaxis = convert_fields_to_si(nonaxis, scales)
    return total, nonaxis


def main(include_nonaxisymmetric=False):
    parser = argparse.ArgumentParser(description=(
        'Plot total and summed n>0 fields on a prescribed grid; SI units by default.'
        if include_nonaxisymmetric else __doc__
    ))
    parser.add_argument('file', nargs='?', type=Path, default=REPO_ROOT / 'jorek000000.h5')
    parser.add_argument('--phi', type=float, default=30.0, help='toroidal angle in degrees')
    parser.add_argument('--grid', nargs=2, type=int, default=(241, 301),
                        metavar=('NR', 'NZ'), help='grid dimensions (default: 241 301)')
    parser.add_argument('--r-range', nargs=2, type=float, default=(1.15, 2.15),
                        metavar=('RMIN', 'RMAX'), help='radial limits in metres')
    parser.add_argument('--z-range', nargs=2, type=float, default=(-0.95, 0.80),
                        metavar=('ZMIN', 'ZMAX'), help='vertical limits in metres')
    parser.add_argument('--variables', default='all', help='comma-separated fields, or all')
    parser.add_argument('--jorek-units', action='store_true',
                        help='use original JOREK normalization instead of SI-derived units')
    parser.add_argument('--output', type=Path, help='save figure(s) instead of opening windows')
    args = parser.parse_args()
    try:
        if min(args.grid) < 2:
            raise ValueError('each grid dimension must be at least 2')
        for limits in (args.r_range, args.z_range):
            if not all(math.isfinite(v) for v in limits) or limits[0] >= limits[1]:
                raise ValueError('grid limits must be finite and strictly increasing')
        if not math.isfinite(args.phi):
            raise ValueError('phi must be finite')
        names = select_variables(args.variables)
        restart = load_restart(args.file)
        R_axis = np.linspace(*args.r_range, args.grid[0])
        Z_axis = np.linspace(*args.z_range, args.grid[1])
        R_grid, Z_grid = np.meshgrid(R_axis, Z_axis, indexing='xy')
        total, nonaxis = interpolate_prescribed_grid(
            restart, R_grid, Z_grid, args.phi, include_nonaxisymmetric,
            si_units=not args.jorek_units,
        )
    except (OSError, KeyError, ValueError) as error:
        parser.error(str(error))

    print(f"{args.file}: {np.isfinite(total['rho']).sum()}/{R_grid.size} points located")
    if restart['si_units']:
        print(f"SI conversion: time scale={restart['time_scale']:.12g} "
              f"({restart['time_scale_source']}); temperatures in eV, time in ms")
    if nonaxis is not None and restart['n_tor'] == 1:
        print('This restart is axisymmetric: all n>0 fields are zero inside the mesh.')
    total_output = (labeled_output_path(args.output, 'total')
                    if include_nonaxisymmetric else args.output)
    figures = [plot_rectangular_fields(
        restart, R_axis, Z_axis, total, names, args.phi, total_output, show=False,
    )]
    if nonaxis is not None:
        figures.append(plot_rectangular_fields(
            restart, R_axis, Z_axis, nonaxis, names, args.phi,
            labeled_output_path(args.output, 'nonaxisymmetric'), show=False,
            component_label='sum of all n > 0 harmonics',
        ))
    if args.output is None:
        plt.show()
    for figure in figures:
        plt.close(figure)


if __name__ == '__main__':
    main()
