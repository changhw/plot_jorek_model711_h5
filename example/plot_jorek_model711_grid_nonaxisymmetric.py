#!/usr/bin/env python3
"""Plot total and summed n>0 fields on the same prescribed grid, in SI by default.

This extends plot_jorek_model711_grid.py, sharing its interpolation and command
line options. Each point is located once and used for both field components.
Use jorek060000_3D.h5 to see nonzero n>0 fields; jorek000000.h5 is axisymmetric.

Example:
    python example/plot_jorek_model711_grid_nonaxisymmetric.py jorek060000_3D.h5 --grid 21 31
"""
if __package__:
    from .plot_jorek_model711_grid import main
else:
    from plot_jorek_model711_grid import main


if __name__ == '__main__':
    main(include_nonaxisymmetric=True)
