#!/usr/bin/env python3
"""Inspect and plot fields from a JOREK model-711 HDF5 restart file.

The script plots both the total Fourier reconstruction and the sum of all
``n>0`` harmonics for the nine model variables, ``psi_eq``, ``Fprof_eq``, and
the derived magnetic field.  ``--point R0 Z0`` locates and evaluates one
physical position.  ``--rectangular-grid NR NZ`` additionally evaluates the
finite-element solution on a uniform physical grid, leaving outside points as
NaN.  Angles passed with ``--phi`` are in degrees.

Example::

    python util/plot_jorek_model711_h5.py jorek000100.h5 --phi 30
    python util/plot_jorek_model711_h5.py jorek000100.h5 --point 2.85 -0.15
    python util/plot_jorek_model711_h5.py jorek000100.h5 --rectangular-grid 200 160
    python util/plot_jorek_model711_h5.py jorek000100.h5 --si-units
    python util/plot_jorek_model711_h5.py jorek000100.h5 --output jorek000100.png
"""

from __future__ import annotations

import argparse
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from numpy.polynomial import Polynomial
from tqdm import tqdm

MODEL_VARIABLES = (
    "A3/psi",
    "AR",
    "AZ",
    "uR",
    "uZ",
    "u_phi",
    "rho",
    "T_i",
    "T_e",
)

# The first nine fields come directly from ``values``.  The equilibrium
# profiles and magnetic components are appended after interpolation.
PLOT_VARIABLES = MODEL_VARIABLES + ("psi_eq", "Fprof_eq", "BR", "BZ", "B_phi")

# These constants match JOREK's models/constants.f90.  In particular, the MHD
# normalization uses the proton mass rather than the atomic mass constant.
MU_ZERO = 4.0e-7 * math.pi
ELEMENTARY_CHARGE = 1.602176565e-19
PROTON_MASS = 1.67262178e-27

SI_BASE_UNITS = {
    "A3/psi": "Wb",
    "AR": "T m",
    "AZ": "T m",
    "uR": "m s^-1",
    "uZ": "m s^-1",
    "u_phi": "m s^-1",
    "rho": "m^-3",
    "T_i": "eV",
    "T_e": "eV",
    "psi_eq": "Wb",
    "Fprof_eq": "T m",
    "BR": "T",
    "BZ": "T",
    "B_phi": "T",
}

# Plain dictionaries keep the data flow visible without requiring users to
# understand classes.  Array and metadata keys are assigned in load_restart().
Restart = dict[str, Any]
Geometry = dict[str, np.ndarray]
PreparedFields = dict[str, Any]

# Starting points and safeguards follow grids/grid_utils/find_RZ.f90.  Several
# starts make Newton's method more reliable for curved finite elements.
NEWTON_STARTS = (
    (0.50, 0.50),
    (0.75, 0.75),
    (0.75, 0.25),
    (0.25, 0.75),
    (0.25, 0.25),
)
NEWTON_MAX_ITERATIONS = 20
NEWTON_MAX_STEP = 0.25


def read_scalar(h5: h5py.File, name: str) -> int | float:
    """Read one scalar metadata value from an open HDF5 file.

    JOREK normally writes scalars as one-item arrays.  Flattening first also
    accepts a true scalar dataset while still rejecting malformed arrays.
    """
    data = np.asarray(h5[name]).reshape(-1)
    if data.size != 1:
        raise ValueError(f"dataset {name!r} is not a scalar")
    return data[0].item()


def load_restart(filename: Path) -> Restart:
    """Load all arrays needed from one model-711 restart file.

    The HDF5 file is closed before this function returns.  The returned
    dictionary contains ordinary NumPy arrays and scalar metadata.
    """
    required = {
        "values",
        "x",
        "size",
        "vertex",
        "psi_eq",
        "Fprof_eq",
        "jorek_model",
        "n_var",
        "n_order",
        "n_tor",
        "n_period",
    }

    with h5py.File(filename, "r") as h5:
        missing = sorted(required - set(h5))
        if missing:
            raise ValueError(f"missing HDF5 datasets: {', '.join(missing)}")

        # h5py exposes reversed Fortran dimensions.  validate_restart() checks
        # the expected Python shapes after all arrays have been loaded.
        restart: Restart = {
            "filename": filename,
            "values": np.asarray(h5["values"], dtype=float),
            "coordinates": np.asarray(h5["x"], dtype=float),
            "size": np.asarray(h5["size"], dtype=float),
            "vertex": np.asarray(h5["vertex"], dtype=np.int64) - 1,
            "psi_eq": np.asarray(h5["psi_eq"], dtype=float),
            "fprof_eq": np.asarray(h5["Fprof_eq"], dtype=float),
            "model": int(read_scalar(h5, "jorek_model")),
            "n_var": int(read_scalar(h5, "n_var")),
            "order": int(read_scalar(h5, "n_order")),
            "n_tor": int(read_scalar(h5, "n_tor")),
            "n_period": int(read_scalar(h5, "n_period")),
            "step": int(read_scalar(h5, "index_now")) if "index_now" in h5 else -1,
            "time": float(read_scalar(h5, "t_now")) if "t_now" in h5 else math.nan,
            "central_density": (
                float(read_scalar(h5, "central_density"))
                if "central_density" in h5
                else math.nan
            ),
            "central_mass": (
                float(read_scalar(h5, "central_mass"))
                if "central_mass" in h5
                else math.nan
            ),
            "si_units": False,
            "time_scale": 1.0,
        }

    validate_restart(restart)
    return restart


def enable_si_units(restart: Restart) -> dict[str, float]:
    """Return model-711 field scales and configure SI time normalization.

    ``central_density`` is stored in units of 1e20 m^-3 and ``central_mass``
    is the central ion mass in proton-mass units.  Vector potentials, magnetic
    flux, equilibrium profiles, and magnetic fields deliberately get no scale.
    """
    central_density = restart["central_density"]
    central_mass = restart["central_mass"]
    if not math.isfinite(central_density) or central_density <= 0.0:
        raise ValueError(
            "--si-units requires a positive central_density dataset in the HDF5 file"
        )
    if not math.isfinite(central_mass) or central_mass <= 0.0:
        raise ValueError(
            "--si-units requires a positive central_mass dataset in the HDF5 file"
        )

    # Match JOREK's data-extractor compatibility rule: current files store a
    # factor in units of 1e20 m^-3, while some older files store absolute m^-3.
    density_scale = (
        central_density if central_density >= 1.0e15 else central_density * 1.0e20
    )
    time_scale = math.sqrt(MU_ZERO * density_scale * central_mass * PROTON_MASS)
    restart["si_units"] = True
    restart["time_scale"] = time_scale
    return {
        "uR": 1.0 / time_scale,
        "uZ": 1.0 / time_scale,
        "u_phi": 1.0 / time_scale,
        "rho": density_scale,
        "T_i": 1.0 / (ELEMENTARY_CHARGE * MU_ZERO * density_scale),
        "T_e": 1.0 / (ELEMENTARY_CHARGE * MU_ZERO * density_scale),
    }


def convert_fields_to_si(
    fields: dict[str, Any], scales: dict[str, float]
) -> dict[str, Any]:
    """Apply SI scales to convertible fields, preserving all other fields."""
    return {name: value * scales.get(name, 1.0) for name, value in fields.items()}


def displayed_time(restart: Restart) -> float:
    """Return ``t_now`` in JOREK units or milliseconds according to CLI mode."""
    factor = 1.0e3 * restart["time_scale"] if restart["si_units"] else 1.0
    return restart["time"] * factor


def field_title(restart: Restart, name: str) -> str:
    """Append converted or already-dimensional units in ``--si-units`` mode."""
    unit = SI_BASE_UNITS.get(name) if restart["si_units"] else None
    return f"{name} [{unit}]" if unit else name


def validate_restart(restart: Restart) -> None:
    """Check model metadata, supported geometry, and all array dimensions."""
    if restart["model"] != 711:
        raise ValueError(f"expected jorek_model=711, found {restart['model']}")
    if restart["n_var"] != len(MODEL_VARIABLES):
        raise ValueError(
            f"model 711 should contain 9 variables, found {restart['n_var']}"
        )
    if restart["order"] not in (3, 5):
        raise ValueError(f"unsupported n_order={restart['order']}; expected 3 or 5")
    if restart["n_tor"] < 1 or restart["n_tor"] % 2 == 0:
        raise ValueError(
            f"invalid n_tor={restart['n_tor']}; expected a positive odd number"
        )
    if restart["n_period"] < 1:
        raise ValueError(f"invalid n_period={restart['n_period']}")

    coordinates = restart["coordinates"]
    if coordinates.ndim == 3:  # Older files omit the n_coord_tor dimension.
        coordinates = coordinates[:, :, np.newaxis, :]
        restart["coordinates"] = coordinates
    if coordinates.ndim != 4 or coordinates.shape[0] < 2:
        raise ValueError(f"unexpected x shape {coordinates.shape}")
    if coordinates.shape[2] != 1:
        raise ValueError("non-axisymmetric grid geometry is not supported")

    # Expected Python shapes after h5py reverses the Fortran dimensions:
    #   values = (variable, degree, harmonic, node)
    #   x      = (coordinate, degree, coordinate harmonic, node)
    #   size   = (degree, vertex, element)
    #   vertex = (vertex, element)
    n_degrees = ((restart["order"] + 1) // 2) ** 2
    n_nodes = coordinates.shape[-1]
    expected_values = (restart["n_var"], n_degrees, restart["n_tor"], n_nodes)
    if restart["values"].shape != expected_values:
        raise ValueError(
            f"unexpected values shape {restart['values'].shape}; "
            f"expected {expected_values}"
        )
    if coordinates.shape[1] != n_degrees:
        raise ValueError(f"unexpected x shape {coordinates.shape}")

    vertex = restart["vertex"]
    if vertex.ndim != 2 or vertex.shape[0] != 4:
        raise ValueError(f"unexpected vertex shape {vertex.shape}")
    if restart["size"].shape != (n_degrees, 4, vertex.shape[1]):
        raise ValueError(f"unexpected size shape {restart['size'].shape}")
    if restart["psi_eq"].shape != (n_degrees, n_nodes):
        raise ValueError(f"unexpected psi_eq shape {restart['psi_eq'].shape}")
    if restart["fprof_eq"].shape != (n_degrees, n_nodes):
        raise ValueError(f"unexpected Fprof_eq shape {restart['fprof_eq'].shape}")
    if np.any(vertex < 0) or np.any(vertex >= n_nodes):
        raise ValueError("vertex contains an invalid node index")


@lru_cache(maxsize=2)
def one_dimensional_basis(order: int) -> list[list[Polynomial]]:
    """Build the one-dimensional JOREK Hermite--Bezier basis.

    ``order`` is the finite-element polynomial order (3 or 5).  The returned
    list is indexed as ``basis[side][degree]``, where side 0/1 denotes the
    left/right endpoint and degree denotes the nodal derivative type.
    """
    if order not in (3, 5):
        raise ValueError(f"unsupported n_order={order}; expected 3 or 5")

    q = Polynomial([0.0, 1.0])
    one_minus_q = Polynomial([1.0, -1.0])

    # Start from B_k^order(q) = C(order,k) q^k (1-q)^(order-k).
    bernstein = [
        math.comb(order, k) * q**k * one_minus_q ** (order - k)
        for k in range(order + 1)
    ]
    endpoint_degrees = (order + 1) // 2
    left: list[Polynomial] = []
    right: list[Polynomial] = []

    # Combine Bernstein control-point functions into the value/derivative
    # endpoint functions used by JOREK's basisfunctions_1D_0.
    for degree in range(endpoint_degrees):
        left.append(
            sum(
                (
                    math.comb(k, degree) * bernstein[k]
                    for k in range(degree, endpoint_degrees)
                ),
                Polynomial([0.0]),
            )
        )
        right.append(
            sum(
                (
                    math.comb(k, degree) * bernstein[order - k]
                    for k in range(degree, endpoint_degrees)
                ),
                Polynomial([0.0]),
            )
        )
    return [left, right]


def degree_pairs(order: int) -> list[tuple[int, int]]:
    """Map each HDF5 degree to its ``(s degree, t degree)`` pair.

    For cubic elements this returns ``(0,0), (1,0), (0,1), (1,1)``.
    Quintic elements add the five pairs containing a second derivative type.
    """
    pairs: list[tuple[int, int]] = []
    for square in range((order + 1) // 2):
        for other in range(square + 1):
            pairs.append((square, other))
            if square != other:
                pairs.append((other, square))
    return pairs


def basis_functions(
    order: int, s: np.ndarray, t: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate the tensor-product basis at local coordinates ``(s,t)``.

    ``s`` and ``t`` are equally sized arrays with values in ``[0,1]``.  Each
    returned array has shape ``(n_degrees, 4, n_points)`` and contains the
    basis value, its s derivative, or its t derivative.
    """
    basis_1d = one_dimensional_basis(order)
    pairs = degree_pairs(order)
    shape = (len(pairs), 4, s.size)
    value = np.empty(shape)
    derivative_s = np.empty(shape)
    derivative_t = np.empty(shape)

    # JOREK vertex order: lower-left, lower-right, upper-right, upper-left.
    vertex_sides = ((0, 0), (1, 0), (1, 1), (0, 1))
    for field_degree, (s_degree, t_degree) in enumerate(pairs):
        for vertex, (s_side, t_side) in enumerate(vertex_sides):
            hs = basis_1d[s_side][s_degree]
            ht = basis_1d[t_side][t_degree]
            # H(s,t) is a tensor product.  Its derivatives therefore require
            # differentiating only the corresponding one-dimensional factor.
            value[field_degree, vertex] = hs(s) * ht(t)
            derivative_s[field_degree, vertex] = hs.deriv()(s) * ht(t)
            derivative_t[field_degree, vertex] = hs(s) * ht.deriv()(t)
    return value, derivative_s, derivative_t


def toroidal_basis(
    n_tor: int, n_period: int, phi: float
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate JOREK's real Fourier basis at ``phi`` (in radians).

    Harmonic index 0 is axisymmetric.  Every subsequent pair stores cosine
    then sine, with physical mode number ``pair * n_period``.  The two returned
    arrays contain the basis value and its analytic phi derivative.
    """
    value = np.zeros(n_tor)
    derivative = np.zeros(n_tor)
    value[0] = 1.0
    for pair in range(1, (n_tor + 1) // 2):
        mode = pair * n_period
        phase = mode * phi
        value[2 * pair - 1] = math.cos(phase)
        value[2 * pair] = math.sin(phase)
        derivative[2 * pair - 1] = -mode * math.sin(phase)
        derivative[2 * pair] = mode * math.cos(phase)
    return value, derivative


def nonaxisymmetric_toroidal_basis(
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Remove ``n=0`` while retaining the sum of every ``n>0`` harmonic."""
    value = fourier_basis[0].copy()
    derivative = fourier_basis[1].copy()
    value[0] = 0.0
    derivative[0] = 0.0
    return value, derivative


def make_sampling_basis(order: int, nsub: int) -> tuple[np.ndarray, ...]:
    """Create the local sample grid and evaluate its poloidal basis arrays."""
    if nsub < 2:
        raise ValueError("--nsub must be at least 2")

    # Include both endpoints, matching JOREK's jorek2vtk sampling convention.
    sample = np.linspace(0.0, 1.0, nsub)
    sample_s, sample_t = np.meshgrid(sample, sample, indexing="ij")
    return basis_functions(order, sample_s.ravel(), sample_t.ravel())


def contract_element_basis(
    nodal_values: np.ndarray,
    vertex: np.ndarray,
    size: np.ndarray,
    basis: np.ndarray,
) -> np.ndarray:
    """Contract nodal coefficients with one poloidal basis array.

    ``nodal_values`` has shape ``(degree,node)`` and ``basis`` has shape
    ``(degree,vertex,point)``.  The result has shape ``(element,point)``.
    """
    # The signed size factors convert nodal derivative degrees of freedom into
    # the correctly scaled and oriented coefficients inside each element.
    element_coefficients = nodal_values[:, vertex] * size
    return np.einsum("dve,dvp->ep", element_coefficients, basis, optimize=True)


def interpolate_field(
    coefficients: np.ndarray,
    vertex: np.ndarray,
    size: np.ndarray,
    poloidal_basis: tuple[np.ndarray, ...],
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate one field at every sampled element point.

    Model variables have shape ``(degree,harmonic,node)``.  Axisymmetric
    equilibrium and geometry arrays have shape ``(degree,node)``.  The result
    is ``(value, d/ds, d/dt, d/dphi)`` with shape ``(element,point)``.
    """
    if coefficients.ndim == 2:
        # Geometry, psi_eq, and Fprof_eq have no toroidal dependence.
        nodal_values = coefficients
        nodal_phi = np.zeros_like(coefficients)
    else:
        # Reconstruct all stored Fourier harmonics at the requested phi before
        # applying the poloidal finite-element interpolation.
        nodal_values = np.einsum(
            "dhn,h->dn", coefficients, fourier_basis[0], optimize=True
        )
        nodal_phi = np.einsum(
            "dhn,h->dn", coefficients, fourier_basis[1], optimize=True
        )

    value = contract_element_basis(nodal_values, vertex, size, poloidal_basis[0])
    derivative_s = contract_element_basis(nodal_values, vertex, size, poloidal_basis[1])
    derivative_t = contract_element_basis(nodal_values, vertex, size, poloidal_basis[2])
    derivative_phi = contract_element_basis(nodal_phi, vertex, size, poloidal_basis[0])
    return value, derivative_s, derivative_t, derivative_phi


def interpolate_prepared_field_in_element(
    nodal_values: np.ndarray,
    nodal_phi: np.ndarray | None,
    element_index: int,
    vertex: np.ndarray,
    size: np.ndarray,
    poloidal_basis: tuple[np.ndarray, ...],
) -> tuple[float, float, float, float]:
    """Interpolate already reconstructed nodal values in one element."""
    nodes = vertex[:, element_index]
    element_size = size[:, :, element_index]

    def contract(values: np.ndarray, basis: np.ndarray) -> float:
        element_coefficients = values[:, nodes] * element_size
        return float(np.einsum("dv,dv->", element_coefficients, basis[:, :, 0]))

    derivative_phi = (
        0.0 if nodal_phi is None else contract(nodal_phi, poloidal_basis[0])
    )
    return (
        contract(nodal_values, poloidal_basis[0]),
        contract(nodal_values, poloidal_basis[1]),
        contract(nodal_values, poloidal_basis[2]),
        derivative_phi,
    )


def reconstruct_nodal_values(
    coefficients: np.ndarray,
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray | None]:
    """Reconstruct nodal values and phi derivatives at one toroidal angle."""
    if coefficients.ndim == 2:
        return coefficients, None
    nodal_values = np.einsum("dhn,h->dn", coefficients, fourier_basis[0], optimize=True)
    nodal_phi = np.einsum("dhn,h->dn", coefficients, fourier_basis[1], optimize=True)
    return nodal_values, nodal_phi


def prepare_point_fields(
    restart: Restart,
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> PreparedFields:
    """Reconstruct all nodal fields once for repeated point evaluations.

    Point location depends only on the axisymmetric mesh geometry.  Once a
    point has been located, its field values depend on ``phi`` through the
    Fourier reconstruction.  A rectangular grid contains many points at the
    same ``phi``, so doing this reconstruction once avoids repeating an
    element-independent calculation at every grid point.

    The returned dictionary stores ``(value, d/dphi)`` nodal arrays for the
    nine model fields and both equilibrium profiles.  It may be passed to
    :func:`evaluate_at_position` or :func:`calculate_values_in_element`.
    """
    return {
        "model": [
            reconstruct_nodal_values(variable, fourier_basis)
            for variable in restart["values"]
        ],
        "psi_eq": reconstruct_nodal_values(restart["psi_eq"], fourier_basis),
        "fprof_eq": reconstruct_nodal_values(restart["fprof_eq"], fourier_basis),
    }


def interpolate_field_in_element(
    coefficients: np.ndarray,
    element_index: int,
    vertex: np.ndarray,
    size: np.ndarray,
    poloidal_basis: tuple[np.ndarray, ...],
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> tuple[float, float, float, float]:
    """Reconstruct and interpolate one field at one local element point.

    This is the pointwise counterpart of :func:`interpolate_field`.  The
    returned tuple is ``(value, d/ds, d/dt, d/dphi)``.
    """
    nodal_values, nodal_phi = reconstruct_nodal_values(coefficients, fourier_basis)
    return interpolate_prepared_field_in_element(
        nodal_values,
        nodal_phi,
        element_index,
        vertex,
        size,
        poloidal_basis,
    )


def calculate_geometry(
    restart: Restart,
    poloidal_basis: tuple[np.ndarray, ...],
    fourier_basis: tuple[np.ndarray, np.ndarray],
) -> Geometry:
    """Interpolate ``R(s,t)``, ``Z(s,t)``, and their geometry Jacobian."""
    vertex = restart["vertex"]
    size = restart["size"]
    coordinates = restart["coordinates"]

    R, R_s, R_t, _ = interpolate_field(
        coordinates[0, :, 0, :], vertex, size, poloidal_basis, fourier_basis
    )
    Z, Z_s, Z_t, _ = interpolate_field(
        coordinates[1, :, 0, :], vertex, size, poloidal_basis, fourier_basis
    )

    # J = d(R,Z)/d(s,t).  Mark numerically singular mappings as NaN so their
    # physical derivatives cannot produce unstable or infinite plot values.
    jacobian = R_s * Z_t - R_t * Z_s
    scale = R_s**2 + R_t**2 + Z_s**2 + Z_t**2
    singular = (~np.isfinite(jacobian)) | (
        np.abs(jacobian) <= 100.0 * np.finfo(float).eps * scale
    )

    return {
        "R": R,
        "Z": Z,
        "R_s": R_s,
        "R_t": R_t,
        "Z_s": Z_s,
        "Z_t": Z_t,
        "jacobian": np.where(singular, np.nan, jacobian),
    }


def geometry_in_element(
    restart: Restart, element_index: int, s: float, t: float
) -> tuple[float, float, float, float, float, float]:
    """Evaluate the coordinate mapping inside one element.

    ``element_index`` is zero-based, while ``s`` and ``t`` lie in ``[0,1]``.
    The result is ``(R, Z, R_s, R_t, Z_s, Z_t)`` at that exact local point.
    """
    local_basis = basis_functions(restart["order"], np.asarray([s]), np.asarray([t]))
    # Geometry is axisymmetric, so these dummy arrays are never used.
    no_fourier_dependence = (np.ones(1), np.zeros(1))
    common_arguments = (
        element_index,
        restart["vertex"],
        restart["size"],
        local_basis,
        no_fourier_dependence,
    )
    R, R_s, R_t, _ = interpolate_field_in_element(
        restart["coordinates"][0, :, 0, :], *common_arguments
    )
    Z, Z_s, Z_t, _ = interpolate_field_in_element(
        restart["coordinates"][1, :, 0, :], *common_arguments
    )
    return R, Z, R_s, R_t, Z_s, Z_t


@lru_cache(maxsize=2)
def coordinate_bernstein_basis(order: int) -> np.ndarray:
    """Convert every 2-D geometry basis function to Bernstein coefficients.

    A polynomial on the unit square lies inside the range of its tensor-product
    Bernstein control coefficients.  Those coefficients therefore provide a
    conservative element box without missing extrema between sample points.
    The result has shape ``(degree, vertex, order+1, order+1)``.
    """
    basis_1d = one_dimensional_basis(order)
    pairs = degree_pairs(order)
    power_basis = np.zeros((len(pairs), 4, order + 1, order + 1))
    vertex_sides = ((0, 0), (1, 0), (1, 1), (0, 1))

    for degree, (s_degree, t_degree) in enumerate(pairs):
        for vertex_index, (s_side, t_side) in enumerate(vertex_sides):
            hs = basis_1d[s_side][s_degree]
            ht = basis_1d[t_side][t_degree]
            power_basis[
                degree,
                vertex_index,
                : hs.coef.size,
                : ht.coef.size,
            ] = np.outer(hs.coef, ht.coef)

    # If p(x) = sum_i a_i*x**i = sum_k b_k*B_k^order(x), then
    # b_k = sum_{i<=k} C(k,i)/C(order,i) * a_i.
    power_to_bernstein = np.zeros((order + 1, order + 1))
    for k in range(order + 1):
        for i in range(k + 1):
            power_to_bernstein[k, i] = math.comb(k, i) / math.comb(order, i)

    return np.einsum(
        "ki,dvij,lj->dvkl",
        power_to_bernstein,
        power_basis,
        power_to_bernstein,
        optimize=True,
    )


def make_element_search_bounds(restart: Restart) -> np.ndarray:
    """Return conservative ``Rmin,Rmax,Zmin,Zmax`` boxes for all elements.

    JOREK stores curved-element bounds in an R-tree.  This standalone script
    obtains safe polynomial bounds from Bernstein control coefficients, then
    passes only boxes containing the requested point to Newton's method.

    The returned array has shape ``(n_elements,4)`` and column order
    ``R_min, R_max, Z_min, Z_max``.  It is computed once and reused for every
    point of a rectangular grid.
    """
    coordinates = restart["coordinates"]
    vertex = restart["vertex"]
    size = restart["size"]
    bernstein_basis = coordinate_bernstein_basis(restart["order"])

    def bounds_for(coordinate: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        element_coefficients = coordinate[:, vertex] * size
        control = np.einsum(
            "dve,dvkl->ekl", element_coefficients, bernstein_basis, optimize=True
        )
        return control.min(axis=(1, 2)), control.max(axis=(1, 2))

    R_min, R_max = bounds_for(coordinates[0, :, 0, :])
    Z_min, Z_max = bounds_for(coordinates[1, :, 0, :])
    bounds = np.column_stack((R_min, R_max, Z_min, Z_max))
    # Include points on shared boundaries despite floating-point roundoff.
    span = max(float(np.ptp(bounds[:, :2])), float(np.ptp(bounds[:, 2:])), 1.0)
    padding = 100.0 * np.finfo(float).eps * span
    bounds[:, (0, 2)] -= padding
    bounds[:, (1, 3)] += padding
    return bounds


def newton_local_coordinates(
    restart: Restart,
    element_index: int,
    R0: float,
    Z0: float,
    tolerance: float,
) -> tuple[float, float] | None:
    """Solve ``R(s,t)=R0`` and ``Z(s,t)=Z0`` in one candidate element.

    The five starts, 20-iteration limit, 0.25 maximum step, and clipping to
    ``[0,1]`` follow JOREK's ``find_RZ_single`` algorithm.

    ``element_index`` is a zero-based candidate selected by its bounding box.
    ``tolerance`` is the accepted physical ``sqrt(dR**2+dZ**2)`` residual.
    The function returns ``(s,t)`` on convergence and ``None`` otherwise.
    """
    best: tuple[float, float, float] | None = None

    for initial_s, initial_t in NEWTON_STARTS:
        s = initial_s
        t = initial_t
        for _ in range(NEWTON_MAX_ITERATIONS):
            R, Z, R_s, R_t, Z_s, Z_t = geometry_in_element(restart, element_index, s, t)
            residual_R = R - R0
            residual_Z = Z - Z0
            residual = math.hypot(residual_R, residual_Z)

            if best is None or residual < best[0]:
                best = (residual, s, t)
            if residual <= tolerance:
                return s, t

            jacobian = R_s * Z_t - R_t * Z_s
            jacobian_scale = R_s**2 + R_t**2 + Z_s**2 + Z_t**2
            if (
                not math.isfinite(jacobian)
                or abs(jacobian) <= 100.0 * np.finfo(float).eps * jacobian_scale
            ):
                break

            # Solve [R_s R_t; Z_s Z_t] [ds;dt] = -[dR;dZ].
            delta_s = (R_t * residual_Z - Z_t * residual_R) / jacobian
            delta_t = (Z_s * residual_R - R_s * residual_Z) / jacobian
            delta_s = float(np.clip(delta_s, -NEWTON_MAX_STEP, NEWTON_MAX_STEP))
            delta_t = float(np.clip(delta_t, -NEWTON_MAX_STEP, NEWTON_MAX_STEP))
            s = float(np.clip(s + delta_s, 0.0, 1.0))
            t = float(np.clip(t + delta_t, 0.0, 1.0))

    if best is not None and best[0] <= tolerance:
        return best[1], best[2]
    return None


def find_element(
    restart: Restart,
    R0: float,
    Z0: float,
    search_bounds: np.ndarray,
    tolerance: float = 1.0e-10,
) -> tuple[int, float, float] | None:
    """Find the zero-based element index and local ``(s,t)`` for ``(R0,Z0)``.

    Candidate bounding boxes play the role of JOREK's element R-tree.  As in
    its standard ``find_RZ`` routine, the first converged candidate is returned.

    The return value is ``(element_index,s,t)`` with a zero-based Python
    element index.  ``None`` means that no candidate contains the point or no
    candidate's inverse mapping converged to the requested physical position.
    """
    if not math.isfinite(R0) or not math.isfinite(Z0):
        return None

    candidate_mask = (
        (search_bounds[:, 0] <= R0)
        & (R0 <= search_bounds[:, 1])
        & (search_bounds[:, 2] <= Z0)
        & (Z0 <= search_bounds[:, 3])
    )
    # A conservative box can produce extra candidates, but never proves that a
    # point is inside an element.  Newton inversion performs that final test.
    candidates = np.flatnonzero(candidate_mask)
    for element_index in candidates:
        local = newton_local_coordinates(restart, int(element_index), R0, Z0, tolerance)
        if local is not None:
            return int(element_index), local[0], local[1]
    return None


def physical_derivatives(
    derivative_s: np.ndarray,
    derivative_t: np.ndarray,
    geometry: Geometry,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform local ``(s,t)`` derivatives to cylindrical ``(R,Z)``."""
    derivative_R = (
        geometry["Z_t"] * derivative_s - geometry["Z_s"] * derivative_t
    ) / geometry["jacobian"]
    derivative_Z = (
        -geometry["R_t"] * derivative_s + geometry["R_s"] * derivative_t
    ) / geometry["jacobian"]
    return derivative_R, derivative_Z


def calculate_values_in_element(
    restart: Restart,
    element_index: int,
    s: float,
    t: float,
    fourier_basis: tuple[np.ndarray, np.ndarray],
    prepared_fields: PreparedFields | None = None,
    include_equilibrium: bool = True,
) -> dict[str, float]:
    """Calculate all 14 plotted quantities at one known ``(element,s,t)``.

    ``fourier_basis`` selects either the full Fourier sum or the ``n>0`` sum.
    Reusing ``prepared_fields`` avoids reconstructing nodal Fourier values for
    every point.  Set ``include_equilibrium=False`` for the non-axisymmetric
    calculation: this makes ``psi_eq`` and ``Fprof_eq`` zero and omits the
    equilibrium ``Fprof_eq/R`` contribution from ``B_phi``.
    """
    local_basis = basis_functions(restart["order"], np.asarray([s]), np.asarray([t]))
    if prepared_fields is None:
        prepared_fields = prepare_point_fields(restart, fourier_basis)

    common_arguments = (
        element_index,
        restart["vertex"],
        restart["size"],
        local_basis,
    )
    interpolated = [
        interpolate_prepared_field_in_element(
            nodal_values, nodal_phi, *common_arguments
        )
        for nodal_values, nodal_phi in prepared_fields["model"]
    ]
    values = {
        name: interpolated[index][0] for index, name in enumerate(MODEL_VARIABLES)
    }

    if include_equilibrium:
        psi_eq = interpolate_prepared_field_in_element(
            *prepared_fields["psi_eq"], *common_arguments
        )[0]
        fprof_eq = interpolate_prepared_field_in_element(
            *prepared_fields["fprof_eq"], *common_arguments
        )[0]
    else:
        # psi_eq and Fprof_eq contain only n=0, so their n>0 sums are zero.
        psi_eq = 0.0
        fprof_eq = 0.0
    values["psi_eq"] = psi_eq
    values["Fprof_eq"] = fprof_eq

    R, _, R_s, R_t, Z_s, Z_t = geometry_in_element(restart, element_index, s, t)
    jacobian = R_s * Z_t - R_t * Z_s

    def to_RZ(derivative_s: float, derivative_t: float) -> tuple[float, float]:
        derivative_R = (Z_t * derivative_s - Z_s * derivative_t) / jacobian
        derivative_Z = (-R_t * derivative_s + R_s * derivative_t) / jacobian
        return derivative_R, derivative_Z

    _, a3_s, a3_t, _ = interpolated[0]
    _, ar_s, ar_t, ar_phi = interpolated[1]
    _, az_s, az_t, az_phi = interpolated[2]

    jacobian_scale = R_s**2 + R_t**2 + Z_s**2 + Z_t**2
    singular = (
        not math.isfinite(jacobian)
        or abs(jacobian) <= 100.0 * np.finfo(float).eps * jacobian_scale
    )
    if singular or R == 0.0:
        values["BR"] = math.nan
        values["BZ"] = math.nan
        values["B_phi"] = math.nan
        return values

    a3_R, a3_Z = to_RZ(a3_s, a3_t)
    _, ar_Z = to_RZ(ar_s, ar_t)
    az_R, _ = to_RZ(az_s, az_t)
    values["BR"] = (a3_Z - az_phi) / R
    values["BZ"] = (ar_phi - a3_R) / R
    values["B_phi"] = az_R - ar_Z + fprof_eq / R
    return values


def evaluate_at_position(
    restart: Restart,
    R0: float,
    Z0: float,
    phi_degrees: float,
    search_bounds: np.ndarray | None = None,
    fourier_basis: tuple[np.ndarray, np.ndarray] | None = None,
    prepared_fields: PreparedFields | None = None,
) -> dict[str, Any] | None:
    """Find and evaluate a physical point ``(R0,Z0)``.

    ``R0`` and ``Z0`` are physical cylindrical coordinates and ``phi_degrees``
    is the toroidal angle in degrees.  The result contains the one-based JOREK
    ``element_rank``, zero-based Python ``element_index``, local coordinates
    ``s`` and ``t``, and a ``values`` dictionary.  ``None`` means that the
    point is outside the mesh or Newton's method did not converge in any
    candidate element.

    The final three arguments are caches for repeated queries.  If supplied,
    they must have been built from the same restart and toroidal angle.  The
    function creates them automatically for a standalone ``--point`` query.
    """
    if not math.isfinite(phi_degrees):
        raise ValueError("phi_degrees must be finite")
    if search_bounds is None:
        search_bounds = make_element_search_bounds(restart)

    # Location uses only the axisymmetric geometry.  Fourier reconstruction is
    # deliberately postponed until an element has actually been found.
    location = find_element(restart, R0, Z0, search_bounds)
    if location is None:
        return None

    if fourier_basis is None:
        fourier_basis = toroidal_basis(
            restart["n_tor"],
            restart["n_period"],
            math.radians(phi_degrees),
        )
    if prepared_fields is None:
        prepared_fields = prepare_point_fields(restart, fourier_basis)

    element_index, s, t = location
    return {
        "element_rank": element_index + 1,
        "element_index": element_index,
        "s": s,
        "t": t,
        "values": calculate_values_in_element(
            restart, element_index, s, t, fourier_basis, prepared_fields
        ),
    }


def interpolate_uniform_grid(
    restart: Restart,
    phi_degrees: float,
    n_R: int,
    n_Z: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    dict[str, np.ndarray],
    dict[str, np.ndarray],
]:
    """Interpolate total and non-axisymmetric fields onto an ``(R,Z)`` grid.

    ``n_R`` and ``n_Z`` are the numbers of uniformly spaced radial and vertical
    coordinates.  The rectangle spans the conservative element bounds, so it
    normally includes points outside the curved plasma mesh.  Such points stay
    ``NaN`` in every field array.

    The return value is ``(R_axis, Z_axis, total_fields, nonaxis_fields)``.
    Both dictionaries contain every name in :data:`PLOT_VARIABLES`; each array
    has shape ``(n_Z,n_R)``.  ``nonaxis_fields`` is the sum of all ``n>0``
    harmonics at the requested angle.  Each physical grid point is located
    only once, then both Fourier components are evaluated at the recovered
    element and local coordinates.
    """
    if n_R < 2 or n_Z < 2:
        raise ValueError("rectangular-grid dimensions must both be at least 2")

    # These conservative limits enclose every curved element.  Their global
    # rectangle also contains space outside a non-rectangular JOREK domain.
    search_bounds = make_element_search_bounds(restart)
    R_axis = np.linspace(search_bounds[:, 0].min(), search_bounds[:, 1].max(), n_R)
    Z_axis = np.linspace(search_bounds[:, 2].min(), search_bounds[:, 3].max(), n_Z)

    # NaN is the default so a failed/outside point needs no special assignment.
    # Array indexing is [Z row, R column], matching
    # pcolormesh(R_axis, Z_axis, field).
    grid_fields = {
        name: np.full((n_Z, n_R), np.nan, dtype=float) for name in PLOT_VARIABLES
    }
    nonaxisymmetric_grid_fields = {
        name: np.full((n_Z, n_R), np.nan, dtype=float) for name in PLOT_VARIABLES
    }
    # Both Fourier reconstructions are reused at every grid point.  The second
    # basis differs only by setting the n=0 coefficient to zero.
    fourier_basis = toroidal_basis(
        restart["n_tor"], restart["n_period"], math.radians(phi_degrees)
    )
    prepared_fields = prepare_point_fields(restart, fourier_basis)
    nonaxisymmetric_basis = nonaxisymmetric_toroidal_basis(fourier_basis)
    nonaxisymmetric_prepared_fields = prepare_point_fields(
        restart, nonaxisymmetric_basis
    )

    rows = tqdm(
        Z_axis,
        desc=f"Interpolating {n_R}x{n_Z} rectangular grid",
        unit="row",
    )
    for row, Z0 in enumerate(rows):
        for column, R0 in enumerate(R_axis):
            # This performs the element search and calculates the total field.
            result = evaluate_at_position(
                restart,
                float(R0),
                float(Z0),
                phi_degrees,
                search_bounds,
                fourier_basis,
                prepared_fields,
            )
            if result is None:
                continue
            for name, value in result["values"].items():
                grid_fields[name][row, column] = value

            # Reuse the recovered element and (s,t) for the n>0 calculation;
            # do not run the expensive physical-point search a second time.
            nonaxisymmetric_values = calculate_values_in_element(
                restart,
                result["element_index"],
                result["s"],
                result["t"],
                nonaxisymmetric_basis,
                nonaxisymmetric_prepared_fields,
                include_equilibrium=False,
            )
            for name, value in nonaxisymmetric_values.items():
                nonaxisymmetric_grid_fields[name][row, column] = value

    return R_axis, Z_axis, grid_fields, nonaxisymmetric_grid_fields


def calculate_fields(
    restart: Restart,
    poloidal_basis: tuple[np.ndarray, ...],
    fourier_basis: tuple[np.ndarray, np.ndarray],
    geometry: Geometry,
    include_equilibrium: bool = True,
) -> dict[str, np.ndarray]:
    """Calculate all stored and derived model-711 fields at plot samples.

    ``fourier_basis`` may represent the total or ``n>0`` Fourier sum.  For the
    latter, pass ``include_equilibrium=False`` so the axisymmetric equilibrium
    profiles and their contribution to ``B_phi`` are excluded consistently.
    """
    vertex = restart["vertex"]
    size = restart["size"]

    # Keep the derivatives because the first three variables are the magnetic
    # potentials needed below.  Each tuple is (value, d/ds, d/dt, d/dphi).
    interpolated = [
        interpolate_field(variable, vertex, size, poloidal_basis, fourier_basis)
        for variable in restart["values"]
    ]
    fields = {
        name: interpolated[index][0].ravel()
        for index, name in enumerate(MODEL_VARIABLES)
    }

    if include_equilibrium:
        psi_eq = interpolate_field(
            restart["psi_eq"], vertex, size, poloidal_basis, fourier_basis
        )[0]
        fprof_eq = interpolate_field(
            restart["fprof_eq"], vertex, size, poloidal_basis, fourier_basis
        )[0]
    else:
        # These two equilibrium datasets have no n>0 Fourier coefficients.
        psi_eq = np.zeros_like(geometry["R"])
        fprof_eq = np.zeros_like(geometry["R"])
    fields["psi_eq"] = psi_eq.ravel()
    fields["Fprof_eq"] = fprof_eq.ravel()

    _, a3_s, a3_t, _ = interpolated[0]
    _, ar_s, ar_t, ar_phi = interpolated[1]
    _, az_s, az_t, az_phi = interpolated[2]
    a3_R, a3_Z = physical_derivatives(a3_s, a3_t, geometry)
    _, ar_Z = physical_derivatives(ar_s, ar_t, geometry)
    az_R, _ = physical_derivatives(az_s, az_t, geometry)

    # Model 711 uses A = AR e_R + AZ e_Z + A3 grad(phi).  The following
    # expressions are curl(A), with the equilibrium F profile added to B_phi.
    radius = geometry["R"]
    with np.errstate(divide="ignore", invalid="ignore"):
        fields["BR"] = ((a3_Z - az_phi) / radius).ravel()
        fields["BZ"] = ((ar_phi - a3_R) / radius).ravel()
        fields["B_phi"] = (az_R - ar_Z + fprof_eq / radius).ravel()
    return fields


def make_triangles(nsub: int, n_elements: int) -> np.ndarray:
    """Split every sampled element patch into two plotting triangles."""
    local_triangles = []
    for i in range(nsub - 1):
        for j in range(nsub - 1):
            lower_left = i * nsub + j
            lower_right = (i + 1) * nsub + j
            upper_right = lower_right + 1
            upper_left = lower_left + 1
            local_triangles.extend(
                (
                    (lower_left, lower_right, upper_right),
                    (lower_left, upper_right, upper_left),
                )
            )

    local_triangles = np.asarray(local_triangles, dtype=np.int64)
    offsets = (np.arange(n_elements, dtype=np.int64) * nsub * nsub)[:, None, None]
    return (local_triangles[None, :, :] + offsets).reshape(-1, 3)


def select_variables(text: str) -> list[str]:
    """Parse ``--variables``, validate names, and remove duplicates."""
    if text.strip().lower() == "all":
        return list(PLOT_VARIABLES)
    names = [name.strip() for name in text.split(",") if name.strip()]
    unknown = [name for name in names if name not in PLOT_VARIABLES]
    if unknown:
        raise ValueError(
            f"unknown variables: {', '.join(unknown)}; "
            f"choose from {', '.join(PLOT_VARIABLES)}"
        )
    if not names:
        raise ValueError("--variables must not be empty")
    return list(dict.fromkeys(names))


def separate_output_paths(output: Path | None) -> tuple[Path | None, Path | None]:
    """Add ``_original`` and ``_rectangular`` to a requested output path."""
    if output is None:
        return None, None
    return (
        output.with_name(f"{output.stem}_original{output.suffix}"),
        output.with_name(f"{output.stem}_rectangular{output.suffix}"),
    )


def labeled_output_path(output: Path | None, label: str) -> Path | None:
    """Insert a component label before an optional output path's extension."""
    if output is None:
        return None
    return output.with_name(f"{output.stem}_{label}{output.suffix}")


def plot_fields(
    restart: Restart,
    triangulation: mtri.Triangulation,
    fields: dict[str, np.ndarray],
    names: list[str],
    phi_degrees: float,
    output: Path | None,
    show: bool = True,
    component_label: str = "total field",
) -> Any:
    """Build one element-local figure with at most four subplot columns.

    With ``show=False``, the figure remains open and is returned without a
    call to ``plt.show()``.  This lets ``main()`` create every total and
    non-axisymmetric figure before displaying all windows together.
    """
    columns = min(4, len(names))
    rows = math.ceil(len(names) / columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(5.0 * columns, 4.0 * rows), squeeze=False
    )

    for axis, name in zip(axes.ravel(), names):
        image = axis.tripcolor(
            triangulation, fields[name], shading="gouraud", cmap="RdBu_r"
        )
        axis.set_title(field_title(restart, name))
        axis.set_xlabel("R [m]" if restart["si_units"] else "R")
        axis.set_ylabel("Z [m]" if restart["si_units"] else "Z")
        axis.set_aspect("equal", adjustable="box")
        figure.colorbar(image, ax=axis, shrink=0.85)
    for axis in axes.ravel()[len(names) :]:
        axis.remove()

    figure.suptitle(
        f"{restart['filename'].name}, {component_label}, element-local grid, "
        f"phi={phi_degrees:g} deg, step={restart['step']}, "
        f"t={displayed_time(restart):.8g}{' [ms]' if restart['si_units'] else ''}"
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=300, bbox_inches="tight")
        print(output)
    elif show:
        plt.show()
    if show:
        plt.close(figure)
    return figure


def plot_rectangular_fields(
    restart: Restart,
    R_axis: np.ndarray,
    Z_axis: np.ndarray,
    fields: dict[str, np.ndarray],
    names: list[str],
    phi_degrees: float,
    output: Path | None,
    show: bool = True,
    component_label: str = "total field",
) -> Any:
    """Build one uniform-grid figure from ``(n_Z,n_R)`` field arrays.

    ``R_axis`` labels array columns and ``Z_axis`` labels rows.  NaNs outside
    the finite-element mesh are masked.  With ``show=False``, the returned
    figure stays open for the shared ``plt.show()`` call in ``main()``.
    """
    columns = min(4, len(names))
    rows = math.ceil(len(names) / columns)
    figure, axes = plt.subplots(
        rows, columns, figsize=(5.0 * columns, 4.0 * rows), squeeze=False
    )

    for axis, name in zip(axes.ravel(), names):
        image = axis.pcolormesh(
            R_axis,
            Z_axis,
            np.ma.masked_invalid(fields[name]),
            shading="auto",
            cmap="RdBu_r",
        )
        axis.set_title(field_title(restart, name))
        axis.set_xlabel("R [m]" if restart["si_units"] else "R")
        axis.set_ylabel("Z [m]" if restart["si_units"] else "Z")
        axis.set_aspect("equal", adjustable="box")
        figure.colorbar(image, ax=axis, shrink=0.85)
    for axis in axes.ravel()[len(names) :]:
        axis.remove()

    figure.suptitle(
        f"{restart['filename'].name}, {component_label}, rectangular grid, "
        f"phi={phi_degrees:g} deg, step={restart['step']}, "
        f"t={displayed_time(restart):.8g}{' [ms]' if restart['si_units'] else ''}"
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=300, bbox_inches="tight")
        print(output)
    elif show:
        plt.show()
    if show:
        plt.close(figure)
    return figure


def print_point_result(
    restart: Restart, R0: float, Z0: float, result: dict[str, Any] | None
) -> None:
    """Print one ``--point`` result, or report that the point was not found."""
    if result is None:
        print(f"Point (R={R0:.12g}, Z={Z0:.12g}) is outside the JOREK mesh")
        return

    print(
        f"Point (R={R0:.12g}, Z={Z0:.12g}): "
        f"element rank={result['element_rank']}, "
        f"s={result['s']:.12g}, t={result['t']:.12g}"
    )
    for name in PLOT_VARIABLES:
        print(f"  {field_title(restart, name)}: {result['values'][name]:.16g}")


def main() -> None:
    """Run the visible load -> interpolate -> calculate -> plot workflow."""
    # Preserve the multi-line feature introduction and examples in --help.
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("file", type=Path, help="model-711 jorek??????.h5 file")
    parser.add_argument(
        "--phi", type=float, default=0.0, help="toroidal angle in degrees"
    )
    parser.add_argument(
        "--nsub", type=int, default=5, help="sample points per element direction"
    )
    parser.add_argument(
        "--variables",
        default="all",
        help="comma-separated field names, or all (default: all)",
    )
    parser.add_argument(
        "--output", type=Path, help="save image(s) instead of opening windows"
    )
    parser.add_argument(
        "--point",
        nargs=2,
        type=float,
        metavar=("R0", "Z0"),
        help=(
            "locate physical (R0,Z0), then print its one-based element rank, "
            "local (s,t), and total field values"
        ),
    )
    parser.add_argument(
        "--rectangular-grid",
        nargs=2,
        type=int,
        metavar=("NR", "NZ"),
        help=(
            "also interpolate total and n>0 fields onto a uniform NR-by-NZ "
            "physical grid (outside points are NaN)"
        ),
    )
    parser.add_argument(
        "--si-units",
        action="store_true",
        help="convert rho, T_i/T_e, uR/uZ/u_phi, and time to SI-derived units",
    )
    args = parser.parse_args()

    try:
        if not math.isfinite(args.phi):
            raise ValueError("--phi must be finite")

        # The complete calculation is deliberately written as a sequence of
        # plain function calls so each processing stage can be read separately.
        names = select_variables(args.variables)
        restart = load_restart(args.file)
        si_scales = enable_si_units(restart) if args.si_units else {}
        n_elements = restart["vertex"].shape[1]

        if args.point:
            R0, Z0 = args.point
            point_result = evaluate_at_position(restart, R0, Z0, args.phi)
            if point_result is not None and si_scales:
                point_result["values"] = convert_fields_to_si(
                    point_result["values"], si_scales
                )
            print_point_result(restart, R0, Z0, point_result)

        # The element-local grid is always calculated.  Its total and summed
        # n>0 fields form two independent figures.  Rectangular mode adds two
        # more figures evaluated at physical points on a uniform (R,Z) grid.
        poloidal_basis = make_sampling_basis(restart["order"], args.nsub)
        phi_radians = math.radians(args.phi)
        fourier_basis = toroidal_basis(
            restart["n_tor"], restart["n_period"], phi_radians
        )
        geometry = calculate_geometry(restart, poloidal_basis, fourier_basis)
        original_fields = calculate_fields(
            restart, poloidal_basis, fourier_basis, geometry
        )
        nonaxisymmetric_basis = nonaxisymmetric_toroidal_basis(fourier_basis)
        nonaxisymmetric_original_fields = calculate_fields(
            restart,
            poloidal_basis,
            nonaxisymmetric_basis,
            geometry,
            include_equilibrium=False,
        )
        if si_scales:
            original_fields = convert_fields_to_si(original_fields, si_scales)
            nonaxisymmetric_original_fields = convert_fields_to_si(
                nonaxisymmetric_original_fields, si_scales
            )
        triangles = make_triangles(args.nsub, n_elements)
        triangulation = mtri.Triangulation(
            geometry["R"].ravel(), geometry["Z"].ravel(), triangles
        )

        if args.rectangular_grid:
            n_R, n_Z = args.rectangular_grid
            (
                R_axis,
                Z_axis,
                rectangular_fields,
                nonaxisymmetric_rectangular_fields,
            ) = interpolate_uniform_grid(restart, args.phi, n_R, n_Z)
            if si_scales:
                rectangular_fields = convert_fields_to_si(
                    rectangular_fields, si_scales
                )
                nonaxisymmetric_rectangular_fields = convert_fields_to_si(
                    nonaxisymmetric_rectangular_fields, si_scales
                )
    except (OSError, KeyError, ValueError) as error:
        parser.error(str(error))

    print(
        f"{args.file}: model=711, variables=9, harmonics={restart['n_tor']}, "
        f"n_period={restart['n_period']}, elements={n_elements}, "
        f"t={displayed_time(restart):.8g}{' [ms]' if restart['si_units'] else ''}"
    )
    figures = []
    if args.rectangular_grid:
        original_output, rectangular_output = separate_output_paths(args.output)
        figures.append(
            plot_fields(
                restart,
                triangulation,
                original_fields,
                names,
                args.phi,
                original_output,
                show=False,
            )
        )
        figures.append(
            plot_rectangular_fields(
                restart,
                R_axis,
                Z_axis,
                rectangular_fields,
                names,
                args.phi,
                rectangular_output,
                show=False,
            )
        )
        figures.append(
            plot_fields(
                restart,
                triangulation,
                nonaxisymmetric_original_fields,
                names,
                args.phi,
                labeled_output_path(args.output, "nonaxisymmetric_original"),
                show=False,
                component_label="sum of all n > 0 harmonics",
            )
        )
        figures.append(
            plot_rectangular_fields(
                restart,
                R_axis,
                Z_axis,
                nonaxisymmetric_rectangular_fields,
                names,
                args.phi,
                labeled_output_path(args.output, "nonaxisymmetric_rectangular"),
                show=False,
                component_label="sum of all n > 0 harmonics",
            )
        )
    else:
        figures.append(
            plot_fields(
                restart,
                triangulation,
                original_fields,
                names,
                args.phi,
                args.output,
                show=False,
            )
        )
        figures.append(
            plot_fields(
                restart,
                triangulation,
                nonaxisymmetric_original_fields,
                names,
                args.phi,
                labeled_output_path(args.output, "nonaxisymmetric"),
                show=False,
                component_label="sum of all n > 0 harmonics",
            )
        )

    # Keeping every figure open until here is important: one plt.show() call
    # displays all two (normal mode) or four (rectangular mode) GUI windows at
    # the same time.  Saved figures are also closed together below.
    if args.output is None:
        plt.show()
    for figure in figures:
        plt.close(figure)


if __name__ == "__main__":
    main()
