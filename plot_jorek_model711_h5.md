# Plotting JOREK model-711 HDF5 restart files

[`plot_jorek_model711_h5.py`](plot_jorek_model711_h5.py) reads a
`jorek??????.h5` restart written by JOREK model 711, reconstructs the stored
finite-element fields at a selected toroidal angle, calculates the cylindrical
magnetic field, and plots the requested quantities in the poloidal `(R,Z)`
plane.

The script is intentionally limited to model 711. Its nine primary variables
have the fixed order

1. `A3/psi`
2. `AR`
3. `AZ`
4. `uR`
5. `uZ`
6. `u_phi`
7. `rho`
8. `T_i`
9. `T_e`

It also plots the axisymmetric `psi_eq` and `Fprof_eq` datasets and derives
`BR`, `BZ`, and `B_phi` from the magnetic potentials. In addition to the total
fields, a separate figure shows the sum of all positive toroidal harmonics.

## Running the script

The examples use the supplied ASDEX Upgrade restart `jorek000000.h5`.
Its conservative element bounds span approximately $R=1.168$–$2.100$ m and
$Z=-0.900$–$0.744$ m. These are bounding limits, not a rectangular plasma
boundary. The point example below is inside this mesh.

Plot all available fields at `phi = 0`:

```console
python plot_jorek_model711_h5.py jorek000000.h5
```

Convert density, temperatures, velocities, and displayed time to SI units:

```console
python plot_jorek_model711_h5.py jorek000000.h5 \
  --si-units --variables uR,uZ,u_phi,rho,T_i,T_e
```

This creates two independent Matplotlib figures: the total fields and the sum
of all `n > 0` harmonics. Both figures are created before one `plt.show()` call.

Plot selected fields at 30 degrees and save both component figures:

```console
python plot_jorek_model711_h5.py jorek000000.h5 \
  --phi 30 \
  --variables A3/psi,rho,T_i,T_e,BR,BZ,B_phi \
  --output jorek000000.png
```

This writes the total fields to `jorek000000.png` and the positive-harmonic sum
to `jorek000000_nonaxisymmetric.png`.

Find the element and local coordinates of one physical point and print all
field values there:

```console
python plot_jorek_model711_h5.py jorek000000.h5 \
  --phi 30 --point 1.8 0.1
```

`--point` is an additional numerical query; it does not suppress the figures.
For the supplied `jorek000000.h5`, the printed result has the form

```text
Point (R=1.8, Z=0.1): element rank=1525, s=0.573974021583, t=0.148416827563
  A3/psi: ...
  ...
  B_phi: ...
```

The element rank is one-based, as in JOREK. The internal Python element index
is zero-based but is not printed by the command-line interface.

Create separate original and uniform-grid figures and show them together:

```console
python plot_jorek_model711_h5.py jorek000000.h5 \
  --phi 30 --rectangular-grid 200 160 \
  --variables rho,T_i,T_e,BR,BZ,B_phi \
  --output jorek000000.png
```

Without `--output`, all figures are created before a single `plt.show()` call,
so their windows appear together. Rectangular mode creates four figures: total
and non-axisymmetric fields on both the original and rectangular grids. With
the output argument shown above, they are saved as
`jorek000000_original.png`, `jorek000000_rectangular.png`,
`jorek000000_nonaxisymmetric_original.png`, and
`jorek000000_nonaxisymmetric_rectangular.png`.

The command-line parameters are:

- `file`: model-711 HDF5 restart file;
- `--phi`: toroidal angle in degrees, with a default of zero;
- `--nsub`: number of sample points along each local element direction,
  including both endpoints; the default is five;
- `--variables`: comma-separated field names, or `all`;
- `--output`: output image path; without it, Matplotlib opens the figure. With
  `--rectangular-grid`, descriptive component and grid labels are inserted
  before the extension to create four output paths. Without rectangular mode,
  the total uses the requested path and `_nonaxisymmetric` is added to the
  second path;
- `--point R0 Z0`: print the containing JOREK element rank, local `(s,t)`
  coordinates, and all 14 total-field values at the requested physical
  position, in addition to producing the requested figures. A failed search
  reports that the point is outside the mesh or point location did not converge;
- `--rectangular-grid NR NZ`: create separate original element-local and
  uniform rectangular-grid figures for the total and `n > 0` components,
  displaying all four with one `plt.show()` call. `nsub` controls the original
  figures but is not used for rectangular-grid interpolation. A terminal
  progress bar reports completed `Z` rows;
- `--si-units`: convert `uR`, `uZ`, and `u_phi` to $\mathrm{m\,s^{-1}}$,
  `rho` to number density in $\mathrm{m}^{-3}$, `T_i` and `T_e` to
  $\mathrm{eV}$, and the displayed and printed time to milliseconds. The point
  query, element-local figures, rectangular-grid figures, and summed $n>0$
  figures all use the same conversion.

## Optional conversion to SI units

Current restarts contain `central_density` in units of
$10^{20}\,\mathrm{m}^{-3}$ and `central_mass` in atomic-mass units. With

$$
n_0=\mathtt{central\_density}\,10^{20}\ \mathrm{m}^{-3},
\qquad
\rho_0=n_0\,\mathtt{central\_mass}\,m_u,
$$

JOREK's time normalization is

$$
t_{\mathrm{norm}}=\sqrt{\mu_0\rho_0}.
$$

Some older files instead store `central_density` directly in
$\mathrm{m}^{-3}$. Matching JOREK's data extractor, a value greater than or
equal to $10^{15}$ is treated as the already dimensional $n_0$ and is not
multiplied by $10^{20}$.

The `--si-units` option applies

$$
\begin{aligned}
n_{\mathrm{SI}} &= \rho\,n_0 &&[\mathrm{m}^{-3}],\\
T_{i,\mathrm{SI}} &=
 \frac{T_i}{e\mu_0n_0} &&[\mathrm{eV}],\\
T_{e,\mathrm{SI}} &=
 \frac{T_e}{e\mu_0n_0} &&[\mathrm{eV}],\\
(u_R,u_Z,u_\phi)_{\mathrm{SI}} &=
 \frac{(u_R,u_Z,u_\phi)}{t_{\mathrm{norm}}}
 &&[\mathrm{m}\,\mathrm{s}^{-1}],\\
t_{\mathrm{display}} &=10^3t_{\mathrm{now}}t_{\mathrm{norm}}
 &&[\mathrm{ms}].
\end{aligned}
$$

The time normalization is selected in this order:

1. A finite positive `t_norm` stored in the restart.
2. A finite positive stored `sqrt_mu0_rho0`.
3. The formula above with $m_u=1.660539040\times10^{-27}$ kg, matching
   current JOREK's atomic mass unit.

Stored normalization preserves older restarts that used proton mass. Invalid
stored values produce a warning and are skipped. If both stored values are
valid but disagree by more than a relative tolerance of $10^{-10}$, conversion
stops with an error. Positive `central_density` is always required;
positive `central_mass` is required only for the calculated fallback.

Only the six normalized model variables above and time are multiplied by
conversion factors. `A3/psi`, `AR`, `AZ`, `psi_eq`, `Fprof_eq`, `BR`, `BZ`,
and `B_phi` are already dimensional in JOREK and remain numerically unchanged.
In SI mode their subplot and point-output labels explicitly show their existing
units:

| Fields | Unit |
| --- | --- |
| `A3/psi`, `psi_eq` | $\mathrm{Wb}=\mathrm{T\,m^2}$ |
| `AR`, `AZ`, `Fprof_eq` | $\mathrm{T\,m}$ |
| `BR`, `BZ`, `B_phi` | $\mathrm{T}$ |

The figure axes similarly use `R [m]` and `Z [m]`. For mass density consistent
with the selected normalization, use
$\rho_{\mathrm{mass,SI}}=\rho_{\mathrm{code}}t_{\mathrm{norm}}^2/\mu_0$.
The plotted `rho` remains number density.

For axisymmetric flux, $\psi=\Psi_{\mathrm{pol}}/(2\pi)$ is flux per radian.
The `Wb` label is dimensionally valid, but does not mean total poloidal flux.
No factor of $2\pi$ is applied to stored coefficients. In a non-axisymmetric
field, `A3` is the toroidal vector-potential coefficient.

## Procedural structure

The implementation uses plain functions and dictionaries rather than a custom
class. This makes the processing order visible in `main()` and allows each
calculation stage to be used independently.

The main sequence is:

```text
load_restart()
  -> make_sampling_basis()
  -> toroidal_basis()
  -> calculate_geometry()
  -> calculate_fields() for total and n>0 components
  -> make_triangles()
  -> create both plot_fields() figures
  -> plt.show()
```

When `--rectangular-grid` is supplied, the final part becomes

```text
interpolate_uniform_grid()
  -> make_element_search_bounds()
  -> prepare_point_fields()
  -> evaluate total and n>0 fields at every rectangular-grid point
  -> create two plot_fields() figures
  -> create two plot_rectangular_fields() figures
  -> plt.show()
```

The functions have the following responsibilities:

| Function | Purpose |
| --- | --- |
| `read_scalar` | Read one scalar metadata dataset |
| `load_restart` | Load all required arrays and call `validate_restart` |
| `enable_si_units` | Validate normalization metadata and calculate SI factors |
| `convert_fields_to_si` | Scale only velocity, density, and temperature fields |
| `displayed_time` | Return normalized time or physical time in milliseconds for output |
| `field_title` | Add units to converted plot and point-output labels |
| `validate_restart` | Check the model number, metadata, and array dimensions |
| `one_dimensional_basis` | Construct endpoint Hermite--Bezier polynomials |
| `degree_pairs` | Reproduce JOREK's derivative-degree ordering |
| `basis_functions` | Evaluate the 2-D basis and its local derivatives |
| `make_sampling_basis` | Create the `nsub * nsub` local sample grid |
| `toroidal_basis` | Evaluate the Fourier basis and its phi derivative |
| `nonaxisymmetric_toroidal_basis` | Remove `n=0` while retaining every `n>0` mode |
| `contract_element_basis` | Apply element sizes and contract nodal coefficients with a basis |
| `interpolate_field` | Reconstruct one field and its local/phi derivatives |
| `reconstruct_nodal_values` | Sum Fourier coefficients at the requested `phi` |
| `prepare_point_fields` | Reconstruct all nodal fields once for repeated queries |
| `interpolate_prepared_field_in_element` | Apply one element's basis to reconstructed nodal values |
| `interpolate_field_in_element` | Evaluate one field at one local point |
| `calculate_geometry` | Interpolate `R`, `Z`, their derivatives, and the Jacobian |
| `geometry_in_element` | Evaluate the coordinate mapping in one element |
| `coordinate_bernstein_basis` | Convert geometry bases to safe bounding coefficients |
| `make_element_search_bounds` | Build element boxes used to select candidates |
| `newton_local_coordinates` | Invert `(R,Z)` to `(s,t)` in one candidate element |
| `find_element` | Return the element index and local coordinates of a point |
| `physical_derivatives` | Transform `(s,t)` derivatives to `(R,Z)` |
| `calculate_fields` | Calculate stored, equilibrium, and magnetic fields |
| `calculate_values_in_element` | Calculate all fields at one known element position |
| `evaluate_at_position` | Locate and evaluate an arbitrary `(R0,Z0)` point |
| `interpolate_uniform_grid` | Locate each rectangular point once and return total and `n>0` grids |
| `make_triangles` | Create plotting triangles for every element |
| `select_variables` | Validate the requested field names |
| `separate_output_paths` | Create distinct paths for two saved figures |
| `labeled_output_path` | Add a component label to a saved-figure path |
| `plot_fields` | Build, display, or save the four-column figure |
| `plot_rectangular_fields` | Plot fields on the uniform rectangular grid |

`load_restart()` returns a `restart` dictionary containing HDF5 arrays and
metadata. Important keys include `values`, `coordinates`, `size`, `vertex`,
`psi_eq`, `fprof_eq`, `order`, `n_tor`, `n_period`, `central_density`, and
`central_mass`. Optional normalization keys are `t_norm` and `sqrt_mu0_rho0`;
optional `axis_node` flags identify nodes on collapsed magnetic-axis edges.

`calculate_geometry()` returns a separate `geometry` dictionary containing
`R`, `Z`, `R_s`, `R_t`, `Z_s`, `Z_t`, and `jacobian`. Keeping these values in
one dictionary avoids a long class definition while giving every array a
descriptive name.

Conceptually, the script performs these operations:

1. Load and validate a model-711 HDF5 restart.
2. Create `nsub * nsub` finite-element sample points.
3. Construct the poloidal and toroidal basis arrays.
4. Interpolate the axisymmetric geometry and calculate its Jacobian.
5. Reconstruct all stored toroidal harmonics at the requested `phi`.
6. Interpolate both the total and summed `n > 0` model variables; include the
   equilibrium profiles only in the total component.
7. Convert local `(s,t)` derivatives to physical `(R,Z)` derivatives.
8. Calculate `BR`, `BZ`, and `B_phi`.
9. Either triangulate the sampled elements or locate every requested
   rectangular-grid point.
10. Create all requested figures, call `plt.show()` once, and close them after
    the shared display returns.

## HDF5 arrays

HDF5 dimensions appear in reverse Fortran order when accessed through
`h5py`. The arrays used by the script therefore have these Python shapes:

| Dataset | Python shape | Meaning |
| --- | --- | --- |
| `values` | `(9, n_degrees, n_tor, n_nodes)` | Model-variable coefficients |
| `x` | `(n_dim, n_degrees, n_coord_tor, n_nodes)` | Grid-coordinate coefficients |
| `size` | `(n_degrees, 4, n_elements)` | Signed element scaling factors |
| `vertex` | `(4, n_elements)` | One-based element node indices |
| `psi_eq` | `(n_degrees, n_nodes)` | Axisymmetric equilibrium flux |
| `Fprof_eq` | `(n_degrees, n_nodes)` | Axisymmetric equilibrium F profile |

Older restart files may store `x` without the `n_coord_tor` dimension. The
script inserts a length-one dimension in that case. A grid with
`n_coord_tor > 1` is rejected because the derivative transformation assumes
axisymmetric geometry.

The one-based node indices in `vertex` are converted to zero-based Python
indices immediately after loading.

Optional `axis_node` flags have shape `(n_nodes,)` and contain `T` or `F`.
They enable the collapsed-axis safeguard during point location. Optional
scalar `t_norm` and `sqrt_mu0_rho0` datasets supply time normalization.

## Sampling inside an element

Every JOREK element uses local coordinates

$$
0 \leq s \leq 1, \qquad 0 \leq t \leq 1.
$$

For each direction, the script creates

$$
s_i=t_i=\frac{i}{n_{\mathrm{sub}}-1},
\qquad i=0,\ldots,n_{\mathrm{sub}}-1.
$$

For example, `nsub=5` gives the coordinates

$$
0,\quad 0.25,\quad 0.5,\quad 0.75,\quad 1.
$$

The Cartesian product of the two coordinate arrays produces 25 sample points
per element. Their field values are evaluated directly from the JOREK
finite-element polynomial; they are not obtained by linearly interpolating the
four vertex values.

Increasing `nsub` changes only the plotting resolution. It does not change the
finite-element solution stored in the restart.

## One-dimensional Hermite--Bezier basis

The script supports the cubic (`n_order=3`) and quintic (`n_order=5`) bases.
It constructs them from Bernstein polynomials

$$
B_k^p(q)=\binom{p}{k}q^k(1-q)^{p-k},
$$

where $p$ is `n_order` and $q$ is either $s$ or $t$. With
$r=(p+1)/2$ derivative-type degrees of freedom at each endpoint, the general
functions used by the Python implementation are

$$
\begin{aligned}
h_{Ld}(q) &= \sum_{k=d}^{r-1}\binom{k}{d}B_k^p(q),\\
h_{Rd}(q) &= \sum_{k=d}^{r-1}\binom{k}{d}B_{p-k}^p(q),
\qquad d=0,\ldots,r-1.
\end{aligned}
$$

For `n_order=3`, the endpoint functions are:

$$
\begin{aligned}
h_{L0}(q) &= (1-q)^2(1+2q)=1-3q^2+2q^3,\\
h_{L1}(q) &= 3q(1-q)^2=3q-6q^2+3q^3,\\
h_{R0}(q) &= q^2(3-2q)=3q^2-2q^3,\\
h_{R1}(q) &= 3q^2(1-q)=3q^2-3q^3.
\end{aligned}
$$

Their first derivatives are

$$
\begin{aligned}
h'_{L0}(q)&=-6q+6q^2, & h'_{L1}(q)&=3-12q+9q^2,\\
h'_{R0}(q)&= 6q-6q^2, & h'_{R1}(q)&=6q-9q^2,
\end{aligned}
$$

and their second derivatives are

$$
\begin{aligned}
h''_{L0}(q)&=-6+12q, & h''_{L1}(q)&=-12+18q,\\
h''_{R0}(q)&= 6-12q, & h''_{R1}(q)&= 6-18q.
\end{aligned}
$$

In particular, $h'_{R1}(1)=-3$; the sign is combined with JOREK's signed
element `size` coefficient. The Python code obtains all these derivatives
analytically from `Polynomial` objects, not by finite differences.

## Cubic two-dimensional basis

For `n_order=3`, every vertex has four JOREK derivative degrees:

| JOREK degree $d$ | Derivative pair $(d_s,d_t)$ | Interpretation |
| ---: | --- | --- |
| 1 | $(0,0)$ | value |
| 2 | $(1,0)$ | $s$-derivative type |
| 3 | $(0,1)$ | $t$-derivative type |
| 4 | $(1,1)$ | mixed-derivative type |

The vertex side ordering is

$$
(L,L),\quad(R,L),\quad(R,R),\quad(L,R).
$$

For a vertex $v$ whose sides are $a,b\in\{L,R\}$, the compact definition is

$$
\begin{aligned}
H_{v1}(s,t)&=h_{a0}(s)h_{b0}(t),\\
H_{v2}(s,t)&=h_{a1}(s)h_{b0}(t),\\
H_{v3}(s,t)&=h_{a0}(s)h_{b1}(t),\\
H_{v4}(s,t)&=h_{a1}(s)h_{b1}(t).
\end{aligned}
$$

The complete set of 16 cubic two-dimensional basis functions is as follows.

For vertex 1, $(L,L)$:

$$
\begin{aligned}
H_{1,1}={}&(1-s)^2(1+2s)(1-t)^2(1+2t),\\
H_{1,2}={}&3s(1-s)^2(1-t)^2(1+2t),\\
H_{1,3}={}&3t(1-s)^2(1+2s)(1-t)^2,\\
H_{1,4}={}&9st(1-s)^2(1-t)^2.
\end{aligned}
$$

For vertex 2, $(R,L)$:

$$
\begin{aligned}
H_{2,1}={}&s^2(3-2s)(1-t)^2(1+2t),\\
H_{2,2}={}&3(1-s)s^2(1-t)^2(1+2t),\\
H_{2,3}={}&3s^2(3-2s)t(1-t)^2,\\
H_{2,4}={}&9(1-s)s^2t(1-t)^2.
\end{aligned}
$$

For vertex 3, $(R,R)$:

$$
\begin{aligned}
H_{3,1}={}&s^2(3-2s)t^2(3-2t),\\
H_{3,2}={}&3(1-s)s^2t^2(3-2t),\\
H_{3,3}={}&3s^2(3-2s)(1-t)t^2,\\
H_{3,4}={}&9(1-s)s^2(1-t)t^2.
\end{aligned}
$$

For vertex 4, $(L,R)$:

$$
\begin{aligned}
H_{4,1}={}&(1-s)^2(1+2s)t^2(3-2t),\\
H_{4,2}={}&3s(1-s)^2t^2(3-2t),\\
H_{4,3}={}&3(1-s)^2(1+2s)(1-t)t^2,\\
H_{4,4}={}&9s(1-s)^2(1-t)t^2.
\end{aligned}
$$

For quintic elements, the script constructs the analogous
tensor products for the nine derivative pairs in JOREK's square-by-square
ordering.

For any tensor product $H=h_s(s)h_t(t)$, its first and second local derivatives
are

$$
\begin{aligned}
H_s    &=h'_s h_t, & H_t    &=h_s h'_t,\\
H_{ss} &=h''_s h_t,& H_{st} &=h'_s h'_t,&
H_{tt} &=h_s h''_t.
\end{aligned}
$$

The model-711 script currently evaluates $H$, $H_s$, and $H_t$. The equations
for $H_{ss}$, $H_{st}$, and $H_{tt}$ are included here to show how second
derivatives would be calculated from the same polynomial basis.

## Poloidal finite-element interpolation

A model variable $P$ inside
element $e$ is

$$
P(s,t,\phi)=
\sum_{v,d,k}
S_{dve}\,P_{dkN(v,e)}\,H_{dv}(s,t)\,T_k(\phi),
$$

where

- $v$ is an element vertex;
- $d$ is the poloidal derivative degree;
- $k$ is the toroidal harmonic index;
- $N(v,e)$ is the node at vertex $v$ of element $e$;
- $S_{dve}$ is the corresponding signed `size` value;
- $P_{dkN(v,e)}$ is the HDF5 coefficient for field $P$;
- $H_{dv}$ is the two-dimensional poloidal basis function;
- $T_k$ is the toroidal sine or cosine basis function.

Multiplication by `size` is essential. The HDF5 values contain nodal
value/derivative degrees of freedom, while `size` supplies the
element-dependent derivative scaling and orientation.

At fixed $\phi$, the local first derivatives are obtained by differentiating
only $H_{dv}$:

$$
\begin{aligned}
P_s&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,s}T_k,\\
P_t&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,t}T_k.
\end{aligned}
$$

The local second derivatives would be

$$
\begin{aligned}
P_{ss}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,ss}T_k,\\
P_{st}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,st}T_k,\\
P_{tt}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,tt}T_k.
\end{aligned}
$$

The same interpolation is applied to `R(s,t)`, `Z(s,t)`, `psi_eq`, and
`Fprof_eq`. Because the equilibrium datasets are axisymmetric, they do not
have a toroidal-harmonic dimension and can be treated as having only $T_0=1$.

## Toroidal reconstruction

The real Fourier basis stored by JOREK is

$$
T_0=1,
$$

followed by a cosine/sine pair for every positive mode. For pair number
$m=1,2,\ldots$, the physical toroidal mode is

$$
n_m=m\,n_{\mathrm{period}}.
$$

The complete field at the requested angle is

$$
P(s,t,\phi)=P_0(s,t)+
\sum_{m\geq1}\left[
P_{m,c}(s,t)\cos(n_m\phi)
+P_{m,s}(s,t)\sin(n_m\phi)
\right].
$$

The analytic toroidal derivatives are:

$$
\begin{aligned}
T'_{m,c}&=-n_m\sin(n_m\phi), &
T''_{m,c}&=-n_m^2\cos(n_m\phi),\\
T'_{m,s}&= n_m\cos(n_m\phi), &
T''_{m,s}&=-n_m^2\sin(n_m\phi).
\end{aligned}
$$

Consequently,

$$
\begin{aligned}
P_\phi&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv}T'_k,\\
P_{\phi\phi}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv}T''_k,\\
P_{s\phi}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,s}T'_k,\\
P_{t\phi}&=\sum_{v,d,k}S_{dve}P_{dkN(v,e)}H_{dv,t}T'_k.
\end{aligned}
$$

The model-711 script uses $P_\phi$ for the magnetic field. The other equations
show how the same basis data would provide second and mixed derivatives.

The separate non-axisymmetric field is

$$
P_{n>0}(s,t,\phi)=
\sum_{m\geq1}\left[
P_{m,c}(s,t)\cos(n_m\phi)
+P_{m,s}(s,t)\sin(n_m\phi)
\right].
$$

This is the sum of every stored positive harmonic, not an isolated magnetic
harmonic. For the magnetic components, the script applies the curl formulas to
the summed `n>0` potentials and omits the axisymmetric equilibrium `Fprof_eq`
term. Thus the result includes the required cosine/sine coupling of toroidal
derivatives. The non-axisymmetric components of `psi_eq` and `Fprof_eq` are
exactly zero because those datasets contain only `n=0`.

## Finding an element from `(R0,Z0)`

The `--point R0 Z0` workflow is:

1. Build one conservative physical bounding box for every finite element.
2. Select only boxes that contain `(R0,Z0)`.
3. Invert the curved coordinate mapping of each candidate from `(R,Z)` to
   local `(s,t)` using Newton iteration.
4. Accept the first candidate whose reconstructed physical position agrees
   with `(R0,Z0)` within the required tolerance.
5. Evaluate all model, equilibrium, and magnetic quantities directly from that
   element's basis functions at the recovered `(s,t)`.

No `nsub` samples or Matplotlib triangles participate in this calculation.
The point value is therefore an evaluation of the stored finite-element
polynomial, rather than an interpolation of pixels or plot vertices.

Source paths below are relative to the root of a JOREK source checkout.

The point-search implementation follows the algorithm in
`grids/grid_utils/find_RZ.f90` and the
candidate-element search in
`tools/mod_element_rtree.f90`. JOREK's
curved-element bound calculation is in
`diagnostics/RZ_minmax.f90`.

JOREK first uses an R-tree containing element bounding boxes to find elements
that might contain the requested point. It then calls `find_RZ_single` for
each candidate. The standalone Python script does not depend on JOREK's C++
R-tree. Instead, `make_element_search_bounds()` converts the polynomial
coordinate mapping to tensor-product Bernstein form. If an order-$p$
coordinate polynomial is written as

$$
P(s,t)=\sum_{i,j=0}^{p}a_{ij}s^it^j
      =\sum_{k,l=0}^{p}b_{kl}B_k^p(s)B_l^p(t),
$$

its Bernstein coefficients are

$$
b_{kl}=\sum_{i=0}^{k}\sum_{j=0}^{l}
\frac{\binom{k}{i}}{\binom{p}{i}}
\frac{\binom{l}{j}}{\binom{p}{j}}a_{ij}.
$$

Because the Bernstein basis is nonnegative and sums to one on the unit square,

$$
\min_{k,l}b_{kl}\leq P(s,t)\leq\max_{k,l}b_{kl}.
$$

Applying this to both `R` and `Z` gives a conservative box

$$
(R_{\min},R_{\max},Z_{\min},Z_{\max}).
$$

Unlike bounds obtained by sampling, this polynomial enclosure cannot miss an
extremum between samples. It may be slightly larger than the tightest possible
box, which only adds harmless extra Newton candidates. Only boxes containing
`(R0,Z0)` are passed to the nonlinear solver. These boxes are an acceleration
filter; successful membership still requires the inverse finite-element
mapping to converge inside the element.

For one candidate element, `newton_local_coordinates()` solves

$$
R(s,t)=R_0,\qquad Z(s,t)=Z_0.
$$

At the current iterate, define

$$
\Delta R=R(s,t)-R_0,\qquad
\Delta Z=Z(s,t)-Z_0.
$$

The Newton correction satisfies

$$
\begin{bmatrix}
R_s & R_t\\
Z_s & Z_t
\end{bmatrix}
\begin{bmatrix}
\Delta s\\
\Delta t
\end{bmatrix}
=-
\begin{bmatrix}
\Delta R\\
\Delta Z
\end{bmatrix}.
$$

Using $J=R_sZ_t-R_tZ_s$, the explicit correction is

$$
\begin{aligned}
\Delta s &= \frac{R_t\Delta Z-Z_t\Delta R}{J},\\
\Delta t &= \frac{Z_s\Delta R-R_s\Delta Z}{J}.
\end{aligned}
$$

As in JOREK, the Python search uses:

- five initial guesses: `(0.5,0.5)`, `(0.75,0.75)`, `(0.75,0.25)`,
  `(0.25,0.75)`, and `(0.25,0.25)`;
- at most 20 Newton updates per initial guess, evaluating the initial point
  and every updated point;
- a maximum absolute step of 0.25 in each local direction;
- clipping of every iterate to `0 <= s,t <= 1`;
- rejection of singular coordinate Jacobians.

When optional `axis_node` metadata is present, an edge with two axis nodes
is treated as collapsed. A clipped update landing there is accepted if its
residual passes; otherwise the affected coordinate is moved halfway back
toward the previous interior iterate. Without these flags, this special
safeguard is unavailable. A stagnant iterate ends the current attempt.

Like current JOREK, the default acceptance criterion uses the L1 residual:

$$
|\Delta R|+|\Delta Z|\leq10^{-13}\max(1,|R_0|,|Z_0|).
$$

An explicit positive `tolerance` argument to `find_element()` overrides this
threshold. The Newton system is scaled by its largest absolute Jacobian entry
before inversion. Unlike JOREK, this script does not implement its relaxed,
topology-dependent approximate fallback. Failure can therefore mean either
an outside point or nonconvergence; it is not proof of being outside.

`find_element()` returns a zero-based Python element index together with
`(s,t)`. The higher-level `evaluate_at_position()` also reports
`element_rank = element_index + 1`, matching JOREK's one-based element
numbering, and evaluates all primary, equilibrium, and magnetic fields using
the basis functions at that exact local position. This evaluation is
independent of `nsub`.

Programmatically, `evaluate_at_position()` returns

```python
{
    "element_rank": element_index + 1,  # one-based JOREK rank
    "element_index": element_index,      # zero-based Python index
    "s": s,
    "t": t,
    "values": {"A3/psi": ..., "AR": ..., "B_phi": ...},  # plus other fields
}
```

It returns `None` if the physical point is non-finite, outside all candidate
boxes, or cannot be recovered by Newton iteration. Optional precomputed search
bounds and Fourier arrays are used by rectangular interpolation to avoid
rebuilding the same data for every query.

## Transformation from `(s,t)` to `(R,Z)`

The interpolated geometry defines

$$
J=R_sZ_t-R_tZ_s.
$$

For any field $P$, the physical derivatives follow from the inverse geometry
Jacobian:

$$
\begin{aligned}
P_R &= \frac{P_sZ_t-P_tZ_s}{J},\\
P_Z &= \frac{P_tR_s-P_sR_t}{J}.
\end{aligned}
$$

### Second derivatives in physical coordinates

Although `plot_jorek_model711_h5.py` currently needs only first derivatives,
the second derivatives can be obtained from the same basis. First evaluate

$$
P_{ss},\ P_{st},\ P_{tt},\qquad
R_{ss},\ R_{st},\ R_{tt},\qquad
Z_{ss},\ Z_{st},\ Z_{tt}
$$

with the local derivative sums given above. The physical derivatives of the
geometry Jacobian are

$$
\begin{aligned}
\frac{\partial J}{\partial R}={1\over J}\big[{}
 &R_{ss}Z_t^2-Z_{ss}R_tZ_t-2R_{st}Z_sZ_t\\
 &+Z_{st}(R_sZ_t+R_tZ_s)+R_{tt}Z_s^2-Z_{tt}R_sZ_s
\big],\\
\frac{\partial J}{\partial Z}={1\over J}\big[{}
 &Z_{tt}R_s^2-R_{tt}Z_sR_s-2Z_{st}R_tR_s\\
 &+R_{st}(Z_tR_s+Z_sR_t)+Z_{ss}R_t^2-R_{ss}Z_tR_t
\big].
\end{aligned}
$$

The cylindrical second derivatives are then

$$
\begin{aligned}
P_{RR}={1\over J^2}\big[{}
 &P_{ss}Z_t^2-2P_{st}Z_sZ_t+P_{tt}Z_s^2\\
 &+P_s(Z_{st}Z_t-Z_{tt}Z_s)
  +P_t(Z_{st}Z_s-Z_{ss}Z_t)
\big]\\
&-{\left(\partial J/\partial R\right)\over J^2}(P_sZ_t-P_tZ_s),\\
P_{RZ}={1\over J^2}\big[{}
 &-P_{ss}Z_tR_t-P_{tt}R_sZ_s
  +P_{st}(Z_sR_t+Z_tR_s)\\
 &-P_s(R_{st}Z_t-R_{tt}Z_s)
  -P_t(R_{st}Z_s-R_{ss}Z_t)
\big]\\
&-{\left(\partial J/\partial R\right)\over J^2}(-P_sR_t+P_tR_s),\\
P_{ZZ}={1\over J^2}\big[{}
 &P_{ss}R_t^2-2P_{st}R_sR_t+P_{tt}R_s^2\\
 &+P_s(R_{st}R_t-R_{tt}R_s)
  +P_t(R_{st}R_s-R_{ss}R_t)
\big]\\
&-{\left(\partial J/\partial Z\right)\over J^2}(-P_sR_t+P_tR_s).
\end{aligned}
$$

Mixed derivatives $P_{R\phi}$ and $P_{Z\phi}$ follow by applying the same
inverse Jacobian to $P_{s\phi}$ and $P_{t\phi}$:

$$
P_{R\phi}=\frac{P_{s\phi}Z_t-P_{t\phi}Z_s}{J},\qquad
P_{Z\phi}=\frac{P_{t\phi}R_s-P_{s\phi}R_t}{J}.
$$

These second-derivative equations are documented for extension and analysis.
They are not evaluated by `plot_jorek_model711_h5.py`: the present script plots
field values and calculates the magnetic field from first derivatives, but it
does not plot current-density profiles. Calculating derivatives of the derived
magnetic components would require the corresponding higher derivatives of the
vector potentials and geometry.

JOREK uses the toroidal-angle convention

$$
x=R\cos\phi,\qquad y=-R\sin\phi,
$$

as seen in
`diagnostics/new_diag/mod_expression.f90`.
Thus positive $\phi$ points opposite to the usual mathematical cylindrical
angle. The ordered triad $(\mathbf e_R,\mathbf e_Z,\mathbf e_\phi)$ is
right-handed. With components $(B_R,B_Z,B_\phi)$, a current-density
extension starting from $\mu_0\mathbf J=\nabla\times\mathbf B$ must therefore
use

$$
\begin{aligned}
\mu_0J_R&=\partial_ZB_\phi-{1\over R}\partial_\phi B_Z,\\
\mu_0J_\phi&=\partial_RB_Z-\partial_ZB_R,\\
\mu_0J_Z&={1\over R}\partial_\phi B_R
          -{1\over R}\partial_R(RB_\phi).
\end{aligned}
$$

These signs agree directly with JOREK's diagnostic expressions
`JpolR = (R*BP_Z-BZ_p)/R`, `Jtor = BZ_R-BR_Z`, and
`JpolZ = (BR_p-R*BP_R-Btor)/R`.

Since each $B$ component already contains first derivatives of $A_3$, $A_R$,
and $A_Z$, these expressions require second derivatives of the potentials
(together with the geometry terms above). They are deliberately documentation
only in this simplified model-711 plotter.

The script identifies numerically singular mappings using the scale-independent
condition

$$
|J|\leq100\,\epsilon_{\mathrm{machine}}
\left(R_s^2+R_t^2+Z_s^2+Z_t^2\right).
$$

Derivatives at those samples are set to `NaN` rather than producing unstable
or infinite plot values.

## Magnetic-field calculation

For model 711, the magnetic potential is represented as

$$
\mathbf A=A_R\mathbf e_R+A_Z\mathbf e_Z+A_3\nabla\phi.
$$

The script follows the full-MHD JOREK equations and calculates

$$
\begin{aligned}
B_R &= \frac{\partial_Z A_3-\partial_\phi A_Z}{R},\\
B_Z &= \frac{\partial_\phi A_R-\partial_R A_3}{R},\\
B_\phi &= \partial_R A_Z-\partial_Z A_R+
          \frac{F_{\mathrm{prof,eq}}}{R}.
\end{aligned}
$$

The poloidal derivatives are obtained from the analytic finite-element
derivatives and the inverse geometry Jacobian. The toroidal derivatives are
obtained from the analytic Fourier derivatives. No finite differences are
used in either calculation. Both pointwise and batched evaluation return
`NaN` for all three magnetic components at singular mappings or `R == 0`,
where these formulas cannot be evaluated directly.

## Uniform rectangular-grid interpolation

`--rectangular-grid NR NZ` adds uniform-grid interpolation without replacing
the original element-local calculation. Its processing stages are:

1. Use all conservative element boxes to define global radial and vertical
   limits.
2. Create `NR` uniformly spaced `R` coordinates and `NZ` uniformly spaced `Z`
   coordinates across that rectangle.
3. Initialize total and non-axisymmetric arrays for every plotted quantity to
   `NaN`, each with shape `(NZ,NR)`.
4. Locate each physical grid point with the same bounding-box and Newton search
   used by `--point`.
5. If found, evaluate both the total Fourier sum and the sum of all `n > 0`
   harmonics using the recovered element and `(s,t)`. If not found, leave every
   value at that grid point as `NaN`.
6. Plot total and non-axisymmetric rectangular grids in figures separate from
   their corresponding element-local figures.

Thus `NR` controls the number of array columns and `NZ` controls the number of
rows. The rectangular grid is a new direct finite-element evaluation grid; it
is not a resampling of the original `nsub` plot.

`interpolate_uniform_grid()` creates uniform coordinate axes spanning all
element search boxes:

$$
R_i=R_{\min}+\frac{i}{N_R-1}(R_{\max}-R_{\min}),\qquad
Z_j=Z_{\min}+\frac{j}{N_Z-1}(Z_{\max}-Z_{\min}).
$$

For every Cartesian product point `(R_i,Z_j)`, it uses
`evaluate_at_position()` to obtain the containing element and exact local
coordinates. It then evaluates the non-axisymmetric fields at the same
location, avoiding a second element search.

Because point searches can take noticeable time on a fine grid, `tqdm` shows
one progress-bar update per completed `Z` row. The total is therefore `NZ`,
while every update represents `NR` newly processed grid points.

All grid points use the same toroidal angle. Therefore,
`prepare_point_fields()` first reconstructs the nodal Fourier series and its
analytic `phi` derivative once for each field. Point evaluation then applies
only the containing element's poloidal basis. This changes no mathematics; it
avoids repeating the same toroidal sum at every rectangular-grid point.

All rectangular field arrays are initialized with `NaN`. If no candidate
element contains a grid point, or Newton's method cannot converge to it, the
corresponding entries remain `NaN`. Consequently, parts of the bounding
rectangle outside the curved JOREK domain are masked in the resulting plot.

The function returns

```python
R_axis, Z_axis, total_fields, nonaxisymmetric_fields
```

Both dictionaries use the names in `PLOT_VARIABLES`. Every dictionary value
has shape

$$
(N_Z,N_R),
$$

so the first index corresponds to `Z_axis` and the second to `R_axis`.

### Interpolating onto a user-defined grid

The command-line option chooses its limits automatically. To use prescribed,
nonuniform, or curvilinear coordinates, construct two arrays of identical shape
containing the desired physical locations

$$
\widehat R_{ij},\qquad \widehat Z_{ij}.
$$

For each pair, first find $(e,s_{ij},t_{ij})$ such that

$$
R_e(s_{ij},t_{ij})=\widehat R_{ij},\qquad
Z_e(s_{ij},t_{ij})=\widehat Z_{ij},
$$

then evaluate

$$
P_{ij}=\sum_{v,d,k}
S_{dve}P_{dkN(v,e)}H_{dv}(s_{ij},t_{ij})T_k(\phi).
$$

Two runnable scripts in `example/` implement the examples below. Run the
commands shown here from the repository root:

- [`plot_jorek_model711_grid.py`](example/plot_jorek_model711_grid.py) plots total fields.
- [`plot_jorek_model711_grid_nonaxisymmetric.py`](example/plot_jorek_model711_grid_nonaxisymmetric.py)
  plots total and summed $n>0$ fields, locating each grid point only once.

Both use SI-derived units by default: velocities in m/s, number density in
m$^{-3}$, temperatures in eV, and displayed time in ms. Use `--jorek-units`
only when the original normalization is wanted. The default input is
`jorek000000.h5` in the repository root, angle is 30 degrees, and grid size is 241 by 301. That full
grid can take considerable time; start with a smaller preview:

```console
python example/plot_jorek_model711_grid.py --grid 21 31 --variables rho,T_i,T_e
python example/plot_jorek_model711_grid_nonaxisymmetric.py jorek060000_3D.h5 \
  --grid 21 31 --variables rho,BR,BZ,B_phi
```

The supplied `jorek000000.h5` is axisymmetric, so its $n>0$ fields are zero.
Use `jorek060000_3D.h5` to display nonzero perturbations. Without `--output`,
the scripts open the figure(s). To save instead:

```console
python example/plot_jorek_model711_grid.py --grid 21 31 --output total_grid.png
python example/plot_jorek_model711_grid_nonaxisymmetric.py jorek060000_3D.h5 \
  --grid 21 31 --output grid.png
```

The second command writes `grid_total.png` and `grid_nonaxisymmetric.png`.
Use `--r-range RMIN RMAX`, `--z-range ZMIN ZMAX`, `--phi DEGREES`, and
`--variables` to customize the plots. The shared `interpolate_prescribed_grid`
function also accepts equal-shaped nonuniform or curvilinear coordinate arrays
and returns total and optional non-axisymmetric field dictionaries.

The following example reuses the point-search bounds and Fourier arrays rather
than rebuilding them for every grid point. The prescribed ranges
$R=1.15$–$2.15$ m and $Z=-0.95$–$0.80$ m enclose the supplied ASDEX Upgrade
mesh with a small margin. Points outside the curved mesh remain `NaN`:

```python
import math
from pathlib import Path

import numpy as np

from plot_jorek_model711_h5 import (
    PLOT_VARIABLES,
    calculate_values_in_element,
    evaluate_at_position,
    enable_si_units,
    convert_fields_to_si,
    load_restart,
    make_element_search_bounds,
    nonaxisymmetric_toroidal_basis,
    prepare_point_fields,
    toroidal_basis,
)

restart = load_restart(Path("jorek000000.h5"))
si_scales = enable_si_units(restart)
phi_degrees = 30.0

# Example prescribed grid. R_grid and Z_grid may instead be any equal-shaped
# arrays of physical coordinates, including a nonuniform or curvilinear grid.
R_values = np.linspace(1.15, 2.15, 241)
Z_values = np.linspace(-0.95, 0.80, 301)
R_grid, Z_grid = np.meshgrid(R_values, Z_values, indexing="xy")

search_bounds = make_element_search_bounds(restart)
fourier_basis = toroidal_basis(
    restart["n_tor"], restart["n_period"], math.radians(phi_degrees)
)
prepared = prepare_point_fields(restart, fourier_basis)

fields = {
    name: np.full(R_grid.shape, np.nan) for name in PLOT_VARIABLES
}

for index in np.ndindex(R_grid.shape):
    result = evaluate_at_position(
        restart,
        float(R_grid[index]),
        float(Z_grid[index]),
        phi_degrees,
        search_bounds,
        fourier_basis,
        prepared,
    )
    if result is None:
        continue  # Outside the mesh or point location did not converge.
    for name, value in result["values"].items():
        fields[name][index] = value

fields = convert_fields_to_si(fields, si_scales)
```

To calculate the summed $n>0$ fields on the same prescribed grid, prepare a
second Fourier basis once:

```python
nonaxis_basis = nonaxisymmetric_toroidal_basis(fourier_basis)
nonaxis_prepared = prepare_point_fields(restart, nonaxis_basis)

# Inside the successful-result branch of the loop above:
nonaxis_values = calculate_values_in_element(
    restart,
    result["element_index"],
    result["s"],
    result["t"],
    nonaxis_basis,
    nonaxis_prepared,
    include_equilibrium=False,
)
nonaxis_values = convert_fields_to_si(nonaxis_values, si_scales)
```

Store `nonaxis_values` in arrays just as for `fields`. The important rule is to
search for the element only once per physical grid point and reuse the returned
element index and $(s,t)$ for every variable and harmonic aggregate. Any point
for which `evaluate_at_position()` returns `None` should remain `NaN`.

## Triangulation and rendering

Adjacent sample points define `(nsub-1) * (nsub-1)` quadrilateral patches per
element. Each patch is split into two triangles for Matplotlib. The exact
finite-element values are attached to the sampled triangle vertices, and
`tripcolor(..., shading="gouraud")` provides smooth color interpolation between
them for display.

The figure contains up to four columns. Every field has its own color scale,
and unused axes in the final row are removed.

Without `--rectangular-grid`, the script creates two independent figures: total
fields and the summed `n > 0` fields on the original `tripcolor` element-local
grid. Rectangular mode adds two `pcolormesh` figures containing the same two
components on the uniform grid, for four figures in total. It keeps every
figure open and then calls `plt.show()` once, allowing all windows to appear at
the same time. Invalid `NaN` grid points are masked, and each figure contains
at most four columns.

## Scope and limitations

- Only `jorek_model=711` with exactly nine model variables is accepted.
- Only polynomial orders three and five are supported.
- Only axisymmetric grid geometry (`n_coord_tor=1`) is supported.
- Fourier modes are consecutive multiples of `n_period`, as in JOREK `det_modes`.
- Total and summed `n > 0` fields are reconstructed at one requested toroidal
  angle; individual positive harmonics are not plotted.
- The script plots values and derived magnetic components, not general field
  derivatives.
- Rectangular-grid interpolation performs an inverse element search for every
  grid point and can therefore be considerably slower than element-local
  plotting on large meshes.
- Plot sampling is element-local, so neighboring elements retain duplicate
  samples along their common edge.

## Regression checks

Run `python -m unittest discover -s tests -v` from this directory. Tests cover
stored and fallback normalization, SI field/time conversion, final-iterate
evaluation, stagnation, outside-point rejection, and failure reporting.
`tests/test_synthetic_model711.py` independently constructs curved cubic and
quintic elements from Bernstein control points with signed element sizes, then
checks magnetic fields against analytic expressions. It also checks
collapsed-axis damping, axis-flag parsing, extreme Jacobian scaling, and
invalid/singular inputs
without any external restart files. Run those portable checks alone with
`python -m unittest discover -s tests -p test_synthetic_model711.py -v`.

When the local restart files are available, the full suite also verifies stored
normalization, 48 near-axis locations,
pointwise/batched field agreement, and zero-radius masking. Restart-based
checks are skipped if the local 3D fixture is absent.
