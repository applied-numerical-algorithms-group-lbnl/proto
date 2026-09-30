# Master / cubedSphereDev integration

This branch combines `origin/master` at `5418a7b9003c16936ba5209a29eb7464cb40d079`
and `origin/cubedSphereDev` at `1ef022043bfaba986df937ff1df35133787b1904`.
Both remote tips were fetched on 2026-09-30. The local `master` and
`cubedSphereDev` refs have not been moved. Nothing has been pushed.

## Interface decisions

- Keep master's AMR, multiblock topology, map-type interfaces, composite operators,
  generic interpolation, tests, and archived examples.
- Port the cubed-sphere MHD geometry, configurable shell radii/stretch, spherical
  boundary interpolation, flux-register accumulation/reflux, `exchange0`,
  `coarsenTo`, checkpoint iteration metadata, and Euler/RK3/RK4 stage hooks.
- General `CubedSphereShell::InterpOp` retains master's requested output ghost
  width. `CubedSphereShell::MHDInterpOp` retains cubedSphereDev's angular stencil
  and its two extra source ghost layers. The distinction is explicit because
  the two algorithms have different caller contracts.
- Keep master's default linear shell (radii 0.9 to 1.0, radial direction 2).
  MHD applications explicitly configure their radii, radial direction 0, and
  radial stretch. A zero stretch selects the linear map without dividing by zero.
- Keep the generic `MBLevelRK4::advance(state, dt, time)` overload alongside the
  MHD overload with radial volumes and temporal order 1, 3, or 4. The generic
  RK4 now uses the intermediate state in stages 3 and 4. A regression test checks
  its result against the RK4 stability polynomial for y'=y.
- Preserve `PR_OPS` / CMake `OPS` as compatibility aliases for LAPACK.
- Skip zero-width multiblock boundary allocations, retain reading checkpoints
  without iteration metadata (iteration defaults to zero), remove stale calls
  to the now-private `MBGraph::fixRotations`, and gate CUDA-only performance
  targets on CUDA builds.

The active `CubedSphereMHDCons` example is updated to the current interfaces and
registered under that name in CMake. Its viscosity stage hook is aligned with
cubedSphereDev's integrator. Master's archived examples remain archived.

## HelioCubed migration

`HelioCubed-CubedSphere.patch` applies to HelioCubed's `CubedSphere` commit
`455a495`. It was checked against a clean export and applied only in a temporary
validation copy. The actual HelioCubed checkout and branches were not modified.

The patch changes:

1. Map template parameters from `template<MemType> typename MAP` to `typename MAP`.
2. Concrete map arguments to `MBMap_CubedSphereShell<HOST>` or `<MEM>`.
3. The MHD interpolation factory call to `CubedSphereShell::MHDInterpOp`.
4. Build flags/includes from OPS to LAPACK, preserving explicit MHD shell parameters.

When ready, apply it on a new HelioCubed branch based on `CubedSphere`:

```sh
git switch -c integrate/proto-master CubedSphere
git apply --check /path/to/proto/integration/HelioCubed-CubedSphere.patch
git apply /path/to/proto/integration/HelioCubed-CubedSphere.patch
```

Then set `PROTO`, compiler, HDF5 and LAPACK locations for the target machine and
rebuild all application objects. The patch is for the CubedSphere solver;
HelioCubed's FullSphere and independent Unstructured solvers were not migrated.

## Reproduce the Proto build

Initialize the existing BLT dependency, then configure a CPU/MPI/HDF5 build:

```sh
git submodule update --init blt
cmake -S . -B build/integration-3d \
  -DDIM=3 -DAMR=ON -DMMB=ON -DLAPACK=ON \
  -DENABLE_MPI=ON -DENABLE_HDF5=ON -DENABLE_EXAMPLES=ON \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-g -O1" \
  -DDISABLE_INIT_CHECKS=ON
cmake --build build/integration-3d --parallel 6
ctest --test-dir build/integration-3d --output-on-failure --timeout 180
```

For 2D, use a separate build directory and `-DDIM=2`. Assertions remain enabled;
`DISABLE_INIT_CHECKS` disables the expensive automatic BoxData value scans.
On this Apple Silicon machine the additional settings were
`HDF5_ROOT=/opt/homebrew/opt/hdf5-mpi`,
`-DLAPACK_ROOT=/opt/homebrew/opt/openblas`,
`-DCMAKE_CXX_FLAGS=-I/opt/homebrew/opt/openblas/include`, and
`-DCMAKE_POLICY_VERSION_MINIMUM=3.5` for the bundled BLT with CMake 4.4.3.

## Validation

CPU builds succeeded in 2D and 3D with AMR, MMB, MPI, HDF5, and LAPACK enabled.
The active MHD example builds. Five added regression cases pass on both one and
two MPI ranks: generic RK4 stages, zero-width boundaries, constant-preserving
angular interpolation, coarsening across different patch sizes, and new/legacy
checkpoint metadata.

| Suite | Result |
| --- | --- |
| 2D, 41 CTest targets | 39 pass; HDF5Tests and MBAMRTests fail |
| 3D, 42 CTest targets | 31 pass; 10 fail; MBInterpOpTests exceeds 240 seconds |

The 2D failures reproduce in an isolated master snapshot. The ten 3D failures
also reproduce on master: OperatorTests (CellAverageMatrixProduct), MBGraphTests,
MBDisjointBoxLayoutTests, MBBoundaryRegisterTests, MBLevelBoxDataTests,
MBMultigridTests, HDF5Tests, MBLevelOpTests, MBLevelFluxRegisterTests, and MBAMRTests.
Master's MBInterpOpTests also exceeded its 180-second comparison timeout; that
long convergence test has **not** been validated to completion. Baseline builds
needed only Clang void-pointer-arithmetic fixes and removal of stale calls to
private `fixRotations`, applied solely in the temporary snapshot. No baseline
numerical algorithms or test assertions were changed.

The initial unoptimized runs also timed out on several targets. Those were
rechecked with `-O1` and assertions enabled; the table reflects the rechecks.
Some upstream tests report themselves disabled, so a passing CTest target does
not imply every advertised upstream case ran.

The patched HelioCubed CubedSphere executable builds against this branch. A
synthetic radial-pulse RK4 step on one rank (six 12x12x12 blocks) and two RK3 steps
on two ranks (six 12x24x24 blocks, 12x12x12 patches), with viscosity and boundary
reflux enabled, completed with finite checkpoint values. Compared with identical
runs using unmodified origin/cubedSphereDev, maximum absolute checkpoint
differences divided by the largest baseline value were 9.67e-17 and 9.63e-17.
A further two-rank, two-step RK3 nonradial magnetic-pulse case (initial condition
5, same grid) also completed with finite values; its checkpoint difference on
the same scale was 1.43e-12. Normalizing each component separately instead
gave a largest difference of 3.56e-11.

These checks used shell radii 0.1 to 2.0 and stretch 1.5. They are short numerical
regressions, not production-scale solar-wind validation.

The companion HelioCubed patch passes `git apply --check` against commit
455a49546dc0c720da544b8dbc5c72c57fd871ce. Logs and temporary reference builds for
this session are in `/private/tmp/proto-integration`; native CTest logs are also
under `build/integration-2d` and `build/integration-3d`.

CUDA/HIP execution, long-time accuracy, and production datasets were not tested.
Resolve or explicitly triage the inherited failures and complete the long
convergence test before promoting this branch to master.
