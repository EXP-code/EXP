#!/usr/bin/env python3
"""Check the portable EXP particle HDF5 schema produced by IC generators."""

import argparse
import sys

import h5py
import numpy as np


def fail(message):
    print(f"HDF5 particle check failed: {message}", file=sys.stderr)
    raise SystemExit(1)


parser = argparse.ArgumentParser()
parser.add_argument("file")
parser.add_argument("--count", type=int, required=True)
parser.add_argument("--float-bytes", type=int, default=4)
parser.add_argument("--mass", type=float, default=1.0)
parser.add_argument("--aux-ints", type=int, default=0)
parser.add_argument("--aux-floats", type=int, default=0)
parser.add_argument("--cube", action="store_true", help="check that particles are in the unit cube")
args = parser.parse_args()

with h5py.File(args.file, "r") as h5:
    for name, expected in (("num_particles", args.count), ("num_aux_ints", args.aux_ints),
                           ("num_aux_floats", args.aux_floats)):
        if name not in h5.attrs or h5.attrs[name] != expected:
            fail(f"attribute {name!r} is not {expected}")

    if "particles" not in h5:
        fail("missing /particles group")
    particles = h5["particles"]
    required = {"m", "x", "y", "z", "u", "v", "w"}
    required.update({f"aux_int_{i}" for i in range(args.aux_ints)})
    required.update({f"aux_float_{i}" for i in range(args.aux_floats)})
    if set(particles) != required:
        fail(f"unexpected datasets: {sorted(particles)}")

    for name in required:
        data = particles[name]
        if data.shape != (args.count,):
            fail(f"/particles/{name} has shape {data.shape}, expected ({args.count},)")
        if name.startswith("aux_int_"):
            if data.dtype.kind != "i":
                fail(f"/particles/{name} has dtype {data.dtype}, expected an integer type")
        else:
            if data.dtype.kind != "f" or data.dtype.itemsize != args.float_bytes:
                fail(f"/particles/{name} has dtype {data.dtype}, expected float{8 * args.float_bytes}")
        if not np.isfinite(data[:]).all():
            fail(f"/particles/{name} contains non-finite values")

    if not np.isclose(particles["m"][:].sum(), args.mass, rtol=2e-4, atol=2e-4):
        fail(f"particle masses do not sum to {args.mass}")
    if args.cube:
        for name in ("x", "y", "z"):
            if (particles[name][:] < 0.0).any() or (particles[name][:] > 1.0).any():
                fail(f"/particles/{name} is outside the unit cube")
