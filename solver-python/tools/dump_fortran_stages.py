#!/usr/bin/env python3
"""Load the raw stage dumps written by the instrumented Fortran reference.

The instrumented driver (``tools/fortran_shims/Callaway_2D_2V_DG_instrumented.f90``)
writes, into ``$PYBTE_DUMP_DIR``,

    manifest.txt   ``<name> <dtype> <ndim> <d1> .. <dn>`` per array
    <name>.bin     the raw payload in Fortran (column-major) element order

This module turns that into numpy arrays that keep the *Fortran* index
order, i.e. ``d["VDF"][L, TRID, J1, J2]`` with 0-based indices.

Usage as a library::

    from dump_fortran_stages import FortranDump
    d = FortranDump("fortran-build/run_gsis_dump/dump")
    d["TRI_ORDER"].shape            # (N_TRIS, NPOLE, NAZIM)
    d.iter_array(1, "vdf_after_sweep")

Usage from the shell::

    python tools/dump_fortran_stages.py DIR [--list] [--show NAME]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

_DTYPE = {"r8": np.dtype("<f8"), "i4": np.dtype("<i4")}


class FortranDump:
    """Lazy, dict-like access to one dump directory."""

    def __init__(self, directory):
        self.dir = Path(directory)
        manifest = self.dir / "manifest.txt"
        if not manifest.exists():
            raise FileNotFoundError(f"no manifest.txt in {self.dir}")
        self.meta: dict[str, tuple[str, tuple[int, ...]]] = {}
        for line in manifest.read_text().splitlines():
            parts = line.split()
            if len(parts) < 3:
                continue
            name, dtype, ndim = parts[0], parts[1], int(parts[2])
            dims = tuple(int(v) for v in parts[3:3 + ndim])
            self.meta[name] = (dtype, dims)
        self._cache: dict[str, np.ndarray] = {}

    # -- mapping protocol ---------------------------------------------
    def __contains__(self, name: str) -> bool:
        return name in self.meta

    def __len__(self) -> int:
        return len(self.meta)

    def keys(self):
        return self.meta.keys()

    def __getitem__(self, name: str) -> np.ndarray:
        if name in self._cache:
            return self._cache[name]
        if name not in self.meta:
            raise KeyError(f"{name!r} not in dump {self.dir}")
        dtype, dims = self.meta[name]
        raw = np.fromfile(self.dir / f"{name}.bin", dtype=_DTYPE[dtype])
        want = int(np.prod(dims)) if dims else 1
        if raw.size != want:
            raise ValueError(f"{name}: file has {raw.size} elements, manifest says {want}")
        arr = raw.reshape(dims, order="F")
        if dims == (1,):
            arr = arr.reshape(())
        arr.flags.writeable = False
        self._cache[name] = arr
        return arr

    def get(self, name, default=None):
        return self[name] if name in self.meta else default

    def scalar(self, name):
        """Scalars are dumped as length-1 arrays under an ``s_`` prefix."""
        key = name if name in self.meta else f"s_{name}"
        v = self[key]
        return v.item()

    # -- convenience ---------------------------------------------------
    def iter_name(self, step: int, tag: str) -> str:
        return f"it{step:06d}_{tag}"

    def iter_array(self, step: int, tag: str) -> np.ndarray:
        return self[self.iter_name(step, tag)]

    def has_iter(self, step: int, tag: str) -> bool:
        return self.iter_name(step, tag) in self.meta

    def steps(self) -> list[int]:
        out = set()
        for k in self.meta:
            if k.startswith("it") and "_" in k:
                head = k[2:k.index("_")]
                if head.isdigit():
                    out.add(int(head))
        return sorted(out)

    # -- derived -------------------------------------------------------
    def csr(self):
        """The assembled HDG global matrix as a scipy CSR matrix."""
        from scipy.sparse import csr_matrix

        data = np.asarray(self["KKAcomp"], dtype=float)
        indices = np.asarray(self["JKKA"], dtype=np.int64) - 1   # 1-based -> 0-based
        indptr = np.asarray(self["IKKA"], dtype=np.int64) - 1
        n = indptr.size - 1
        return csr_matrix((data, indices, indptr), shape=(n, n))

    def __repr__(self) -> str:
        return f"FortranDump({str(self.dir)!r}, {len(self.meta)} arrays)"


def _main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory")
    p.add_argument("--list", action="store_true")
    p.add_argument("--show", action="append", default=[])
    args = p.parse_args(argv)

    d = FortranDump(args.directory)
    if args.list or not args.show:
        for name in sorted(d.keys()):
            dtype, dims = d.meta[name]
            print(f"{name:44s} {dtype}  {dims}")
    for name in args.show:
        a = d[name]
        print(f"\n=== {name} shape={a.shape} dtype={a.dtype} ===")
        print(f"min={np.min(a)!r} max={np.max(a)!r}")
        print(a)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
