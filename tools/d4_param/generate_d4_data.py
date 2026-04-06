#!/usr/bin/env python3
"""
Generate CP2K DFT-D4 data modules from the dftd4 library source.

Reads the dftd4 Fortran source tree and produces two CP2K-compatible modules:
  - qs_dispersion_d4_data.F : element property tables + model constants
  - qs_dispersion_d4_ref.F  : reference system data (replaces reference.inc)

Usage:
    python generate_d4_data.py /path/to/dftd4/src/dftd4 [output_dir]

The dftd4 source can be obtained from: https://github.com/dftd4/dftd4
"""

import argparse
import math
import os
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
MAX_ELEM = 118
MAX_REF = 7
N_FREQ = 23
N_SEC = 17

# aatoau = 1 / bohr(Angstrom).  CODATA 2014 value used by mctc-lib.
AATOAU = 1.0 / 0.52917721067

CP2K_HEADER = """\
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!
"""


# ===========================================================================
# Fortran formatting helpers
# ===========================================================================
def fmt_real(vals, per_line=4, fmt="{:20.14f}_dp"):
    """Format a flat list of reals as a Fortran DATA value block."""
    lines = []
    for i in range(0, len(vals), per_line):
        chunk = vals[i : i + per_line]
        s = ", ".join(fmt.format(v).strip() for v in chunk)
        s += ", &" if i + per_line < len(vals) else "/"
        pfx = "      /" if i == 0 else "      &      "
        lines.append(pfx + s)
    return "\n".join(lines)


def fmt_int(vals, per_line=14):
    """Format a flat list of ints as a Fortran DATA value block."""
    lines = []
    for i in range(0, len(vals), per_line):
        chunk = vals[i : i + per_line]
        s = ", ".join(str(int(v)) for v in chunk)
        s += ", &" if i + per_line < len(vals) else "/"
        pfx = "      /" if i == 0 else "      &      "
        lines.append(pfx + s)
    return "\n".join(lines)


def fmt_param(vals, arrname, fmt="{:20.14f}_dp", max_len=100):
    """Format the full PARAMETER array declaration line with greedy fprettify-matching logic."""
    prefix = f"      REAL(KIND=dp), DIMENSION(d4_max_elem), PARAMETER :: {arrname} = ["
    lines = []
    line_items = []
    current_len = len(prefix)
    indent = "         "
    for i, v in enumerate(vals):
        val_str = fmt.format(v).strip()
        is_last = i == len(vals) - 1
        tail_len = 1 if is_last else 2
        item_prefix_len = 0 if not line_items else 2
        if not line_items:
            line_items.append(val_str)
            current_len += len(val_str)
        elif current_len + item_prefix_len + len(val_str) + tail_len <= max_len:
            line_items.append(val_str)
            current_len += item_prefix_len + len(val_str)
        else:
            line_str = (prefix if not lines else indent) + ", ".join(line_items)
            if len(line_str) + 3 <= max_len:
                lines.append(line_str + ", &")
            else:
                lines.append(line_str + ",&")
            line_items = [val_str]
            current_len = len(indent) + len(val_str)
    if line_items:
        line_str = (prefix if not lines else indent) + ", ".join(line_items)
        lines.append(line_str + "]")
    return "\n".join(lines)


def flat2(arr, n_fast, n_slow):
    """Flatten arr[n_slow][n_fast] to Fortran column-major order."""
    return [arr[j][i] for j in range(n_slow) for i in range(n_fast)]


def flat3(arr, n1, n2, n3):
    """Flatten arr[n3][n2][n1] to Fortran column-major order."""
    return [arr[k][j][i] for k in range(n3) for j in range(n2) for i in range(n1)]


# ===========================================================================
# Parsers for dftd4 source files
# ===========================================================================
def extract_fortran_array(text, array_name):
    """Extract numeric values from a Fortran PARAMETER array declaration.

    Handles both `real(wp), parameter :: name(N) = [ ... ]` and
    multi-line continuations with `&`.
    """
    # Find the array declaration
    pattern = (
        rf"{array_name}\s*\([^)]*\)\s*=\s*"
        r"(?:aatoau\s*\*\s*)?"  # optional aatoau multiplication
        r"\[\s*(.*?)\s*\]"
    )
    m = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
    if not m:
        raise ValueError(f"Array '{array_name}' not found")
    raw = m.group(1)
    # Check if aatoau factor is present
    has_aatoau = bool(
        re.search(rf"{array_name}\s*\([^)]*\)\s*=\s*aatoau", text, re.IGNORECASE)
    )
    vals = [float(v) for v in re.findall(r"([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)", raw)]
    return vals, has_aatoau


def parse_element_data(dftd4_dir):
    """Parse the five element-property data files from dftd4/data/."""
    data_dir = dftd4_dir / "data"

    # --- Covalent radii ---
    text = (data_dir / "covrad.f90").read_text()
    covrad_aa, _ = extract_fortran_array(text, "covalent_rad_2009")
    # D3-type covalent radii = 4/3 * aatoau * Pyykko-Atsumi values
    covrad = [4.0 / 3.0 * AATOAU * v for v in covrad_aa]

    # --- Electronegativities ---
    text = (data_dir / "en.f90").read_text()
    en, _ = extract_fortran_array(text, "pauling_en")

    # --- Chemical hardness ---
    text = (data_dir / "hardness.f90").read_text()
    eta, _ = extract_fortran_array(text, "chemical_hardness")

    # --- r4/r2 expectation values ---
    text = (data_dir / "r4r2.f90").read_text()
    r4r2_raw, _ = extract_fortran_array(text, "r4_over_r2")
    # Precompute sqrt(0.5 * r4r2 * sqrt(Z))
    sqzr4r2 = [math.sqrt(0.5 * r4r2_raw[i] * math.sqrt(i + 1)) for i in range(MAX_ELEM)]

    # --- Effective nuclear charges ---
    text = (data_dir / "zeff.f90").read_text()
    zeff_raw, _ = extract_fortran_array(text, "effective_nuclear_charge")
    zeff = [float(v) for v in zeff_raw]

    # Validate lengths
    for name, arr in [
        ("covrad", covrad),
        ("en", en),
        ("eta", eta),
        ("sqzr4r2", sqzr4r2),
        ("zeff", zeff),
    ]:
        assert (
            len(arr) == MAX_ELEM
        ), f"{name}: expected {MAX_ELEM} values, got {len(arr)}"

    return covrad, en, eta, sqzr4r2, zeff


def parse_reference_inc(dftd4_dir):
    """Parse reference.inc and extract all arrays needed by the D4 model."""
    inc_path = dftd4_dir / "reference.inc"
    raw = inc_path.read_text()

    # Preprocess: strip comments, join continuations
    code = []
    for line in raw.split("\n"):
        s = line.strip()
        if s.startswith("!") or s == "":
            continue
        if "!" in s:
            s = s[: s.index("!")]
        code.append(s)
    text = "\n".join(code)
    text = re.sub(r"&\s*\n\s*&\s*", " ", text)
    text = re.sub(r"&\s*\n\s*", " ", text)

    # Split into individual DATA statements
    stmts = []
    for line in text.split("\n"):
        for part in line.split(";"):
            part = part.strip()
            if part.lower().startswith("data "):
                stmts.append(part)

    # Initialise arrays
    refn = [0] * MAX_ELEM
    refcovcn = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    refcn = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    refsys = [[0] * MAX_REF for _ in range(MAX_ELEM)]
    clsq = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    clsh = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    eeqbcq = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    eeqbch = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    hcount = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    ascale = [[0.0] * MAX_REF for _ in range(MAX_ELEM)]
    alphaiw = [[[0.0] * N_FREQ for _ in range(MAX_REF)] for _ in range(MAX_ELEM)]
    secq = [0.0] * N_SEC
    sscale = [0.0] * N_SEC
    seccn = [0.0] * N_SEC
    secaiw = [[0.0] * N_FREQ for _ in range(N_SEC)]

    # Dispatcher: array_name -> (target, rank, is_int)
    scalar_2d = {
        "refcovcn": refcovcn,
        "refcn": refcn,
        "clsq": clsq,
        "clsh": clsh,
        "eeqbcq": eeqbcq,
        "eeqbch": eeqbch,
        "hcount": hcount,
        "ascale": ascale,
    }
    scalar_2d_int = {"refsys": refsys}
    vector_3d = {"alphaiw": alphaiw}
    scalar_1d_sec = {"secq": secq, "sscale": sscale, "seccn": seccn}
    vector_2d_sec = {"secaiw": secaiw}

    for stmt in stmts:
        m = re.match(r"data\s+(\w+)\s*\(([^)]+)\)\s*/\s*(.*?)\s*/", stmt, re.IGNORECASE)
        if not m:
            continue
        name = m.group(1).lower()
        idx_raw = m.group(2).strip()
        val_raw = m.group(3).strip()

        idx_parts = [x.strip() for x in idx_raw.split(",")]
        indices = [None if p == ":" else int(p) - 1 for p in idx_parts]

        vals = []
        for v in re.split(r",\s*", val_raw):
            v = v.strip().replace("_wp", "").replace("_dp", "").rstrip("/")
            if v:
                vals.append(float(v))

        if name == "refn" and len(indices) == 1:
            refn[indices[0]] = int(vals[0])
        elif name in scalar_2d and len(indices) == 2:
            scalar_2d[name][indices[1]][indices[0]] = vals[0]
        elif name in scalar_2d_int and len(indices) == 2:
            scalar_2d_int[name][indices[1]][indices[0]] = int(vals[0])
        elif name in vector_3d and len(indices) == 3 and indices[0] is None:
            ie, ir = indices[2], indices[1]
            for k, v in enumerate(vals[:N_FREQ]):
                vector_3d[name][ie][ir][k] = v
        elif name in scalar_1d_sec and len(indices) == 1:
            scalar_1d_sec[name][indices[0]] = vals[0]
        elif name in vector_2d_sec and len(indices) == 2 and indices[0] is None:
            isec = indices[1]
            for k, v in enumerate(vals[:N_FREQ]):
                vector_2d_sec[name][isec][k] = v

    total_refs = sum(refn)
    n_with_refs = sum(1 for x in refn if x > 0)
    print(
        f"  reference.inc: {len(stmts)} statements, "
        f"{n_with_refs} elements with refs, {total_refs} total refs"
    )

    return {
        "refn": refn,
        "refcovcn": refcovcn,
        "refcn": refcn,
        "refsys": refsys,
        "clsq": clsq,
        "clsh": clsh,
        "eeqbcq": eeqbcq,
        "eeqbch": eeqbch,
        "hcount": hcount,
        "ascale": ascale,
        "alphaiw": alphaiw,
        "secq": secq,
        "sscale": sscale,
        "seccn": seccn,
        "secaiw": secaiw,
    }


# ===========================================================================
# File generators
# ===========================================================================
def generate_d4_data(covrad, en, eta, sqzr4r2, zeff):
    """Generate qs_dispersion_d4_data.F content."""

    def _func(fname, retname, arrname, vals, desc, fmt, per_line, retdoc):
        """Generate one ELEMENTAL FUNCTION with embedded PARAMETER array."""
        body = f"""\
! **************************************************************************************************
!> \\brief {desc}
!> \\param num atomic number (1-{MAX_ELEM})
!> \\return {retdoc}
! **************************************************************************************************
   ELEMENTAL FUNCTION {fname}(num) RESULT({retname})
      INTEGER, INTENT(IN)                                :: num
      REAL(KIND=dp)                                      :: {retname}

{fmt_param(vals, arrname, fmt=fmt)}

      IF (num > 0 .AND. num <= d4_max_elem) THEN
         {retname} = {arrname}(num)
      ELSE
         {retname} = 0.0_dp
      END IF

   END FUNCTION {fname}
"""
        return body

    out = CP2K_HEADER + f"""\

! **************************************************************************************************
!> \\brief Element-specific data and model constants for the DFT-D4 dispersion model.
!>        Data from dftd4 library (github.com/dftd4/dftd4):
!>          E. Caldeweyher et al, JCP 150: 154122 (2019)
!>          E. Caldeweyher et al, PCCP 22: 8499 (2020)
!>        Covalent radii: Pyykko & Atsumi, Chem. Eur. J. 15: 188 (2009)
!>        Electronegativities: Pauling scale
!>        r4/r2 values: PBE0/def2-QZVP (Grimme 2010, Mewes 2018)
!>        Effective nuclear charges: from def2-ECPs
!>        Chemical hardnesses: element-specific for D4 charge scaling
!>
!>        All values are precomputed for {MAX_ELEM} elements (H-Og).
!>        This module is auto-generated — DO NOT EDIT MANUALLY.
!>        Re-generate with: tools/d4_param/generate_d4_data.py
! **************************************************************************************************
MODULE qs_dispersion_d4_data

   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d4_data'

   INTEGER, PARAMETER, PUBLIC :: d4_max_elem = {MAX_ELEM}

   !> Default charge scaling height for partial charge extrapolation
   REAL(KIND=dp), PARAMETER, PUBLIC :: d4_ga = 3.0_dp
   !> Default charge scaling steepness for partial charge extrapolation
   REAL(KIND=dp), PARAMETER, PUBLIC :: d4_gc = 2.0_dp
   !> Default weighting factor for coordination number interpolation
   REAL(KIND=dp), PARAMETER, PUBLIC :: d4_wf = 6.0_dp
   !> Maximum number of reference systems per element
   INTEGER, PARAMETER, PUBLIC :: d4_max_ref = {MAX_REF}

   PUBLIC :: d4_get_covalent_rad, d4_get_electronegativity, d4_get_hardness, &
             d4_get_r4r2_val, d4_get_effective_charge

CONTAINS

"""
    out += _func(
        "d4_get_covalent_rad",
        "rad",
        "covrad",
        covrad,
        "D3-type covalent radii for coordination number computation [Bohr].\n"
        "!>        = 4/3 * aatoau * Pyykko-Atsumi-2009 covalent radii (metals -10%).",
        "{:20.14f}_dp",
        4,
        "covalent radius in Bohr",
    )

    out += _func(
        "d4_get_electronegativity",
        "en",
        "pauling_en",
        en,
        "Pauling electronegativities for the covalent coordination number.",
        "{:10.2f}_dp",
        6,
        "Pauling electronegativity",
    )

    out += _func(
        "d4_get_hardness",
        "eta",
        "chemical_hardness",
        eta,
        "Element-specific chemical hardness for the D4 charge scaling function.",
        "{:16.8f}_dp",
        5,
        "chemical hardness",
    )

    out += _func(
        "d4_get_r4r2_val",
        "val",
        "sqrt_z_r4_over_r2",
        sqzr4r2,
        "Expectation value sqrt(Z * <r4>/<r2>) / 2 for C8/C6 extrapolation.\n"
        "!>        Precomputed as sqrt(0.5 * r4_over_r2(Z) * sqrt(Z)).\n"
        "!>        r4/r2 from PBE0/def2-QZVP (S. Grimme 2010, J. Mewes 2018).",
        "{:20.14f}_dp",
        4,
        "sqrt(Z*r4/r2)/2",
    )

    # Effective charge: integer values stored as dp
    out += _func(
        "d4_get_effective_charge",
        "zeff",
        "effective_nuclear_charge",
        zeff,
        "Effective nuclear charges from def2-ECPs for reference polarizabilities.",
        "{:6.1f}_dp",
        10,
        "effective nuclear charge",
    )

    out += "END MODULE qs_dispersion_d4_data\n"
    return out


def generate_d4_ref(ref):
    """Generate qs_dispersion_d4_ref.F content."""
    out = CP2K_HEADER + f"""\

! **************************************************************************************************
!> \\brief Reference system data for the DFT-D4 dispersion model.
!>        Machine-generated from dftd4 library reference.inc — DO NOT EDIT MANUALLY.
!>        Re-generate with: tools/d4_param/generate_d4_data.py
!>        Source: github.com/dftd4/dftd4
!>        E. Caldeweyher et al, JCP 150: 154122 (2019); PCCP 22: 8499 (2020)
!>
!>        Contains reference coordination numbers, EEQ/EEQBC charges,
!>        dynamic polarizabilities ({N_FREQ} frequency points), and
!>        secondary reference system (SEC) data for {MAX_ELEM} elements.
! **************************************************************************************************
MODULE qs_dispersion_d4_ref

   USE kinds,                           ONLY: dp
   USE qs_dispersion_d4_data,           ONLY: d4_max_elem,&
                                              d4_max_ref
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d4_ref'

   INTEGER, PARAMETER, PUBLIC :: d4_n_freq = {N_FREQ}
   INTEGER, PARAMETER, PUBLIC :: d4_n_sec = {N_SEC}

   INTEGER, PUBLIC, SAVE :: d4_refn(d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refcovcn(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refcn(d4_max_ref, d4_max_elem)
   INTEGER, PUBLIC, SAVE :: d4_refsys(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refq_eeq(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refh_eeq(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refq_eeqbc(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_refh_eeqbc(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_hcount(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_ascale(d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_ref_alphaiw(d4_n_freq, d4_max_ref, d4_max_elem)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_secq(d4_n_sec)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_sscale(d4_n_sec)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_seccn(d4_n_sec)
   REAL(KIND=dp), PUBLIC, SAVE :: d4_secaiw(d4_n_freq, d4_n_sec)

   ! ============================================================================
   ! DATA blocks — machine-generated, do not edit
   ! ============================================================================
"""
    out += f"\n   DATA d4_refn &\n{fmt_int(ref['refn'])}\n"

    for vname, key, desc in [
        ("d4_refcovcn", "refcovcn", "Ref covalent CN"),
        ("d4_refcn", "refcn", "Ref CN"),
        ("d4_refq_eeq", "clsq", "EEQ ref charges"),
        ("d4_refh_eeq", "clsh", "EEQ ref scaling"),
        ("d4_refq_eeqbc", "eeqbcq", "EEQBC ref charges"),
        ("d4_refh_eeqbc", "eeqbch", "EEQBC ref scaling"),
        ("d4_hcount", "hcount", "H count"),
        ("d4_ascale", "ascale", "Alpha scale"),
    ]:
        out += f"\n   ! {desc}\n"
        out += f"   DATA {vname} &\n{fmt_real(flat2(ref[key], MAX_REF, MAX_ELEM))}\n"

    out += f"\n   ! Ref system indices\n"
    out += f"   DATA d4_refsys &\n{fmt_int(flat2(ref['refsys'], MAX_REF, MAX_ELEM))}\n"

    out += f"\n   ! Reference polarizabilities alphaiw(freq,ref,elem)\n"
    out += f"   DATA d4_ref_alphaiw &\n"
    out += f"{fmt_real(flat3(ref['alphaiw'], N_FREQ, MAX_REF, MAX_ELEM))}\n"

    out += f"\n   DATA d4_secq &\n{fmt_real(ref['secq'])}\n"
    out += f"\n   DATA d4_sscale &\n{fmt_real(ref['sscale'])}\n"
    out += f"\n   DATA d4_seccn &\n{fmt_real(ref['seccn'])}\n"
    out += f"\n   ! SEC polarizabilities\n"
    out += f"   DATA d4_secaiw &\n{fmt_real(flat2(ref['secaiw'], N_FREQ, N_SEC))}\n"

    out += "\nEND MODULE qs_dispersion_d4_ref\n"
    return out


# ===========================================================================
# Verification
# ===========================================================================
def verify_d4_data(text, covrad, en, eta, sqzr4r2, zeff):
    """Round-trip verify generated qs_dispersion_d4_data.F."""
    errors = 0

    def _check_param(name, expected, tol=1e-10, is_int=False):
        nonlocal errors
        pattern = rf"PARAMETER :: {name} = \[\s*(.*?)\s*\]"
        m = re.search(pattern, text, re.DOTALL)
        if not m:
            print(f"    FAIL: {name} not found")
            errors += 1
            return
        if is_int:
            vals = [float(x) for x in re.findall(r"([-\d.]+)_dp", m.group(1))]
        else:
            vals = [float(x) for x in re.findall(r"([-+\d.eE]+)_dp", m.group(1))]
        if len(vals) != len(expected):
            print(f"    FAIL: {name} length {len(vals)} != {len(expected)}")
            errors += 1
            return
        maxdiff = max(abs(a - b) for a, b in zip(vals, expected))
        if maxdiff > tol:
            errors += 1
            print(f"    FAIL: {name} max diff = {maxdiff:.2e}")
        else:
            print(f"    OK: {name} ({len(vals)} values, max diff {maxdiff:.2e})")

    _check_param("covrad", covrad)
    _check_param("pauling_en", en, tol=1e-6)
    _check_param("chemical_hardness", eta)
    _check_param("sqrt_z_r4_over_r2", sqzr4r2)
    _check_param("effective_nuclear_charge", zeff, tol=0.2, is_int=True)
    return errors


def verify_d4_ref(text, ref):
    """Round-trip verify generated qs_dispersion_d4_ref.F."""
    errors = 0

    def _extract(name, is_int=False):
        pattern = rf"DATA {name}\s+&\s*\n\s*/\s*(.*?)\s*/"
        m = re.search(pattern, text, re.DOTALL)
        if not m:
            raise ValueError(f"DATA {name} not found")
        raw = m.group(1)
        if is_int:
            return [int(x) for x in re.findall(r"(-?\d+)", raw)]
        return [float(x) for x in re.findall(r"([-+\d.eE]+)_dp", raw)]

    # refn
    vals = _extract("d4_refn", is_int=True)
    if vals != ref["refn"]:
        print("    FAIL: d4_refn mismatch")
        errors += 1
    else:
        print(f"    OK: d4_refn ({sum(vals)} total refs)")

    # 2D real arrays
    for vname, key in [
        ("d4_refcovcn", "refcovcn"),
        ("d4_refcn", "refcn"),
        ("d4_refq_eeq", "clsq"),
        ("d4_refh_eeq", "clsh"),
        ("d4_refq_eeqbc", "eeqbcq"),
        ("d4_refh_eeqbc", "eeqbch"),
        ("d4_hcount", "hcount"),
        ("d4_ascale", "ascale"),
    ]:
        vals = _extract(vname)
        expected = flat2(ref[key], MAX_REF, MAX_ELEM)
        maxdiff = max(abs(a - b) for a, b in zip(vals, expected))
        if maxdiff > 1e-6:
            print(f"    FAIL: {vname} max diff = {maxdiff:.2e}")
            errors += 1
        else:
            print(f"    OK: {vname}")

    # refsys
    vals = _extract("d4_refsys", is_int=True)
    expected = flat2(ref["refsys"], MAX_REF, MAX_ELEM)
    if vals != expected:
        print("    FAIL: d4_refsys")
        errors += 1
    else:
        print("    OK: d4_refsys")

    # alphaiw
    vals = _extract("d4_ref_alphaiw")
    expected = flat3(ref["alphaiw"], N_FREQ, MAX_REF, MAX_ELEM)
    maxdiff = max(abs(a - b) for a, b in zip(vals, expected))
    if maxdiff > 1e-6:
        print(f"    FAIL: d4_ref_alphaiw max diff = {maxdiff:.2e}")
        errors += 1
    else:
        print(f"    OK: d4_ref_alphaiw ({len(vals)} values)")

    # SEC
    for vname, key in [
        ("d4_secq", "secq"),
        ("d4_sscale", "sscale"),
        ("d4_seccn", "seccn"),
    ]:
        vals = _extract(vname)
        maxdiff = max(abs(a - b) for a, b in zip(vals, ref[key]))
        if maxdiff > 1e-6:
            errors += 1
        print(f"    OK: {vname}" if maxdiff <= 1e-6 else f"    FAIL: {vname}")

    vals = _extract("d4_secaiw")
    expected = flat2(ref["secaiw"], N_FREQ, N_SEC)
    maxdiff = max(abs(a - b) for a, b in zip(vals, expected))
    print(f"    OK: d4_secaiw" if maxdiff <= 1e-6 else f"    FAIL: d4_secaiw")
    if maxdiff > 1e-6:
        errors += 1

    return errors


# ===========================================================================
# Main
# ===========================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Generate CP2K DFT-D4 data modules from dftd4 source."
    )
    parser.add_argument(
        "dftd4_dir", type=Path, help="Path to dftd4/src/dftd4 directory"
    )
    parser.add_argument(
        "output_dir",
        type=Path,
        nargs="?",
        default=Path("."),
        help="Output directory (default: current dir)",
    )
    parser.add_argument(
        "--no-verify", action="store_true", help="Skip round-trip verification"
    )
    args = parser.parse_args()

    dftd4_dir = args.dftd4_dir
    out_dir = args.output_dir

    # Validate input
    if not (dftd4_dir / "data" / "covrad.f90").exists():
        sys.exit(
            f"Error: {dftd4_dir}/data/covrad.f90 not found.\n"
            f"Provide the path to the dftd4/src/dftd4 directory."
        )
    if not (dftd4_dir / "reference.inc").exists():
        sys.exit(f"Error: {dftd4_dir}/reference.inc not found.")

    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Parse ---
    print("Parsing element data...")
    covrad, en, eta, sqzr4r2, zeff = parse_element_data(dftd4_dir)
    print(f"  5 tables × {MAX_ELEM} elements OK")

    print("Parsing reference.inc...")
    ref = parse_reference_inc(dftd4_dir)

    # --- Generate ---
    print("Generating qs_dispersion_d4_data.F...")
    data_text = generate_d4_data(covrad, en, eta, sqzr4r2, zeff)
    data_path = out_dir / "qs_dispersion_d4_data.F"
    data_path.write_text(data_text)
    nlines = data_text.count("\n")
    print(f"  Written: {data_path} ({nlines} lines, {len(data_text) // 1024} KB)")

    print("Generating qs_dispersion_d4_ref.F...")
    ref_text = generate_d4_ref(ref)
    ref_path = out_dir / "qs_dispersion_d4_ref.F"
    ref_path.write_text(ref_text)
    nlines = ref_text.count("\n")
    print(f"  Written: {ref_path} ({nlines} lines, {len(ref_text) // 1024} KB)")

    # --- Verify ---
    if not args.no_verify:
        print("\nVerifying qs_dispersion_d4_data.F...")
        e1 = verify_d4_data(data_text, covrad, en, eta, sqzr4r2, zeff)

        print("Verifying qs_dispersion_d4_ref.F...")
        e2 = verify_d4_ref(ref_text, ref)

        total_errors = e1 + e2
        if total_errors == 0:
            print(f"\nAll verifications passed.")
        else:
            print(f"\n{total_errors} verification(s) FAILED.")
            sys.exit(1)


if __name__ == "__main__":
    main()
