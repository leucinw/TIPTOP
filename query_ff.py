#!/usr/bin/env python3
"""
Query AMOEBA/AMOEBA+ force field parameters from Poltype2 output files.

Usage:
    python query_ff.py --atoms <spec> [--terms <spec>] [--xyz final.xyz] [--key final.key]

<spec> for --atoms : comma-separated atom indices (1-based, from XYZ) or atom types.
                     Auto-detected: numbers matching valid XYZ indices are treated as
                     atom indices; others are treated as atom types.

<spec> for --terms : comma-separated energy-term names (case-insensitive).
                     Omit to show all terms.

Type vs. class:
    In AMOEBA/AMOEBA+, 'atom' lines define both atom type (col 1) and atom class (col 2).
    multipole and polarize parameters are indexed by atom TYPE.
    All other parameters (bond, angle, torsion, vdw, etc.) are indexed by atom CLASS.
    This script handles the distinction automatically.

Supported term names:
    atom        atom type/class definitions
    bond        bond stretching
    angle       in-plane + out-of-plane angle bending  (anglep + angle)
    anglep      in-plane angle only
    strbnd      stretch-bend coupling
    opbend      out-of-plane bending
    torsion     torsional (dihedral)
    vdw         van der Waals
    multipole   permanent multipoles          (alias: mpole)
    polarize    atomic polarizability          (alias: polar)
    chgpen      charge penetration
    chgtrn      charge transfer
    bndcflux    bond charge flux
    angcflux    angle charge flux
    bonded      all bonded terms (bond, angle, anglep, strbnd, opbend, torsion)
    electrostatic  all electrostatic terms (multipole, polarize, chgpen, chgtrn,
                                            bndcflux, angcflux)
    all         every term (default)

Examples:
    python query_ff.py --atoms 1,2
    python query_ff.py --atoms 401,402 --terms bond,angle,torsion
    python query_ff.py --atoms 9      --terms vdw,polarize,multipole
    python query_ff.py --atoms 401    --terms all
"""

import argparse
import sys

# ── Number of atom identifier columns each keyword carries ───────────────────
KEYWORD_N_ATOMS = {
    "atom":      1,   # col 1 = atom type  (col 2 = atom class, handled separately)
    "vdw":       1,
    "polarize":  1,
    "chgpen":    1,
    "chgtrn":    1,
    "multipole": 1,   # only the first field is the subject; the rest define the frame
    "bond":      2,
    "bndcflux":  2,
    "angle":     3,
    "anglep":    3,
    "strbnd":    3,
    "angcflux":  3,
    "opbend":    4,
    "torsion":   4,
}

# Keywords whose identifiers are atom TYPES (not classes).
# Everything else uses atom CLASS.
TYPE_BASED_KEYWORDS = {"multipole", "polarize", "atom"}

# Number of *extra* lines that follow the header line for multi-line blocks
MULTILINE_EXTRA = {
    "multipole": 4,
}

# Friendly aliases → list of actual keywords
TERM_ALIASES = {
    "atom":          ["atom"],
    "bond":          ["bond"],
    "angle":         ["angle", "anglep"],
    "anglep":        ["anglep"],
    "strbnd":        ["strbnd"],
    "opbend":        ["opbend"],
    "torsion":       ["torsion"],
    "vdw":           ["vdw"],
    "multipole":     ["multipole"],
    "mpole":         ["multipole"],
    "polarize":      ["polarize"],
    "polar":         ["polarize"],
    "chgpen":        ["chgpen"],
    "chgtrn":        ["chgtrn"],
    "bndcflux":      ["bndcflux"],
    "angcflux":      ["angcflux"],
    "bonded":        ["bond", "angle", "anglep", "strbnd", "opbend", "torsion"],
    "electrostatic": ["multipole", "polarize", "chgpen", "chgtrn",
                      "bndcflux", "angcflux"],
    "all":           list(KEYWORD_N_ATOMS.keys()),
}


# ── File parsers ──────────────────────────────────────────────────────────────

def parse_xyz(filename):
    """Return (idx_to_type, type_to_idxs) from a Tinker XYZ file."""
    idx_to_type = {}
    type_to_idxs = {}
    with open(filename) as f:
        lines = f.readlines()
    n_atoms = int(lines[0].split()[0])
    for line in lines[1: n_atoms + 1]:
        parts = line.split()
        if len(parts) < 6:
            continue
        idx   = int(parts[0])
        atype = int(parts[5])
        idx_to_type[idx] = atype
        type_to_idxs.setdefault(atype, []).append(idx)
    return idx_to_type, type_to_idxs


def build_type_to_class(key_filename):
    """
    Scan 'atom' lines in the KEY file and return a {atom_type: atom_class} dict.
    atom line format: atom  <type>  <class>  <element>  "<name>"  ...
    """
    type_to_class = {}
    with open(key_filename) as f:
        for line in f:
            parts = line.split()
            if not parts or parts[0].lower() != 'atom':
                continue
            if len(parts) < 3:
                continue
            try:
                atype  = int(parts[1])
                aclass = int(parts[2])
                type_to_class[atype] = aclass
            except ValueError:
                continue
    return type_to_class


def parse_key(filename):
    """
    Parse a Tinker KEY file.

    Returns a list of dicts:
        {'keyword': str, 'types': list[int], 'lines': list[str]}

    'types' holds the atom identifier integers from the parameter line
    (absolute values; multipole frame signs are stripped).
    For type-based keywords these are atom types; for class-based keywords
    these are atom classes.
    Comment lines are stored with keyword='#'.
    """
    entries = []
    with open(filename) as f:
        raw = f.readlines()

    i = 0
    while i < len(raw):
        line = raw[i]
        stripped = line.strip()

        if not stripped:
            i += 1
            continue

        if stripped.startswith('#'):
            entries.append({'keyword': '#', 'types': [], 'lines': [line.rstrip('\n')]})
            i += 1
            continue

        parts = stripped.split()
        kw = parts[0].lower()

        if kw not in KEYWORD_N_ATOMS:
            i += 1
            continue

        n_id_cols = KEYWORD_N_ATOMS[kw]
        ids = []
        for j in range(1, n_id_cols + 1):
            if j < len(parts):
                try:
                    ids.append(abs(int(parts[j])))
                except ValueError:
                    break

        entry_lines = [line.rstrip('\n')]
        for _ in range(MULTILINE_EXTRA.get(kw, 0)):
            i += 1
            if i < len(raw):
                entry_lines.append(raw[i].rstrip('\n'))

        entries.append({'keyword': kw, 'types': ids, 'lines': entry_lines})
        i += 1

    return entries


# ── Argument resolution helpers ───────────────────────────────────────────────

def resolve_target_types(atom_specs, idx_to_type, type_to_idxs):
    """
    Map user atom specs → set of atom types.

    A number is first tried as an XYZ atom index; if not found there,
    it is tried as an atom type.  Warns and skips unknown values.
    """
    target = set()
    for spec in atom_specs:
        spec = spec.strip()
        if not spec:
            continue
        try:
            n = int(spec)
        except ValueError:
            print(f"Warning: '{spec}' is not an integer, skipping.", file=sys.stderr)
            continue
        if n in idx_to_type:
            target.add(idx_to_type[n])
        elif n in type_to_idxs:
            target.add(n)
        else:
            print(f"Warning: '{n}' not found as an atom index or atom type, skipping.",
                  file=sys.stderr)
    return target


def resolve_terms(term_specs):
    """Expand user term specs → set of keyword strings."""
    keywords = set()
    for t in term_specs:
        t = t.strip().lower()
        if not t:
            continue
        if t in TERM_ALIASES:
            keywords.update(TERM_ALIASES[t])
        elif t in KEYWORD_N_ATOMS:
            keywords.add(t)
        else:
            print(f"Warning: unknown energy term '{t}', skipping.", file=sys.stderr)
    return keywords


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Query AMOEBA/AMOEBA+ force field parameters from Poltype2 output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--atoms", required=True,
        help="Comma-separated atom indices (1-based) or atom types, e.g. '1,2' or '401,402'.",
    )
    parser.add_argument(
        "--terms", default="all",
        help="Comma-separated energy term names (default: all). See --help for the full list.",
    )
    parser.add_argument("--xyz", default="final.xyz",
                        help="Tinker XYZ file (default: final.xyz)")
    parser.add_argument("--key", default="final.key",
                        help="Tinker KEY file (default: final.key)")
    args = parser.parse_args()

    # ── Build mappings ────────────────────────────────────────────────────────
    idx_to_type, type_to_idxs = parse_xyz(args.xyz)
    type_to_class = build_type_to_class(args.key)

    # ── Resolve target atom types and their classes ───────────────────────────
    atom_specs   = args.atoms.split(',')
    target_types = resolve_target_types(atom_specs, idx_to_type, type_to_idxs)
    if not target_types:
        sys.exit("Error: no valid atom types resolved from --atoms.")

    # For class-based parameters, map each queried type to its class.
    # Fall back to the type number itself if the atom line is absent.
    target_classes = {type_to_class.get(t, t) for t in target_types}

    term_specs      = args.terms.split(',')
    target_keywords = resolve_terms(term_specs)
    if not target_keywords:
        sys.exit("Error: no valid energy terms resolved from --terms.")

    entries = parse_key(args.key)

    # ── Print header ──────────────────────────────────────────────────────────
    print(f"# Atom types queried : {sorted(target_types)}")
    detail = "  |  ".join(
        f"type {t} (class {type_to_class.get(t, '?')}) → XYZ atoms {type_to_idxs.get(t, [])}"
        for t in sorted(target_types)
    )
    print(f"# ({detail})")
    print(f"# Energy terms       : {sorted(target_keywords)}")
    print(f"# Note: multipole/polarize matched by atom TYPE; "
          f"all other terms matched by atom CLASS.")

    # ── Filter and print ──────────────────────────────────────────────────────
    prev_kw = None
    matched = 0
    for entry in entries:
        kw = entry['keyword']
        if kw == '#':
            continue
        if kw not in target_keywords:
            continue

        # Choose the correct identifier set for this keyword
        if kw in TYPE_BASED_KEYWORDS:
            match_set = target_types
        else:
            match_set = target_classes

        if not any(ident in match_set for ident in entry['types']):
            continue

        if kw != prev_kw:
            print(f"\n{'─' * 60}")
            label = f"  {kw.upper()}"
            if kw in TYPE_BASED_KEYWORDS:
                label += "  [matched by atom type]"
            else:
                label += "  [matched by atom class]"
            print(label)
            print(f"{'─' * 60}")
            prev_kw = kw

        for line in entry['lines']:
            print(line)
        matched += 1

    print(f"\n# {matched} parameter block(s) printed.")


if __name__ == "__main__":
    main()
