"""Check that regenerated mechanism files match the committed ones.

Run after regenerating mechanisms in a git checkout. Every file that git
reports as modified is compared with its committed version:

- floating-point literals may differ by a small relative tolerance, because the
  last printed digit of fitted coefficients (e.g. transport fits computed
  with numpy/LAPACK) can depend on the CPU the generator runs on;
- the ``cantera-version``, ``git-commit`` and ``date`` header lines that
  Cantera writes into ``qssa.yaml`` are ignored;
- everything else (code structure, order of operations, names) must be
  identical.

Untracked files under the checked paths are reported as failures.
"""

import argparse
import pathlib
import re
import subprocess
import sys

NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")
VOLATILE_YAML_KEYS = ("cantera-version:", "git-commit:", "date:")


def git(*args, cwd=None):
    """Run a git command and return its stdout."""
    return subprocess.run(
        ["git", *args], check=True, capture_output=True, text=True, cwd=cwd
    ).stdout


def split_numbers(line):
    """Split a line into its non-numeric skeleton and its numeric literals."""
    return NUMBER.sub("#", line), NUMBER.findall(line)


def is_integer(literal):
    """Return True for integer literals (indices, counts, ids)."""
    return not any(c in literal for c in ".eE")


def numbers_close(a, b, rtol):
    """Compare two numeric literals: integers exactly, floats with rtol."""
    if a == b:
        return True
    if is_integer(a) or is_integer(b):
        return False
    x, y = float(a), float(b)
    return abs(x - y) <= rtol * max(abs(x), abs(y))


def compare(path, committed, regenerated, rtol):
    """Return a list of differences that are not tolerated."""
    old = committed.splitlines()
    new = regenerated.splitlines()
    if path.endswith(".yaml"):
        old = [x for x in old if not x.startswith(VOLATILE_YAML_KEYS)]
        new = [x for x in new if not x.startswith(VOLATILE_YAML_KEYS)]
    if len(old) != len(new):
        return [f"{path}: {len(old)} lines committed, {len(new)} regenerated"]
    problems = []
    for lineno, (a, b) in enumerate(zip(old, new, strict=True), start=1):
        if a == b:
            continue
        skel_a, num_a = split_numbers(a)
        skel_b, num_b = split_numbers(b)
        if (
            skel_a != skel_b
            or len(num_a) != len(num_b)
            or not all(
                numbers_close(x, y, rtol) for x, y in zip(num_a, num_b, strict=True)
            )
        ):
            problems.append(f"{path}:{lineno}:\n  - {a.strip()}\n  + {b.strip()}")
    return problems


def main():
    """Compare the regenerated files with the committed ones."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "paths",
        nargs="*",
        default=["Mechanisms"],
        help="paths relative to the repository root (default: Mechanisms)",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=1e-6,
        help="relative tolerance for numeric literals (default: 1e-6)",
    )
    parser.add_argument(
        "--max-report", type=int, default=20, help="differences to print per file"
    )
    args = parser.parse_args()

    root = pathlib.Path(git("rev-parse", "--show-toplevel").strip())
    untracked = git(
        "ls-files", "--others", "--exclude-standard", "--", *args.paths, cwd=root
    ).split()
    modified = git("diff", "--name-only", "HEAD", "--", *args.paths, cwd=root).split()

    failed = False
    tolerated = 0
    for path in untracked:
        print(f"untracked file after regeneration: {path}")
        failed = True
    for path in modified:
        problems = compare(
            path,
            git("show", f"HEAD:{path}", cwd=root),
            (root / path).read_text(),
            args.rtol,
        )
        if problems:
            failed = True
            print(f"{path}: {len(problems)} difference(s) beyond tolerance")
            for p in problems[: args.max_report]:
                print(p)
        else:
            tolerated += 1
    print(
        f"{len(modified)} modified file(s), {tolerated} within tolerance, "
        f"{len(untracked)} untracked"
    )
    if failed:
        print(
            "Committed mechanisms do not match the regenerated ones. Regenerate"
            " them with the locked versions (uv sync --locked) and commit the"
            " result, or apply the regenerated-mechanisms patch of the CI run"
            " (see Docs/sphinx/Ceptr.rst)."
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
