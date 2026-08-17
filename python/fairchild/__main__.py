"""The `fairchild` command, as installed by pip.

Two ways in, one implementation: the console script declared in pyproject.toml
(`fairchild …`) and `python -m fairchild …` both land in `main()`, which hands
argv to the same Rust CLI the standalone binary runs.

The Rust side returns an exit code rather than calling `exit()`, so a bad
netlist here raises nothing and kills nothing — it prints and returns non-zero,
exactly as it would from a shell.
"""

import sys

from .fairchild import _cli_main


def main() -> int:
    return _cli_main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
