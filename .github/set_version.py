"""Write the version of a release tag (e.g. ``v0.1.1``) into ``src/twowayfeweights/__init__.py``.

Used by .github/workflows/publish.yml before building, so the version on PyPI is the tag's.

    python .github/set_version.py v0.1.1
"""

import re
import sys
from pathlib import Path

INIT = Path(__file__).resolve().parents[1] / "src" / "twowayfeweights" / "__init__.py"

tag = sys.argv[1].strip()
m = re.fullmatch(r"v?(\d+\.\d+\.\d+(?:(?:a|b|rc)\d+)?(?:\.post\d+)?(?:\.dev\d+)?)", tag)
if m is None:
    sys.exit(f"Tag {tag!r} is not a version: use v<major>.<minor>.<patch>, e.g. v0.1.1")
version = m.group(1)

text = INIT.read_text(encoding="utf-8")
new, n = re.subn(r'^__version__ = "[^"]*"$', f'__version__ = "{version}"', text, flags=re.M)
if n != 1:
    sys.exit(f"__version__ not found in {INIT}")
INIT.write_text(new, encoding="utf-8")
print(f"version set to {version}")
