"""Where this repository and its sibling repositories live.

The org keeps every repository side by side (<org root>/<repo>/).  The scripts
used to hard-code the server layout (/store/shuvam/E-motioner-X-SBS/...) and
several read from and wrote to the superseded `datasets/co-evolution/`
snapshot.  They now resolve everything from here:

* REPO_ROOT     this repository (the canonical copy of every input and result)
* ORG_ROOT      $KMAP_ORG_ROOT if set, else the parent of REPO_ROOT
* NARY_SRC      <org root>/n-ary-kmap/src           (nkmap encoders)
* KMAP_SBM_SRC  <org root>/kmap-sbm-validation/src  (kmap_sbm: Gray code, QM)

On the original server the default resolves to the same paths as before.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
ORG_ROOT = Path(os.environ.get("KMAP_ORG_ROOT") or REPO_ROOT.parent)
NARY_SRC = ORG_ROOT / "n-ary-kmap" / "src"
KMAP_SBM_SRC = ORG_ROOT / "kmap-sbm-validation" / "src"
