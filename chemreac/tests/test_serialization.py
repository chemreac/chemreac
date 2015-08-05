import os
import tempfile

from chemreac import ReactionDiffusion
from chemreac.serialization import dump, load


def test_dump_load():
    rd1 = ReactionDiffusion(2, [[0]], [[1]], [1])
    f = tempfile.NamedTemporaryFile(delete=False)
    try:
        f.close()
        dump(rd1, f.name)
        rd2 = load(f.name)
        assert rd1.n == rd2.n
        assert rd1.k == rd2.k
        assert rd1.stoich_active == rd2.stoich_active
        assert rd1.stoich_prod == rd2.stoich_prod
    finally:
        os.unlink(f.name)
