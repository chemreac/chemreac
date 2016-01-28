# -*- coding: utf-8 -*-

from operator import attrgetter

from chemreac import ReactionDiffusion
from chemreac.chemistry import (
    Substance, Reaction, ReactionSystem,
    mk_sn_dict_from_names
)


def test_ReactionSystem__to_ReactionDiffusion():
    sbstncs = mk_sn_dict_from_names('AB')
    r1 = Reaction({'A': 2}, {'B': 1}, k=3.0)
    rsys = ReactionSystem([r1], sbstncs)
    rd = ReactionDiffusion.from_ReactionSystem(rsys)
    assert rd.stoich_active == [[0, 0]]
    assert rd.stoich_prod == [[1]]
    assert rd.k == [3.0]


def test_ReactionSystem__from_ReactionDiffusion():
    rd = ReactionDiffusion(2, [[0]], [[1]], [1])
    rsys = rd.to_ReactionSystem('AB')
    assert len(rsys.rxns) == 1
