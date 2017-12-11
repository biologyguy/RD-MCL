#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import random
import pandas as pd
from collections import OrderedDict
from types import SimpleNamespace
from buddysuite import buddy_resources as br

from .. import mcmcmc


def test_variable_init():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    assert var.name == "foo"
    assert var.min == 2
    assert var.max == 5
    assert round(var.rand_gen.random(), 12) == 0.847433736937
    assert var.current_value == 2.403092732337
    assert var.draw_value == 2.403092732337
    assert var.history == OrderedDict([('draws', [2.403092732337]), ('accepts', [])])


def test_variable_draw_new_value():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    var.draw_new_value(heat=0.1)
    assert var.draw_value == 2.695965617119
    assert var.history == OrderedDict([('draws', [2.403092732337, 2.695965617119]), ('accepts', [])])

    var.current_value = 1.9
    var.draw_new_value(heat=0.1)
    assert var.draw_value == 2.517084979257
    assert var.history == OrderedDict([('draws', [2.403092732337, 2.695965617119, 2.517084979257]), ('accepts', [])])

    var.current_value = 5.1
    var.draw_new_value(heat=0.1)
    assert var.draw_value == 4.911174134664
    assert var.history == OrderedDict([('draws', [2.403092732337, 2.695965617119, 2.517084979257, 4.911174134664]),
                                       ('accepts', [])])

    var.current_value = 100000
    with pytest.raises(RuntimeError) as err:
        var.draw_new_value(0.1)

    assert "Popped safety valve in Variable.draw_new_value() (draw val = 99700.3507217117)" in str(err)


def test_variable_draw_random():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    var.draw_random()
    assert var.current_value == 4.542301210812
    assert var.draw_value == 4.542301210812


def test_variable_set_value():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    var.set_value(1234)
    assert var.current_value == 1234
    assert var.draw_value == 1234


def test_variable_accept_draw():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    var.draw_value = 1234
    var.accept_draw()
    assert var.current_value == 1234
    assert var.history == OrderedDict([('draws', [2.403092732337]), ('accepts', [1234])])


def test_variable_str():
    var = mcmcmc.Variable("foo", min_val=2, max_val=5, r_seed=1)
    assert str(var) == """
Name: foo
Min: 2
Max: 5
Current value: 2.403092732337
Draw value: 2.403092732337
History: [('draws', [2.403092732337]), ('accepts', [])]
"""


def test_walker_init(capsys):
    rand_gen = random.Random(1)
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 4, heat=0.1, params=["foo", "bar"],
                            quiet=False, r_seed=1, lava=True, ice=False, min_max=(1, 5))

    assert walker.variables == [foo_var, bar_var]
    assert walker.function() == 4
    assert walker.params == ["foo", "bar"]
    assert walker.lava is True
    assert walker.ice is False
    assert walker.heat == 0.1
    assert walker.current_score is None
    assert walker.proposed_score is None
    assert walker.score_history == [1, 5]
    assert round(walker.rand_gen.random(), 12) == 0.83576510392
    assert walker.name == "iK2ZWeqhFWCEPyYngFb5"

    out, err = capsys.readouterr()
    assert "User-defined initial chain parameters: 1, 5\n" in out

    with pytest.raises(AttributeError) as err:
        mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 4, heat=0.1, lava=True, ice=True)
    assert "The _Walker.lava and .ice parameters" in str(err)

    with pytest.raises(ValueError) as err:
        mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 4, heat=0.1, min_max=True)
    assert "min_max value passed into _Walker is not of type [num, different num]" in str(err)

    # Find starting params
    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: rand_gen.random(), heat=0.1,
                            params=["foo", "bar"], quiet=False, r_seed=1)
    assert walker.score_history == [0.134364244112, 0.847433736937]
    out, err = capsys.readouterr()
    assert out == """\
Setting initial chain parameters:
\tStep 1: foo = 0.15, bar = 0.51, Score = 0.13436424411240122
\tStep 2: foo = 0.15, bar = 0.51, Score = 0.8474337369372327
"""

    with pytest.raises(RuntimeError) as err:
        mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, quiet=False)

    assert "Popped the safety valve while initializing chain." in str(err)
    out, err = capsys.readouterr()
    assert "Step 30:" in out
    assert "Step 31:" not in out

    mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5), quiet=True)
    out, err = capsys.readouterr()
    assert not out + err


def test_walker_accpet(capsys):
    foo_var = SimpleNamespace(draw_random=lambda: True, accept_draw=lambda: print("foo accepted"), draw_value=0.1,
                              name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, accept_draw=lambda: print("bar accepted"), draw_value=0.5,
                              name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5), quiet=True)
    walker.proposed_score = 0.1234
    walker.accept()
    assert walker.proposed_score == 0.1234
    assert walker.current_score == 0.1234

    out, err = capsys.readouterr()
    assert out == "foo accepted\nbar accepted\n"


def test_walker_set_heat():
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5), quiet=True)
    assert walker.heat == 0.1
    walker.set_heat(0.75)
    assert walker.heat == 0.75

    with pytest.raises(ValueError) as err:
        walker.set_heat(-1)
    assert "heat values must be positive, between 0.000001 and 1.0." in str(err)

    with pytest.raises(ValueError) as err:
        walker.set_heat(0)
    assert "heat values must be positive, between 0.000001 and 1.0." in str(err)

    with pytest.raises(ValueError) as err:
        walker.set_heat(1.01)
    assert "heat values must be positive, between 0.000001 and 1.0." in str(err)


def test_walker_dump_obj():
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5),
                            r_seed=1, quiet=True)
    dump = walker._dump_obj()
    assert dump["vars"][0].name == "foo"
    assert dump["lava"] is False
    assert dump["ice"] is False
    assert dump["heat"] == 0.1
    assert dump["cur_score"] is None
    assert dump["prop_score"] is None
    assert dump["score_hist"] == [1.0, 5.0]
    assert dump["name"] == "iK2ZWeqhFWCEPyYngFb5"


def test_walker_apply_dump(monkeypatch):
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5),
                            r_seed=1, quiet=True)

    new_vars = {"vars": [foo_var, bar_var], "lava": True, "ice": False, "heat": 0.75, "cur_score": 2.5,
                "prop_score": 3.1, "score_hist": [1.12, 3.42], "name": "SomEOtheRnAme"}

    monkeypatch.setattr(mcmcmc._Walker, "set_heat", lambda self, heat: setattr(self, "heat", heat))
    walker._apply_dump(new_vars)
    assert walker.variables == [foo_var, bar_var]
    assert walker.lava is True
    assert walker.ice is False
    assert walker.heat == 0.75
    assert walker.current_score == 2.5
    assert walker.proposed_score == 3.1
    assert walker.score_history == [1.12, 3.42]
    assert walker.name == "SomEOtheRnAme"

    new_vars["ice"] = True
    with pytest.raises(AttributeError) as err:
        walker._apply_dump(new_vars)

    assert "The _Walker.lava and .ice parameters both set to True. This is not valid." in str(err)


def test_walker_str():
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)

    walker = mcmcmc._Walker(variables=[foo_var, bar_var], func=lambda *_: 1, heat=0.1, min_max=(1, 5),
                            r_seed=1, quiet=True)

    assert str(walker) == """\
Walker iK2ZWeqhFWCEPyYngFb5
\tfoo:\t0.15
\tbar:\t0.51
\tScore:\tNone"""


def test_chain_init():
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)
    walker1 = SimpleNamespace(variables=[foo_var, bar_var])
    walker2 = SimpleNamespace(variables=[foo_var, bar_var])

    tmp_file = br.TempFile()

    chain = mcmcmc._Chain(walkers=[walker1, walker2], outfile=tmp_file.path, cold_heat=0.01, hot_heat=0.2)
    assert chain.walkers == [walker1, walker2]
    assert chain.outfile == tmp_file.path
    assert chain.cold_heat == 0.01
    assert chain.hot_heat == 0.2
    assert chain.step_counter == 0
    assert chain.best_score_ever_seen == 0
    assert tmp_file.read() == """\
Gen\tfoo\tbar\tresult
"""


def test_chain_swap_hot_cold(monkeypatch, capsys):
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)
    lava_foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.222)
    lava_bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="bar", current_value=0.999)
    ice_foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.123)
    ice_bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="bar", current_value=0.321)

    walker1 = SimpleNamespace(variables=[foo_var, bar_var], lava=False, ice=False, current_score=35,
                              set_heat=lambda heat: print("Setting walker1 heat = %s" % heat))
    walker2 = SimpleNamespace(variables=[foo_var, bar_var], lava=False, ice=False, current_score=15,
                              set_heat=lambda heat: print("Setting walker2 heat = %s" % heat))
    lavawalker = SimpleNamespace(variables=[lava_foo_var, lava_bar_var], lava=True, ice=False, current_score=45,
                                 set_heat=lambda heat: print("Changing lava_walker heat! Oh Nos!"))
    ice_walker = SimpleNamespace(variables=[ice_foo_var, ice_bar_var], lava=False, ice=True, current_score=10,
                                 set_heat=lambda heat: print("Changing ice_walker heat! Oh Nos!"))

    tmp_file = br.TempFile()

    monkeypatch.setattr(mcmcmc._Chain, "get_best_walker", lambda *_: walker1)
    monkeypatch.setattr(mcmcmc._Chain, "get_cold_walker", lambda *_: walker2)
    monkeypatch.setattr(mcmcmc._Chain, "get_ice_walker", lambda *_: False)

    chain = mcmcmc._Chain(walkers=[walker1, walker2], outfile=tmp_file.path, cold_heat=0.01, hot_heat=0.2)
    chain.swap_hot_cold()
    out, err = capsys.readouterr()
    assert "Setting walker1 heat = 0.01" in out
    assert "Setting walker2 heat = 0.2" in out
    assert chain.best_score_ever_seen == 35

    monkeypatch.setattr(mcmcmc._Chain, "get_best_walker", lambda *_: lavawalker)
    monkeypatch.setattr(mcmcmc._Chain, "get_cold_walker", lambda *_: walker1)

    chain.walkers.append(lavawalker)
    chain.swap_hot_cold()
    out, err = capsys.readouterr()
    assert not out
    assert chain.best_score_ever_seen == 45
    assert foo_var.current_value == 0.222
    assert bar_var.current_value == 0.999

    monkeypatch.setattr(mcmcmc._Chain, "get_ice_walker", lambda *_: ice_walker)

    lavawalker.current_score = 55
    chain.walkers.append(ice_walker)
    chain.swap_hot_cold()
    out, err = capsys.readouterr()
    assert not out
    assert chain.best_score_ever_seen == 55
    assert ice_foo_var.current_value == 0.222
    assert ice_bar_var.current_value == 0.999

    # Ice chain returned as best, but is lower than best ever, so do not copy values
    monkeypatch.setattr(mcmcmc._Chain, "get_best_walker", lambda *_: ice_walker)

    ice_foo_var.current_value = 0.01
    ice_bar_var.current_value = 0.10101

    chain.swap_hot_cold()
    out, err = capsys.readouterr()
    assert not out
    assert chain.best_score_ever_seen == 55
    assert foo_var.current_value == 0.222
    assert bar_var.current_value == 0.999

    # Now give ice walker the best score ever
    monkeypatch.setattr(mcmcmc._Chain, "get_best_walker", lambda *_: ice_walker)

    ice_walker.current_score = 100

    chain.swap_hot_cold()
    out, err = capsys.readouterr()
    assert not out
    assert chain.best_score_ever_seen == 100
    assert foo_var.current_value == 0.01
    assert bar_var.current_value == 0.10101


def test_chain_get_best_walker():
    walker1 = SimpleNamespace(current_score=35)
    walker2 = SimpleNamespace(current_score=15)

    chain = SimpleNamespace(walkers=[walker1, walker2], get_best_walker=mcmcmc._Chain.get_best_walker)
    assert chain.get_best_walker(chain) == walker1

    chain.walkers = [walker2, walker1]
    assert chain.get_best_walker(chain) == walker1


def test_chain_get_cold_walker():
    walker1 = SimpleNamespace(heat=0.1)
    walker2 = SimpleNamespace(heat=0.4)

    chain = SimpleNamespace(walkers=[walker1, walker2], get_cold_walker=mcmcmc._Chain.get_cold_walker, cold_heat=0.1)
    assert chain.get_cold_walker(chain) == walker1

    chain.walkers = [walker2, walker1]
    assert chain.get_cold_walker(chain) == walker1


def test_chain_get_ice_walker():
    walker1 = SimpleNamespace(ice=False)
    walker2 = SimpleNamespace(ice=False)
    ice_walker = SimpleNamespace(ice=True)

    chain = SimpleNamespace(walkers=[walker1, walker2], get_ice_walker=mcmcmc._Chain.get_ice_walker)
    assert chain.get_ice_walker(chain) is False

    chain.walkers = [walker2, ice_walker, walker1]
    assert chain.get_ice_walker(chain) == ice_walker


def test_chain_get_results():
    tmp_file = br.TempFile()
    tmp_file.write("""rec_id1,rec_id2,r_square
BOL-PanxαB,Bab-PanxαA,0.016894041431
BOL-PanxαB,Bch-PanxαA,0.087311057754
BOL-PanxαB,Bfo-PanxαE,0.274041115357""")

    chain = SimpleNamespace(outfile=tmp_file.path, get_results=mcmcmc._Chain.get_results)
    assert type(chain.get_results(chain)) == pd.DataFrame
    assert str(chain.get_results(chain)) == """\
      rec_id1     rec_id2        r_square
0  BOL-PanxαB  Bab-PanxαA  0.016894041431
1  BOL-PanxαB  Bch-PanxαA  0.087311057754
2  BOL-PanxαB  Bfo-PanxαE  0.274041115357""", print(chain.get_results(chain))


def test_chain_write_sample():
    foo_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.1, name="foo", current_value=0.15)
    bar_var = SimpleNamespace(draw_random=lambda: True, draw_value=0.5, name="bar", current_value=0.51)
    walker1 = SimpleNamespace(variables=[foo_var, bar_var], lava=False, ice=False, current_score=35,
                              heat=0.1)

    tmp_file = br.TempFile()
    chain = SimpleNamespace(step_counter=2, get_cold_walker=lambda *_: walker1, outfile=tmp_file.path,
                            write_sample=mcmcmc._Chain.write_sample)

    chain.write_sample(chain)
    assert tmp_file.read() == "2\t0.15\t0.51\t35\n", print(tmp_file.read())


def test_chain_dump_obj():
    walker1 = SimpleNamespace(_dump_obj=lambda *_: "walker1")
    walker2 = SimpleNamespace(_dump_obj=lambda *_: "walker2")
    tmp_file = br.TempFile()
    tmp_file.write("outfile results")

    chain = SimpleNamespace(walkers=[walker1, walker2], outfile=tmp_file.path, cold_heat=0.1, hot_heat=0.2,
                            step_counter=20, best_score_ever_seen=100, _dump_obj=mcmcmc._Chain._dump_obj)

    dump = chain._dump_obj(chain)
    assert dump["walkers"] == ["walker1", "walker2"]
    assert dump["cold_heat"] == 0.1
    assert dump["hot_heat"] == 0.2
    assert dump["step_count"] == 20
    assert dump["best_score"] == 100
    assert dump["results"] == "outfile results"


def test_chain_apply_dump(capsys):
    walker1 = SimpleNamespace(_apply_dump=lambda *_: print("Applying dump to walker1"))
    walker2 = SimpleNamespace(_apply_dump=lambda *_: print("Applying dump to walker2"))

    tmp_file = br.TempFile()
    chain = SimpleNamespace(walkers=[walker1, walker2], outfile=tmp_file.path, cold_heat=None, hot_heat=None,
                            step_counter=None, best_score_ever_seen=None, _apply_dump=mcmcmc._Chain._apply_dump)

    var_dict = {"walkers": [None, None], "cold_heat": 0.1, "hot_heat": 0.2,
                "step_count": 20, "best_score": 100, "results": "Some results"}
    chain._apply_dump(chain, var_dict)
    assert chain.walkers == [walker1, walker2]
    out, err = capsys.readouterr()
    assert out == "Applying dump to walker1\nApplying dump to walker2\n"
    assert chain.cold_heat == 0.1
    assert chain.hot_heat == 0.2
    assert chain.step_counter == 20
    assert chain.best_score_ever_seen == 100
    assert tmp_file.read() == "Some results"
