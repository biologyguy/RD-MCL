#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import random
from collections import OrderedDict
from types import SimpleNamespace
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
