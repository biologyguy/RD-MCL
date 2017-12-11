#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
from collections import OrderedDict
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
