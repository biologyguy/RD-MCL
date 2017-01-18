#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 25 2015

"""
Use the Metropolis-Hastings algorithm to estimate optimal parameters for a given function.
ToDo: Auto-detect when to end the run (when the stationary distribution is more or less established)
"""

import os
import sys
import random
import string
import numpy as np
from scipy.stats import norm
from buddysuite.buddy_resources import DynamicPrint, TempDir, usable_cpu_count
from copy import deepcopy
from multiprocessing import Process


class Variable:
    def __init__(self, _name, _min, _max, variance_covariate=0.05, r_seed=None):
        self.name = _name
        self.min = _min
        self.max = _max
        _range = _max - _min
        self.variance = _range * variance_covariate

        # select a random start value
        self.rand_gen = random.Random() if not r_seed else random.Random(r_seed)
        self.current_value = self.rand_gen.random() * _range + _min
        self.draw_value = self.current_value
        self.history = {"draws": [self.draw_value], "accepts": []}

    def draw_new_value(self):
        #  NOTE: Might need to tune _variance if the acceptance rate is to low or high.
        draw_val = np.random.normal(self.current_value, self.variance)
        safety_check = 100
        while draw_val < self.min or draw_val > self.max:
            if draw_val < self.min:
                draw_val = self.min + abs(self.min - draw_val)

            elif draw_val > self.max:
                draw_val = self.max - abs(self.max - draw_val)

            else:
                sys.exit("Error: %s" % draw_val)

            if safety_check < 0:
                import traceback
                raise RuntimeError("Popped safety valve on step()\n%s\n%s" % (traceback.print_stack(), draw_val))
            safety_check -= 1
        self.draw_value = draw_val
        self.history["draws"].append(self.draw_value)
        return

    def draw_random(self):
        self.current_value = self.rand_gen.random() * (self.max - self.min) + self.min
        self.draw_value = self.current_value
        return

    def set_value(self, value):
        self.current_value = value
        self.draw_value = value
        return

    def accept_draw(self):
        self.current_value = self.draw_value
        self.history["accepts"].append(self.current_value)
        return


class _Chain:
    def __init__(self, variables, function, params=None, quiet=False, r_seed=None):
        self.variables = variables
        self.function = function
        self.params = params
        self.heat = 0.32
        self.current_raw_score = 0.
        self.proposed_raw_score = 0.
        self.current_score = 0.
        self.proposed_score = 0.
        self.score_history = []
        self.rand_gen = random.Random() if not r_seed else random.Random(r_seed)
        self.name = "".join([self.rand_gen.choice(string.ascii_letters + string.digits) for _ in range(20)])

        # Sample `function` for starting min/max scores
        self.raw_min = 0.
        self.raw_max = 0.
        valve = 0
        if not quiet:
            print("Setting initial chain parameters:")
        while len(self.score_history) < 2 or min(self.score_history) == max(self.score_history):
            if valve == 30:
                raise RuntimeError("Popped the safety valve while initializing chain.")
            valve += 1

            output = "\tStep %s:" % valve
            func_args = []
            for variable in self.variables:
                variable.draw_random()
                func_args.append(variable.draw_value)
                output += " %s = %s," % (variable.name, variable.current_value)

            func_args.append(self.rand_gen.randint(1, 999999999999999))  # Always add a new seed for the target function
            score = self.function(func_args) if not self.params else self.function(func_args, self.params)
            self.score_history.append(score)
            output += " Score = %s" % score

            if not quiet:
                print(output)

        for variable in self.variables:
            variable.draw_random()

        # To adapt to any range of min/max, transform the set of possible scores to {0, (max - min)}
        self.raw_max = max(self.score_history)
        self.raw_min = min(self.score_history)
        trans_max = self.raw_max - self.raw_min

        self.gaussian_pdf = norm(trans_max, trans_max * self.heat)

    def step(self):
        func_args = []
        for variable in self.variables:
            variable.draw_new_value()
            func_args.append(variable.draw_value)

        func_args.append(self.rand_gen.randint(1, 999999999999999))  # Always add a new seed for the target function
        self.proposed_raw_score = self.function(func_args) if not self.params else self.function(func_args, self.params)
        if len(self.score_history) >= 1000:
            self.score_history.pop(0)

        self.score_history.append(self.proposed_raw_score)
        self.raw_max = max(self.score_history)
        self.raw_min = min(self.score_history)
        self.set_gaussian()

        self.proposed_score = self.proposed_raw_score - self.raw_min

        if self.proposed_score >= self.current_score:
            self.accept()

        else:
            rand_check_val = self.rand_gen.random()
            accept_check = self.gaussian_pdf.pdf(self.proposed_score) / self.gaussian_pdf.pdf(self.current_score)

            if accept_check > rand_check_val:
                self.accept()

        return

    def accept(self):
        for variable in self.variables:
            variable.accept_draw()

        self.current_raw_score = self.proposed_raw_score
        self.current_score = self.proposed_score
        return

    def set_gaussian(self):
        trans_max = self.raw_max - self.raw_min
        self.gaussian_pdf = norm(trans_max, trans_max * self.heat)
        return

    def set_heat(self, heat):
        # Higher heat = more likely to accept step.
        # Default is 0.32, try 0.1 for low heat
        if heat <= 0. or heat > 1.:
            raise ValueError("heat values must be positive, between 0.000001 and 1.0.")

        self.heat = heat
        self.set_gaussian()
        return

    def __str__(self):
        output = "Chain %s" % self.name
        for variable in self.variables:
            output += "\n\t%s:\t%s" % (variable.name, variable.current_value)
        output += "\n\tScore:\t%s" % self.current_raw_score
        return output


class MCMCMC:
    def __init__(self, variables, function, params=None, steps=10000, sample_rate=1, num_chains=3,
                 outfile='./mcmcmc_out.csv', burn_in=100, quiet=False, r_seed=None):

        self.global_variables = variables
        self.steps = steps
        self.sample_rate = sample_rate
        self.output = ""
        self.outfile = os.path.abspath(outfile)
        with open(self.outfile, "w") as ofile:
            heading = "Gen\t"
            for var in self.global_variables:
                heading += "%s\t" % var.name
            heading += "result\n"
            ofile.write(heading)
        self.rand_gen = random.Random() if not r_seed else random.Random(r_seed)
        self.chains = [_Chain(deepcopy(self.global_variables), function, params=params,
                              quiet=quiet, r_seed=self.rand_gen.randint(1, 999999999999999)) for _ in range(num_chains)]
        self.best = {"score": 0., "variables": {x.name: 0. for x in variables}}

        # Set a cold chain. The cold chain should always be set at index 0, even if a chain swap occurs
        self.cold_heat = 0.1
        self.hot_heat = 0.32
        self.chains[0].set_heat(self.cold_heat)
        self.burn_in = burn_in
        self.quiet = quiet

        self.samples = {"vars": [], "score": []}
        self.printer = DynamicPrint("stderr", quiet=quiet)

    def run(self):
        """
        NOTE: Gibbs sampling is a way of selecting variables one at a time instead of all at once. This is beneficial in
        high dimensional variable space because it will increase the probability of accepting a new sample. It isn't
        implemented here, but good to keep in mind.
        """
        def mc_step_run(_chain, args):
            _func_args, out_path = args
            score = _chain.function(func_args) if not _chain.params else _chain.function(func_args, _chain.params)
            with open(out_path, "w") as ofile:
                ofile.write(str(score))
            return

        def step_parse(_chain):
            with open("%s/%s" % (temp_dir.path, _chain.name), "r") as ifile:
                _chain.proposed_raw_score = float(ifile.read())

            if len(_chain.score_history) >= 1000:
                _chain.score_history.pop(0)

            _chain.score_history.append(_chain.proposed_raw_score)
            _chain.raw_max = max(_chain.score_history)
            _chain.raw_min = min(_chain.score_history)
            _chain.set_gaussian()

            _chain.proposed_score = _chain.proposed_raw_score - _chain.raw_min

            if _chain.proposed_score >= _chain.current_score:
                _chain.accept()

            else:
                rand_check_val = self.rand_gen.random()
                prop_gaus = _chain.gaussian_pdf.pdf(_chain.proposed_score)
                cur_gaus = _chain.gaussian_pdf.pdf(_chain.current_score)
                accept_check = prop_gaus / cur_gaus

                if accept_check > rand_check_val:
                    _chain.accept()
            return

        temp_dir = TempDir()
        counter = 0
        while counter <= self.steps:
            running_processes = 0
            child_list = {}
            for chain in self.chains:  # Note that this will spin off as many new processes as there are chains
                # Start new process
                func_args = []
                for variable in chain.variables:
                    variable.draw_new_value()
                    func_args.append(variable.draw_value)

                # Always add a new seed for the target function
                func_args.append(self.rand_gen.randint(1, 999999999999999))
                outfile = "%s/%s" % (temp_dir.path, chain.name)
                p = Process(target=mc_step_run, args=(chain, [func_args, outfile]))
                p.start()
                child_list[chain.name] = p
                running_processes += 1

            # wait for remaining processes to complete
            finished_processes = []
            while len(child_list) > 0:
                for i in child_list:
                    if child_list[i].is_alive():
                        continue
                    else:
                        finished_processes.append(child_list[i])
                        del child_list[i]
                        break

            # By sorting on PID (ie. name), a specific order is ensured so r_seed works
            finished_processes = sorted(finished_processes, key=lambda l: l.name)
            for finished_chain in self.chains:
                for fin_process in finished_processes:
                    if finished_chain.name == fin_process:
                        step_parse(finished_chain)

            for chain in self.chains:
                if chain.current_raw_score > self.best["score"]:
                    self.best["score"] = chain.current_raw_score
                    self.best["variables"] = {x.name: x.current_value for x in chain.variables}

                # Pseudo burn in, replace
                if counter == self.burn_in:
                    chain.score_history = [chain.raw_max - np.std(chain.score_history), chain.raw_max]

            # Swap any hot chain into the cold chain position if the hot chain score is better than the cold chain score
            best_score = self.chains[0].current_score
            best_chain_index = 0
            for i in range(1, len(self.chains[1:])):
                if self.chains[i].current_score > best_score:
                    best_score = self.chains[i].current_score
                    best_chain_index = i
            if best_chain_index != 0:
                self.chains[0].set_heat(self.hot_heat)
                best_chain = self.chains.pop(best_chain_index)
                self.chains.insert(0, best_chain)

            # Send output to file
            if counter % self.sample_rate == 0:
                self.samples["vars"].append(self.chains[0].variables)
                self.samples["score"].append(self.chains[0].current_raw_score)

                self.output += "%s\t" % counter
                for var in self.chains[0].variables:
                    self.output += "%s\t" % var.current_value
                self.output += "%s\n" % self.chains[0].current_raw_score
                self.write()

            printer_vars = ""
            for x in self.chains[0].variables:
                printer_vars += "%s: %6.3f\t" % (x.name, round(x.current_value, 3))

            self.printer.write("%5.0f: %8.3f\t(%s)" % (counter, round(self.chains[0].current_raw_score, 3), printer_vars.strip()))
            counter += 1

        self.printer.new_line()
        return

    def reset_params(self, params):
        """
        :param params: list of new input parameters pushed to chains
        :return: None
        """
        if len(params) != len(self.chains[0].params):
            raise AttributeError("Incorrect number of params supplied in reset_params().\n%s expected\n%s supplied\n%s"
                                 % (len(self.chains[0].params), len(params), str(params)))

        for _chain in self.chains:
            _chain.params = params
        return

    def write(self):
        with open(self.outfile, "a") as ofile:
            ofile.write(self.output)
        self.output = ''
        return

    def stdout(self):
        print(self.output)
        self.output = ''
        return

if __name__ == '__main__':
    # A couple of examples
    def parabola(var):
        x = var[0]
        output = (-6 * (x ** 2)) + (2 * x) + 4
        return output

    def paraboloid(_vars):
        x, y = _vars
        output = (-6 * (x ** 2)) + (-2 * (y ** 2)) + 2
        # print(output)
        return output

    parabola_variables = [Variable("x", -4, 4)]
    paraboloid_variables = [Variable("x", -100, 100, 0.01), Variable("y", -100, 100, 0.01)]

    mcmcmc = MCMCMC(parabola_variables, parabola, steps=3000, sample_rate=1)
    mcmcmc.run()
    print(mcmcmc.best)
