#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Mar 25 2015

"""
Use the Metropolis-Hastings algorithm to estimate optimal parameters for a given function.
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
from collections import OrderedDict
import pandas as pd


class Variable:
    def __init__(self, _name, _min, _max, variance_covariate=0.05, r_seed=None):
        self.name = _name
        self.min = _min
        self.max = _max
        _range = _max - _min
        self.variance = _range * variance_covariate

        # select a random start value
        self.rand_gen = random.Random(r_seed)
        self.current_value = self.rand_gen.random() * _range + _min
        self.draw_value = self.current_value
        self.history = OrderedDict([("draws", [self.draw_value]), ("accepts", [])])

    def draw_new_value(self):
        #  NOTE: Might need to tune _variance if the acceptance rate is to low or high.
        draw_val = self.rand_gen.gauss(self.current_value, self.variance)
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

    def __str__(self):
        return """
Name: {}
Min: {}
Max: {}
Variance: {}
Random: {}
Current value: {}
Draw value: {}
History: {}
""".format(self.name, self.min, self.max, self.variance, self.rand_gen.getstate(),
           self.current_value, self.draw_value, self.history)


class _Walker:
    def __init__(self, variables, func, params=None, quiet=False, r_seed=None, lava=False):
        self.variables = variables
        self.function = func
        self.params = params
        self.heat = 0.32
        self.lava = lava  # This will cause the walker to draw brand new parameters every time
        self.current_score = None
        self.proposed_score = None
        self.score_history = []
        self.rand_gen = random.Random(r_seed)
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
        self.gaussian_pdf = norm(self.raw_max, self.raw_max * self.heat)

    def accept(self):
        for variable in self.variables:
            variable.accept_draw()
        self.current_score = self.proposed_score
        return

    def set_gaussian(self):
        self.gaussian_pdf = norm(self.raw_max, self.raw_max * self.heat)
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
        output += "\n\tScore:\t%s" % self.current_score
        return output


class _Chain(object):
    def __init__(self, walkers, outfile):
        self.walkers = walkers
        self.outfile = outfile

        with open(self.outfile, "w") as ofile:
            heading = "Gen\t"
            for var in walkers[0].variables:
                heading += "%s\t" % var.name
            heading += "result\n"
            ofile.write(heading)

        self.cold_heat = 0.1
        self.hot_heat = 0.32  # Default given to each new _Walker on instantiation
        # Set a cold walker
        self.walkers[0].set_heat(self.cold_heat)  # Set a cold walker
        self.step_counter = 0

    def swap_hot_cold(self):
        # Swap any hot chain into the cold chain position if the hot chain score is better than the cold chain score
        best_walker = self.get_best_walker()
        cold_walker = self.get_cold_walker()
        if best_walker.lava:
            for cold_var, lava_var in zip(cold_walker.variables, best_walker.variables):
                cold_var.current_value = float(lava_var.current_value)
        else:
            cold_walker.set_heat(self.hot_heat)
            best_walker.set_heat(self.cold_heat)
        return

    def get_best_walker(self):
        best_walker = sorted(self.walkers, key=lambda walker: walker.current_score)[-1]
        return best_walker

    def get_cold_walker(self):
        cold_walker = sorted(self.walkers, key=lambda walker: walker.heat)[0]
        return cold_walker

    def get_results(self):
        results = pd.read_csv(self.outfile)
        return results

    def write_sample(self):
        output = "%s\t" % self.step_counter
        best_walker = self.get_best_walker()
        for var in best_walker.variables:
            output += "%s\t" % var.current_value
        output += "%s\n" % best_walker.current_score
        with open(self.outfile, "a") as ofile:
            ofile.write(output)
        return


class MCMCMC:
    """
    Sets up the infrastructure to run a Metropolis Hasting random walk
    """
    def __init__(self, variables, func, params=None, steps=0, sample_rate=1, num_walkers=3, num_chains=3,
                 include_lava=False, outfiles='./chain', burn_in=100, quiet=False, r_seed=None, convergence=1.05):
        self.global_variables = variables
        assert steps >= 100 or steps == 0
        self.steps = steps
        self.sample_rate = sample_rate
        self.outfile = os.path.abspath(outfiles)
        self.rand_gen = random.Random(r_seed)
        self.chains = []
        # The deepcopy below duplicates r_seed, so walkers need to be updated with a new one
        for i in range(num_chains):
            walkers = []
            for j in range(num_walkers):
                walker = _Walker(deepcopy(self.global_variables), func, params=params,
                                 quiet=quiet, r_seed=self.rand_gen.randint(1, 999999999999999))
                for variable in walker.variables:
                    variable.rand_gen.seed(self.rand_gen.randint(1, 999999999999999))
                walkers.append(walker)
            if include_lava:
                walker = _Walker(deepcopy(self.global_variables), func, params=params, lava=True,
                                 quiet=quiet, r_seed=self.rand_gen.randint(1, 999999999999999))
                for variable in walker.variables:
                    variable.rand_gen.seed(self.rand_gen.randint(1, 999999999999999))
                walkers.append(walker)
            chain = _Chain(walkers, "%s_%s.csv" % (self.outfile, i + 1))
            self.chains.append(chain)
        self.best = OrderedDict([("score", None), ("variables", OrderedDict([(x.name, None) for x in variables]))])
        self.burn_in = burn_in
        if convergence > 1.0:
            self.convergence = convergence
        else:
            raise ValueError("Gelman-Rubin convergence ratio must be greater than 1")
        self.quiet = quiet

    def run(self):
        """
        NOTE: Gibbs sampling is a way of selecting variables one at a time instead of all at once. This is beneficial in
        high dimensional variable space because it will increase the probability of accepting a new sample. It isn't
        implemented here, but good to keep in mind.
        """
        def mc_step_run(_walker, args):
            _func_args, out_path = args
            score = _walker.function(func_args) if not _walker.params else _walker.function(func_args, _walker.params)
            with open(out_path, "w") as ofile:
                ofile.write(str(score))
            return

        def step_parse(_chain):  # Implements Metropolis-Hastings
            with open(os.path.join(temp_dir.path, _chain.name), "r") as ifile:
                _chain.proposed_score = float(ifile.read())

            if len(_chain.score_history) >= 1000:
                _chain.score_history.pop(0)

            _chain.score_history.append(_chain.proposed_score)
            _chain.raw_max = max(_chain.score_history)
            _chain.raw_min = min(_chain.score_history)
            _chain.set_gaussian()

            # If the score hasn't been set or the new score is better, the step is accepted
            if _chain.current_score is None or _chain.proposed_score >= _chain.current_score:
                _chain.accept()

            # Even if the score is worse, there's a chance of accepting it relative to how much worse it is
            else:
                rand_check_val = _chain.rand_gen.random()
                # Calculate acceptance ratio as α=f(x')/f(xt), then compare to rand_check_val
                prop_gaus = _chain.gaussian_pdf.pdf(_chain.proposed_score)
                cur_gaus = _chain.gaussian_pdf.pdf(_chain.current_score)
                cur_gaus = 2.2250738585072014e-308 if cur_gaus == 0 else cur_gaus  # Set '0' values to smallest float
                accept_check = prop_gaus / cur_gaus

                if accept_check > rand_check_val:
                    _chain.accept()
            return

        temp_dir = TempDir()
        counter = 0
        while not self._check_convergence() and (counter <= self.steps or self.steps == 0):
            counter += 1
            child_list = OrderedDict()
            for chain in self.chains:  # Note that this will spin off as many new processes as there are walkers
                for walker in chain.walkers:
                    # Start new process
                    func_args = []
                    for variable in walker.variables:
                        if walker.lava:
                            variable.draw_random()
                        else:
                            variable.draw_new_value()
                        func_args.append(variable.draw_value)

                    # Always add a new seed for the target function
                    func_args.append(self.rand_gen.randint(1, 999999999999999))
                    outfile = os.path.join(temp_dir.path, walker.name)
                    p = Process(target=mc_step_run, args=(walker, [func_args, outfile]))
                    p.start()
                    child_list[walker.name] = p

            # wait for remaining processes to complete
            while len(child_list) > 0:
                for _name, child in child_list.items():
                    if child.is_alive():
                        continue
                    else:
                        del child_list[_name]
                        break

            for chain in self.chains:
                for walker in chain.walkers:
                    step_parse(walker)
                    if self.best["score"] is None or walker.current_score > self.best["score"]:
                        self.best["score"] = walker.current_score
                        self.best["variables"] = OrderedDict([(x.name, x.current_value) for x in walker.variables])

            for chain in self.chains:
                chain.swap_hot_cold()

                # Send output to files
                if counter % self.sample_rate == 0:
                    chain.step_counter += 1
                    chain.write_sample()
        return

    def _check_convergence(self):
        """
        Implements the Gelman-Rubin statistic:
            1) Discard burn in (10%)
            2) Calculate within-chain variance for each scalar
            3) Calculate between-chain variance for each scalar
            4) Compute 'potential scale reduction factor' (PSRF) for each scalar --> aiming for < 1.1 for everything

        :return: True or False
        """
        if self.chains[0].step_counter < 100:
            return False

        data = [pd.read_csv(chain.outfile, sep="\t") for chain in self.chains]
        # 1) Strip out 10% burn in
        start_indx = round(len(data[0]) * 0.1)
        data = [x.iloc[start_indx:, :] for x in data]

        # 2) Calculate within-chain variance for each scalar
        within_variance = OrderedDict()
        for df in data:
            scalars = df.drop(["Gen", "result"], axis=1)
            for name, series in scalars.items():
                mean = series.mean()
                series_mean_sqr = (series - mean) ** 2
                variance = series_mean_sqr.sum()/(len(series_mean_sqr) - 1)
                within_variance.setdefault(name, 0.0000001)
                within_variance[name] += variance / len(scalars.columns)

        for name, variance in within_variance.items():
            within_variance[name] /= len(data)

        # 3) Calculate between-chain variance for each scalar
        between_variance = OrderedDict()
        counter = 0
        for i, chain_i in enumerate(data[:-1]):
            scalars_i = chain_i.drop(["Gen", "result"], axis=1)
            for j, chain_j in enumerate(data[i+1:]):
                counter += 2  # Between-variances will be added in both directions
                scalars_j = chain_j.drop(["Gen", "result"], axis=1)
                for name, series_i in scalars_i.items():
                    mean_i = series_i.mean()
                    mean_j = scalars_j[name].mean()

                    series_mean_sqr_i = (series_i - mean_j) ** 2
                    variance_i = series_mean_sqr_i.sum()/(len(series_mean_sqr_i) - 1)
                    between_variance.setdefault(name, 0)
                    between_variance[name] += variance_i / len(scalars_i.columns)

                    series_mean_sqr_j = (scalars_j[name] - mean_i) ** 2
                    variance_j = series_mean_sqr_j.sum()/(len(series_mean_sqr_j) - 1)
                    between_variance[name] += variance_j / len(scalars_j.columns)

        for name, variance in between_variance.items():
            between_variance[name] /= counter

        # 4) Compute 'potential scale reduction factor' (PSRF) for each scalar
        n = len(data[0].columns)
        for name, b in between_variance.items():
            w = within_variance[name]
            psrf = (((n - 1) / n) + (b / (n * w))) ** (1/2)
            if psrf >= self.convergence:
                return False
        return True  # If all PSRFs are below 1.1, time to call it quits

    def reset_params(self, params):
        """
        :param params: list of new input parameters pushed to chains
        :return: None
        """
        if len(params) != len(self.chains[0].walkers[0].params):
            raise AttributeError("Incorrect number of params supplied in reset_params().\n%s expected\n%s supplied\n%s"
                                 % (len(self.chains[0].params), len(params), str(params)))
        for _chain in self.chains:
            for walker in _chain.walkers:
                walker.params = params
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

    mcmcmc = MCMCMC(parabola_variables, parabola, steps=300, sample_rate=1)
    mcmcmc.run()
    print(mcmcmc.best)
