#!/usr/bin/python
# -*- coding: utf-8 -*-

# Actions of cyclic groups on an algebraic curve.
# Copyright © 2012 Stefano Maggiolo <s.maggiolo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import argparse


### PRINTING FUNCTIONS ###

def to_text(g, branch, results):
    string = ""
    string += "Faithful actions of non-trivial cyclic groups " \
              "on a curve of genus %d" % g
    if sum(branch) != 0:
        string += " with:\n"
        last_branch = max((i for i in range(len(branch)) if branch[i] != 0))
        for i, f in enumerate(branch):
            if f != 0:
                string += "  %d points with counterimage " \
                          "with %d points" % (f, i + 1)
                if i != last_branch:
                    string += ",\n"
    string += ".\n\n"
    for idx_r, r in enumerate(sorted(results)):
        for idx_h, h in enumerate(sorted(results[r], reverse=True)):
            for idx_x, x in enumerate(results[r][h][1]):
                ci_string = ""
                for idx_y, y in enumerate(x):
                    counter_images = (r - results[r][h][0][idx_y])
                    if y == 1:
                        ci_string += "(%s)" % ("*" * counter_images)
                    elif y > 1:
                        ci_string += "(%s)^%d" % ("*" * counter_images, y)
                string += "Z_%d, h = %d: %s\n" % (r, h, ci_string)
        if idx_r != len(results) - 1:
            string +=  "\n"
    return string


def to_latex(g, branch, results):
    multirow = "\\multirow{%d}{*}{$%d$} "
    string = ""
    string += """\
\\begin{table}
  \\centering
  \\begin{tabular}{lll}
    \\toprule
    $r$ & $h$ & Additional ramification\\\\
    \\midrule[1pt]
"""
    for idx_r, r in enumerate(sorted(results)):
        size_r = sum((len(results[r][x][1]) for x in results[r]))
        for idx_h, h in enumerate(sorted(results[r], reverse=True)):
            size_h = len(results[r][h][1])
            for idx_x, x in enumerate(results[r][h][1]):
                ci_string = ""
                for idx_y, y in enumerate(x):
                    counter_images = (r - results[r][h][0][idx_y])
                    if y == 1:
                        ci_string += "(%s)" % ("\\bullet" * counter_images)
                    elif y > 1:
                        ci_string += "(%s)^%d" % ("\\bullet" * counter_images, y)
                string += "    "
                if idx_h == 0 and idx_x == 0:
                    string += multirow % (size_r, r)
                string += "& "
                if idx_x == 0:
                    string +=  multirow % (size_h, h)
                string += "& "
                string += "$%s$\\\\\n" % ci_string
            if idx_h != len(results[r]) - 1:
                string +=  "    \\cmidrule{2-3}\n"
        if idx_r != len(results) - 1:
            string +=  "    \\midrule\n"
    caption = "Cyclic groups acting on a curve of genus $%d$" % g
    if sum(branch) != 0:
        caption += " with"
        last_branch = max((i for i in range(len(branch)) if branch[i] != 0))
        for i, f in enumerate(branch):
            if f != 0:
                caption += " $%d$ points with counterimage " \
                          "consisting of $%d$ points" % (f, i + 1)
                if i != last_branch:
                    caption += ","
    caption += "."
    string += """\
    \\bottomrule
  \\end{tabular}
  \\caption{%s}
  \\label{tab:cyclic_group_actions}
\\end{table}
""" % caption
    return string


### UTILITY FUNCTIONS ###

def gcd(a, b):
    """Return the gcd of a and b.

    a, b (int): integers.

    return (int): the gcd of a and b.

    """
    _gcd, tmp = a, b
    while tmp != 0:
        _gcd, tmp = tmp, _gcd % tmp
    return _gcd


def lcm(a, b):
    """Return the lcm of a and b.

    a, b (int): integers.

    return (int): the lcm of a and b.

    """
    return a * b / gcd(a, b)


### MAIN FUNCTIONS ###

def discrepancies(r):
    """Return a list of all possible contribution to the discrepancy
    coming from a branch point in the quotient curve.

    r (int): the order of the cyclic group.

    return ([int]): all possible contributions.

    """
    ret = []
    for i in xrange(2, r + 1):
        if r % i == 0:
            ret.append((i - 1) * r / i)
    return ret


def all_discrepancies(discrepancies, Q, start=0):
    """Given the list of possible contributions to the discrepancy and
    the discrepancy Q to reach, return a list of all coefficients of
    the linear combinations of elements of discrepancies that sum up
    to Q.

    discrepancies ([int]): the list of possible contribution as
                           returned by discrepancies.
    Q (int): the total discrepancy to reach.
    start (int): used internally: we assume to use only discrepancies
                                  starting from the index start.

    return ([[int]]): the coefficients of the linear combination of
                      discrepancies summing up to Q.

    """
    elements_to_use = len(discrepancies) - start
    if Q == 0:
        return [[0] * elements_to_use]
    elif elements_to_use <= 0:
        return []
    ret = []
    curr = discrepancies[start]
    for i in xrange(Q / curr, -1, -1):
        r = all_discrepancies(discrepancies, Q - i * curr, start + 1)
        for x in r:
            ret.append([i] + x)
    return ret


def _ac_check(r, a, start=0, total=0):
    """Helper function for ac_check. We can safely assign n_1 = 1
    because in any case we can apply an automorphism of C_r sending
    n_1 to 1.

    r (int): the size of the cyclic group.
    a ([int]): the a_i's in the equation.
    start (int): the current position to assign.
    total (int): the current total degree in the right of the equation.

    return (bool): True if the degree equation has a solution.

    """
    if start >= len(a):
        return total % r == 0
    if start == 0:
        return _ac_check(r, a, 1, r // a[0])

    for i in xrange(1, a[start]):
        if gcd(i, a[start]) == 1:
            result = _ac_check(r, a, start + 1, total + i * r / a[start])
            if result:
                return True
    return False


def ac_check(r, br, disc, ret):
    """Check a consequence of the abelian cover condition, that states
    that if there is an cyclic cover with point P_i with preimage
    stabilized by the subgroup of order of a_1, ..., a_k, then there
    is a line bundle L on the quotient curve and coefficients n_1,
    ..., n_k such that rL = sum n_i * r * P_i / a_i and gcd(n_i, a_i)
    = 1. This function check that this equation has a solution at
    least when we pass to the degrees of the line bundles.

    r (int): the size of the cyclic group.
    br ([int]): the prescribed branch points.
    disc ([int]): the possible discrepancies;
    ret ([int]): the additional branch points relative to disc.

    return (bool): True if the degree equation has a solution.

    """
    a = []
    for idx_x, x in enumerate(br):
        if r != idx_x + 1:
            a += ([r / (idx_x + 1)] * x)
    for idx_x, x in enumerate(ret):
        if r != idx_x + 1:
            a += ([r / (r - disc[idx_x])] * x)
    return _ac_check(r, sorted(a, reverse=True))


def run(g, branch):
    """Main algorithm.

    g (int): genus of the curve above.
    branch ([int]): the i-th element is the number of points of the
                    quotient with (i+1) points in the counterimage.
    latex (bool): whether to provide output in LaTeX format.

    """
    lcm_ = 1
    N = 0
    c = 0

    for i, f in enumerate(branch):
        if f != 0:
            lcm_ = lcm(lcm_, i + 1)
            N += f
            c += (i + 1) * f
    cp = 2 * g - 2 + c

    upper_limit = 2 * (2 * g + 1)  # Using Wiman

    if (2 - N < 0):
        upper_limit = min(upper_limit, cp // (N - 2))
    elif (2 - N == 0):
        upper_limit = min(upper_limit, 2 * cp)

    results = {}
    for r in range(max(lcm_, 2), upper_limit + 1, lcm_):
        h_den = 2 * r
        h_num = (2 - N) * r + cp

        disc = discrepancies(r)

        for h in xrange(h_num / h_den, -1, -1):
            Q = h_num - h_den * h
            rets = all_discrepancies(disc, Q)
            for ret in rets:
                if not ac_check(r, branch, disc, ret):
                    continue
                if r not in results:
                    results[r] = {}
                if h not in results[r]:
                    results[r][h] = [disc, []]
                results[r][h][1].append(ret)

    return results


def main():
    """Analyze command line arguments and call the main function.

    """

    parser = argparse.ArgumentParser(description="Compute all possible "
                                     "actions of a cyclic group on a curve.")
    parser.add_argument("genus", type=int,
                        help="genus of the starting curve")
    parser.add_argument("branch", type=int, nargs="*",
                        help="the i-th number is the number of branch points "
                        "with i+1 points in the counterimage")
    parser.add_argument("-l", "--latex", action="store_true",
                        help="output also in LaTeX format")

    args = parser.parse_args()
    if (args.genus < 0 or any([x < 0 for x in args.branch])):
        parser.usage()
        return
    results = run(args.genus, args.branch)
    print_function = to_text
    if args.latex:
        print_function = to_latex
    print print_function(args.genus, args.branch, results)


if __name__ == "__main__":
    main()
