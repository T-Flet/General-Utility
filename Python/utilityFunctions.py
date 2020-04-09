from itertools import chain
from functools import reduce


# Note: versions of functions with the '_h' suffix only work and are optimal for collections of hashable elements


def flatten(list_of_lists): return list(chain.from_iterable(list_of_lists))


def partition(p, xs): # Haskell's partition function: partition p xs == (filter p xs, filter (not . p) xs)
    def select(acc, x):
        acc[not p(x)].append(x)
        return acc
    return reduce(select, xs, ([],[]))


def unique(xs):
    seen = [] # Note: 'in' tests x is z or x == z, hence it works with __eq__ overloading
    return [x for x in xs if x not in seen and not seen.append(x)] # Neat short-circuit 'and' trick
def unique_h(xs): return type(xs)(set(xs))


def eq(xs, ys):
    cys = list(ys) # make a mutable copy
    try:
        for x in xs: cys.remove(x)
    except ValueError: return False
    return not cys
def eq_h(xs, ys): return xs == ys # Pointless, but here for completeness


def diff(xs, ys):
    cxs = list(xs) # make a mutable copy
    try:
        for y in ys: cxs.remove(y)
    except ValueError: pass
    return cxs
def diff_h(xs, ys): return type(xs)(set(xs) - set(ys))


def chunk(xs, n): return (xs[i:i + n] for i in range(0, len(xs), n))


def interval_overlap(a, b): # Two interval tuples
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


# Update a dictionary's entries with those of another using a given function, e.g. appending (operator.add is ideal for this)
# NOTE: This modifies d0, so might want to give it a deepcopy
def update_dict_with(d0, d1, f):
    for k, v in d1.items(): d0[k] = f(d0[k], v) if k in d0 else v
    return d0


