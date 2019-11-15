from itertools import chain


def flatten(list_of_lists):
    return list(chain.from_iterable(list_of_lists))


def unique(xs):
    seen = [] # Note: 'in' tests x is z or x == z, hence it works with __eq__ overloading
    return [x for x in xs if x not in seen and not seen.append(x)] # Neat short-circuit 'and' trick


def chunk(xs, n):
    return (xs[i:i + n] for i in range(0, len(xs), n))