def lsame(a, b):
    return a.lower() == b.lower()


def sign(x):
    return 1 if x >= 0 else -1


def slice_(N, inc):
    if inc > 0:
        return slice(None, N * inc, inc)
    else:
        return slice(-(N - 1) * inc, None, inc)
