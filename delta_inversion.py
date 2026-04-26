import numpy as np
from enum import Enum


class OptionType(Enum):
    CALL = 'CALL'
    PUT = 'PUT'


tol = 1e-8
sqrt_2pi = 2.50662827463


def norm_cdf_poly(xabs: float):
    """Hart's approximation for (scaled) Mills ratio"""
    assert xabs >= 0
    norm_poly = 1.0
    if xabs > 37.0:
        norm_poly = 0

    if xabs < 7.07106781186547:
        build = 0.0352624965998911 * xabs + 0.700383064443688
        build = build * xabs + 6.37396220353165
        build = build * xabs + 33.912866078383
        build = build * xabs + 112.079291497871
        build = build * xabs + 221.213596169931
        build = build * xabs + 220.206867912376
        norm_poly = norm_poly * build
        build = 0.0883883476483184 * xabs + 1.75566716318264
        build = build * xabs + 16.064177579207
        build = build * xabs + 86.7807322029461
        build = build * xabs + 296.564248779674
        build = build * xabs + 637.333633378831
        build = build * xabs + 793.826512519948
        build = build * xabs + 440.413735824752
        norm_poly = norm_poly / build
    else:
        build = xabs + 0.65
        build = xabs + 4 / build
        build = xabs + 3 / build
        build = xabs + 2 / build
        build = xabs + 1 / build
        norm_poly = norm_poly / (build * 2.506628274631)

    return norm_poly


def norm_cdf(x: float) -> float:
    """Normal CDF"""
    if np.isnan(x):
        return np.nan
    if x <= 0:
        return norm_cdf_poly(abs(x)) * np.exp(-0.5 * x ** 2)
    else:
        return 1.0 - norm_cdf_poly(abs(x)) * np.exp(-0.5 * x ** 2)


def inv_normal_cdf(p: float):
    """Inverse Normal CDF"""
    sup = p > 0.5
    up = 1.0 - p if sup else p

    a0 = 2.50662823884
    a1 = -18.61500062529
    a2 = 41.39119773534
    a3 = -25.44106049637

    b0 = -8.47351093090
    b1 = 23.08336743743
    b2 = -21.06224101826
    b3 = 3.13082909833

    c0 = 0.3374754822726147
    c1 = 0.9761690190917186
    c2 = 0.1607979714918209
    c3 = 0.0276438810333863
    c4 = 0.0038405729373609
    c5 = 0.0003951896511919
    c6 = 0.0000321767881768
    c7 = 0.0000002888167364
    c8 = 0.0000003960315187

    x = up - 0.5
    r = np.nan

    if abs(x) < 0.42:
        r = x * x
        r = x * (((a3 * r + a2) * r + a1) * r + a0) / ((((b3 * r + b2) * r + b1) * r + b0) * r + 1.0)
        return -r if sup else r

    r = up
    r = np.log(-np.log(r))
    r = c0 + r * (c1 + r * (c2 + r * (c3 + r * (c4 + r * (c5 + r * (c6 + r * (c7 + r * c8)))))))
    # return r if sup else -r

    return r if sup else -r


def R(x):
    """Mills ratio"""
    if x >= 0:
        return sqrt_2pi * norm_cdf_poly(x)
    else:
        return sqrt_2pi * (np.exp(0.5 * x ** 2) - norm_cdf_poly(-x))


def f1(z: float, total_vol: float) -> float:
    return z / total_vol - 0.5 * total_vol


def f2(z: float, total_vol: float) -> float:
    return z / total_vol + 0.5 * total_vol


def bracket_old(y: float) -> tuple[float, float]:
    sqrt_half_pi = np.sqrt(0.5 * np.pi)
    xx = 1.0 / y
    if 0 < y <= sqrt_half_pi:
        if 0 < y <= 0.5:
            x2 = 0.5 * (1 + np.sqrt(1 - 4 * y * y)) / y
            return (x2, xx)
        else:  # y > 0.5:
            return (0.0, xx)
    elif y >= sqrt_half_pi:
        return (sqrt_half_pi - y, 0.0)
    else:
        raise NotImplementedError

def bracket_new(y: float) -> tuple[float, float]:
    sqrt_half_pi = np.sqrt(0.5 * np.pi)
    two_over_pi = 2.0 / np.pi
    if 0 < y <= sqrt_half_pi:
        ll = 1.0 / y * (np.pi - 2.0 * y ** 2) / (np.pi - 1.0 + np.sqrt(1 + 2.0 * (np.pi - 2.0) * y ** 2))
        rr = 1.0 / y - 2.0 * y / np.pi
        return (ll, rr)
    elif y >= sqrt_half_pi:
        return (-np.sqrt(2.0 * np.log(y / np.sqrt(2.0 * np.pi) + 0.5)),
                         np.sqrt(two_over_pi) - np.sqrt(two_over_pi + 2.0 * np.log(np.sqrt(two_over_pi) * y)))
    else:
        raise NotImplementedError


def bisect(func, y0: float, a: float, b: float,
           tol: float = tol, is_bounds_to_nan: bool = True) -> float:
    """
    find root of f(x)=y0 over x\in[a,b] via bisection
    """
    assert a < b
    f = func(a) - y0
    fmid = func(b) - y0

    if f == 0.0:
        return a
    if fmid == 0.0:
        return b

    if f * fmid < 0.0:
        if f < 0.0:
            rtb = a
            dx = b - a
        else:
            rtb = b
            dx = a - b
        xmid = rtb
        for j in range(0, 40):
            dx = dx * 0.5
            xmid = rtb + dx
            fmid = func(xmid) - y0
            if fmid <= 0.0:
                rtb = xmid
            if np.abs(fmid) < tol:
                break
        x_out = xmid

    else:
        if f < 0:
            x_out = a
        else:
            x_out = b

    if is_bounds_to_nan:  # in case vol was inferred it will return nan
        if np.abs(x_out - a) < tol or np.abs(x_out - b) < tol:
            x_out = np.nan
    return x_out


def R_inv(y0: float):
    """Calculates solution x0 of the equation x0=R^{-1}(y0), where R is Mills ratio using bisection"""
    if y0 <= 0.0 or not np.isfinite(y0):
        raise ValueError(f"y0 must be positive and finite, y0 = {y0}")
    lb, ub = bracket_new(y0)
    x0 = bisect(R, y0, lb, ub, is_bounds_to_nan=False)
    return x0


def delta_prem_adj(total_vol: float, z, opt_type: OptionType) -> float:
    w = 1 if opt_type == OptionType.CALL else -1
    delta_fwd = w * np.exp(z) * norm_cdf(-w * f2(z, total_vol))
    return delta_fwd

def z_prem_adj(total_vol: float, delta_fwd: float, opt_type: OptionType) -> float:
    """resolves log-moneyness z for premium-adjusted forward delta"""
    w = 1 if opt_type == OptionType.CALL else -1
    assert 0.0 < delta_fwd * w < 1.0
    fo = R_inv(1.0 / total_vol)
    z_R = -w * total_vol * inv_normal_cdf(w * delta_fwd) + 0.5 * total_vol ** 2
    if opt_type == OptionType.PUT:
        z_L = np.log(np.fabs(delta_fwd))
    else:
        z_L = total_vol * (fo - 0.5 * total_vol)
        delta_max = np.exp(z_L) * norm_cdf(-f2(z_L, total_vol))
        if delta_fwd >= delta_max:
            return np.nan

    # premium adjusted delta
    func = lambda z: delta_prem_adj(total_vol, z, opt_type)
    z_sol = bisect(func, delta_fwd, z_L, z_R)

    err = np.fabs(func(z_sol) - delta_fwd)
    # print(f"vol={sigma}, error={err}")
    assert err < tol
    return z_sol


