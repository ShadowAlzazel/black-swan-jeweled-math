#!/usr/bin/env python3
"""Superstable period-doubling centers for f_c(z)=z*z+c on the real axis.

Numerical experiment, not an interval-certified proof. mpmath is the only
dependency. Output is written next to this script. The selected sign-changing
brackets stay on the principal real period-doubling cascade, and explicit
proper-divisor checks reject lower-period roots.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp

DECIMAL_PRECISION = 80
REPORTED_SIGNIFICANT_DIGITS = 72
MAX_N = 9
BISECTION_STEPS = 250


def orbit(c, steps):
    z = mp.mpf(0)
    values = []
    for _ in range(steps):
        z = z * z + c
        values.append(z)
    return values


def residual(c, period):
    z = mp.mpf(0)
    for _ in range(period):
        z = z * z + c
    return z


def parameter_derivative(c, period):
    z, derivative = mp.mpf(0), mp.mpf(0)
    for _ in range(period):
        derivative = 2 * z * derivative + 1
        z = z * z + c
    return derivative


def digits(value):
    return mp.nstr(value, REPORTED_SIGNIFICANT_DIGITS)


def find_center(n, centers):
    if n == 2:
        left, right = mp.mpf('-1.4'), mp.mpf('-1.2')
    else:
        previous_gap = centers[-2] - centers[-1]
        left = centers[-1] - previous_gap / 4
        right = centers[-1] - previous_gap / 5
    original = [left, right]
    period = 2 ** n
    left_value, right_value = residual(left, period), residual(right, period)
    assert left_value * right_value < 0
    # This scan is a numerical sanity check, not a root-counting certificate.
    grid = [left + (right - left) * k / 64 for k in range(65)]
    sampled_signs = [mp.sign(residual(x, period)) for x in grid]
    sampled_crossings = sum(a * b < 0 for a, b in zip(sampled_signs, sampled_signs[1:]))
    assert sampled_crossings == 1
    for _ in range(BISECTION_STEPS):
        midpoint = (left + right) / 2
        middle_value = residual(midpoint, period)
        if left_value * middle_value <= 0:
            right, right_value = midpoint, middle_value
        else:
            left, left_value = midpoint, middle_value
    return (left + right) / 2, {
        'initial_bracket': [digits(x) for x in original],
        'final_bracket': [digits(left), digits(right)],
        'sampled_sign_crossings': sampled_crossings,
    }


def calculate():
    mp.mp.dps = DECIMAL_PRECISION
    centers = [mp.mpf(0), mp.mpf(-1)]
    brackets = [None, None]
    for n in range(2, MAX_N + 1):
        c, bracket = find_center(n, centers)
        centers.append(c)
        brackets.append(bracket)
    rows, half_returns = [], [None]
    for n, c in enumerate(centers):
        period = 2 ** n
        values = orbit(c, period)
        return_residual = abs(values[-1])
        assert return_residual < mp.mpf('1e-65')
        divisors = [{'period': 2 ** k, 'absolute_residual': digits(abs(values[2 ** k - 1]))}
                    for k in range(n)]
        assert all(mp.mpf(item['absolute_residual']) > mp.mpf('1e-10') for item in divisors)
        row = {
            'n': n, 'period': period, 'center': digits(c),
            'absolute_return_residual': digits(return_residual),
            'parameter_derivative': digits(parameter_derivative(c, period)),
            'proper_power_of_two_divisor_checks': divisors,
            'bracketing': brackets[n],
        }
        if n >= 1:
            half_return = values[2 ** (n - 1) - 1]
            half_returns.append(half_return)
            closest_index = min(range(period - 1), key=lambda k: abs(values[k])) + 1
            assert closest_index == 2 ** (n - 1)
            row['half_period_critical_return'] = digits(half_return)
            row['closest_nonzero_orbit_step'] = closest_index
        if n >= 2:
            delta = (centers[n - 1] - centers[n - 2]) / (centers[n] - centers[n - 1])
            alpha = half_returns[n - 1] / half_returns[n]
            row['delta_estimate'] = digits(delta)
            row['signed_alpha_estimate'] = digits(alpha)
        rows.append(row)
    return {
        'map': 'f_c(z)=z^2+c',
        'cascade': 'Principal real period-doubling cascade from c=0 toward c≈-1.401155189.',
        'working_precision_decimal_digits': DECIMAL_PRECISION,
        'reported_significant_digits': REPORTED_SIGNIFICANT_DIGITS,
        'mpmath_version': mp.__version__,
        'method': '250 bisections inside sign-changing brackets; 65-point sign scan; proper-divisor and closest-critical-return checks.',
        'delta_definition': '(c_(n-1)-c_(n-2))/(c_n-c_(n-1)); parameter-interval scaling.',
        'alpha_definition': 'd_(n-1)/d_n where d_n=f_(c_n)^(2^(n-1))(0); signed spatial scaling near the critical point. Some conventions quote its positive magnitude.',
        'interpretation': 'These are superstable-center ratios, not the exact period-doubling bifurcation parameter values. Their limiting ratios are the Feigenbaum constants for the quadratic universality class.',
        'certification': 'High-precision numerical observations; the brackets and scans are not an interval-arithmetic proof of uniqueness or universality.',
        'rows': rows,
        'tests': ['Return residual <1e-65 for every reported center.',
                  'Every proper power-of-two divisor has return residual >1e-10.',
                  'The half-period return is the closest nonzero point to the critical point in each computed cycle.',
                  'Each selected nontrivial bracket changes sign and has exactly one sign crossing on the 65-point scan.'],
    }


def main():
    result = calculate()
    target = Path(__file__).with_name('feigenbaum_centers.json')
    target.write_text(json.dumps(result, indent=2) + '\n')
    for row in result['rows']:
        print(row['period'], row['center'][:38], row.get('delta_estimate', '')[:23],
              row.get('signed_alpha_estimate', '')[:23])
    print(target)


if __name__ == '__main__':
    main()
