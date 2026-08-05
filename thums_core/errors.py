"""Exception types for THUMS.

Design rule: the solver never substitutes a value on failure. Every degraded path
raises. `runner.run_study` is the ONLY place permitted to catch these; it records
`converged=False` and continues, so a failed run appears in the results table
instead of vanishing into a NaN that `pandas.mean()` silently skips.
"""


class ThumsError(Exception):
    """Base class for all THUMS errors."""


class ThumsPropertyError(ThumsError):
    """A thermophysical property could not be evaluated.

    Raised by `thums_core.props`. Replaces the legacy pattern of returning
    `(1e-9, 0, 0, 0)`, which made a well segment silently adiabatic.
    """


class ThumsConvergenceError(ThumsError):
    """An iterative solve did not converge.

    Replaces the legacy patterns of returning `0.0` from the melt-front solve
    (which routes into the zero-resistance branch and reports maximum heat
    transfer) and of returning a bracket bound from a bisection.
    """

    def __init__(self, what, iterations=None, residual=None, last=None):
        self.what = what
        self.iterations = iterations
        self.residual = residual
        self.last = last
        bits = [f"{what} did not converge"]
        if iterations is not None:
            bits.append(f"after {iterations} iterations")
        if residual is not None:
            bits.append(f"residual={residual:.3e}")
        if last is not None:
            bits.append(f"last iterate={last!r}")
        super().__init__(", ".join(bits))


class ThumsGeometryError(ThumsError):
    """A geometric constraint was violated.

    Examples: the melt front left the domain (delta > delta_max), melt fronts of
    adjacent tubes merged, or V_melt exceeded V_well.
    """


class ThumsValidityError(ThumsError):
    """A modelling assumption was pushed outside its stated validity range.

    Example: the Stefan number exceeded the threshold at which the quasi-steady
    melt-layer approximation is defensible.
    """
