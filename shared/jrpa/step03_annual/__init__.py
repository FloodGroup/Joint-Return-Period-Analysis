"""03. Estimate univariate design rainfall (S1) and coastal water levels (S2)."""


def __getattr__(name):
    if name == "step03_annual":
        from .run import step03_annual

        return step03_annual
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
