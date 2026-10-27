"""04. Fit event margins and copula dependence."""


def __getattr__(name):
    if name == "step04_copula":
        from .run import step04_copula

        return step04_copula
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
