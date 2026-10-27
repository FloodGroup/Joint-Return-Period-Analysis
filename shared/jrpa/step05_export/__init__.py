"""05. Calculate Kendall design points and export tables and figures."""


def __getattr__(name):
    if name == "step05_export":
        from .run import step05_export

        return step05_export
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
