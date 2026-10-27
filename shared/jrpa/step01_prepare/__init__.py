"""01. Read observations, apply quality rules and vertical references."""


def __getattr__(name):
    if name == "step01_prepare":
        from .run import step01_prepare

        return step01_prepare
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
