"""02. Merge rainfall sources and select independent event pairs."""


def __getattr__(name):
    if name == "step02_events":
        from .run import step02_events

        return step02_events
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
