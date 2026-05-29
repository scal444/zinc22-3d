"""Optional NVTX range helpers for profiling.

Set ZINC_NVTX=1 to emit ranges. Without that variable, or without an NVTX
backend available, the context manager is a no-op.
"""

from contextlib import contextmanager
import os


_TRUE_VALUES = {"1", "true", "yes", "on"}
_BACKEND = None
_BACKEND_LOADED = False


def enabled():
    return os.environ.get("ZINC_NVTX", "").strip().lower() in _TRUE_VALUES


def _load_backend():
    global _BACKEND
    global _BACKEND_LOADED

    if not enabled():
        return None
    if _BACKEND_LOADED:
        return _BACKEND

    _BACKEND_LOADED = True
    try:
        import nvtx

        _BACKEND = ("nvtx", nvtx)
        return _BACKEND
    except Exception:
        pass

    try:
        from torch.cuda import nvtx as torch_nvtx

        _BACKEND = ("torch", torch_nvtx)
        return _BACKEND
    except Exception:
        _BACKEND = None
        return None


@contextmanager
def nvtx_range(message):
    backend = _load_backend()
    if backend is None:
        yield
        return

    kind, module = backend
    if kind == "nvtx":
        with module.annotate(message=message):
            yield
        return

    pushed = False
    try:
        module.range_push(message)
        pushed = True
    except Exception:
        pass

    try:
        yield
    finally:
        if pushed:
            try:
                module.range_pop()
            except Exception:
                pass
