"""Lightweight nested stage timer used by both the reference and fast pipelines."""
import time
from contextlib import contextmanager


class StageTimer:
    def __init__(self):
        self.records = []  # (depth, name, seconds)
        self._depth = 0

    @contextmanager
    def __call__(self, name):
        idx = len(self.records)
        self.records.append([self._depth, name, None])
        self._depth += 1
        t0 = time.perf_counter()
        try:
            yield
        finally:
            self.records[idx][2] = time.perf_counter() - t0
            self._depth -= 1

    def report(self, title="Stage timings"):
        lines = [title, "-" * 64]
        for depth, name, secs in self.records:
            secs = float("nan") if secs is None else secs
            lines.append(f"{'  ' * depth}{name:<{48 - 2 * depth}} {secs:9.2f} s")
        return "\n".join(lines)

    def as_dict(self):
        return [{"depth": d, "name": n, "seconds": s} for d, n, s in self.records]


# a module-level timer that stage code can use without plumbing
TIMER = StageTimer()
