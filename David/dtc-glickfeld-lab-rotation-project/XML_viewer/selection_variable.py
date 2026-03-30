"""Random-without-replacement selection variable matching MWEL semantics."""

import random


class SelectionVariable:
    """Implements MWorks selection variable behavior.

    Each call to reset() shuffles the pool. accept() commits the current value
    and advances. next() moves to the next value in the pool.

    Usage matches MWEL:
        tStimulusNumber = svStimNumber    # read current_value
        ...end of trial...
        accept_selections(sv)             # commit
        if tNAccepted >= 80:
            reset_selection(sv)
        else:
            next_selection(sv)
    """

    def __init__(self, name: str, values: list[int] | None = None, nsamples: int = 80):
        self.name = name
        self.values = list(values) if values else list(range(80))
        self.nsamples = nsamples
        self._pool: list[int] = []
        self._index: int = 0
        self._counter: int = 0  # number of accepted samples
        self.reset()

    @property
    def current_value(self) -> int:
        """The current drawn value (what MWEL reads as svStimNumber)."""
        if not self._pool:
            return 0
        return self._pool[self._index]

    def reset(self):
        """Shuffle pool and start over (matches reset_selection)."""
        self._pool = list(self.values)
        random.shuffle(self._pool)
        self._index = 0
        self._counter = 0

    def accept(self):
        """Commit the current value (matches accept_selections)."""
        self._counter += 1

    def next(self):
        """Advance to the next value in pool (matches next_selection)."""
        self._index += 1
        if self._index >= len(self._pool):
            # Wrap around (shouldn't happen if reset is called at 80)
            self._index = 0

    def needs_reset(self) -> bool:
        """True if counter >= nsamples (time to reset_selection)."""
        return self._counter >= self.nsamples

    @property
    def counter(self) -> int:
        """Number of accepted samples since last reset."""
        return self._counter
