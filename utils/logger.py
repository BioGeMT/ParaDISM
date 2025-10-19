"""Logging helpers for the homologous-region mapper pipeline."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Iterable


class PipelineLogger:
    """Simple timestamped logger that writes to a single file."""

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.touch(exist_ok=True)

    def section(self, message: str) -> None:
        """Write a bannered section entry."""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        with self.path.open('a', encoding='utf-8') as handle:
            handle.write(f"\n{'=' * 48}\n")
            handle.write(f"[{timestamp}] {message}\n")
            handle.write(f"{'=' * 48}\n")

    def write(self, text: str) -> None:
        """Append arbitrary text to the log."""
        if not text:
            return
        with self.path.open('a', encoding='utf-8') as handle:
            handle.write(text)

    def tail(self, lines: int = 20) -> Iterable[str]:
        """Yield the last ``lines`` lines from the log file."""
        if lines <= 0:
            return []
        with self.path.open('r', encoding='utf-8') as handle:
            content = handle.readlines()
        return content[-lines:]
