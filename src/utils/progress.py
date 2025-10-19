"""Spinner/progress helpers for running shell commands in the pipeline."""

from __future__ import annotations

import codecs
import os
import select
import subprocess
import sys
import threading
import time
from typing import Iterable, Sequence

from utils.logger import PipelineLogger


class ProgressRunner:
    """Utility to execute commands while rendering rich terminal feedback."""

    def __init__(self, logger: PipelineLogger, spinner_chars: str = '⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏') -> None:
        self.logger = logger
        self.spinner_chars = spinner_chars

    # ------------------------------------------------------------------ #
    # Internal helpers
    # ------------------------------------------------------------------ #
    def _tail_log(self) -> None:
        print("Last 20 lines:", file=sys.stderr)
        for line in self.logger.tail(20):
            print("  " + line.rstrip(), file=sys.stderr)

    # ------------------------------------------------------------------ #
    # Public API
    # ------------------------------------------------------------------ #
    def run_with_spinner(self, cmd: Sequence[str] | str, message: str, *, shell: bool = False) -> subprocess.CompletedProcess:
        """Run a short-lived command, showing a spinner until completion."""

        self.logger.section(message)

        stop_event = threading.Event()

        def spin() -> None:
            index = 0
            while not stop_event.is_set():
                char = self.spinner_chars[index % len(self.spinner_chars)]
                print(f"\r{char} {message}", end="", file=sys.stderr)
                sys.stderr.flush()
                index += 1
                time.sleep(0.1)

        thread = threading.Thread(target=spin, daemon=True)
        thread.start()

        try:
            result = subprocess.run(
                cmd,
                shell=shell,
                capture_output=True,
                text=True,
                check=True,
            )
            self.logger.write(result.stdout or "")
            self.logger.write(result.stderr or "")
            stop_event.set()
            thread.join()
            print(f"\r\033[0;36m✓\033[0m \033[0;36m{message}\033[0m", file=sys.stderr)
            return result
        except subprocess.CalledProcessError:
            stop_event.set()
            thread.join()
            print(f"\r\033[0;31m✗\033[0m \033[0;31m{message}\033[0m", file=sys.stderr)
            print("Error occurred. Check log:", self.logger.path, file=sys.stderr)
            self._tail_log()
            sys.exit(1)

    def run_with_progress(self, cmd: Sequence[str], message: str) -> None:
        """Run a command that emits progress to stdout, preserving live updates."""

        self.logger.section(message)

        spinner_interval = 0.1
        spinner_index = 0
        last_refresh = 0.0
        display_progress = ""
        progress_fragment = ""

        def render_frame(spinner_char: str, progress_line: str) -> None:
            formatted_progress = f"  {progress_line}" if progress_line else ""
            sys.stderr.write("\033[2A\r")
            sys.stderr.write(f"{spinner_char} {message}\033[K\n")
            if formatted_progress:
                sys.stderr.write(f"{formatted_progress}\033[K\n")
            else:
                sys.stderr.write("\033[K\n")
            sys.stderr.flush()

        sys.stderr.write("\033[?25l")
        sys.stderr.write(f"{self.spinner_chars[0]} {message}\n\n")
        sys.stderr.flush()

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        stdout_fd = process.stdout.fileno() if process.stdout else None
        stderr_fd = process.stderr.fileno() if process.stderr else None
        stdout_done = stdout_fd is None
        stderr_done = stderr_fd is None

        stdout_decoder = codecs.getincrementaldecoder("utf-8")() if not stdout_done else None
        stderr_decoder = codecs.getincrementaldecoder("utf-8")() if not stderr_done else None

        exit_code = None

        try:
            with open(self.logger.path, 'a', encoding='utf-8') as log_handle:
                while True:
                    read_fds: list[int] = []
                    if not stdout_done:
                        read_fds.append(stdout_fd)
                    if not stderr_done:
                        read_fds.append(stderr_fd)

                    if read_fds:
                        ready, _, _ = select.select(read_fds, [], [], spinner_interval)
                    else:
                        ready = []
                        time.sleep(spinner_interval)

                    now = time.time()
                    if now - last_refresh >= spinner_interval:
                        spinner_char = self.spinner_chars[spinner_index % len(self.spinner_chars)]
                        spinner_index += 1
                        render_frame(spinner_char, display_progress)
                        last_refresh = now

                    for fd in ready:
                        if fd == stdout_fd and not stdout_done:
                            chunk = os.read(stdout_fd, 4096)
                            if not chunk:
                                stdout_done = True
                                remaining = stdout_decoder.decode(b"", final=True) if stdout_decoder else ""
                                if remaining:
                                    log_handle.write(remaining)
                                    for ch in remaining:
                                        if ch in ("\r", "\n"):
                                            if progress_fragment:
                                                display_progress = progress_fragment
                                            progress_fragment = ""
                                        else:
                                            progress_fragment += ch
                                    if progress_fragment:
                                        display_progress = progress_fragment
                                        progress_fragment = ""
                            else:
                                text = stdout_decoder.decode(chunk)
                                if text:
                                    log_handle.write(text)
                                    for ch in text:
                                        if ch in ("\r", "\n"):
                                            if progress_fragment:
                                                display_progress = progress_fragment
                                            progress_fragment = ""
                                        else:
                                            progress_fragment += ch
                                    if progress_fragment:
                                        display_progress = progress_fragment
                        elif fd == stderr_fd and not stderr_done:
                            chunk = os.read(stderr_fd, 4096)
                            if not chunk:
                                stderr_done = True
                                remaining = stderr_decoder.decode(b"", final=True) if stderr_decoder else ""
                                if remaining:
                                    log_handle.write(remaining)
                            else:
                                text = stderr_decoder.decode(chunk)
                                if text:
                                    log_handle.write(text)

                    if process.poll() is not None and stdout_done and stderr_done:
                        break

            if progress_fragment:
                display_progress = progress_fragment

            exit_code = process.wait()

            status_line = (
                f"\033[0;36m✓\033[0m \033[0;36m{message}\033[0m"
                if exit_code == 0
                else f"\033[0;31m✗\033[0m \033[0;31m{message}\033[0m"
            )

            formatted_progress_line = f"  {display_progress}" if exit_code != 0 and display_progress else ""

            sys.stderr.write("\033[2A\r")
            sys.stderr.write(f"{status_line}\033[K\n")
            if formatted_progress_line:
                sys.stderr.write(f"{formatted_progress_line}\033[K\n")
            else:
                sys.stderr.write("\033[K")
                sys.stderr.write("\r")
            sys.stderr.flush()

        finally:
            sys.stderr.write("\033[?25h")
            sys.stderr.flush()

        if exit_code != 0:
            print("Error occurred. Check log:", self.logger.path, file=sys.stderr)
            self._tail_log()
            sys.exit(1)

