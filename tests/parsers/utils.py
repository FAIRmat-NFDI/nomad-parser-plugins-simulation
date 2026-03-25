"""Shared test utilities for parser tests."""

import pytest


def approx(value, abs=0, rel=1e-6):
    """Helper for floating-point comparisons."""
    return pytest.approx(value, abs=abs, rel=rel)
