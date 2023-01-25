"""
Unit and regression test for the pipeline package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pipeline


def test_pipeline_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pipeline" in sys.modules
