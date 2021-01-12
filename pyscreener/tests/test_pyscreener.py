"""
Unit and regression test for the pyscreener package.
"""

# Import package, test suite, and other packages as needed
import pyscreener
import pytest
import sys

def test_pyscreener_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "pyscreener" in sys.modules
