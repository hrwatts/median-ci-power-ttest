from __future__ import annotations

import math
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import median_power_ci as core


def test_chi_square_endpoint_order_and_finiteness() -> None:
    spec = core.interval_spec(n=8, left_tail=0.01)
    assert 0.0 <= spec.A < spec.B
    assert math.isfinite(spec.A)
    assert math.isfinite(spec.B)


def test_one_sided_power_in_unit_interval() -> None:
    value = core.power_one_sample_t(n=8, alpha=0.05, kappa=0.5)
    assert 0.0 <= float(value) <= 1.0
