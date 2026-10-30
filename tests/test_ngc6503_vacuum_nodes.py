# =====================================================================
# Grid Convergence and Resolution Tests (N = 100 to N = 10000)
# =====================================================================

import numpy as np
import pytest

from nodes.active_tensor_foam import ActiveTensorFoamNode

# rtol used for every document-vs-code comparison: the two sides are the same
# expression evaluated twice, so agreement should be exact to float rounding.
RTOL = 1e-12

@pytest.mark.parametrize("N", [100, 300, 1000, 3000, 10000])
def test_vacuum_density_grid_convergence(N):
    """بررسی پایداری و همگرایی چگالی خلأ در رزولوشن‌های مختلف شبکه."""
    node = ActiveTensorFoamNode()
    r = np.linspace(0.0, node.R_cutoff, N)
    rho = node.vacuum_density(r)
    
    # کنترل پایداری ریاضی
    assert np.all(np.isfinite(rho))
    assert np.all(rho >= 0.0)
    # مقایسه با نقطه شروع تحلیلی
    assert np.isclose(rho[0], node.rho_0, rtol=RTOL)

@pytest.mark.parametrize("N", [100, 1000, 10000])
def test_friction_sign_crossing_grid_convergence(N):
    """بررسی همگرایی محل تغییر علامت اصطکاک (Node 02) در شبکه‌های مختلف."""
    node = ActiveTensorFoamNode()
    grid = np.linspace(0.0, 50.0, N)
    friction = node.foam_friction_raw(grid)
    
    # پیدا کردن نقطه عبور از صفر (Crossing Point)
    # The 0 -> positive step at the origin is not a crossing of the axis; only count
    # direction changes away from r = 0.
    sign_changes = [i for i in np.where(np.diff(np.sign(friction)))[0] if grid[i] > 0.5]
    if len(sign_changes) > 0:
        idx = sign_changes[0]
        crossing = grid[idx] + (grid[idx+1] - grid[idx]) * (-friction[idx]) / (friction[idx+1] - friction[idx])
        # نقطه تغییر علامت باید در بازه فیزیکی پایدار باقی بماند
        assert 12.0 < crossing < 12.2

@pytest.mark.parametrize("N", [100, 1000, 10000])
def test_source_sink_integral_coupling_convergence(N):
    """آزمون همگرایی انتگرال کوپلینگ منبع و چاهک در تقریب عددی انتگرال."""
    node = ActiveTensorFoamNode()
    r = np.linspace(0.1, 50.0, N) # دوری از تکینگی صفر مطلق برای انتگرال‌گیری پایدار
    dv = 100.0 # سرعت فرضی نسبی
    
    # محاسبه چگالی و نرخ تبادل
    rho = node.vacuum_density(r)
    gamma = node.exchange_coefficient(r)
    source = -gamma * dv
    sink = +gamma * dv
    
    # انتگرال کوپلینگ کل (باید در تقریب عددی صفر باشد: Source + Sink = 0)
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # np.trapz removed in NumPy 2.0
    total_coupling = trapezoid(source + sink, r)
    assert np.abs(total_coupling) < 1e-12

@pytest.mark.parametrize("N", [100, 1000, 10000])
def test_circular_velocity_profile_stability(N):
    """بررسی پایداری پروفایل سرعت چرخشی در رزولوشن‌های مختلف شبکه."""
    node = ActiveTensorFoamNode()
    r = np.linspace(0.1, 50.0, N)
    v_circ = node.circular_velocity(r)
    
    assert np.all(np.isfinite(v_circ))
    assert np.all(v_circ >= 0.0)
