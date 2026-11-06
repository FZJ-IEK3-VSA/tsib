import pandas as pd
import pyomo.environ as p
import pytest

from tsib.thermal.model5R1C import Building5R1C


def thermal_model(tilt=0, thermal_class="very light"):
    cfg = {
        "A_ref": 100,
        "thermalClass": thermal_class,
        "A_Door_1": 0,
        "U_Door_1": 0,
        "comfortT_lb": 20,
        "comfortT_ub": 26,
        "roofTilt": tilt,
        "roofOrientation": 0,
        "latitude": 50,
        "longitude": 8,
        "F_sh_vert": 1,
        "F_sh_hor": 1,
        "F_f": 0,
        "F_w": 1,
        "c_m": 80,
    }
    for direction in ["North", "East", "South", "West", "Horizontal"]:
        cfg["A_Window_" + direction] = 0
    cfg["weather"] = pd.DataFrame(
        {"DNI": [0.0], "GHI": [0.0], "DHI": [0.0]},
        index=pd.date_range("2020-06-21 12:00", periods=1, tz="UTC"),
    )
    model = Building5R1C(cfg)
    state = p.ConcreteModel()
    state.bInsul = {"Windows": ["base"], "Walls": ["base"], "Roof": ["base"]}
    state.bU = {"Windows": {"base": 1}}
    state.bH = {"Walls": {"base": 0.2}, "Roof": {"base": 0.2}}
    state.bg_gl = {"base": 0.5}
    state.profiles = {}
    model._init5R1C(state)
    return model, state


@pytest.mark.parametrize("tilt, factor", [(0, 1), (60, 0.75), (90, 0.5)])
def test_dark_sky_profiles_respect_surface_hemisphere(tilt, factor):
    _, state = thermal_model(tilt)
    # Same opaque conductance and zero shortwave for walls/roof. Full hemisphere
    # reference is a horizontal surface; a vertical one sees half the sky.
    assert state.profiles["bQ_sol_Wallsbase"][0] == pytest.approx(-0.198)
    assert state.profiles["bQ_sol_Roofbase"][0] == pytest.approx(-0.396 * factor)
    assert state.profiles["bQ_sol_Windowsbase"][0] == 0
