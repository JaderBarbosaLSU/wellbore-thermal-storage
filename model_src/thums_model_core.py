"""THUMS model -- flat, self-contained.

Assembled from the repository source by docs/assemble.py. The functions in the
`cycles`, `network`, `front_closed_form`, `march_legacy` and `hydraulics`
sections are the ORIGINAL v0.1 notebook functions, extracted verbatim; only
warning prints and duplicated imports were removed. Nothing was rewritten, so
this module reproduces the conference-paper numbers exactly.

The `front_energy_balance`, `march`, `sizing` and `run` sections are the v0.2
corrections and are new code.
"""
import numpy as np
import CoolProp.CoolProp as CP
from CoolProp.CoolProp import PropsSI
from math import log, pi
import pandas as pd



# ==========================================================================
# properties
# ==========================================================================


def make_cp_state(fluid: str):
    # HEOS is usually fine; if you use REFPROP, swap backend.
    return CP.AbstractState("HEOS", fluid)


def get_props_state(state, T: float, P: float):
    # Returns rho [kg/m3], mu [Pa.s], cp [J/kg.K], k [W/m.K]
    state.update(CP.PT_INPUTS, P, T)
    rho = state.rhomass()
    mu  = state.viscosity()
    cp  = state.cpmass()
    k   = state.conductivity()
    return rho, mu, cp, k


def compute_h_i_from_state(state, r_i, m_dot_d_well, T, P):
    # Same logic as your compute_h_i_coolprop, but using a pre-built state
    try:
        rho, mu, cp, k = get_props_state(state, float(T), float(P))
    except Exception:
        return 1e-9, 0.0, 0.0, 0.0

    D_h = 2.0 * r_i
    A = pi * r_i**2
    u = m_dot_d_well / (rho * A) if (rho * A) > 1e-30 else 0.0

    Re = rho * u * D_h / mu if mu > 1e-30 else 0.0
    Pr = cp * mu / k if (k > 1e-30 and mu > 1e-30) else 0.0

    if Re < 2300.0:
        Nu = 3.66
    else:
        # Your turbulent correlation
        try:
            f = (0.79 * np.log(Re) - 1.64)**(-2)
            if f > 1e-30:
                Nu = (f / 8.0) * (Re - 1000.0) * Pr / (1.0 + 12.7 * np.sqrt(f / 8.0) * (Pr**(2/3) - 1.0))
            else:
                Nu = 3.66
        except Exception:
            Nu = 3.66

    h_i = Nu * k / D_h if (D_h > 1e-30 and k > 1e-30) else 1e-9
    return h_i, Re, Pr, cp

########################################
########################################
########################################



# ==========================================================================
# cycles
# ==========================================================================


def double_stage_rankine(fluid, T_1e, T_3e):
  """
  Calculates parameters for a two-stage Rankine cycle with regeneration/reheating

  Args:
    fluid: The name of the working fluid.
    T_1e: Temperature of state 1e in Kelvin.
    T_3e: Temperature of state 3e in Kelvin.

  Returns:
    A dictionary containing calculated Rankine cycle parameters,
    including enthalpy, entropy, temperature, and pressure for each state,
    the Rankine efficiency, and saturation pressures.
  """
  # State 3e: Saturated vapor at the 1st stage (turbine) inlet
  p_evape = CP.PropsSI('P', 'T', T_3e, 'Q', 1, fluid) # Condenser pressure (same as 3e)
  p_3e = p_evape
  h_3e = CP.PropsSI('H', 'T', T_3e, 'Q', 1, fluid)/1000.
  s_3e = CP.PropsSI('S', 'T', T_3e, 'Q', 1, fluid)

  # State 1e: Saturated liquid at condenser outlet
  p_conde = CP.PropsSI('P', 'T', T_1e, 'Q', 0, fluid) # Evaporator pressure
  p_1e = p_conde
  h_1e = CP.PropsSI('H', 'T', T_1e, 'Q', 0, fluid)/1000.
  s_1e = CP.PropsSI('S', 'T', T_1e, 'Q', 0, fluid)

  # Intermediate pressure
#  p_inte = (p_evape * p_conde)**(1./2.)
  p_inte = (p_evape + p_conde)*(1./2.)

  # State 5e: 1st stage turbine outlet
  p_5e = p_inte
  s_5e = s_3e
  h_5e = CP.PropsSI('H', 'S', s_5e, 'P', p_5e, fluid)/1000.
  T_5e = CP.PropsSI('T', 'S', s_5e, 'P', p_5e, fluid)

  # State 6e: 2nd stage turbine inlet
  p_6e = p_5e
  T_6e = T_3e #ideal reheat
  h_6e = CP.PropsSI('H', 'T', T_6e, 'P', p_6e, fluid)/1000.
  s_6e = CP.PropsSI('S', 'T', T_6e, 'P', p_6e, fluid)

  # State 4e: 2nd stage turbine outlet
  p_4e = p_conde
  s_4e = s_6e
  h_4e = CP.PropsSI('H', 'S', s_4e, 'P', p_4e, fluid)/1000.
  T_4e = CP.PropsSI('T', 'S', s_4e, 'P', p_4e, fluid)

  # State 2e: Pump 1 outlet
  p_2e = p_evape
  s_2e = s_1e
  h_2e = CP.PropsSI('H', 'S', s_2e, 'P', p_2e, fluid)/1000.
  T_2e = CP.PropsSI('T', 'S', s_2e, 'P', p_2e, fluid)

  # State 7e: Pump 2 inlet
  p_7e = p_inte
  h_7e = CP.PropsSI('H', 'P', p_7e, 'Q', 0, fluid)/1000.
  s_7e = CP.PropsSI('S', 'P', p_7e, 'Q', 0, fluid)
  T_7e = CP.PropsSI('T', 'P', p_7e, 'Q', 0, fluid)

  # State 8e: Pump 2 outlet
  p_8e = p_evape
  s_8e = s_7e
  h_8e = CP.PropsSI('H', 'P', p_8e, 'S', s_8e, fluid)/1000.
  T_8e = CP.PropsSI('T', 'P', p_8e, 'S', s_8e, fluid)

  # State 9e:
  p_9e = p_evape
  T_9e = T_7e
  h_9e = CP.PropsSI('H', 'T', T_9e, 'P', p_9e, fluid)/1000.
  s_9e = CP.PropsSI('S', 'T', T_9e, 'P', p_9e, fluid)

  # Energy balance in the regenerator
  y_frac = (h_5e - h_7e) / ((h_9e - h_2e) + (h_5e - h_7e))

  # State 10e:
  p_10e = p_evape
  h_10e = y_frac * h_9e + (1. - y_frac) * h_8e
  s_10e = CP.PropsSI('S', 'H', h_10e*1000, 'P', p_10e, fluid)
  T_10e = CP.PropsSI('T', 'H', h_10e*1000, 'P', p_10e, fluid)

  # cycle thermal efficiency
  numerator_eff = (h_3e - h_10e) + y_frac * (h_6e - h_5e) - y_frac * (h_4e - h_1e)
  denominator_eff = (h_3e - h_10e) + y_frac * (h_6e - h_5e)
  if abs(denominator_eff) > 1e-9:
      rank_eff = numerator_eff / denominator_eff
  else:
      rank_eff = float('inf') # Assign large value if denominator is zero or near zero

  return {
      'T_1e': T_1e, 'p_1e': p_1e, 'h_1e': h_1e, 's_1e': s_1e,
      'T_2e': T_2e, 'p_2e': p_2e, 'h_2e': h_2e, 's_2e': s_2e,
      'T_3e': T_3e, 'p_3e': p_3e, 'h_3e': h_3e, 's_3e': s_3e,
      'T_4e': T_4e, 'p_4e': p_4e, 'h_4e': h_4e, 's_4e': s_4e,
      'T_5e': T_5e, 'p_5e': p_5e, 'h_5e': h_5e, 's_5e': s_5e,
      'T_6e': T_6e, 'p_6e': p_6e, 'h_6e': h_6e, 's_6e': s_6e,
      'T_7e': T_7e, 'p_7e': p_7e, 'h_7e': h_7e, 's_7e': s_7e,
      'T_8e': T_8e, 'p_8e': p_8e, 'h_8e': h_8e, 's_8e': s_8e,
      'T_9e': T_9e, 'p_9e': p_9e, 'h_9e': h_9e, 's_9e': s_9e,
      'T_10e': T_10e, 'p_10e': p_10e, 'h_10e': h_10e, 's_10e': s_10e,
      'rank_eff': rank_eff,
      'p_evape': p_evape,
      'p_conde': p_conde,
      'p_inte': p_inte,
      'y_frac': y_frac # Added y_frac to the return dictionary
  }

# @title


def two_stage_htheatpump_2regs(refrig, T_13h, T_2h, DT_sub):
  """
  Calculates parameters for a two-stage high-temp heat pump with a liquid separator and 2 regenerators

  Args:
    refrig: The name of the refrigerant.
    T_13h: Temperature of state 13h in Kelvin.
    T_2h: Temperature of state 2h in Kelvin.
    DT_sub: Refrigerant subcooling degree at condenser outlet.

  Returns:
    A dictionary containing calculated Rankine cycle parameters.
  """
  # State 2h: Saturated Vapor at T_2h (Condenser Outlet)
  p_condh = CP.PropsSI('P', 'T', T_2h, 'Q', 0, refrig) # Condensing pressure
  p_2h = p_condh
  h_2h = CP.PropsSI('H', 'T', T_2h, 'Q', 1, refrig)/1000.
  s_2h = CP.PropsSI('S', 'T', T_2h, 'Q', 1, refrig)

  # State 3h: Subcooled Liquid at (T_2h - DT_sub) and p_condh (Condenser Outlet with subcooling)
  T_3h = T_2h - DT_sub
  p_3h = p_condh
  h_3h = CP.PropsSI('H', 'T', T_3h, 'P', p_condh, refrig)/1000.
  s_3h = CP.PropsSI('S', 'T', T_3h, 'P', p_condh, refrig)

  # State 13h: Saturated vapor at evaporating pressure
  p_evaph = CP.PropsSI('P', 'T', T_13h, 'Q', 1, refrig) # Evaporating pressure
  p_13h = p_evaph
  h_13h = CP.PropsSI('H', 'T', T_13h, 'Q', 1, refrig)/1000.
  s_13h = CP.PropsSI('S', 'T', T_13h, 'Q', 1, refrig)


  # Intermediate pressure
  p_inth = (p_evaph * p_condh)**(1./2.)
#  p_inth = (p_evaph + p_condh)*(1./2.)

  # State 1h: At T_1h and p_evaph (Evaporator Outlet)
  s_1h = s_2h # Isentropic compression to 2h
  p_1h = p_evaph
  T_1h = CP.PropsSI('T', 'P', p_1h, 'S', s_1h, refrig)
  h_1h = CP.PropsSI('H', 'T', T_1h, 'P', p_evaph, refrig)/1000.

  # State 11h: Low pressure compressor outlet (p_inth, isentropic from 1h)
  s_11h = s_1h # Isentropic compression
  p_11h = p_inth
  T_11h = CP.PropsSI('T', 'S', s_11h, 'P', p_11h, refrig)
  h_11h = CP.PropsSI('H', 'S', s_11h, 'P', p_11h, refrig)/1000.

  # State 5h: Low Pressure Compressor outlet
  T_5h = T_11h
  p_5h = p_11h
  h_5h = h_11h
  s_5h = s_11h

  # State 6h: High Pressure Compressor inlet
  T_6h = T_11h
  p_6h = p_11h
  h_6h = h_11h
  s_6h = s_11h

  # State 9h
  p_9h = p_inth
  T_9h = CP.PropsSI('T', 'P', p_9h, 'Q', 0, refrig) # Saturated liquid temp at p_inth
  h_9h = CP.PropsSI('H', 'P', p_9h, 'Q', 0, refrig)/1000.
  s_9h = CP.PropsSI('S', 'P', p_9h, 'Q', 0, refrig)

  # State 10h
  p_10h = p_inth
  T_10h = CP.PropsSI('T', 'P', p_10h, 'Q', 1, refrig) # Saturated vapor temp at p_inth
  h_10h = CP.PropsSI('H', 'P', p_10h, 'Q', 1, refrig)/1000.
  s_10h = CP.PropsSI('S', 'P', p_10h, 'Q', 1, refrig)

  # State 14h
  p_14h = p_evaph
  T_14h = CP.PropsSI('T', 'P', p_14h, 'Q', 0, refrig) # Saturated liquid temp at p_evaph
  h_14h = CP.PropsSI('H', 'P', p_14h, 'Q', 0, refrig)/1000.
  s_14h = CP.PropsSI('S', 'P', p_14h, 'Q', 0, refrig)

  # IHX-1 effectiveness (based on T_3h and T_10h)
  # Ensure denominator is non-zero
  denominator_epsilon1 = (T_3h - T_10h)
  if abs(denominator_epsilon1) > 1e-9:
      epsilon_IHX_1 = (T_11h - T_10h) / denominator_epsilon1
  else:
      epsilon_IHX_1 = 1.0 # Assume 100% effectiveness if temperature difference is zero

  # CP_ratio_1:
  try:
      cpf_1 = CP.PropsSI('C', 'P', p_condh, 'Q', 0, refrig)
      cpv_1 = CP.PropsSI('C', 'P', p_inth, 'Q', 1, refrig)
      # Avoid division by zero
      if abs(cpf_1) > 1e-9:
           cpr_1 = cpv_1 / cpf_1
      else:
           cpr_1 = 1.0 # Assume ratio is 1 if cpf is zero
  except ValueError:
      cpr_1 = 1.0

  # --- the IHX network, closed by an exact enthalpy balance ---------------
  # The two regenerators are closed by ENERGY, not by effectiveness-times-
  # cp-ratio correlations. Which stream passes through each one matters:
  #
  #   IHX-1   liquid 3h -> 12h  (1 kg)   heats the FLASHED VAPOUR, x_8h,
  #                                      from 10h to 11h
  #   IHX-2   liquid 12h -> 7h  (1 kg)   heats the SUCTION stream, which is
  #                                      what left the separator as LIQUID,
  #                                      (1 - x_8h), from 13h to 1h
  #
  # x_8h appears on both sides, so the three relations are solved as a fixed
  # point. The map is linear and strongly contracting here, so it converges
  # in a handful of steps.
  #
  # Until v0.9 IHX-2 used x_4h -- the vapour QUALITY at the evaporator inlet
  # -- where the flow SPLIT (1 - x_8h) belongs. A quality is not a flow
  # fraction, and the cycle did not close: IHX-2 took 9.19 kJ/kg out of the
  # liquid and put 21.76 kJ/kg into the vapour, and the condenser reported
  # 271.43 kJ/kg against 258.33 kJ/kg of work plus evaporator heat. That is
  # a 4.8 % creation of energy, and it inflated the COP.
  p_7h = p_condh
  dh_IHX1 = h_11h - h_10h            # per kg of flashed vapour
  dh_IHX2 = h_1h - h_13h             # per kg of suction (low-stage) flow
  x_8h = 0.3
  for _n_ihx in range(200):
      h_12h = h_3h - x_8h * dh_IHX1
      h_7h = h_12h - (1.0 - x_8h) * dh_IHX2
      x_new = (h_7h - h_9h) / (h_10h - h_9h)          # h_8h = h_7h
      converged = abs(x_new - x_8h) < 1e-12
      x_8h = x_new
      if converged:
          break
  else:
      raise RuntimeError(
          "two_stage_htheatpump_2regs: the IHX enthalpy balance did not "
          f"converge in 200 iterations (x_8h = {x_8h!r}). The pressure "
          "levels or the subcooling are probably inconsistent.")
  if not np.isfinite(x_8h) or not (0.0 <= x_8h <= 1.0):
      raise RuntimeError(
          f"two_stage_htheatpump_2regs: separator vapour fraction "
          f"x_8h = {x_8h:.6g} is outside [0, 1], so there is no physical "
          "flash split. Check T_2h, T_13h and DT_sub.")
  h_12h = h_3h - x_8h * dh_IHX1
  h_7h = h_12h - (1.0 - x_8h) * dh_IHX2
  h_8h = h_7h                        # isenthalpic expansion into the separator
  h_4h = h_9h                        # isenthalpic expansion into the evaporator
  T_7h = CP.PropsSI('T', 'H', h_7h * 1000., 'P', p_condh, refrig)
  T_12h = CP.PropsSI('T', 'H', h_12h * 1000., 'P', p_condh, refrig)

  # Both of these are now DIAGNOSTICS rather than inputs. State 1h is fixed
  # by the isentropic compression back from 2h, so the effectiveness IHX-2
  # would need in order to deliver it is an OUTPUT of the cycle, and worth
  # reading: a value above 1 would mean the regenerator cannot exist.
  x_4h = (h_4h - h_14h) / (h_13h - h_14h)
  epsilon_IHX_2 = (T_1h - T_13h) / (T_12h - T_13h)
  cpr_2 = (CP.PropsSI('C', 'P', p_evaph, 'Q', 1, refrig)
           / CP.PropsSI('C', 'P', p_condh, 'Q', 0, refrig))


  # State 7h (using the converged T_7h)
  # h_7h and s_7h were calculated inside the loop, re-calculate if needed outside
  try:
       h_7h = CP.PropsSI('H', 'T', T_7h, 'P', p_7h, refrig)/1000.
       s_7h = CP.PropsSI('S', 'T', T_7h, 'P', p_7h, refrig)
  except ValueError:
       h_7h, s_7h = np.nan, np.nan


  # State 12h (using the converged T_12h)
  # T_12h was calculated inside the loop, re-calculate if needed outside
  try:
       p_12h = p_condh
       h_12h = CP.PropsSI('H', 'T', T_12h, 'P', p_12h, refrig)/1000.
       s_12h = CP.PropsSI('S', 'T', T_12h, 'P', p_12h, refrig)
  except ValueError:
       h_12h, s_12h = np.nan, np.nan


  # State 8h (using the converged h_8h=h_7h and p_8h=p_inth)
  # h_8h was calculated inside the loop, re-calculate if needed outside
  try:
       p_8h = p_inth
       T_8h = CP.PropsSI('T', 'P', p_8h, 'H', h_8h*1000., refrig)
       s_8h = CP.PropsSI('S', 'P', p_8h, 'H', h_8h*1000., refrig)
  except ValueError:
       T_8h, s_8h = np.nan, np.nan


  # State 4h (using the converged h_4h=h_9h and p_4h=p_evaph)
  # h_4h was calculated inside the loop, re-calculate if needed outside
  try:
       p_4h = p_evaph
       T_4h = CP.PropsSI('T', 'P', p_4h, 'H', h_4h*1000., refrig)
       s_4h = CP.PropsSI('S', 'P', p_4h, 'H', h_4h*1000., refrig)
  except ValueError:
       T_4h, s_4h = np.nan, np.nan


  # HP COP
  # Ensure denominator is non-zero before division and check for NaN values
  # Need to use the final converged x_8h and h_5h, h_1h
  denominator_cop = (h_2h - h_1h) - x_8h * (h_5h - h_1h)
  if abs(denominator_cop) > 1e-9 and not np.isnan(h_2h) and not np.isnan(h_1h) and not np.isnan(x_8h) and not np.isnan(h_5h):
      hp_cop = (h_2h - h_3h) / denominator_cop
  else:
      hp_cop = float('inf') # Assign a large value if denominator is zero or near zero or NaN


  return {
      'T_1h': T_1h, 'p_1h': p_1h, 'h_1h': h_1h, 's_1h': s_1h,
      'T_2h': T_2h, 'p_2h': p_2h, 'h_2h': h_2h, 's_2h': s_2h,
      'T_3h': T_3h, 'p_3h': p_3h, 'h_3h': h_3h, 's_3h': s_3h,
      'T_4h': T_4h, 'p_4h': p_4h, 'h_4h': h_4h, 's_4h': s_4h,
      'T_5h': T_5h, 'p_5h': p_5h, 'h_5h': h_5h, 's_5h': s_5h,
      'T_6h': T_6h, 'p_6h': p_6h, 'h_6h': h_6h, 's_6h': s_6h,
      'T_7h': T_7h, 'p_7h': p_7h, 'h_7h': h_7h, 's_7h': s_7h,
      'T_8h': T_8h, 'p_8h': p_8h, 'h_8h': h_8h, 's_8h': s_8h,
      'T_9h': T_9h, 'p_9h': p_9h, 'h_9h': h_9h, 's_9h': s_9h,
      'T_10h': T_10h, 'p_10h': p_10h, 'h_10h': h_10h, 's_10h': s_10h,
      'T_11h': T_11h, 'p_11h': p_11h, 'h_11h': h_11h, 's_11h': s_11h,
      'T_12h': T_12h, 'p_12h': p_12h, 'h_12h': h_12h, 's_12h': s_12h,
      'T_13h': T_13h, 'p_13h': p_13h, 'h_13h': h_13h, 's_13h': s_13h,
      'hp_cop': hp_cop,
      'p_evaph': p_evaph,
      'p_condh': p_condh, # Changed from p_evaph to p_condh
      'p_inth': p_inth,
      'epsilon_IHX_1': epsilon_IHX_1,
      'epsilon_IHX_2': epsilon_IHX_2,
      'x_8h': x_8h,
  }

# @title
#import numpy as np
#import CoolProp.CoolProp as CP
#import matplotlib.pyplot as plt



# ==========================================================================
# network
# ==========================================================================


def compute_T_re(T0, T1, T_m, r_i, r_e, k_w, h_i, Rf_i_prime, k_m, delta, L_tube, fin_t, fin_L, num_fins):
    """
    Computes the temperature at the external radius of the inner tube (T_re)
    considering the thermal resistances, including the melt layer around finned tubes.

    Args:
        T0: Fluid temperature at the start of the segment [K].
        T1: Fluid temperature at the end of the segment [K].
        T_m: Melting temperature of the PCM [K].
        r_i: Internal radius of the tube [m].
        r_e: External radius of the tube [m].
        k_w: Wall thermal conductivity (also used for fin conductivity) [W/m.K]. # Updated docstring
        h_i: Internal convective heat transfer coefficient [W/m.K].
        Rf_i_prime: Internal fouling resistance [m²·K/W].
        k_m: PCM thermal conductivity [W/m.K].
        delta: Melt layer thickness [m].
        L_tube: Length of the tube [m].
        fin_t: Thickness of each fin [m].
        fin_L: Length of each fin extending radially outwards [m].
        num_fins: The number of fins.

    Returns:
        Temperature at the external radius of the inner tube [K].
    """
    T_avg = 0.5 * (T0 + T1)

    # R1 remains unchanged
    R1 = Rf_i_prime + 1/h_i + (r_i / k_w) * np.log(r_e / r_i)

    # Calculate h_e (heat transfer coefficient for the outer layer)
    if delta > 1e-9 and r_e > 1e-9: # Check for non-zero delta and r_e
        log_term = np.log(1 + delta / r_e)
        if log_term > 1e-9: # Check for non-zero log term
            h_e = k_m / r_e / log_term
        else:
             h_e = k_m / delta # Use the previous definition for very small delta
    else:
        h_e = 1e9 # Assign a small value, leading to large resistance

    # Calculate fin parameters for efficiency (needed for eta_o)
    fin_Lc = fin_L + fin_t / 2.0
    # Ensure h_e is not zero or negative before taking sqrt, use k_w for fin conductivity
    if h_e > 1e-9 and k_w > 1e-9 and fin_t > 1e-9:
        fin_m = np.sqrt(2 * h_e / k_w / fin_t) # Used k_w instead of KW
    else:
        fin_m = 1e9 # Assign a large value if denominator is near zero

    # Calculate fin efficiency
    fin_m_Lc = fin_m * fin_Lc
    if fin_m_Lc != 0:
        fin_eff = np.tanh(fin_m_Lc) / fin_m_Lc
    else:
        fin_eff = 1.0 # Efficiency is 100% if fin_m * fin_Lc is zero

    # Calculate fin surface area (needed for eta_o)
    A_f = L_tube * (2 * fin_L + fin_t)

    # Calculate total outer surface area (needed for eta_o)
    A_t = 2 * L_tube * (np.pi * r_e + num_fins * fin_L)

    # Calculate overall surface efficiency
    if A_t > 0:
        eta_o = 1 - (num_fins * A_f / A_t) * (1 - fin_eff)
    else:
        eta_o = 1.0

    # Calculate perim_T
    perim_T = 2 * (np.pi * r_e + num_fins * fin_L)

    # Calculate R2
    # Ensure perim_T, h_e, and eta_o are effectively non-zero
    if perim_T > 1e-9 and h_e > 1e-9 and eta_o > 1e-9:
         R2 = (2 * np.pi * r_i) / (perim_T * h_e * eta_o)
    else:
         R2 = float('inf') # Assign infinite resistance if denominator is zero or near zero

    # Avoid division by zero or infinity when calculating T_re
    denominator = R1 + R2
    if denominator == 0 or np.isinf(denominator):
         # If R1 + R2 is zero or infinite, T_re cannot be calculated in this way.
         # This might indicate an issue with resistances (e.g., R2 is inf).
         # In this case, T_re is likely T_avg or T_m depending on the context.
         # Given that R2 being inf means no heat transfer from the outer layer,
         # the temperature at the outer radius should approach the fluid temperature T_avg.
         T_re = T_avg
    else:
         numerator = R1 * (T_avg - T_m)
         T_re = T_avg - numerator / denominator

    return T_re

########################################
########################################
########################################


def compute_U_i(h_i, r_i, r_e, k_w, Rf_i_prime, k_m, delta, L_tube, fin_t, fin_L, num_fins):
    """
    Computes the overall heat transfer coefficient based on internal surface area,
    including thermal resistances for internal convection, wall conduction, fouling,
    and the outer finned surface with a melt layer and fin efficiency.

    Args:
        h_i: Internal convective heat transfer coefficient [W/m²·K].
        r_i: Internal radius of the tube [m].
        r_e: External radius of the tube [m].
        k_w: Casing wall thermal conductivity (also used for fin conductivity) [W/m·K]. # Updated docstring
        Rf_i_prime: Internal fouling resistance [m²·K/W].
        k_m: PCM thermal conductivity [W/m.K].
        delta: Melt layer thickness [m].
        L_tube: Length of the tube [m].
        fin_t: Thickness of each fin [m].
        fin_L: Length of each fin extending radially outwards [m].
        num_fins: The number of fins.

    Returns:
        Overall heat transfer coefficient based on internal surface area [W/m²·K].
    """
    R_conv_internal = 1.0 / h_i
    R_wall = (r_i / k_w) * log(r_e / r_i)
    R_fouling = Rf_i_prime

    # Calculate h_e (heat transfer coefficient for the outer layer)
    if delta > 1e-9 and r_e > 1e-9: # Check for non-zero delta and r_e
        log_term = np.log(1 + delta / r_e)
        if log_term > 1e-9: # Check for non-zero log term
            h_e = k_m / r_e / log_term
        else:
             # If log_term is near zero, delta/r_e is very small, meaning very thin layer.
             # This would imply a large h_e, but setting it too large might cause issues.
             # Assign a large, but finite value. Or, consider the limit as delta approaches 0.
             # As delta -> 0, log(1 + delta/r_e) ~ delta/r_e. h_e ~ k_m / r_e / (delta/r_e) = k_m/delta.
             # So the previous definition is consistent with the limit.
             # Let's revert to the previous definition for small delta, or handle it as an infinite h_e (zero resistance).
             # Given the error handling downstream, a large finite number for h_e is safer than infinity here.
             h_e = k_m / delta # Use the previous definition for very small delta
    else:
        # If delta is zero, there is no melt layer, so resistance is infinite
        h_e = 1e9 # Assign a small value, leading to large resistance

    # Calculate fin parameters for efficiency
    fin_Lc = fin_L + fin_t / 2.0
    # Ensure h_e is non-zero or negative before taking sqrt, use k_w for fin conductivity
    if h_e > 1e-9 and k_w > 1e-9 and fin_t > 1e-9:
         fin_m = np.sqrt(2 * h_e / k_w / fin_t) # Used k_w instead of KW
    else:
         fin_m = 1e9 # Assign a large value if denominator is near zero


    # Calculate fin efficiency
    fin_m_Lc = fin_m * fin_Lc
    if fin_m_Lc != 0:
        fin_eff = np.tanh(fin_m_Lc) / fin_m_Lc
    else:
        fin_eff = 1.0 # Efficiency is 100% if fin_m * fin_Lc is zero

    # Calculate fin surface area
    A_f = L_tube * (2 * fin_L + fin_t)

    # Calculate total outer surface area
    A_t = 2 * L_tube * (np.pi * r_e + num_fins * fin_L)

    # Calculate overall surface efficiency
    # Ensure A_t is not zero to avoid division error
    if A_t > 0:
        eta_o = 1 - (num_fins * A_f / A_t) * (1 - fin_eff)
    else:
        eta_o = 1.0 # If total area is zero, efficiency is irrelevant, assume 1

    # Calculate the outer perimeter of the finned tube (used in R_outer_layer calculation)
    perim_T = 2 * pi * r_e + 2 * num_fins * fin_L

    if eta_o > 1e-9 and h_e > 1e-9: # Check if eta_o and h_e are effectively non-zero
         R_outer_layer = (2 * pi * r_i) / (perim_T * eta_o * h_e)
    else:
        R_outer_layer = float('inf') # Assign infinite resistance if heat transfer is blocked

    R_total = R_conv_internal + R_wall + R_fouling + R_outer_layer

    if R_total != 0 and not np.isinf(R_total):
        U_i = 1.0 / R_total
    else:
        U_i = 0.0 # If total resistance is infinite or zero, U_i is 0.0 (no heat transfer)

    return U_i


########################################
########################################
########################################



# ==========================================================================
# front_closed_form
# ==========================================================================


def compute_delta2_fast(r_e, k_m, cp_m, rho_m, h_m, T_m, T_re, t,tol=1e-8, max_iter=30, delta_max=1.0):
    # Physical/degenerate guards
    if r_e <= 0.0 or t <= 0.0:
        return 0.0
    if abs(T_re - T_m) < 1e-9 or T_re <= T_m:
        return 0.0

    alpha_m = k_m / (rho_m * cp_m)
    Fo = alpha_m * t / (r_e * r_e)
    Ph = abs(h_m / (cp_m * (T_re - T_m)))  # >0 since T_re>T_m

    # term(x) = 0.5 x^2 ln x - 0.25 x^2 + 0.25, x=1+delta/r_e
    def F(delta):
        x = 1.0 + delta / r_e
        if x <= 0.0:
            return np.nan
        term = 0.5 * x * x * np.log(x) - 0.25 * x * x + 0.25
        return Fo - Ph * term

    def dF(delta):
        x = 1.0 + delta / r_e
        if x <= 0.0:
            return np.nan
        # dF/dδ = -Ph*(x ln x)/r_e
        return -(Ph / r_e) * (x * np.log(x))

    # Bracket on [0, delta_max]
    a, b = 0.0, float(delta_max)
    fa, fb = F(a), F(b)

    # If no root is bracketed, fall back to 0 (consistent with your physical logic)
    if not np.isfinite(fa) or not np.isfinite(fb) or fa * fb > 0:
        # if fa is already ~0, return 0
        if np.isfinite(fa) and abs(fa) < tol:
            return 0.0
        return 0.0

    # Initial guess: small-delta asymptotic, clipped to bracket
    delta = r_e * np.sqrt(max(0.0, 2.0 * Fo / Ph))
    delta = float(np.clip(delta, a, b))

    # Safeguarded Newton (Newton when it stays in bracket; otherwise bisection)
    for _ in range(max_iter):
        f = F(delta)
        if not np.isfinite(f):
            delta = 0.5 * (a + b)
            continue
        if abs(f) < tol:
            return delta

        df = dF(delta)
        if np.isfinite(df) and abs(df) > 1e-14:
            delta_new = delta - f / df
        else:
            delta_new = np.nan

        # Keep bracket; if Newton goes out, bisect
        if (not np.isfinite(delta_new)) or (delta_new <= a) or (delta_new >= b):
            delta_new = 0.5 * (a + b)

        f_new = F(delta_new)
        if not np.isfinite(f_new):
            delta_new = 0.5 * (a + b)
            f_new = F(delta_new)

        if abs(f_new) < tol:
            return delta_new

        # Update bracket
        if fa * f_new < 0:
            b, fb = delta_new, f_new
        else:
            a, fa = delta_new, f_new

        delta = delta_new

    # best available (midpoint of final bracket)
    return 0.5 * (a + b)

#####################################
#####################################
#####################################



# ==========================================================================
# march_legacy
# ==========================================================================


def temperature_profile_melt(geom_par_vector, T_inlet, N_lay, T_m_lay, k_w, Rf_i_prime, m_dot_d_well, P, fluid, t,
    k_m, cp_m, rho_m, h_m, n_segments, delta_max=1.0,
    delta_tol=1e-5, delta_maxiter=20
):
    """
    Calculates the temperature and melt layer profiles along the well at a given time.

    Args:
        geom_par_vector: An array or list containing [L_well, L_tube, D_well, r_i, r_e, D_i_tube, D_e_tube, fin_t, fin_L, num_fins, num_tubes].
        T_inlet: Inlet temperature of the fluid in Kelvin.
        N_lay: Number of layers in the PCM.
        T_m_lay: Vector (size N_lay + 1) with the melting temperatures of the PCM layers in Kelvin.
        k_w: Wall thermal conductivity (also used for fin conductivity) [W/m.K]. # Updated docstring
        Rf_i_prime: Internal fouling resistance [m²·K/W].
        m_dot_d_well: Mass flow rate of the fluid per well [kg/s].
        P: Pressure of the fluid [Pa].
        fluid: Name of the working fluid.
        t: Time elapsed [s].
        k_m: PCM thermal conductivity [W/m.K].
        cp_m: PCM specific heat capacity [J/kg.K].
        rho_m: PCM density [kg/m³].
        h_m: Latent heat of fusion of the PCM [J/kg].
        n_segments: Number of segments to divide the well length.
        delta_max: Maximum possible melt layer thickness for convergence [m].
        delta_tol: Tolerance for melt layer thickness convergence.
        delta_maxiter: Maximum iterations for melt layer thickness calculation.

    Returns:
        A tuple containing:
            - df: DataFrame with the temperature and other profiles along the well.
            - Q: Heat transfer rate [W].
            - V_melt: Melt volume [m³].
    """
    # Unpack the geometric parameters
    L_well, L_tube, D_well, r_i, r_e, D_i_tube, D_e_tube, fin_t, fin_L, num_fins, num_tubes = geom_par_vector
    L = L_tube # Use L_tube as the length for the profile calculation

    state = make_cp_state(fluid)

    def compute_h_i_coolprop(r_i, m_dot_d_well, fluid, T, P):
        try:
            rho = PropsSI("D", "T", T, "P", P, fluid)
            mu = PropsSI("V", "T", T, "P", P, fluid)
            cp = PropsSI("C", "T", T, "P", P, fluid)
            k = PropsSI("L", "T", T, "P", P, fluid)
        except ValueError:
            # Handle cases where T or P are outside fluid range
            return 1e-9, 0, 0, 0 # Return a small h_i and zero other properties

        D_h = 2 * r_i
        A = pi * r_i**2
        if rho * A > 1e-9:
             u = m_dot_d_well / (rho * A)
        else:
             u = 0 # Avoid division by zero

        if mu > 1e-9: # Avoid division by zero
             Re = rho * u * D_h / mu
        else:
             Re = 0 # If viscosity is zero or near zero, Re is effectively infinite, but 0 avoids division

        if k > 1e-9 and cp > 1e-9 and mu > 1e-9: # Avoid division by zero
             Pr = cp * mu / k
        else:
             Pr = 0 # If k, cp, or mu is zero or near zero, Pr is effectively infinite, but 0 avoids division


        if Re < 2300:
            Nu = 3.66
        else:
            try:
                f = (0.79 * np.log(Re) - 1.64)**(-2)
                # Ensure f is non-negative and non-zero before taking sqrt
                if f > 1e-9:
                    Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * np.sqrt(f / 8) * (Pr**(2/3) - 1))
                else:
                    Nu = 3.66 # Fallback to laminar Nu if f is problematic
            except (ValueError, RuntimeWarning):
                 Nu = 3.66 # Fallback to laminar Nu in case of log(Re) issues for small Re


        if D_h > 1e-9 and k > 1e-9: # Avoid division by zero
            h_i = Nu * k / D_h
        else:
             h_i = 1e-9 # Assign a small value if Dh or k is zero

        return h_i, Re, Pr, cp


    dz = L / n_segments
    z = np.linspace(0, L, n_segments + 1)
    z_lay = np.linspace(0, L, N_lay + 1)    # positions of the PCM layer transitions
    T = [T_inlet]
    NTUs = []
    U_is = []
    deltas = []
    Res = []
    Prs = []
    T_res = []


    # Assign T_m to each z based on coarse grid
    T_m = np.zeros_like(z)

    for j in range(N_lay):
        # interval [z_lay[j], z_lay[j+1])
        mask = (z >= z_lay[j]) & (z < z_lay[j+1])
        T_m[mask] = T_m_lay[j]

    # Ensure the last point gets the last layer temperature
    T_m[-1] = T_m_lay[-1]




    for i in range(n_segments):
        z0, z1 = z[i], z[i + 1]
        T0 = T[-1]

        #h_i, Re, Pr, cp_d = compute_h_i_coolprop(r_i, m_dot_d_well, fluid, T0, P)
        h_i, Re, Pr, cp_d = compute_h_i_from_state(state, r_i, m_dot_d_well, T0, P)


        # Initial guess for delta and T_re
        #delta_guess = 0.001 # Start with a small melt layer guess
        #T1_guess = T0
        #T_re_guess = T0

        # Initial guess for delta and T_re (warm-start)
        if i == 0:
            delta_guess = 0.001
            T_re_guess = T0
        else:
            # warm start from previous segment
            delta_guess = deltas[-1]
            T_re_guess = T_res[-1]

        delta_guess = float(np.clip(delta_guess, 0.0, delta_max))
        T1_guess = T0  # keep simple; you already have an NTU-based T1 update below






        # Predict NTU and T1 with base guess for delta
        # Call compute_U_i without KW
        U_i_guess = compute_U_i(h_i, r_i, r_e, k_w, Rf_i_prime, k_m, delta_guess, L_tube, fin_t, fin_L, num_fins)

        # Avoid division by zero if cp_d is zero or near zero
        if m_dot_d_well * cp_d > 1e-9:
             NTU_guess = (2 * pi * r_i * U_i_guess * dz) / (m_dot_d_well * cp_d)
        else:
             NTU_guess = 0 # If denominator is zero, NTU is zero, no temperature change

        # Ensure exp argument is not too large or too small
        if NTU_guess > 50: # Arbitrary large value to prevent overflow
             exp_term = 0
        elif NTU_guess < -50: # Arbitrary small value
             exp_term = np.inf # Or handle appropriately, likely indicates an issue
        else:
             exp_term = np.exp(-NTU_guess)

        T1_guess = T_m[i] + (T0 - T_m[i]) * exp_term


        for _ in range(delta_maxiter):
            # Call compute_T_re without KW
            T_re_new = compute_T_re(
                T0, T1_guess, T_m[i], r_i, r_e, k_w, h_i, Rf_i_prime, k_m, delta_guess, L_tube, fin_t, fin_L, num_fins
            )
            try:
                delta_new = compute_delta2_fast(
                    r_e, k_m, cp_m, rho_m, h_m, T_m[i], T_re_new, t,
                    tol=delta_tol, max_iter=50, delta_max=delta_max
                )
            except Exception as e:
                delta_new = 0.0 # Default to zero melt layer on error

            # Call compute_U_i without KW
            U_i_new = compute_U_i(h_i, r_i, r_e, k_w, Rf_i_prime, k_m, delta_new, L_tube, fin_t, fin_L, num_fins)

             # Avoid division by zero if cp_d is zero or near zero
            if m_dot_d_well * cp_d > 1e-9:
                 NTU_new = (2 * pi * r_i * U_i_new * dz) / (m_dot_d_well * cp_d)
            else:
                 NTU_new = 0 # If denominator is zero, NTU is zero, no temperature change

             # Ensure exp argument is not too large or too small
            if NTU_new > 50: # Arbitrary large value to prevent overflow
                 exp_term_new = 0
            elif NTU_new < -50: # Arbitrary small value
                 exp_term_new = np.inf # Or handle appropriately
            else:
                 exp_term_new = np.exp(-NTU_new)


            T1_new = T_m[i] + (T0 - T_m[i]) * exp_term_new


            # Check for convergence of delta and T_re
            if np.abs(delta_new - delta_guess) < delta_tol and np.abs(T_re_new - T_re_guess) < 1e-3:
                break
            delta_guess = delta_new
            T1_guess = T1_new
            T_re_guess = T_re_new


        # After convergence loop, store the results for this segment
        deltas.append(delta_new)
        U_is.append(U_i_new)
        NTUs.append(NTU_new)
        T.append(T1_new)
        Res.append(Re)
        Prs.append(Pr)
        T_res.append(T_re_new)

    # Compute enthalpies and heat transfer rate using the final temperatures
    try:
        h_in = PropsSI("H", "T", T[0], "P", P, fluid)
        h_out = PropsSI("H", "T", T[-1], "P", P, fluid)
        Q = m_dot_d_well * (h_in - h_out)  # [W]
    except ValueError:
        Q = 0.0 # Assign zero heat transfer if enthalpy calculation fails


    # Compute the melt cross-section area at each segment and total melt volume
    # deltas has length n_segments, z has length n_segments + 1
    # A_melt should correspond to the segment centers or ends. Using deltas at segment ends (i=1 to n_segments)
    # or segment average delta. Let's use the delta calculated for each segment.
    deltas_array = np.array(deltas) # This has length n_segments
    # Ensure r_e + deltas_array is non-negative before squaring
    radii_plus_delta_sq = (r_e + np.maximum(0, deltas_array))**2
    A_melt = np.pi * (radii_plus_delta_sq - r_e**2)  # [m²], len=n_segments
    A_melt = np.maximum(0, A_melt) # Ensure melt area is non-negative

    V_melt = np.sum(A_melt * dz)  # [m³], simple Riemann sum using segment length dz


    df = pd.DataFrame({
        "z [m]": z[1:], # Use segment end points for plotting profile
        "T [K]": T[1:],
        "NTU [-]": NTUs,
        "U_i [W/m²·K]": U_is,
        "delta [m]": deltas,
        "T_re [K]": T_res,
        "Re [-]": Res,
        "Pr [-]": Prs,
        "A_melt [m2]": A_melt
    })

    return df, Q, V_melt


########################################
########################################
########################################


def time_profiles_melt(times, geom_par_vector, T_inlet, N_lay, T_m_lay, k_w, Rf_i_prime, m_dot_d_well, P, fluid, k_m, cp_m, rho_m, h_m,
    n_segments, delta_max=1.0, delta_tol=1e-5, delta_maxiter=20
):
    """
    Computes the temperature and melt layer profiles along the well at several times,
    returns the heat transfer rate, melt volume at each time,
    and the cumulative heat transferred (in Joules).

    Parameters:
    -----------
    times : array-like
        List or array of times [s] at which to compute the profiles.
    geom_par_vector: An array or list containing [L_well, L_tube, D_well, r_i, r_e, D_i_tube, D_e_tube, fin_t, fin_L, num_fins, num_tubes].

    Returns:
    --------
    profiles : dict
        Dictionary: {time: DataFrame}, one DataFrame per time step
    Qs : dict
        Dictionary: {time: Q}, heat transfer rate [W] at each time step
    Vmelts : dict
        Dictionary: {time: V_melt}, melt volume [m³] at each time step
    Q_cumulative : float
        Cumulative heat transferred to the melt [J]
    """
    profiles = {}
    Qs = {}
    Vmelts = {}
    for t in times:
        df, Q, V_melt = temperature_profile_melt(geom_par_vector, T_inlet, N_lay, T_m_lay, k_w, Rf_i_prime,m_dot_d_well, P, fluid, t,
            k_m, cp_m, rho_m, h_m, n_segments, delta_max, delta_tol, delta_maxiter
        )
        profiles[t] = df
        Qs[t] = Q
        Vmelts[t] = V_melt

    # Convert Qs and times to arrays for integration
    Q_values = np.array(list(Qs.values()))
    times_arr = np.array(list(Qs.keys()))
    Q_cumulative = np.trapezoid(Q_values, times_arr)/1000.  # [kJ]

    return profiles, Qs, Vmelts, Q_cumulative



# ==========================================================================
# hydraulics
# ==========================================================================


def calculate_pressure_drop(geom_par_vector, m_dot_well_tube, fluid, T_inlet, T_outlet, P_inlet):
    """
    Calculates the pressure drop and pumping power in a single tube of the hairpin heat exchanger.

    Args:
        geom_par_vector: An array or list containing [L_well, L_tube, D_well, r_i, r_e, D_i_tube, D_e_tube, fin_t, fin_L, num_fins, num_tubes].
                         We will use L_tube and D_i_tube (2*r_i) from this vector.
        m_dot_well_tube: Mass flow rate of the fluid through a single tube [kg/s].
        fluid: Name of the working fluid.
        T_inlet: Inlet temperature of the fluid [K].
        T_outlet: Outlet temperature of the fluid [K].
        P_inlet: Inlet pressure of the fluid [Pa].


    Returns:
        A tuple containing:
            - Pressure drop in the tube [Pa].
            - Pumping power [W]. # Added pumping power to return values
    """
    L_well, L_tube, D_well, r_i, r_e, D_i_tube, D_e_tube, fin_t, fin_L, num_fins, num_tubes = geom_par_vector
    L = L_tube # Length of the tube
    D_i = D_i_tube # Inner diameter of the tube
    r_i_tube = D_i / 2.0 # Inner radius

    # Average temperature and pressure for fluid properties
    T_avg = 0.5 * (T_inlet + T_outlet)
    # Assume average pressure is close to inlet pressure for property calculation for now
    # A more rigorous approach might iterate, but this is a good starting point.
    P_avg = P_inlet # Approximation


    try:
        # Get fluid properties at average conditions
        rho = CP.PropsSI("D", "T", T_avg, "P", P_avg, fluid)
        mu = CP.PropsSI("V", "T", T_avg, "P", P_avg, fluid) # Dynamic viscosity
        # You might also need Specific Heat (Cp) and Thermal Conductivity (k) for friction factor correlations,
        # but for basic pressure drop (Darcy-Weisbach), rho and mu are key.

    except ValueError:
        print(f"Warning: Could not retrieve fluid properties for {fluid} at T={T_avg:.2f}K, P={P_avg:.2f}Pa.")
        return float('nan'), float('nan') # Return NaN for both if properties fail


    # Calculate flow velocity
    A_flow = np.pi * r_i_tube**2
    if rho * A_flow > 1e-9:
        u = m_dot_well_tube / (rho * A_flow)
    else:
        return 0.0, 0.0 # No flow, no pressure drop and no pumping power

    # Calculate Reynolds number
    if mu > 1e-9 and D_i > 1e-9:
        Re = rho * u * D_i / mu
    else:
        return 0.0, 0.0

    # Calculate friction factor (Darcy friction factor)
    f = 0.0
    if Re < 2000: # Laminar flow
        f = 64.0 / Re
    else: # Turbulent flow (using Haaland correlation as an example)
        # More accurate correlations like Colebrook-White could be used if roughness is known.
        # For simplicity, using an approximation that doesn't require iteration for now.
        # If pipe roughness (epsilon) is available, a better correlation can be used.
        # Assuming smooth pipe for now (epsilon = 0).
        # The Haaland correlation: 1/sqrt(f) = -1.8 * log10((epsilon/D_i)/3.7 + (6.9/Re))
        # With epsilon=0: 1/sqrt(f) = -1.8 * log10(6.9/Re)
        # f = (1 / (-1.8 * log10(6.9/Re)))**2

        # Using a simpler turbulent correlation like Blasius for smooth pipes (Re < 100,000)
        # or a transition correlation for wider Re range.
        # Let's use a common approximation for smooth pipes for simplicity:
        if Re <= 1e5: # Blasius correlation
             f = 0.3164 / (Re**0.25)
        else: # For higher Re, use a different correlation, e.g., from Moody chart or a more general formula
             # Swamee–Jain equation (explicit approximation of Colebrook–White)
             # Requires pipe roughness (epsilon). Let's assume a typical value for smooth drawn tubing.
             # Example epsilon for drawn tubing might be 0.0015 mm = 1.5e-6 m
             # If epsilon is not a parameter, we have to stick to correlations that don't need it,
             # or make an assumption. Let's assume smooth pipe and use Blasius/extend it.
             # A more robust approach would take epsilon as an input parameter.
             # For now, extending Blasius or using a simple fit:
             # Let's use a piecewise approach based on Re.
             # For very high Re, f approaches a constant for smooth pipes.
             # A simple fit across a wider range might be better.
             # Let's stick to the Blasius for Re <= 100,000 and
             # a rough approximation for higher Re if needed, or add roughness as parameter.

             # For now, let's refine the Blasius range and consider a higher Re approximation if necessary.
             # The Blasius correlation is typically for Re up to 10^5.
             # For 10^5 < Re < 10^6, a common approximation is f = 0.0032 + 0.221 / Re^0.237.
             # For Re >= 10^6, f ~ 0.018 (fully turbulent smooth pipe).
             # Let's add a simple transition:
             if Re <= 1e6:
                  f = 0.0032 + 0.221 / Re**0.237
             else:
                  f = 0.018 # Approximation for very high Re, smooth pipe


    # Frictional pressure drop (Darcy-Weisbach equation)
    # Ensure D_i is non-zero before division
    if D_i > 1e-9:
        delta_P_friction = f * (L / D_i) * (0.5 * rho * u**2)
    else:
         delta_P_friction = 0.0


    # Minor losses (e.g., bends, entrance/exit)
    # These depend on fittings and geometry. Assuming a simple hairpin with two 180-degree bends
    # and entrance/exit losses. Need minor loss coefficients (K).
    # For simplicity, assuming K_bend for 180 deg return bend ~ 2.2 (for flanged or welded)
    # K_entrance ~ 0.5, K_exit ~ 1.0
    # Total minor loss coefficient K_total = K_entrance + 2 * K_bend + K_exit
    # K_total = 0.5 + 2 * 2.2 + 1.0 = 5.9 (example values)
    # Minor loss pressure drop: delta_P_minor = K_total * 0.5 * rho * u**2

    # Let's include minor losses as a parameter or use assumed values.
    # For this function, let's assume typical minor loss coefficients for a hairpin.
    K_entrance = 0.5 # Sharp-edged entrance
    K_bend_180 = 2.2 # 180 degree standard elbow/bend
    K_exit = 1.0 # Exit into a large reservoir/pipe

    K_total = K_entrance + 2 * K_bend_180 + K_exit # For one hairpin (two tubes and two bends)

    delta_P_minor = K_total * 0.5 * rho * u**2

    # Total pressure drop
    delta_P_total = delta_P_friction + delta_P_minor

    # Calculate pumping power [W]
    # Ensure rho is non-zero before division
    if rho > 1e-9:
        pumping_power = delta_P_total * m_dot_well_tube / rho
    else:
        pumping_power = 0.0


    return delta_P_total, pumping_power # Return both pressure drop and pumping power

# @title

# Function to evaluate the energy ratio for a given value of N_wells (Charging)



# ==========================================================================
# solvers
# ==========================================================================


def evaluate_Q_ratio_ch(m_dot_c, N_wells_val, num_tubes, geom_par_vector, times_ch,
                        T_4c_kelvin, N_lay, T_m_lay, k_m_l, Rf_i_prime, P, fluid2, cp_m_l, rho_m_l, h_m,
                        n_segments, delta_max, D_E_out_HP):
    """
    Evaluates Q_ratio_ch for a given m_dot_c_well1.
    """
    m_dot_c_well = m_dot_c / N_wells_val
    m_dot_c_well1 = m_dot_c_well / num_tubes

    profsl, Qliq, Vliq_well, Q_cum_l_well1 = time_profiles_melt(times_ch, geom_par_vector,
        T_4c_kelvin, N_lay, T_m_lay, k_m_l, Rf_i_prime, m_dot_c_well1, P, fluid2, k_m_l, cp_m_l, rho_m_l, h_m,
        n_segments, delta_max)

    Q_cum_l_well = num_tubes * Q_cum_l_well1
    Q_cum_l = N_wells_val * Q_cum_l_well

    if abs(Q_cum_l) > 1e-9:
        return D_E_out_HP / Q_cum_l
    else:
        return float('inf') # Assign a large value if denominator is zero

#########################################
#########################################
#########################################

# Function to evaluate the energy ratio for a given value of m_dot_d_well1 (Discharging)


def evaluate_Q_ratio_dc(N_wells, m_dot_d_well1_val, num_tubes, geom_par_vector, times_dc,
                        T_3d_kelvin, N_lay, T_m_lay_dc, k_w, Rf_i_prime, P, fluid2, k_m_s, cp_m_s, rho_m_s, h_m,
                        n_segments, delta_max, D_E_in_ORC):

    profss, Qsol, Vsol_well, Q_cum_s_well1 = time_profiles_melt(times_dc, geom_par_vector,
        T_3d_kelvin, N_lay, T_m_lay_dc, k_w, Rf_i_prime, m_dot_d_well1_val, P, fluid2, k_m_s, cp_m_s, rho_m_s, h_m,
        n_segments, delta_max)

    Q_cum_s_well = -num_tubes * Q_cum_s_well1
    Q_cum_s = N_wells * Q_cum_s_well

    if abs(Q_cum_s) > 1e-9:
        return D_E_in_ORC / Q_cum_s
    else:
        return float('inf') # Assign a large value if denominator is zero


#########################################
#########################################
#########################################


def find_N_wells_for_Q_ratio_ch_fast(target_Q_ratio, tolerance,N_wells_low, N_wells_high,
    m_dot_c, num_tubes, geom_par_vector, times_ch,
    T_4c_kelvin, N_lay, T_m_lay, k_m_l, Rf_i_prime, P, fluid2, cp_m_l, rho_m_l, h_m,
    n_segments, delta_max, D_E_out_HP,
    N_hint=None,
    max_iterations=50,
    cache_round_ndigits=6,
    verbose=False
):
    """
    Faster bracketed solver for N_wells such that Q_ratio_ch ~ target_Q_ratio.

    Improvements included:
      (1) Memoization of expensive function evaluations
      (2) Illinois regula falsi (bracket-preserving, faster than bisection in practice)
      (3) Smarter bracketing using an optional hint N_hint (e.g., previous optimum in sweeps)

    Parameters
    ----------
    N_hint : float or None
        If provided, attempts to bracket the solution around N_hint within [N_wells_low, N_wells_high].
        Use this in sweeps: pass the previous N_wells solution as N_hint to reduce iterations.

    Notes
    -----
    - Assumes Q_ratio_ch is reasonably monotone with N_wells (common here since m_dot per well ~ 1/N).
    - Returns a float N_wells (you can round/ceil outside if you need an integer).
    """

    # --- Small helpers ---
    def _clamp_positive(x):
        return max(1.0, float(x))

    # Memoization cache (keyed by rounded N)
    cache = {}

    def func(N_wells_val):
        # clamp to physical domain
        N_wells_val = _clamp_positive(N_wells_val)

        key = round(N_wells_val, cache_round_ndigits)
        if key in cache:
            return cache[key]

        # quick reject: if flow per tube is ~0, Q_ratio gets pathological
        m_dot_c_well1_val = (m_dot_c / N_wells_val) / num_tubes
        if m_dot_c_well1_val <= 1e-12:
            val = 1e9
            cache[key] = val
            return val

        try:
            q_ratio = evaluate_Q_ratio_ch(
                m_dot_c, N_wells_val, num_tubes, geom_par_vector, times_ch,
                T_4c_kelvin, N_lay, T_m_lay, k_m_l, Rf_i_prime, P, fluid2, cp_m_l, rho_m_l, h_m,
                n_segments, delta_max, D_E_out_HP
            )
            val = q_ratio - target_Q_ratio
        except Exception:
            val = 1e9

        cache[key] = val
        return val

    # --- Step (3): smarter bracketing around N_hint, if provided ---
    lo = _clamp_positive(N_wells_low)
    hi = max(lo + 1e-12, float(N_wells_high))  # avoid hi==lo

    if N_hint is not None and np.isfinite(N_hint):
        # Try a tight bracket around the hint first, clipped to [lo, hi]
        center = float(N_hint)
        a_try = max(lo, 0.6 * center)
        b_try = min(hi, 1.4 * center)

        # If the bracket collapses (hint near bounds), widen gracefully
        if b_try <= a_try * (1.0 + 1e-12):
            a_try = lo
            b_try = hi

        fa_try = func(a_try)
        fb_try = func(b_try)

        # If it brackets, use it; else fall back to original full bracket
        if np.isfinite(fa_try) and np.isfinite(fb_try) and fa_try * fb_try < 0:
            lo, hi = a_try, b_try

    # Evaluate endpoints
    f_lo = func(lo)
    f_hi = func(hi)

    # If not bracketed, attempt automatic expansion around midpoint (best-effort)
    if not (np.isfinite(f_lo) and np.isfinite(f_hi)) or f_lo * f_hi > 0:
        if verbose:
            print(f"[warn] Target not bracketed initially in [{lo}, {hi}]. Attempting expansion.")

        # Expand in a few steps, but stay positive and within some sane upper bound
        mid = 0.5 * (lo + hi)
        a, b = lo, hi
        fa, fb = f_lo, f_hi

        # Try expand up to 6 times (geometric)
        for _ in range(6):
            # expand outward
            a_new = max(1.0, a / 1.8)
            b_new = b * 1.8
            fa_new = func(a_new)
            fb_new = func(b_new)

            if np.isfinite(fa_new) and np.isfinite(fb_new) and fa_new * fb_new < 0:
                lo, hi = a_new, b_new
                f_lo, f_hi = fa_new, fb_new
                break

            a, b, fa, fb = a_new, b_new, fa_new, fb_new

        # If still not bracketed: return whichever endpoint is closer in residual
        if f_lo * f_hi > 0:
            if verbose:
                print("[warn] Could not bracket root after expansion. Returning closest bound.")
            return lo if abs(f_lo) < abs(f_hi) else hi

    # --- Step (2): Illinois regula falsi loop (faster than bisection, still bracketed) ---
    a, b = float(lo), float(hi)
    fa, fb = float(f_lo), float(f_hi)

    if verbose:
        print(f"Starting Illinois regula falsi on [{a:.6g}, {b:.6g}] with target={target_Q_ratio:.4g}")

    for it in range(1, max_iterations + 1):
        denom = (fb - fa)
        if abs(denom) < 1e-30:
            # function almost flat between endpoints; return midpoint
            c = 0.5 * (a + b)
            return _clamp_positive(c)

        # regula falsi point
        c = (a * fb - b * fa) / denom
        c = _clamp_positive(c)

        fc = float(func(c))

        if verbose:
            print(f"  it={it:02d} a={a:.6g} b={b:.6g} c={c:.6g} fa={fa:.3e} fb={fb:.3e} fc={fc:.3e}")

        # Convergence check on residual
        if abs(fc) < tolerance:
            return c

        # Update bracket with Illinois modification to avoid endpoint stalling
        if fa * fc < 0:
            b, fb = c, fc
            fa *= 0.5
        else:
            a, fa = c, fc
            fb *= 0.5

    # If max iterations reached, return best available (midpoint of final bracket)
    return 0.5 * (a + b)


#########################################
#########################################
#########################################


def find_m_dot_d_well1_for_Q_ratio_dc_fast(
    target_Q_ratio, tolerance,
    ratio_low, ratio_high,
    m_dot_c_well1, N_wells, num_tubes, geom_par_vector, times_dc,
    T_3d_kelvin, N_lay, T_m_lay_dc, k_w, Rf_i_prime, P, fluid2, k_m_s, cp_m_s, rho_m_s, h_m,
    n_segments, delta_max, D_E_in_ORC,
    ratio_hint=None,
    max_iterations=50,
    cache_round_ndigits=6,
    verbose=False
):
    """
    Faster bracketed solver for discharging m_dot_d_well1 using a root-find on ratio:
        ratio = m_dot_d_well1 / m_dot_c_well1
        find ratio such that Q_ratio_dc(ratio) ~ target_Q_ratio

    Improvements included:
      (1) Memoization of expensive function evaluations (ratio -> residual)
      (2) Illinois regula falsi (bracket-preserving, faster than bisection in practice)
      (3) Smarter bracketing using optional ratio_hint (e.g., previous optimum in sweeps)

    Returns
    -------
    m_dot_d_well1_opt : float
        Optimal discharging mass flow rate per tube [kg/s]
    """

    # ---- helpers ----
    def _clamp_ratio(x):
        # ratio must be positive
        return max(1e-12, float(x))

    # Memoization cache: key is rounded ratio
    cache = {}

    def func(ratio_val):
        ratio_val = _clamp_ratio(ratio_val)
        key = round(ratio_val, cache_round_ndigits)
        if key in cache:
            return cache[key]

        m_dot_d_well1_val = ratio_val * m_dot_c_well1

        # quick reject
        if m_dot_d_well1_val <= 1e-12:
            val = 1e9
            cache[key] = val
            return val

        try:
            q_ratio = evaluate_Q_ratio_dc(
                N_wells, m_dot_d_well1_val, num_tubes, geom_par_vector, times_dc,
                T_3d_kelvin, N_lay, T_m_lay_dc, k_w, Rf_i_prime, P, fluid2,
                k_m_s, cp_m_s, rho_m_s, h_m,
                n_segments, delta_max, D_E_in_ORC
            )
            val = q_ratio - target_Q_ratio
        except Exception:
            val = 1e9

        cache[key] = val
        return val

    # ---- Step (3): smarter bracketing around ratio_hint ----
    a = _clamp_ratio(ratio_low)
    b = max(a * (1.0 + 1e-12), float(ratio_high))

    if ratio_hint is not None and np.isfinite(ratio_hint):
        center = float(ratio_hint)
        a_try = max(a, 0.6 * center)
        b_try = min(b, 1.4 * center)

        if b_try <= a_try * (1.0 + 1e-12):
            a_try, b_try = a, b

        fa_try, fb_try = func(a_try), func(b_try)

        if np.isfinite(fa_try) and np.isfinite(fb_try) and fa_try * fb_try < 0:
            a, b = a_try, b_try

    fa = func(a)
    fb = func(b)

    # If not bracketed, best-effort expansion
    if not (np.isfinite(fa) and np.isfinite(fb)) or fa * fb > 0:
        if verbose:
            print(f"[warn] Target not bracketed initially in ratio [{a}, {b}]. Attempting expansion.")

        a0, b0 = a, b
        fa0, fb0 = fa, fb

        for _ in range(6):
            a_new = max(1e-12, a0 / 1.8)
            b_new = b0 * 1.8
            fa_new = func(a_new)
            fb_new = func(b_new)
            if np.isfinite(fa_new) and np.isfinite(fb_new) and fa_new * fb_new < 0:
                a, b, fa, fb = a_new, b_new, fa_new, fb_new
                break
            a0, b0, fa0, fb0 = a_new, b_new, fa_new, fb_new

        # Still not bracketed → return closer bound
        if fa * fb > 0:
            if verbose:
                print("[warn] Could not bracket root after expansion. Returning closest bound.")
            ratio_best = a if abs(fa) < abs(fb) else b
            return ratio_best * m_dot_c_well1

    # ---- Step (2): Illinois regula falsi loop ----
    if verbose:
        print(f"Starting Illinois regula falsi on ratio [{a:.6g}, {b:.6g}] target={target_Q_ratio:.4g}")

    for it in range(1, max_iterations + 1):
        denom = (fb - fa)
        if abs(denom) < 1e-30:
            ratio_mid = 0.5 * (a + b)
            return ratio_mid * m_dot_c_well1

        c = (a * fb - b * fa) / denom
        c = _clamp_ratio(c)

        fc = float(func(c))

        if verbose:
            print(f"  it={it:02d} a={a:.6g} b={b:.6g} c={c:.6g} fa={fa:.3e} fb={fb:.3e} fc={fc:.3e}")

        if abs(fc) < tolerance:
            return c * m_dot_c_well1

        if fa * fc < 0:
            b, fb = c, fc
            fa *= 0.5  # Illinois damping
        else:
            a, fa = c, fc
            fb *= 0.5

    # fallback: midpoint
    return 0.5 * (a + b) * m_dot_c_well1

#############################################
#############################################
#############################################
