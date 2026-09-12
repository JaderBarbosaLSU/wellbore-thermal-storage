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

  # Convergence loop for T_7h
  T_7h = T_3h - 4.  # initial guess for T_7h
  tol = 1e-4  # Tolerance for convergence
  max_iter = 100 # Maximum iterations
  p_7h = p_condh # State 7h is at condenser pressure

  for _ in range(max_iter):
      try:
          # Calculate h_7h based on current T_7h and p_7h
          h_7h_calc = CP.PropsSI('H', 'T', T_7h, 'P', p_7h, refrig)/1000.

          # State 8h: Flash Gas two-phase mixture (p_inth, h_8h=h_7h_calc)
          p_8h = p_inth
          h_8h_calc = h_7h_calc # Isenthalpic expansion from 7h
          # Ensure denominator is non-zero before calculating x_8h
          if abs(h_10h - h_9h) > 1e-9:
               x_8h_calc = (h_8h_calc - h_9h) / (h_10h - h_9h)
          else:
               x_8h_calc = 0.0 # If denominator is zero, assume quality is zero

          # Ensure x_8h_calc is within [0, 1] bounds
          x_8h_calc = max(0.0, min(1.0, x_8h_calc))

          # Calculate T_12h based on IHX-1 effectiveness
          T_12h_calc = T_3h - epsilon_IHX_1 * x_8h_calc * cpr_1 * (T_3h - T_10h)

          # IHX-2 effectiveness (based on T_1h and T_13h, and T_12h and T_13h)
          # Ensure denominator is non-zero
          denominator_epsilon2 = (T_12h_calc - T_13h)
          if abs(denominator_epsilon2) > 1e-9:
               epsilon_IHX_2_calc = (T_1h - T_13h) / denominator_epsilon2
          else:
               epsilon_IHX_2_calc = 1.0 # Assume 100% effectiveness if temperature difference is zero


          # CP_ratio_2:
          try:
              cpf_2 = CP.PropsSI('C', 'P', p_condh, 'Q', 0, refrig)
              cpv_2 = CP.PropsSI('C', 'P', p_evaph, 'Q', 1, refrig)
              # Avoid division by zero
              if abs(cpf_2) > 1e-9:
                   cpr_2_calc = cpv_2 / cpf_2
              else:
                   cpr_2_calc = 1.0 # Assume ratio is 1 if cpf is zero
          except ValueError:
              cpr_2_calc = 1.0

          # State 4h: two-phase mixture (p_evaph, h_4h=h_9h)
          p_4h = p_evaph
          h_4h_calc = h_9h # Isenthalpic expansion from 9h
          # Ensure denominator is non-zero before calculating x_4h
          if abs(h_13h - h_14h) > 1e-9:
              x_4h_calc = (h_4h_calc - h_14h) / (h_13h - h_14h)
          else:
               x_4h_calc = 0.0 # If denominator is zero, assume quality is zero

          # Ensure x_4h_calc is within [0, 1] bounds
          x_4h_calc = max(0.0, min(1.0, x_4h_calc))

          # Calculate T_7h_new based on IHX-2 effectiveness
          # Ensure denominator is non-zero before calculation
          if abs(cpr_2_calc * (T_12h_calc - T_13h)) > 1e-9:
               T_7h_new = T_12h_calc - epsilon_IHX_2_calc * x_4h_calc * cpr_2_calc * (T_12h_calc - T_13h)
          else:
               T_7h_new = T_7h # No change if denominator is zero

          # Check for convergence of T_7h
          if abs(T_7h_new - T_7h) < tol:
              T_7h = T_7h_new # Update T_7h to the converged value
              # After convergence, update all dependent state properties with the converged T_7h
              h_7h = h_7h_calc
              h_8h = h_8h_calc
              x_8h = x_8h_calc
              T_12h = T_12h_calc
              epsilon_IHX_2 = epsilon_IHX_2_calc
              cpr_2 = cpr_2_calc
              h_4h = h_4h_calc
              x_4h = x_4h_calc
              break # Exit the loop if converged

          T_7h = T_7h_new # Update guess for the next iteration

      except ValueError as e:
          print(f"Warning: CoolProp calculation failed in T_7h convergence loop: {e}. Breaking loop.")
          # Assign NaN to dependent properties on error and break
          h_7h, x_8h, T_12h, epsilon_IHX_2, cpr_2, h_4h, x_4h = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
          break # Exit loop on CoolProp error

  else:
      # After max iterations without convergence, assign the values from the last iteration
      # Need to recalculate dependent variables one last time based on the final T_7h
      try:
          h_7h = CP.PropsSI('H', 'T', T_7h, 'P', p_7h, refrig)/1000.
          if abs(h_10h - h_9h) > 1e-9:
               x_8h = (h_7h - h_9h) / (h_10h - h_9h)
          else:
               x_8h = 0.0
          x_8h = max(0.0, min(1.0, x_8h))
          T_12h = T_3h - epsilon_IHX_1 * x_8h * cpr_1 * (T_3h - T_10h)
          denominator_epsilon2 = (T_12h - T_13h)
          if abs(denominator_epsilon2) > 1e-9:
               epsilon_IHX_2 = (T_1h - T_13h) / denominator_epsilon2
          else:
               epsilon_IHX_2 = 1.0
          try:
               cpf_2 = CP.PropsSI('C', 'P', p_condh, 'Q', 0, refrig)
               cpv_2 = CP.PropsSI('C', 'P', p_evaph, 'Q', 1, refrig)
               if abs(cpf_2) > 1e-9:
                    cpr_2 = cpv_2 / cpf_2
               else:
                    cpr_2 = 1.0
          except ValueError:
               cpr_2 = 1.0
          h_4h = h_9h
          if abs(h_10h - h_9h) > 1e-9:
              x_4h = (h_4h - h_9h) / (h_10h - h_9h)
          else:
               x_4h = 0.0
          x_4h = max(0.0, min(1.0, x_4h))

      except ValueError as e:
          print(f"Warning: Final CoolProp calculation failed after max iterations: {e}.")
          h_7h, x_8h, T_12h, epsilon_IHX_2, cpr_2, h_4h, x_4h = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


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


#####################################
#####################################
#####################################



# ==========================================================================
# march_legacy
# ==========================================================================


########################################
########################################
########################################


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


#########################################
#########################################
#########################################

# Function to evaluate the energy ratio for a given value of m_dot_d_well1 (Discharging)


#########################################
#########################################
#########################################


#########################################
#########################################
#########################################


#############################################
#############################################
#############################################
