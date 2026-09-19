"""THUMS -- latent heat storage in a repurposed wellbore.

The live model, v0.8. One formulation (enthalpy, "Formulation C"), one sizing
framing (specify the hardware and march to cyclic steady state), no root
finding anywhere.

    plant      cycle_state_points   ->  state points and the two efficiencies
               energy_budget        ->  how much energy must move per cycle

    borehole   pcm_capacities       ->  E'_lat, C_s, C_l  per unit tube length
               pcm_state            ->  E'  ->  (A_melt, T_pcm)
               advance_segment      ->  branchwise-exact step of E'
               march_h              ->  one half-cycle over the whole exchanger

    cycle      simulate_css         ->  repeat until the state returns to itself
               simulate_css_corrected  ->  ... with one ORC correction pass
               css_report           ->  print it

Everything is driven by one frozen dataclass, `Case`. To study a parameter,
copy `CASE` with the field changed -- `CASE.with_(N_lay=12)` -- and re-run.
`validate_case` checks a case for the inconsistencies that are easy to create
and hard to see.

Superseded material is not here. The first-cycle sizing chain, the closed-form
Stefan front and the area-based melt front were removed; the reasoning that
removed them is in DESIGN_NOTES.md and in the project report, and the code is
in version control.
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


# ==========================================================================
# config
# ==========================================================================
# One flat Case object. v0.1 passed sixteen positional arguments, which is how
# the PCM conductivity ended up in the slot meant for the tube wall for months
# without anything contradicting it. Named fields make that mistake visible.

from dataclasses import dataclass, replace

FT = 0.3048
IN = 0.0254


@dataclass(frozen=True)
class Case:
    # ---- PCM (Saher 2024 trimodal) ----
    T_m_C: float = 150.0          # melting temperature            [C]
    rho_s: float = 1550.0         # solid density                  [kg/m3]
    rho_l: float = 1450.0         # liquid density                 [kg/m3]
    cp_s: float = 1280.0          # solid specific heat            [J/kg/K]
    cp_l: float = 1800.0          # liquid specific heat           [J/kg/K]
    k_s: float = 0.60             # solid conductivity             [W/m/K]
    k_l: float = 0.45             # liquid conductivity            [W/m/K]
    h_m: float = 380_000.0        # latent heat of fusion          [J/kg]
    cost_per_kWh: float = 20.0

    # ---- borehole and tube geometry ----
    L_ft: float = 5000.0          # well depth                     [ft]
    D_in: float = 7.0             # borehole diameter              [in]
    r_i: float = 0.01623          # tube inner radius              [m]
    r_e: float = 0.02108          # tube outer radius              [m]
    fin_t: float = 0.0015         # fin thickness                  [m]
    fin_L: float = 0.0075         # fin radial length              [m]
    num_fins: int = 24            # fins per leg
    num_tubes: int = 2            # HAIRPINS per borehole (2 legs each)
    k_wall: float = 45.0          # steel conductivity             [W/m/K]
    Rf_i: float = 1e-3            # internal fouling resistance    [m2K/W]

    # ---- power block and operating conditions ----
    fluid: str = "cyclopentane"   # ORC working fluid
    fluid2: str = "Water"         # borehole secondary fluid
    refrig: str = "cyclopentane"  # heat-pump refrigerant
    P: float = 1e6                # secondary-fluid pressure       [Pa]
    W_dot_el_out: float = 1000.0  # ORC net electrical output      [kW]
    Turb_eff: float = 0.85
    ElG_eff: float = 0.95
    Comp_eff: float = 0.85
    ElH_eff: float = 0.95
    T_sink_C: float = 20.0
    DT_E_sink: float = 5.0
    DT_2D_3E: float = 12.0
    # Discharge-inlet subcooling below the COLDEST cascade layer, state 1d in
    # the plant diagram; the symmetric partner of DT_4C_M. See DN-8.
    DT_M_1D: float = 10.0
    T_source_C: float = 60.0
    DT_4C_M: float = 10.0
    DT_3A_4A: float = 10.0
    DT_3A_13H: float = 10.0
    DT_2H_3C: float = 10.0
    DT_sub: float = 2.0
    DT_3C_2C: float = 55.0        # secondary-fluid glide          [K]
    t_ch: float = 10.0            # charging duration              [h]
    t_dc: float = 10.0            # discharging duration           [h]
    loss_surplus: float = 0.05    # assumed storage loss, lambda
    N_lay: int = 9                # PCM layers along the well

    # ---- numerics ----
    n_segments: int = 100
    n_times: int = 40
    delta_max: float = 0.5
    delta_tol: float = 1e-5
    delta_maxiter: int = 20
    N_wells_bracket: tuple = (3.0, 120.0)
    ratio_bracket: tuple = (0.2, 2.0)
    tol_Q_ratio: float = 1e-3
    max_iterations: int = 40

    # ---- switches ----
    # The ONLY remaining switch. True carries PCM enthalpy (latent + sensible in
    # both phases), so a fully melted segment superheats and a fully solid one
    # subcools. False pins the PCM at T_m, which is retained not as a model but
    # as the control case for the self-levelling comparison: it isolates what
    # the sensible branches do, at a fixed well count.
    sensible_heat: bool = True

    # Which material lies between the tube and the phase front, and therefore
    # which shell the heat has to cross. See `conduction_shell` and DN-13.
    #   'directional'  melting -> liquid shell (k_l); freezing -> frozen
    #                  shell (k_s). The physical arrangement in both directions.
    #   'melt_side'    the melted thickness always, with k_l whenever any melt
    #                  exists. Correct while melting, INVERTED while freezing.
    #                  Retained only to reproduce results published before the
    #                  correction; it is not a model.
    front_geometry: str = "directional"

    # ---- derived ----
    @property
    def T_m(self):
        return self.T_m_C + 273.15

    @property
    def L_well(self):
        return self.L_ft * FT

    @property
    def L_tube(self):
        """Developed length of ONE hairpin: down and back."""
        return 2.0 * self.L_well

    @property
    def D_well(self):
        return self.D_in * IN

    @property
    def V_borehole(self):
        return np.pi * self.D_well ** 2 / 4.0 * self.L_well

    @property
    def V_well(self):
        """PCM volume in one borehole: the hole, less tubes and fins."""
        n_legs = 2 * self.num_tubes
        A_tube = np.pi * self.r_e ** 2
        A_fins = self.num_fins * self.fin_t * self.fin_L
        return self.V_borehole - n_legs * (A_tube + A_fins) * self.L_well

    @property
    def rho_latent(self):
        """Density for the latent inventory, the SAME in both directions.

        The front is carried as a solid-equivalent area, so V_melt/V_well is a
        melted MASS fraction and the energy banked on charging equals the energy
        available on discharging. Using rho_l to melt and rho_s to freeze -- the
        directional choice that is correct for the conductivity k_m -- would
        return rho_s/rho_l = 1.069 times more energy than went in.
        """
        return self.rho_s

    @property
    def r_cell(self):
        """Radius of the equivalent PCM cell owned by ONE tube leg.

        The 2*num_tubes legs share the borehole cross-section equally, so each
        owns pi R^2 / (2 n_t) of it; writing that as pi r_cell^2 gives

            r_cell = R / sqrt(2 n_t),      R = D_well / 2.

        This is the radius at which a leg's melt front meets its neighbour's.
        It is NOT the borehole wall: with four legs in a 7 in hole r_cell is
        44.45 mm while the wall is at 88.90 mm, so a limit imposed at the wall
        allows delta = 67.82 mm where the real limit is 23.37 mm -- a factor
        2.9 too permissive, and therefore silent at every design point.
        """
        return (self.D_well / 2.0) / np.sqrt(2.0 * self.num_tubes)

    @property
    def delta_merge(self):
        """Melt-layer thickness at which neighbouring fronts touch.

        Hypothesis H3 (concentric annulus, no interaction between legs) is a
        statement about delta < delta_merge. At and beyond this thickness the
        remaining solid sits in the corners between legs, reached through a
        constricted path that the annular resistance ln(1+delta/r_e)/(2 pi k)
        does not represent.

        Under Formulation C this is not a bound that can be VIOLATED --
        area_from_delta(delta_merge) equals A_avail identically, and the
        enthalpy cap holds A_melt <= A_avail -- so delta <= delta_merge always.
        It is a bound the solution can SIT ON, which is what it does here.
        """
        return self.r_cell - self.r_e

    def geom_vector(self):
        """The eleven positional geometry arguments the geometry helpers expect."""
        return [self.L_well, self.L_tube, self.D_well, self.r_i, self.r_e,
                2 * self.r_i, 2 * self.r_e, self.fin_t, self.fin_L,
                self.num_fins, self.num_tubes]

    def stefan_number(self, dT):
        """Ste = cp dT / h_m. The quasi-steady melt layer degrades as this grows."""
        return self.cp_l * dT / self.h_m

    def with_(self, **kw):
        return replace(self, **kw)


CASE = Case()          # the design point


# ==========================================================================
# front_energy_balance   (the v0.2 correction)
# ==========================================================================


def delta_from_area(A_melt, r_e, n_f=0, t_f=0.0, L_f=0.0):
    """Invert `area_from_delta`. Closed form; the fins only shift the quadratic.

    For delta <= L_f the fin metal grows with the front:

        pi d^2 + (2 pi r_e - n_f t_f) d - A = 0

    and for delta > L_f it is constant at n_f t_f L_f:

        pi d^2 + 2 pi r_e d - (A + n_f t_f L_f) = 0

    With n_f = 0 both collapse to  d = sqrt(r_e^2 + A/pi) - r_e,  the bare
    concentric annulus used before this correction.
    """
    A = np.maximum(np.asarray(A_melt, dtype=float), 0.0)
    # branch 1: front still inside the fins
    b = 2.0 * np.pi * r_e - n_f * t_f
    d_in = (-b + np.sqrt(b * b + 4.0 * np.pi * A)) / (2.0 * np.pi)
    # branch 2: front beyond the fin tips
    d_out = np.sqrt(r_e ** 2 + (A + n_f * t_f * L_f) / np.pi) - r_e
    return np.where(d_in <= L_f, d_in, d_out)


# ==========================================================================
# march   (segment-by-segment integration down the tube)
# ==========================================================================

_PROP_CACHE = {}


def fluid_props(fluid, T, P):
    """Density, viscosity, conductivity, cp. Cached on 0.05 K bins."""
    key = (fluid, round(T * 20.0), round(P))
    hit = _PROP_CACHE.get(key)
    if hit is not None:
        return hit
    rho = CP.PropsSI("D", "T", T, "P", P, fluid)
    mu = CP.PropsSI("V", "T", T, "P", P, fluid)
    k = CP.PropsSI("L", "T", T, "P", P, fluid)
    cp = CP.PropsSI("C", "T", T, "P", P, fluid)
    _PROP_CACHE[key] = (rho, mu, k, cp)
    return rho, mu, k, cp


def h_internal(fluid, T, P, r_i, m_dot):
    """Internal convective coefficient. Gnielinski if turbulent, else laminar.

        Re = 4 m_dot / (pi D mu)
        f  = (0.79 ln Re - 1.64)^-2
        Nu = (f/8)(Re-1000) Pr / [1 + 12.7 sqrt(f/8) (Pr^(2/3) - 1)]
    """
    rho, mu, k, cp = fluid_props(fluid, T, P)
    D = 2 * r_i
    Re = 4.0 * m_dot / (np.pi * D * mu)
    Pr = cp * mu / k
    if Re < 2300:
        Nu = 3.66                      # fully developed, isothermal wall
    else:
        f = (0.79 * np.log(Re) - 1.64) ** -2
        Nu = (f / 8.0) * (Re - 1000) * Pr / (
            1 + 12.7 * np.sqrt(f / 8.0) * (Pr ** (2.0 / 3.0) - 1))
    return Nu * k / D, cp, Re, Pr


def layer_map(T_m_lay, n_segments, N_lay):
    """Melting temperature per segment, on the v0.1 node convention."""
    z = np.linspace(0.0, 1.0, n_segments + 1)
    z_lay = np.linspace(0.0, 1.0, N_lay + 1)
    T_m = np.empty(n_segments + 1)
    for j in range(N_lay):
        T_m[(z >= z_lay[j]) & (z < z_lay[j + 1])] = T_m_lay[j]
    T_m[-1] = T_m_lay[-1]
    return T_m[1:]


# ==========================================================================
# enthalpy state   (v0.4: sensible heat in both phases)
# ==========================================================================
# Formulation B tracks A_melt and nothing else, so a segment that runs out of
# PCM has nowhere to put further heat: the march clips it and discards the
# remainder. At the design point the discarded heat on discharge was 148% of
# the heat delivered -- i.e. most of what the network asked for.
#
# The fix is to carry ENTHALPY per unit tube length instead of melted area, with
# the melted area recovered from it. One scalar per segment, measured from the
# fully-solid-at-T_m datum:
#
#     E' < 0            solid, SUBCOOLED     T = T_m + E'/C_s
#     0 <= E' <= E'_lat two-phase            T = T_m,  A_melt = E'/(rho h_m)
#     E' > E'_lat       liquid, SUPERHEATED  T = T_m + (E'-E'_lat)/C_l
#
# Three things follow, and the third is the reason to do it:
#
#   1. no clipping anywhere -- heat always has somewhere to go;
#   2. the melt fraction is bounded in [0,1] BY CONSTRUCTION, so the
#      eps_local > 1 that the latent-only model produced cannot occur;
#   3. desuperheating and subcooling during discharge are recovered, which is
#      the energy the latent-only model threw away.
#
# delta still follows from A_melt exactly as before, so the heat-transfer
# coefficient is unchanged in form -- only its driving temperature is now
# T_pcm rather than T_m.


def pcm_capacities(case):
    """Per unit tube length: latent capacity and the two sensible capacities."""
    A_avail = case.V_well / (case.num_tubes * case.L_tube)   # PCM area per tube [m2]
    E_lat = case.rho_latent * case.h_m * A_avail             # J/m to melt it all
    C_s = case.rho_latent * case.cp_s * A_avail              # J/m/K, solid
    C_l = case.rho_latent * case.cp_l * A_avail              # J/m/K, liquid
    return A_avail, E_lat, C_s, C_l


def pcm_state(E, T_m_seg, case):
    """Enthalpy per unit tube length -> (melted area, PCM temperature).

    The inverse of the enthalpy curve above. `A_melt` saturates at the segment's
    own PCM content, which is what makes eps_local <= 1 structural.
    """
    A_avail, E_lat, C_s, C_l = pcm_capacities(case)
    T_m_seg = np.asarray(T_m_seg, dtype=float)
    A = np.clip(E, 0.0, E_lat) / (case.rho_latent * case.h_m)
    if not case.sensible_heat:                    # v0.3 behaviour: PCM pinned at T_m
        return A, T_m_seg.copy()
    T = np.where(E < 0.0, T_m_seg + E / C_s,
                 np.where(E > E_lat, T_m_seg + (E - E_lat) / C_l, T_m_seg))
    return A, T


def advance_segment(E, T_m, T0, K, dt, E_lat, C_s, C_l, sensible):
    """Integrate ONE segment's enthalpy over dt, exactly, branch by branch.

    The segment obeys   dE'/dt = K (T0 - T_pcm(E'))   with T0 held over the step.
    That is linear inside each branch of the enthalpy curve, so it can be solved
    in closed form rather than stepped:

      plateau  T_pcm = T_m constant           -> dE'/dt constant, E' linear in t
      solid    T_pcm = T_m + E'/C_s           -> exponential relaxation to
      liquid   T_pcm = T_m + (E'-E_lat)/C_l      E'_eq = C (T0 - T_m) [+E_lat]

    Explicit Euler is NOT usable here. On the latent plateau the segment has
    effectively infinite heat capacity, so any step is stable; the moment
    sensible heat is added the capacity becomes finite and the stability limit
    is dt < C/K, about 620 s at the design point. The time grid's last step is
    8491 s, fourteen times that, and an explicit update diverges -- it drove the
    secondary fluid to 261 K on the first attempt.

    Returns (E_new, q_avg) with q_avg = (E_new - E)/dt, so the heat taken from
    the fluid is exactly the enthalpy the PCM gained: closure stays structural.
    """
    if K <= 0.0 or dt <= 0.0:
        return E, 0.0
    E0, t_left = E, dt
    for _ in range(6):   # at most a few branch crossings
        if t_left <= 0.0:
            break
        if not sensible:                     # v0.3: PCM pinned at T_m
            E = E + K * (T0 - T_m) * t_left
            break
        if 0.0 <= E <= E_lat:                # ---- latent plateau: linear ----
            rate = K * (T0 - T_m)
            if rate == 0.0:
                break
            t_edge = ((E_lat - E) / rate) if rate > 0 else (E / -rate)
            if t_edge >= t_left:
                E = E + rate * t_left
                break
            # Land on the boundary and step just PAST it. Without the nudge the
            # next iteration re-enters the plateau at zero rate and the segment
            # sticks at E_lat for ever -- which showed up as every segment
            # reporting a melt fraction of exactly 1.000 and no superheat.
            eps = 1e-9 * max(E_lat, 1.0)
            E = (E_lat + eps) if rate > 0 else -eps
            t_left -= t_edge
        else:                                # ---- sensible: exponential ----
            if E < 0.0:
                C, E_ref, lo, hi = C_s, 0.0, -np.inf, 0.0
            else:
                C, E_ref, lo, hi = C_l, E_lat, E_lat, np.inf
            E_eq = E_ref + C * (T0 - T_m)
            E_new = E_eq + (E - E_eq) * np.exp(-K * t_left / C)
            if lo <= E_new <= hi:
                E = E_new
                break
            edge = hi if E_new > hi else lo   # crosses back onto the plateau
            num, den = E_eq - edge, E_eq - E
            if den == 0.0 or num / den <= 0.0:
                E = edge
                break
            t_edge = -C / K * np.log(num / den)
            if t_edge >= t_left or t_edge <= 0.0:
                E = E_new
                break
            E = edge
            t_left -= t_edge
    return E, (E - E0) / dt


def segment_profile(case, T_inlet, T_m_seg, m_dot, k_wall, E, n_segments=None):
    """Instantaneous fluid and conductance profile for a frozen PCM state.

    The inner loop of `march_h` without the time advance: it evaluates the
    resistance network, the segment conductance and the fluid temperature
    profile for the state `E` as it stands, and returns the instantaneous
    q' = K (T_f - T_pcm).

    It exists so that the last recorded frame of a march can carry a fluid
    profile belonging to the *end* state rather than to the state one step
    earlier. Everything it computes is already computed inside `march_h`;
    this is one extra sweep per recorded march, not one per time level.
    """
    n_segments = n_segments or len(E)
    dz = case.L_tube / n_segments
    A, T_pcm = pcm_state(E, T_m_seg, case)
    delta = delta_from_area(A, case.r_e, case.num_fins, case.fin_t, case.fin_L)
    A_avail, E_lat = pcm_capacities(case)[:2]
    d_melt, k_melt = conduction_shell(case, A, A_avail, True, E, E_lat)
    d_froz, k_froz = conduction_shell(case, A, A_avail, False, E, E_lat)

    T_prof = np.empty(n_segments + 1)
    T_prof[0] = T_inlet
    NTU_prof = np.empty(n_segments)
    U_prof = np.empty(n_segments)
    q_prime = np.empty(n_segments)

    T0 = T_inlet
    for i in range(n_segments):
        melting = T0 > T_pcm[i]
        d_path = float(d_melt[i] if melting else d_froz[i])
        k_m = float(k_melt[i] if melting else k_froz[i])
        h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, case.r_i, m_dot)
        U_i = compute_U_i(h_i, case.r_i, case.r_e, k_wall, case.Rf_i, k_m,
                          d_path, case.L_tube, case.fin_t,
                          case.fin_L, case.num_fins)
        NTU = float(np.clip((2 * np.pi * case.r_i * U_i * dz) / (m_dot * cp_d),
                            -50.0, 50.0))
        K = m_dot * cp_d * (1.0 - np.exp(-NTU)) / dz
        q_prime[i] = K * (T0 - T_pcm[i])
        T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
        T_prof[i + 1] = T0
        NTU_prof[i] = NTU
        U_prof[i] = U_i

    return dict(T_fluid=T_prof, NTU=NTU_prof, U_i=U_prof, q_prime=q_prime,
                A_melt=A, T_pcm=T_pcm, delta=delta, E=np.array(E, float))


def bulk_shape_factor(case):
    """Dimensionless resistance of a SINGLE-PHASE cell, mean to tube surface.

    A cell that is entirely one phase has no front. Heat entering at the tube is
    stored throughout its volume, so the temperature varies with radius and the
    lumped state T_pcm is the VOLUME MEAN. The resistance that belongs with that
    mean is not the surface-to-surface annulus value.

    Take the cell as an annulus r_e <= r <= r_o, adiabatic at r_o (its neighbour
    is identical, so the boundary is a symmetry plane), storing sensible heat at
    a uniform volumetric rate s = rho c dT/dt. Quasi-steady, the heat crossing
    radius r is what the material beyond it stores,

        -k 2 pi r dT/dr = s pi (r_o^2 - r^2),

    which integrates to

        T(r_e) - T(r) = (s/2k) [ r_o^2 ln(r/r_e) - (r^2 - r_e^2)/2 ].

    Averaging over the volume and dividing by the total Q' = s pi (r_o^2-r_e^2)
    gives R'_bulk = S/(2 pi k) with, writing beta = r_o/r_e,

        S = [ beta^4 ln(beta) - beta^4/2 + beta^2/2 - (beta^2-1)^2/4 ]
            / (beta^2 - 1)^2 .

    S depends on GEOMETRY ALONE -- k cancels -- so one shape factor serves both
    phases and only the conductivity changes between them. For the default cell,
    beta = 2.1086 and S = 0.34672, against ln(beta) = 0.74604 for the
    surface-to-surface annulus: the mean-to-surface resistance is smaller by a
    factor 2.152, the cylindrical analogue of the familiar 1/3 for a slab with
    uniform generation.

    Verified against numerical quadrature of the same integral to six decimals.
    """
    beta = case.r_cell / case.r_e
    b2 = beta * beta
    return ((b2 * b2 * np.log(beta) - b2 * b2 / 2.0 + b2 / 2.0
             - (b2 - 1.0) ** 2 / 4.0) / (b2 - 1.0) ** 2)


def bulk_equivalent_delta(case):
    """The shell thickness that reproduces the bulk resistance.

    The network in `compute_U_i` takes a thickness, not a resistance, and wraps
    it in the fin efficiency and the finned-area ratio. Feeding it an equivalent
    thickness therefore keeps the single-phase branches inside exactly the same
    machinery as the two-phase ones, rather than bolting a second path onto it:

        ln(1 + delta_eq/r_e) = S    =>    delta_eq = r_e (e^S - 1).

    Like S it is a property of the geometry alone -- the same delta_eq serves
    subcooled solid and superheated liquid, and only k changes. For the default
    cell it is 8.74 mm, against 0 and 23.37 mm, which are the two values the
    model used before this was worked out.
    """
    return case.r_e * (np.exp(bulk_shape_factor(case)) - 1.0)


def conduction_shell(case, A_melt, A_avail, melting, E=None, E_lat=None):
    """Thickness and conductivity of the shell between the tube and the front.

    The PCM-side resistance in `compute_U_i` is that of an annulus growing
    OUTWARD FROM THE TUBE WALL,

        R' = ln(1 + delta/r_e) / (2 pi k),

    so the question this function answers is: which material is that annulus
    made of, and how thick is it?

    It depends on the direction of the phase change, because the front always
    grows away from the tube -- the tube is the driven boundary in both
    half-cycles.

      MELTING (fluid hotter than the PCM). Melting begins at the tube wall and
      the front moves outward, so the shell is LIQUID and its thickness is that
      of the melted area. Heat crosses the melt to reach the remaining solid.

      FREEZING (fluid colder). Solidification also begins at the tube wall, so
      the shell is the FROZEN material and its thickness is that of the solid
      area, A_avail - A_melt. The liquid is displaced outward, beyond the front,
      and is no longer in the conduction path at all.

    All four limits come out right without special-casing:

      fully solid, melting     A_melt = 0            -> delta = 0, no resistance
      fully melted, melting    A_melt = A_avail      -> the full liquid annulus
      fully melted, freezing   A_avail - A_melt = 0  -> delta = 0, no resistance
      fully solid, freezing    A_avail - A_melt = A  -> the full frozen annulus

    The last of those is the one that matters. Until v0.7 the melted thickness
    was used in BOTH directions, which put the conduction path on the wrong
    side of the front during discharge: a segment that had frozen solid was
    given ZERO PCM-side resistance, exactly where the physical resistance is
    largest. The consequence was that the modelled exchanger IMPROVED as it
    froze, when it should degrade. See DN-13 for what that was worth.

    LIMIT OF ANY LUMPED FRONT. One thickness can describe one front. After the
    first half-cycle from a fully solid store the geometry is generally
    three-region -- at cyclic steady state the charge begins with residual
    liquid left in the OUTER part of the cell, so melting produces liquid at the
    tube, solid in the middle and liquid outside. Neither this function nor the
    retired one can represent that; this one is the better approximation, not a
    correct treatment. It is hypothesis H5, stated honestly.
    """
    if case.front_geometry not in ("directional", "melt_side"):
        raise ValueError(f"front_geometry must be 'directional' or "
                         f"'melt_side', not {case.front_geometry!r}")

    A_melt = np.asarray(A_melt, float)
    if case.front_geometry == "melt_side":          # retired, see above
        return (delta_from_area(A_melt, case.r_e, case.num_fins, case.fin_t,
                                case.fin_L),
                np.full(A_melt.shape, case.k_l))

    # ---- two-phase: the shell between the tube and the front --------------
    A_path = A_melt if melting else (A_avail - A_melt)
    delta = delta_from_area(A_path, case.r_e, case.num_fins, case.fin_t,
                            case.fin_L)
    k_path = np.full(A_melt.shape, case.k_l if melting else case.k_s)

    # ---- single phase: there is no front, so there is no side ------------
    # `E` is optional only so that callers which do not have it (none, now)
    # keep the two-phase behaviour; when it is given, the single-phase cells
    # get the bulk resistance of `bulk_shape_factor`, the SAME in both
    # directions, because a one-phase cell has no front and cannot care which
    # way the heat is going.
    if E is not None:
        E = np.asarray(E, float)
        d_bulk = bulk_equivalent_delta(case)
        solid = E < 0.0                              # subcooled solid
        liquid = E > E_lat                           # superheated liquid
        delta = np.where(solid | liquid, d_bulk, delta)
        k_path = np.where(solid, case.k_s,
                          np.where(liquid, case.k_l, k_path))
    return delta, k_path


def mixed_mean_outlet(res):
    """Flow-mixed mean outlet temperature of a march [K].

    The temperature at which the stream returns to the plant, averaged over the
    half-cycle. At constant mass flow and c_p the mixing average reduces to the
    time average, which on a logarithmic grid still has to be weighted by the
    step -- the plain mean of `T_out` is dominated by the early levels and is
    wrong by several kelvin.

    This is the quantity the ORC actually sees, and since v0.6 it is a RESULT
    rather than a prescribed state point (`cycle_state_points`).
    """
    t = np.asarray(res["t"], float)
    dt = np.diff(np.concatenate(([0.0], t)))
    return float(np.sum(np.asarray(res["T_out"], float) * dt) / np.sum(dt))


def unmirror_march(res):
    """Re-index a reversed-flow march result into charge (depth) coordinates.

    Discharge marches from the opposite end of the tube, so its segment index
    runs backwards relative to the charge. `run_cycle` and `simulate_css` hand
    it a mirrored cascade and a mirrored initial state and un-mirror the result
    when they close the loop; anything that *plots* a discharge result has to
    do the same or it draws the well upside down.

    Returns a shallow copy with every per-segment array reversed, including the
    recorded history. `z`, `t` and `t_hist` are coordinates, not fields, and
    are left alone. Scalars are left alone. The cascade in depth coordinates is
    the same for both half-cycles, which is the point: after this call a
    discharge result is plotted against `T_m_seg`, not against `T_m_seg[::-1]`.
    """
    n = len(res["E"])
    out = dict(res)
    for k in ("E", "A_melt", "T_pcm", "delta", "eps_local", "merge_proximity"):
        if k in out and out[k] is not None:
            out[k] = np.asarray(out[k])[::-1].copy()
    h = res.get("history")
    if h:
        hh = {}
        for k, v in h.items():
            v = np.asarray(v)
            # T_fluid has n+1 nodes (inlet plus one per segment); the rest have n
            hh[k] = v[:, ::-1].copy() if v.ndim == 2 and v.shape[1] in (n, n + 1) else v
        out["history"] = hh
    return out


def march_h(case, T_inlet, T_m_seg, m_dot, k_wall, times,
            n_segments=None, E0=None, record=False):
    """Segment march on the ENTHALPY state. Replaces `march` from v0.4.

    Identical to `march` in every respect except what is integrated: the state
    is E' [J per m of tube] rather than A_melt [m2], so no heat is ever clipped
    and the PCM temperature is free to leave T_m.

    Recording convention. Every history frame is stamped, in `t_hist`, with the
    instant at which its state arrays (`E`, `A_melt`, `T_pcm`, `delta`) are
    exact. The profile arrays (`T_fluid`, `NTU`, `U_i`, `q_prime`) in the same
    frame are the ones evaluated from that state, i.e. the profile that drives
    the step *beginning* there. There are `len(t)+1` frames: one at t = 0 and
    one at the end of each accepted step, so the final state is in the record.

    Earlier builds recorded the pre-step state but stamped it with the
    post-step time, which on a logarithmic grid put the last frame a quarter of
    a charge behind its label, and appended `E` after the update while the
    other state arrays were from before it. Neither affected any marched or
    reported quantity -- `Q_cum`, `closure`, `eps_local` and the returned state
    all come from the march itself -- but every recorded plot was lagged.
    """
    n_segments = n_segments or case.n_segments
    dz = case.L_tube / n_segments
    r_i, r_e = case.r_i, case.r_e
    A_avail, E_lat, C_s, C_l = pcm_capacities(case)

    E = np.zeros(n_segments) if E0 is None else np.array(E0, float)
    E_start = E.copy()
    t_prev, Q_cum = 0.0, 0.0
    ts, Qs, T_outs = [], [], []
    hist = {k: [] for k in ("T_fluid", "T_pcm", "A_melt", "delta", "NTU",
                            "U_i", "q_prime", "E")} if record else None
    t_hist = []

    for t in np.asarray(times, float):
        dt = t - t_prev
        if dt <= 0:
            t_prev = t
            continue

        A, T_pcm = pcm_state(E, T_m_seg, case)
        # `delta` stays the MELT thickness: it is what the H3 merge diagnostic
        # is about. The conduction path is a separate question -- which side of
        # the front the heat has to cross -- and depends on the direction.
        delta = delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L)
        d_melt, k_melt = conduction_shell(case, A, A_avail, True, E, E_lat)
        d_froz, k_froz = conduction_shell(case, A, A_avail, False, E, E_lat)
        T0 = T_inlet
        q_prime = np.zeros(n_segments)
        if record:
            E_before = E.copy()          # the state this frame is stamped at
            T_prof = np.empty(n_segments + 1); T_prof[0] = T_inlet
            NTU_prof = np.empty(n_segments); U_prof = np.empty(n_segments)

        for i in range(n_segments):
            # Which shell the heat crosses depends on which way the front is
            # moving; the sign of the driving difference settles it before K is
            # known, and it has the same sign as q'.
            melting = T0 > T_pcm[i]
            d_path = float(d_melt[i] if melting else d_froz[i])
            k_m = float(k_melt[i] if melting else k_froz[i])
            h_i, cp_d, _, _ = h_internal(case.fluid2, T0, case.P, r_i, m_dot)
            U_i = compute_U_i(h_i, r_i, r_e, k_wall, case.Rf_i, k_m,
                              d_path, case.L_tube, case.fin_t,
                              case.fin_L, case.num_fins)
            NTU = float(np.clip((2 * np.pi * r_i * U_i * dz) / (m_dot * cp_d),
                                -50.0, 50.0))
            # conductance of the segment to the fluid, per unit tube length:
            # q' = K (T0 - T_pcm), from the NTU effectiveness
            K = m_dot * cp_d * (1.0 - np.exp(-NTU)) / dz
            E[i], q_prime[i] = advance_segment(E[i], T_m_seg[i], T0, K, dt,
                                               E_lat, C_s, C_l,
                                               case.sensible_heat)
            # the fluid gives up exactly what the PCM took
            T0 = T0 - q_prime[i] * dz / (m_dot * cp_d)
            if record:
                T_prof[i + 1] = T0; NTU_prof[i] = NTU; U_prof[i] = U_i

        if record:
            # stamped at t_prev, where E_before / A / T_pcm / delta are exact
            t_hist.append(t_prev)
            hist["T_fluid"].append(T_prof.copy()); hist["T_pcm"].append(T_pcm.copy())
            hist["A_melt"].append(A.copy());       hist["delta"].append(delta.copy())
            hist["NTU"].append(NTU_prof.copy());   hist["U_i"].append(U_prof.copy())
            hist["q_prime"].append(q_prime.copy()); hist["E"].append(E_before)

        Q_cum += float(np.sum(q_prime) * dz) * dt
        ts.append(t); Qs.append(float(np.sum(q_prime) * dz)); T_outs.append(T0)
        t_prev = t

    if record:
        # the end state, which the loop above never gets to record
        end = segment_profile(case, T_inlet, T_m_seg, m_dot, k_wall, E,
                              n_segments=n_segments)
        t_hist.append(t_prev)
        for k in hist:
            hist[k].append(end[k])

    A, T_pcm = pcm_state(E, T_m_seg, case)
    dE = float(np.sum(E - E_start) * dz)
    delta_end = delta_from_area(A, r_e, case.num_fins, case.fin_t, case.fin_L)

    # ---- H3 proximity ----------------------------------------------------
    # delta can never EXCEED delta_merge here (see Case.delta_merge), so the
    # useful diagnostic is not a violation flag but how close to the limit the
    # solution runs, and over how much of the well. The annular conduction
    # resistance is least trustworthy where this approaches 1.
    d_merge = case.delta_merge
    prox = delta_end / d_merge if d_merge > 0 else np.full_like(delta_end, np.nan)

    return {
        "E": E, "A_melt": A, "T_pcm": T_pcm,
        "delta": delta_end,
        "Q_cum_J": Q_cum,
        "closure": abs(Q_cum - dE) / abs(Q_cum) if Q_cum else np.nan,
        "t": np.array(ts), "Q": np.array(Qs), "T_out": np.array(T_outs),
        "V_melt_tube": float(np.sum(A * dz)),
        "eps_local": A / A_avail,
        "delta_merge": d_merge,
        "merge_proximity": prox,
        "merge_proximity_max": float(np.max(prox)),
        "f_near_merge": float(np.mean(prox > 0.90)),
        "E_sensible_J": float(np.sum(np.where(E > E_lat, E - E_lat,
                                              np.where(E < 0, E, 0.0))) * dz),
        "z": np.linspace(0.0, case.L_tube, n_segments),
        "history": {k: np.array(v) for k, v in hist.items()} if record else None,
        "t_hist": np.array(t_hist) if record else None,
    }


# ==========================================================================
# cycles + energy budget
# ==========================================================================

def T_m_bottom(case):
    """Melting temperature of the COLDEST cascade layer [K].

    The cascade is graded so that each layer sits DT_4C_M below the charging
    fluid at its own leading face, under a linear glide. The layers are spaced
    DT_3C_2C / N_lay apart, so the span of the melting temperatures is one
    layer width SHORT of the glide:

        T_m,bottom = T_m,top - DT_3C_2C * (N_lay - 1) / N_lay

    For the default case that is 150 - 55*8/9 = 101.111 C against a 55 K glide.
    The discharge enters the exchanger at this end, which is why this is the
    reference for DT_M_1D and not T_m,top.
    """
    return case.T_m - case.DT_3C_2C * (case.N_lay - 1) / case.N_lay


def cycle_state_points(case, T_2d=None):
    """ORC and heat-pump state points, and the two cycle efficiencies.

    Note the expansions and compressions carry no isentropic efficiency, so
    eta_ORC and COP are idealised -- only the electrical and mechanical
    efficiencies appear.

    THE DISCHARGE CLOSURE (changed in v0.6). Only the two exchanger INLETS are
    boundary conditions on the march; both outlets are results of the heat
    transfer and cannot be prescribed. Until v0.5a the discharge was closed the
    other way round -- the outlet was pinned at T_2d = T_m,top - DT_m_2D with
    DT_m_2D = 0, and the inlet derived by subtracting the glide. That is not
    merely arbitrary, it is contradicted by the model's own solution: at the
    start of a discharge at CSS the water leaves at 157.93 C, nearly 8 K above
    the top layer's melting point, because the PCM is
    superheated. It is also silently tied to N_lay -- the implied inlet
    approach was DT_3C_2C / N_lay, one layer width, chosen by nobody.

    The inlets are now prescribed symmetrically:

        charge:    T_4c = T_m,top    + DT_4C_M     (superheat, drives melting)
        discharge: T_3d = T_m,bottom - DT_M_1D     (subcooling, drives freezing)

    and T_2d is DEMOTED to a provisional estimate, the value it would take if
    the water achieved the full glide. It is provisional because the ORC needs
    an evaporating temperature before the store has been marched: T_3e is tied
    to T_2d. Pass the realised outlet back in as `T_2d` for the correction pass
    (see `simulate_css`); the estimate is good enough that one pass suffices.

    Setting DT_M_1D = DT_3C_2C / N_lay recovers the retired closure exactly,
    which is how the change is verified to be a no-op.
    """
    T_3d_C = T_m_bottom(case) - 273.15 - case.DT_M_1D
    # provisional unless the caller supplies the realised value
    T_2d_C = (T_3d_C + case.DT_3C_2C) if T_2d is None else (T_2d - 273.15)
    T = dict(
        T_1e=case.T_sink_C + case.DT_E_sink + 273.15,
        T_3e=T_2d_C - case.DT_2D_3E + 273.15,
        T_2d=T_2d_C + 273.15,
        T_3d=T_3d_C + 273.15,                     # the borehole inlet, state 1d
        T_4c=case.T_m_C + case.DT_4C_M + 273.15,
        T_4a=case.T_source_C - case.DT_3A_4A + 273.15,
        T_13h=case.T_source_C - case.DT_3A_13H + 273.15,
    )
    T["T_3c"] = T["T_4c"]
    T["T_2c"] = T["T_3c"] - case.DT_3C_2C
    T["T_2h"] = T["T_3c"] + case.DT_2H_3C

    rank = double_stage_rankine(case.fluid, T["T_1e"], T["T_3e"])
    hp = two_stage_htheatpump_2regs(case.refrig, T["T_13h"], T["T_2h"],
                                    case.DT_sub)
    for name, v in (("rank_eff", rank.get("rank_eff")),
                    ("hp_cop", hp.get("hp_cop"))):
        if v is None or not np.isfinite(v):
            raise RuntimeError(f"{name} is not finite ({v})")
    return rank, hp, T


def energy_budget(case, rank_eff, hp_cop, T):
    """Work backwards from the specified electrical output to the charging duty."""
    W_dot_T = case.W_dot_el_out / case.Turb_eff / case.ElG_eff
    Q_dot_in_ORC = W_dot_T / rank_eff
    D_E_in_ORC = Q_dot_in_ORC * case.t_dc * 3600.0                 # kJ
    D_E_out_HP = D_E_in_ORC * (1.0 + case.loss_surplus)
    Q_dot_out_HP = D_E_out_HP / case.t_ch / 3600.0
    W_dot_C_HP = Q_dot_out_HP / hp_cop
    W_dot_el_in = W_dot_C_HP / case.Comp_eff / case.ElH_eff

    cp_w = CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P",
                      case.P, case.fluid2) / 1000.0
    m_dot_w_ch = Q_dot_out_HP / cp_w / case.DT_3C_2C
    return dict(Q_dot_in_ORC=Q_dot_in_ORC, D_E_in_ORC=D_E_in_ORC,
                D_E_out_HP=D_E_out_HP, Q_dot_out_HP=Q_dot_out_HP,
                W_dot_el_in=W_dot_el_in, m_dot_w_ch=m_dot_w_ch)


def melting_temperatures(case):
    """Layer melting temperatures for charging and for discharging."""
    dTm = case.DT_3C_2C / case.N_lay
    i = np.arange(case.N_lay, dtype=float)
    return ((case.T_m - i * dTm).tolist(),
            (case.T_m - case.DT_3C_2C + (i + 1.0) * dTm).tolist())


# ==========================================================================
# sizing
# ==========================================================================


# ==========================================================================
# run   (charge, size, discharge, KPIs)
# ==========================================================================


# ==========================================================================
# cyclic steady state   (v0.5)
# ==========================================================================
# run_cycle() above answers "how many wells does one charge from a COLD START
# need?". That is a first-cycle question, and the cycle does not return to the
# state it started from, so the answer is not what a plant repeating a 24 h
# schedule would see.
#
# WHY SIZING CANNOT SIMPLY BE MOVED TO CSS. At cyclic steady state the state
# returns to itself, so the enthalpy change over a cycle is zero. This model has
# no loss path -- H2 no axial conduction, H6 adiabatic wall, H9 no azimuthal
# exchange -- so charge energy EQUALS discharge energy identically:
#
#       eta_storage(CSS) == 1     exactly, for every N and every flow.
#
# Two consequences follow, and both are structural rather than numerical.
#
#   1. The assumed loss surplus lambda is UNATTAINABLE at CSS. The budget
#      requires D_E_out_HP = (1+lambda) D_E_in_ORC, so the ratio of CSS charge
#      to requirement pins at 1/(1+lambda) for EVERY N. There is no root.
#      On the first cycle the surplus is not lost but PARKED in the store as
#      residual melt; at CSS there is nowhere left to park it.
#
#   2. Even with lambda = 0 the rate criterion goes DEGENERATE. Charge and
#      discharge are then the same number, so one equation is left for two
#      unknowns: it fixes the FLOW RATIO and says nothing about N.
#
# THE REFORMULATION. Stop asking the energy balance to determine N. Specify the
# hardware, march to CSS, and report what comes out:
#
#       N        from latent heat alone, D_E_out_HP / (rho V_well h_m)
#       m_ch     pinned by the charging glide
#       m_dc     pinned by the discharging glide
#
# Nothing is solved -- there is no root-find anywhere in this routine -- so the
# question of convergence does not arise. What was an unsatisfiable constraint
# becomes a reported output: the deviation of the delivered energy, and hence of
# the net electrical output, from its target.
#
# This also dissolves an older objection. size_well_field's docstring argues
# that the discharge flow cannot be pinned to its glide because the delivered
# energy then saturates below the requirement for any well count, leaving no
# root. Under this framing that saturation is not a failure; it is the answer.


def simulate_css(case, N=None, n_cycles=80, tol=1e-9, record=False, T_2d=None):
    """March a specified well field to cyclic steady state.

    One pass. `T_2d` overrides the provisional exchanger outlet used to set the
    ORC evaporating temperature; `simulate_css_corrected` is the wrapper that
    supplies the realised value. Called directly, this is the uncorrected
    result.

    No root-finding. `N` defaults to the latent-heat-only well count and both
    flows are pinned by their glides, so the routine is a pure simulation.

    `lambda` is forced to zero regardless of `case.loss_surplus`, because a
    non-zero surplus cannot be absorbed at CSS in a model with no loss path
    (see the note above). The override is reported in the returned dict as
    `lambda_overridden` so it can never happen silently.

    Returns the converged half-cycle states, the convergence history, and the
    deviation of the delivered energy from the target.
    """
    lambda_overridden = case.loss_surplus != 0.0
    case = case.with_(loss_surplus=0.0) if lambda_overridden else case

    rank, hp, T = cycle_state_points(case, T_2d=T_2d)
    Eb = energy_budget(case, rank["rank_eff"], hp["hp_cop"], T)
    n = case.n_segments
    T_m_lay, T_m_lay_dc = melting_temperatures(case)
    T_m_seg = layer_map(T_m_lay, n, case.N_lay)
    T_m_seg_dc = T_m_seg[::-1].copy()     # flow reverses; see run_cycle
    times_ch = np.logspace(0.0, np.log10(case.t_ch * 3600.0), case.n_times)
    times_dc = np.logspace(0.0, np.log10(case.t_dc * 3600.0), case.n_times)

    # ---- hardware and flows: all SPECIFIED, none solved -------------------
    if N is None:
        N = Eb["D_E_out_HP"] * 1000.0 / (case.rho_latent * case.V_well * case.h_m)
    cp_ch = CP.PropsSI("C", "T", 0.5 * (T["T_3c"] + T["T_2c"]), "P",
                       case.P, case.fluid2) / 1000.0
    cp_dc = CP.PropsSI("C", "T", 0.5 * (T["T_2d"] + T["T_3d"]), "P",
                       case.P, case.fluid2) / 1000.0
    m1_ch = Eb["Q_dot_out_HP"] / cp_ch / case.DT_3C_2C / (N * case.num_tubes)
    m1_dc = Eb["Q_dot_in_ORC"] / cp_dc / case.DT_3C_2C / (N * case.num_tubes)

    # ---- march until the state repeats -----------------------------------
    def cycle(E0, rec=False):
        ch = march_h(case, T["T_4c"], T_m_seg, m1_ch, case.k_wall, times_ch,
                     n_segments=n, E0=E0, record=rec)
        dc = march_h(case, T["T_3d"], T_m_seg_dc, m1_dc, case.k_wall, times_dc,
                     n_segments=n, E0=ch["E"][::-1], record=rec)
        return ch, dc

    E = np.zeros(n)
    history = []
    for k in range(n_cycles):
        ch, dc = cycle(E)
        Q_ch = ch["Q_cum_J"] * case.num_tubes * N / 1000.0          # kJ, field
        Q_dc = abs(dc["Q_cum_J"]) * case.num_tubes * N / 1000.0
        drift = float(np.max(np.abs(dc["E"][::-1] - E)))
        history.append(dict(cycle=k + 1, Q_charge_kJ=Q_ch, Q_discharge_kJ=Q_dc,
                            drift=drift,
                            eps_end_charge=float(ch["eps_local"].mean()),
                            eps_end_discharge=float(dc["eps_local"].mean())))
        E = dc["E"][::-1]                 # back to charge-march indexing
        if drift < tol:
            break

    if record:
        # One extra cycle from the converged state, this time recording. The
        # state repeats to `tol`, so this reproduces the last cycle; doing it
        # here rather than inside the loop means the profile sweep that closes
        # each recorded march is paid once instead of `cycles` times.
        ch, dc = cycle(E, rec=True)

    # The discharge marches from the far end, so its segment index runs
    # backwards. Everything returned to the caller is in CHARGE (depth)
    # indexing, which is why there is one cascade, `T_m_seg`, and not two.
    dc = unmirror_march(dc)

    req = Eb["D_E_in_ORC"]
    deviation = Q_dc / req - 1.0

    # ---- the outlet temperatures, which are results and not inputs ---------
    T_2d_realised = mixed_mean_outlet(dc)      # what the ORC evaporator sees
    T_2c_realised = mixed_mean_outlet(ch)      # what returns to the HP condenser
    glide_dc = T_2d_realised - T["T_3d"]
    glide_ch = T["T_4c"] - T_2c_realised
    return dict(
        N_wells=N, m1_ch=m1_ch, m1_dc=m1_dc, flow_ratio=m1_dc / m1_ch,
        cycles=len(history), converged=history[-1]["drift"] < tol, tol=tol,
        history=history, charge=ch, discharge=dc, E_css=E,
        Q_charge_kJ=Q_ch, Q_discharge_kJ=Q_dc, required_kJ=req,
        deviation=deviation,
        W_el_out_implied=case.W_dot_el_out * (1.0 + deviation),
        eta_storage=Q_dc / Q_ch, lambda_overridden=lambda_overridden,
        T_2d_assumed=T["T_2d"], T_2d_realised=T_2d_realised,
        T_2c_assumed=T["T_2c"], T_2c_realised=T_2c_realised,
        glide_dc=glide_dc, glide_ch=glide_ch,
        glide_ratio_dc=glide_dc / case.DT_3C_2C,
        T_2d_was_provisional=T_2d is None,
        T=T, budget=Eb, T_m_seg=T_m_seg)
        # no `T_m_seg_dc`: the discharge is returned in depth indexing, where
        # the cascade is the same physical object as on charge. The mirrored
        # copy is an internal detail of the march and does not leave here.


def performance_indices(case, r):
    """The performance indicators of the factorial study, evaluated at CSS.

    These were defined for the first-cycle framing, where the field delivered
    the target energy by construction. At cyclic steady state it delivers
    whatever it delivers, so every one of them is computed from the ENERGY THE
    FIELD ACTUALLY MOVES rather than from the target. That matters: charging
    with the target input while crediting the realised output would flatter the
    round-trip efficiency by exactly the deviation.

        eta_T       thermal -> net electric conversion on the discharge side,
                    W_el_out / Q_in_ORC. Identically eta_ORC * eta_turb * eta_gen,
                    so it is a property of the power block alone.
        eta_T_eff   the same with the discharge pumping parasitic charged
                    against it. This is the "effective discharge efficiency".
        eta_RTE     round-trip: net electricity out over electricity in, both
                    including their pumping parasitics.
        eps_RTE     eta_RTE weighted by the fraction of the store that actually
                    CYCLES, eps(end of charge) - eps(end of discharge).
        dE_therm    thermal energy delivered per well per cycle   [kWh]
        dE_elec     net electric energy delivered per well,
                    dE_therm * eta_T_eff                          [kWh]

    A note on eps_RTE, because it changed meaning and the change matters when
    reading a factorial table. The first-cycle version weighted by eps_m, the
    melted fraction at end of charge. Under the enthalpy formulation that
    SATURATES -- once a segment has melted all its PCM the further energy goes
    into superheat, which eps_m cannot see -- and at this design point it is
    1.0000 exactly, so eps_RTE was identical to eta_RTE and carried no
    information at all. The cycled fraction does not saturate: a store that
    fills completely but only half empties returns 0.5, which is the useful
    statement.

    Note also that eta_RTE is INDEPENDENT of the deviation, exactly. A field
    that delivers 1 % more than target also drew 1 % more in, because the charge
    flow is pinned by its glide in the same way. The deviation is a statement
    about plant size, not about efficiency.
    """
    T, Eb = r["T"], r["budget"]
    N = r["N_wells"]
    t_ch_s, t_dc_s = case.t_ch * 3600.0, case.t_dc * 3600.0

    # ---- what the field actually moved, per cycle ------------------------
    Q_dot_in_ORC = r["Q_discharge_kJ"] / t_dc_s          # kW, field
    Q_dot_out_HP = r["Q_charge_kJ"] / t_ch_s             # kW, field
    rank_eff = Eb["Q_dot_in_ORC"] and (case.W_dot_el_out / case.Turb_eff
                                       / case.ElG_eff / Eb["Q_dot_in_ORC"])
    W_el_out = Q_dot_in_ORC * rank_eff * case.Turb_eff * case.ElG_eff
    cop = Eb["Q_dot_out_HP"] / (Eb["W_dot_el_in"] * case.Comp_eff
                                * case.ElH_eff)
    W_el_in = Q_dot_out_HP / cop / case.Comp_eff / case.ElH_eff

    # ---- parasitics, at the realised temperatures ------------------------
    gv = case.geom_vector()
    _, pp_dc = calculate_pressure_drop(gv, r["m1_dc"], case.fluid2,
                                       r["T_2d_realised"], T["T_3d"], case.P)
    _, pp_ch = calculate_pressure_drop(gv, r["m1_ch"], case.fluid2,
                                       T["T_4c"], r["T_2c_realised"], case.P)
    pumping_dc = pp_dc * case.num_tubes * N / 1000.0     # kW, field
    pumping_ch = pp_ch * case.num_tubes * N / 1000.0

    # ---- the indicators --------------------------------------------------
    eta_T = W_el_out / Q_dot_in_ORC
    eta_T_eff = eta_T - pumping_dc / Q_dot_in_ORC
    eta_RTE = ((W_el_out - pumping_dc) * case.t_dc
               / ((W_el_in + pumping_ch) * case.t_ch))
    eps_cycled = float(r["charge"]["eps_local"].mean()
                       - r["discharge"]["eps_local"].mean())
    dE_therm = r["Q_discharge_kJ"] / 3600.0 / N          # kWh per well
    return dict(
        eta_T=eta_T, eta_T_eff=eta_T_eff,
        eta_RTE=eta_RTE, eps_RTE=eta_RTE * eps_cycled,
        eps_cycled=eps_cycled,
        dE_therm_kWh=dE_therm, dE_elec_kWh=dE_therm * eta_T_eff,
        pumping_ch_kW=pumping_ch, pumping_dc_kW=pumping_dc,
        f_pump=(pumping_ch + pumping_dc) / W_el_out,
        W_el_out_kW=W_el_out, W_el_in_kW=W_el_in,
        rho_E_kWh_m3=dE_therm / case.V_well)


def kpi_report(case, r, k=None):
    """Print the factorial-study indicators."""
    k = k or performance_indices(case, r)
    print("PERFORMANCE INDICATORS   (at cyclic steady state, per cycle)")
    print("=" * 66)
    print(f"  eta_T       discharge thermal -> electric   {k['eta_T']:9.4f}")
    print(f"  eta_T,eff   the same, less discharge pumping{k['eta_T_eff']:9.4f}"
          f"   <- effective")
    print(f"  eta_RTE     round trip, both parasitics     {k['eta_RTE']:9.4f}")
    print(f"  eps_RTE     eta_RTE x cycled fraction       {k['eps_RTE']:9.4f}")
    print()
    print(f"  dE_therm    thermal energy per well      {k['dE_therm_kWh']:12.1f} kWh")
    print(f"  dE_elec     net electric energy per well {k['dE_elec_kWh']:12.1f} kWh")
    print(f"  rho_E       thermal energy density       {k['rho_E_kWh_m3']:12.2f} kWh/m3")
    print()
    print(f"  cycled fraction of the store   {k['eps_cycled']:8.4f}"
          f"   (melted {r['charge']['eps_local'].mean():.4f},"
          f" residual {r['discharge']['eps_local'].mean():.4f})")
    print(f"  pumping, charge / discharge    {k['pumping_ch_kW']:8.2f}"
          f" / {k['pumping_dc_kW']:.2f} kW"
          f"   = {100*k['f_pump']:.2f} % of gross output")
    print()
    print("  eta_RTE is independent of the deviation: a field that delivers")
    print("  1 % above target also drew 1 % more in. The deviation is a")
    print("  statement about plant SIZE, not about efficiency.")


def simulate_css_corrected(case, N=None, record=False, **kw):
    """CSS with a single ORC correction pass. The v0.6 default path.

    The ORC needs an evaporating temperature before the store has been marched,
    because T_3e is tied to T_2d. Since v0.6 T_2d is a result, so the first pass
    uses the provisional estimate -- the outlet the water would reach if it
    achieved the full glide -- and the second re-evaluates the ORC, the energy
    budget, the well count and both flows at the outlet the first pass actually
    produced.

    ONE pass is enough, and that is a property of the problem rather than a
    hopeful choice: the outlet is pinned by the store, not by the plant. Over a
    sweep of the discharge subcooling from 6 to 15 K -- which moves the inlet by
    9 K and the deviation by twelve percentage points -- the realised outlet
    moves by less than 0.3 K. `dT_2d_pass` is returned so the assumption is
    visible rather than assumed; if it is ever large, iterate.

    This keeps the plant level decoupled from the store in the sense that
    matters for Figure "procedure": there is no outer loop, only a second
    evaluation.
    """
    first = simulate_css(case, N=N, record=False, **kw)
    second = simulate_css(case, N=N, record=record,
                          T_2d=first["T_2d_realised"], **kw)
    second["pass1"] = first
    second["dT_2d_pass"] = second["T_2d_realised"] - first["T_2d_realised"]
    second["corrected"] = True
    return second


def css_report(case, r):
    """Print the cyclic-steady-state result."""
    print("CYCLIC STEADY STATE")
    print("=" * 66)
    if r["lambda_overridden"]:
        print("  NOTE  lambda forced to 0. A non-zero loss surplus cannot be")
        print("        absorbed at CSS: the state returns to itself and this")
        print("        model has no loss path, so charge == discharge exactly.")
        print()
    print(f"  N_wells          {r['N_wells']:10.4f}   from latent heat alone, not solved")
    print(f"  m_dot charge     {r['m1_ch']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  m_dot discharge  {r['m1_dc']:10.4f} kg/s per leg-pair, pinned by the glide")
    print(f"  flow ratio       {r['flow_ratio']:10.4f}   a consequence, not a solve")
    print()
    if r["converged"]:
        print(f"  converged in {r['cycles']} cycles"
              f"   (drift {r['history'][-1]['drift']:.1e})")
    else:
        # Loud, because the result is still printed and looks perfectly normal.
        # Slow convergence is not rare: it scales with the well count, and a
        # case with low latent heat or many wells can need three times the
        # default. Re-run with a larger n_cycles before believing anything
        # below this line.
        print(f"  *** NOT CONVERGED *** stopped at the n_cycles limit of"
              f" {r['cycles']} with drift {r['history'][-1]['drift']:.2e},"
              f" tolerance {r['tol']:.0e}")
        print("      Re-run with a larger n_cycles. Everything below assumes")
        print("      the state has returned to itself, and it has not.")
    print(f"  eta_storage      {r['eta_storage']:10.6f}   must be 1: adiabatic, closed cycle")
    print()
    print(f"  required   D_E_in_ORC  {r['required_kJ']:12.5e} kJ")
    print(f"  delivered  at CSS      {r['Q_discharge_kJ']:12.5e} kJ")
    print(f"  DEVIATION              {100*r['deviation']:+12.3f} %")
    print()
    print(f"  implied net output     {r['W_el_out_implied']/1000:10.4f} MWe"
          f"   (target {case.W_dot_el_out/1000:.4f})")
    print()
    print("  the deviation IS the glide shortfall")
    print(f"    realised discharge glide / assumed   {r['glide_ratio_dc']:10.6f}")
    print(f"    delivered energy / required          "
          f"{r['Q_discharge_kJ']/r['required_kJ']:10.6f}")
    print("    The flow is pinned by the ASSUMED glide, so the energy")
    print("    delivered is just the ratio of realised to assumed glide.")
    print()
    print("  exchanger outlets are RESULTS, not inputs (v0.6)")
    print(f"    discharge, to the ORC       assumed {r['T_2d_assumed']-273.15:8.3f} C"
          f"   realised {r['T_2d_realised']-273.15:8.3f} C")
    print(f"    charge, to the HP condenser assumed {r['T_2c_assumed']-273.15:8.3f} C"
          f"   realised {r['T_2c_realised']-273.15:8.3f} C")
    if r.get("corrected"):
        print(f"    ORC corrected at the realised outlet; second pass moved it"
              f" {r['dT_2d_pass']:+.4f} K")
    elif r["T_2d_was_provisional"]:
        print("    UNCORRECTED: the ORC still runs on the provisional outlet.")
        print("    Use simulate_css_corrected for the reported result.")
    ch, dc = r["charge"], r["discharge"]
    print()
    print("  melted fraction at CSS")
    print(f"    end of charge     mean {ch['eps_local'].mean():.4f}"
          f"   min {ch['eps_local'].min():.4f}   max {ch['eps_local'].max():.4f}")
    print(f"    end of discharge  mean {dc['eps_local'].mean():.4f}"
          f"   min {dc['eps_local'].min():.4f}   max {dc['eps_local'].max():.4f}")
    if dc["eps_local"].mean() > 0.05 and ch["eps_local"].mean() > 0.99:
        print()
        print("    The store melts completely but does not refreeze completely,")
        print("    so the shortfall is a DISCHARGE RATE limit, not an inventory")
        print("    limit. More PCM per well is not the lever.")


def validate_case(case, verbose=True):
    """Check a Case for the inconsistencies that are easy to create.

    This is not a check that the model is right -- that lives in the
    verification notebook. It is a check that YOUR PARAMETERS are consistent,
    and it fires on the traps that have actually caught people:

      * `n_segments` not a multiple of `N_lay`, which makes `layer_map` place
        a layer boundary mid-segment;
      * a discharge inlet that is not below the coldest layer, or a charge
        inlet not above the hottest, so the half-cycle cannot run;
      * a glide so large that the coldest layer falls below the sink, or so
        small that the cascade collapses;
      * a time grid too coarse for the fastest segment time constant.

    Returns a list of (level, message). `level` is 'error' for a case that
    cannot be run and 'warn' for one that can but probably should not be.
    """
    msgs = []
    def err(m): msgs.append(('error', m))
    def warn(m): msgs.append(('warn', m))

    # --- the cascade -----------------------------------------------------
    if case.N_lay < 1:
        err(f'N_lay = {case.N_lay}; must be at least 1')
    if case.n_segments % case.N_lay:
        warn(f'n_segments = {case.n_segments} is not a multiple of '
             f'N_lay = {case.N_lay}: layer_map puts a layer boundary inside a '
             f'segment, so one segment carries the wrong melting temperature '
             f'by up to {case.DT_3C_2C / case.N_lay:.3f} K. '
             f'Nearest clean values: '
             f'{case.N_lay * (case.n_segments // case.N_lay)} or '
             f'{case.N_lay * (case.n_segments // case.N_lay + 1)}')

    T_bot = T_m_bottom(case) - 273.15
    spacing = case.DT_3C_2C / case.N_lay

    # --- the two inlets --------------------------------------------------
    T_in_ch = case.T_m_C + case.DT_4C_M
    T_in_dc = T_bot - case.DT_M_1D
    if case.DT_4C_M <= 0:
        err(f'DT_4C_M = {case.DT_4C_M} K; the charge inlet must be ABOVE the '
            f'hottest layer or nothing melts')
    if case.DT_M_1D <= 0:
        err(f'DT_M_1D = {case.DT_M_1D} K; the discharge inlet must be BELOW '
            f'the coldest layer or nothing freezes')
    if 0 < case.DT_4C_M < spacing:
        warn(f'DT_4C_M = {case.DT_4C_M:.3f} K is less than one layer width '
             f'({spacing:.3f} K), so the charging fluid drops below a layer\'s '
             f'melting point before leaving it')
    if 0 < case.DT_M_1D < spacing:
        warn(f'DT_M_1D = {case.DT_M_1D:.3f} K is less than one layer width '
             f'({spacing:.3f} K); the discharge pinches at each layer exit')

    # --- glide and the sink ----------------------------------------------
    if case.DT_3C_2C <= 0:
        err(f'DT_3C_2C = {case.DT_3C_2C} K; the glide sets both mass flows')
    if T_in_dc <= case.T_sink_C:
        err(f'discharge inlet {T_in_dc:.1f} C is at or below the sink '
            f'{case.T_sink_C:.1f} C')
    if T_in_ch <= T_bot:
        err(f'charge inlet {T_in_ch:.1f} C is below the coldest layer '
            f'{T_bot:.1f} C')

    # --- properties ------------------------------------------------------
    for f in ('h_m', 'cp_s', 'cp_l', 'k_s', 'k_l', 'rho_s', 'rho_l'):
        v = getattr(case, f, None)
        if v is not None and v <= 0:
            err(f'{f} = {v}; must be positive')

    # --- the time grid ---------------------------------------------------
    if case.n_times < 8:
        warn(f'n_times = {case.n_times}; the logarithmic grid needs enough '
             f'levels that the last step does not span most of the window')

    if verbose:
        if not msgs:
            print('validate_case: OK')
            print(f'  cascade   {case.T_m_C:.3f} C down to {T_bot:.3f} C '
                  f'in {case.N_lay} layers of {spacing:.3f} K')
            print(f'  inlets    charge {T_in_ch:.3f} C, '
                  f'discharge {T_in_dc:.3f} C')
        for level, m in msgs:
            print(f'{level.upper():5s} {m}')
    return msgs
