"""Build thums.py -- the LIVE model, for reading and for adapting.

The engineering build (`thums_model.py`) carries the first-cycle sizing chain
alongside the cyclic-steady-state path. That chain is not merely unused: the
model now argues that it does not apply to a repeating duty (report S7.2), and
it is the single most confusing thing in the sources -- two well counts, two
criteria, a root find, and a reporting class whose job is to make retired key
names fail loudly.

This script removes it, by REACHABILITY rather than by judgement. Anything not
called, directly or transitively, from the entry points below is dropped.

The result is the module the student imports and the module the notebooks run.
Nothing here is retyped: every function is the same text as in the engineering
build, so the two cannot say different things.
"""
import ast
import pathlib
import re

SRC = pathlib.Path('/tmp/build/thums_model.py')
OUT = pathlib.Path('/tmp/build/thums.py')

# The public surface of the v0.6 procedure. `calculate_pressure_drop` is an
# entry point in its own right because the parasitic power is reported but is
# not called by the CSS march.
ENTRY = ('Case', 'simulate_css_corrected', 'simulate_css', 'css_report',
         'performance_indices', 'kpi_report',
         'cycle_state_points', 'energy_budget', 'melting_temperatures',
         'T_m_bottom', 'layer_map', 'march_h', 'pcm_state', 'pcm_capacities',
         'segment_profile', 'unmirror_march', 'mixed_mean_outlet',
         'orc_pinch', 'feasible_rankine',
         'conduction_shell', 'bulk_shape_factor', 'bulk_equivalent_delta',
         'delta_from_area', 'calculate_pressure_drop')

HEADER = '''"""THUMS -- latent heat storage in a repurposed wellbore.

The live model, v0.9. One formulation (enthalpy, "Formulation C"), one sizing
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
'''

VALIDATE = '''

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
             f'({spacing:.3f} K), so the charging fluid drops below a layer\\'s '
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

    # --- the working fluid is liquid water, and must stay that way --------
    T_sat_w = CP.PropsSI('T', 'P', case.P, 'Q', 0, 'Water') - 273.15
    if T_in_ch >= T_sat_w:
        err(f'the charging inlet T_4c = {T_in_ch:.2f} C is at or above the '
            f'saturation temperature of water at P = {case.P/1e5:.2f} bar, '
            f'which is {T_sat_w:.2f} C. The model assumes single-phase '
            f'liquid throughout (H7) and carries no boiling correlation, so '
            f'it would report nonsense rather than fail. Raise case.P or '
            f'lower T_m_C + DT_4C_M.')
    elif T_sat_w - T_in_ch < 10.0:
        warn(f'the charging inlet is {T_sat_w - T_in_ch:.2f} K below the '
             f'saturation temperature of water at {case.P/1e5:.2f} bar; '
             f'there is very little margin against boiling')

    # --- approaches that are zero by construction -------------------------
    # None of these is an error. Each is a place where the model quietly
    # assumes an infinite exchanger, and a parametric study that leans on it
    # will report an efficiency no hardware can reach.
    if case.DT_3A_4A == case.DT_3A_13H:
        warn(f'DT_3A_4A == DT_3A_13H == {case.DT_3A_4A:.3f} K, so the '
             f'heat-pump source leaves the evaporator at exactly the '
             f'evaporating temperature: a zero approach at that end, which '
             f'needs infinite area. Make DT_3A_4A the smaller of the two for '
             f'a finite evaporator.')
    warn('the ORC regenerator has a zero approach by construction '
         '(T_9e = T_7e in double_stage_rankine), so its duty is an upper '
         'bound rather than a design value. This is structural, not a '
         'setting you can change from Case.')

    # --- the ORC evaporator has to obey the second law --------------------
    try:
        rank, _hp, Tst = cycle_state_points(case)
        got = rank.get('pinch')
        old = rank.get('T_3e_uncorrected')
        if old is not None and rank['T_3e'] < old - 1e-6:
            rank_old = double_stage_rankine(case.fluid, Tst['T_1e'], old)
            warn(f'the ORC pinch forced the boiling temperature DOWN, from '
                 f'{old - 273.15:.2f} C to {rank["T_3e"] - 273.15:.2f} C, and '
                 f'with it eta_ORC from {rank_old["rank_eff"]:.4f} to '
                 f'{rank["rank_eff"]:.4f}. The hot-end rule DT_2D_3E = '
                 f'{case.DT_2D_3E:.1f} K is not what sets the evaporating '
                 f'temperature here; the pinch at DT_pinch_ORC = '
                 f'{case.DT_pinch_ORC:.1f} K is. This is the honest cost of a '
                 f'{case.DT_3C_2C:.0f} K water glide against a fluid that '
                 f'boils at one temperature.')
        elif got is not None:
            msgs.append(('info', f'ORC evaporator pinch {got:+.2f} K, '
                                 f'requirement {case.DT_pinch_ORC:.1f} K'))
    except Exception as e:                       # never block on a diagnostic
        warn(f'could not evaluate the ORC pinch: {e}')

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
'''


def build():
    src = SRC.read_text()
    tree = ast.parse(src)
    defs = {n.name: n for n in tree.body
            if isinstance(n, (ast.FunctionDef, ast.ClassDef))}

    calls = {}
    for name, node in defs.items():
        s = set()
        for sub in ast.walk(node):
            if isinstance(sub, ast.Call):
                f = sub.func
                if isinstance(f, ast.Name):
                    s.add(f.id)
                elif isinstance(f, ast.Attribute):
                    s.add(f.attr)
        calls[name] = s & set(defs)

    live, stack = set(), list(ENTRY)
    while stack:
        x = stack.pop()
        if x in live or x not in defs:
            continue
        live.add(x)
        stack += list(calls[x])
    dropped = sorted(set(defs) - live)

    lines = src.splitlines(True)
    keep = []
    # module-level statements that are not dropped definitions, in file order
    spans = sorted((defs[n].lineno, defs[n].end_lineno, n) for n in defs)
    drop_ranges = []
    for a, b, n in spans:
        if n in dropped:
            # swallow any decorator lines and the comment block above
            start = a - 1
            while start > 0 and (lines[start - 1].lstrip().startswith('@')):
                start -= 1
            while start > 0 and lines[start - 1].lstrip().startswith('#'):
                start -= 1
            drop_ranges.append((start, b))
    drop = set()
    for a, b in drop_ranges:
        drop.update(range(a, b))

    body = ''.join(l for i, l in enumerate(lines) if i not in drop)
    # the original module docstring goes; ours replaces it
    body = re.sub(r'\A""".*?"""\n', '', body, count=1, flags=re.S)
    body = re.sub(r'\n{4,}', '\n\n\n', body)

    out = HEADER + body.rstrip() + '\n' + VALIDATE
    OUT.write_text(out)

    n_all = len(src.splitlines())
    n_out = len(out.splitlines())
    print(f'wrote {OUT}')
    print(f'  kept {len(live)} of {len(defs)} objects; '
          f'{n_all} -> {n_out} lines')
    print(f'  dropped: {", ".join(dropped)}')
    return dropped


if __name__ == '__main__':
    build()
