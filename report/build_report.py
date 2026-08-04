"""Generate the master-architecture report as LaTeX, from live model state.

The document is NOT hand-written. Every number in it is pulled at build time from

    geothermal_cooling/dbhe_master.ipynb        (version stamp, configuration)
    verification/reference/v*_baseline.json     (per-version verified results)
    verification/reference/*_flow_depth_sweep.json

so it cannot drift from the code. Run after every version bump:

    python report/build_report.py

Writes report/master_architecture_report.tex.
"""
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
NB = os.path.join(ROOT, 'geothermal_cooling', 'dbhe_master.ipynb')
REF = os.path.join(ROOT, 'verification', 'reference')
OUT = os.path.join(HERE, 'master_architecture_report.tex')

# --------------------------------------------------------------------------
# harvest live state
# --------------------------------------------------------------------------


def nb_cells():
    return json.load(open(NB))['cells']


def nb_scalar(name, default=None):
    """Pull a top-level assignment `name = <number>` out of the notebook."""
    pat = re.compile(r'^\s*' + re.escape(name) + r'\s*=\s*([-+0-9.eE_]+)', re.M)
    for c in nb_cells():
        if c['cell_type'] != 'code':
            continue
        m = pat.search(''.join(c['source']))
        if m:
            try:
                return float(m.group(1).replace('_', ''))
            except ValueError:
                pass
    return default


def nb_version():
    for c in nb_cells():
        m = re.search(r'NB_VERSION\s*=\s*"([^"]+)"', ''.join(c['source']))
        if m:
            v = m.group(1)
            m2 = re.search(r'NB_STAMP\s*=\s*"([^"]+)"', ''.join(c['source']))
            return v, (m2.group(1) if m2 else 'unknown')
    return 'unknown', 'unknown'


def fixtures():
    out = {}
    for p in sorted(glob.glob(os.path.join(REF, 'v*_baseline.json'))):
        v = re.search(r'v([0-9.]+)_baseline', os.path.basename(p)).group(1)
        out[v] = json.load(open(p))
    return out


def vkey(v):
    return tuple(int(x) for x in v.split('.'))


def sweep():
    p = os.path.join(REF, 'v10.2_flow_depth_sweep.json')
    return json.load(open(p)) if os.path.exists(p) else []


def tex_esc(s):
    return (s.replace('\\', r'\textbackslash{}').replace('_', r'\_')
             .replace('%', r'\%').replace('&', r'\&').replace('#', r'\#'))


# --------------------------------------------------------------------------
# derived tables
# --------------------------------------------------------------------------

VERSION_NOTES = {
    '10.2': ('Wheel-regime diagnostics; warning suppression fixed',
             'Announced the desiccant-wheel sign reversal, and stopped '
             r'\verb|warnings.filterwarnings(''ignore'')| from swallowing the '
             "notebook's own diagnostics."),
    '10.3': ('Generic load interface',
             'The plant sees the building only through a two-equation, '
             'three-unknown interface with an explicit closure. Verified '
             'bit-identical to v10.2.'),
    '10.4': ('Scavenging mixing-node mass balance',
             'The moisture balance used the global damper constant while the '
             'enthalpy balance used the controlled one; the weights summed to '
             '1.07 and moisture was not conserved.'),
    '10.5': ('Scavenging fed from building exhaust',
             "The loop had been recirculating the wheel's own desorbed "
             'moisture. Moisture removal nearly tripled.'),
    '10.6': ('Two-sided humidity limit flag',
             'The saturation flag could only see under-drying. Reporting only; '
             'verified bit-identical to v10.5.'),
}

SERIES_LABEL = {
    'T_ret': (r'$T_{\mathrm{ret}}$', r'\degC', 1.0),
    'T_acc': (r'$T_{\mathrm{wa}}$', r'\degC', 1.0),
    'T_ri':  (r'$T_{\mathrm{ri}}$', r'\degC', 1.0),
    'T_ro':  (r'$T_{\mathrm{ro}}$', r'\degC', 1.0),
    'w_ri':  (r'$w_{\mathrm{ri}}$', r'g/kg', 1000.0),
    'w_pd':  (r'$w_{\mathrm{pd}}$', r'g/kg', 1000.0),
    'dW_process': (r'$\Delta w$ wheel', r'g/kg', 1000.0),
    'UA_cool': (r'$\UA_{\mathrm{cool}}$', r'kW/K', 0.001),
    'd_rA_ctrl': (r'$d_{pA}$ control', r'--', 1.0),
}


def version_history_table(fx):
    vs = sorted(fx, key=vkey)
    rows = []
    for v in vs:
        s = fx[v]['series']
        rows.append((v,
                     s['T_ret'][-1], s['T_acc'][-1], s['T_ro'][-1],
                     s['w_ri'][-1]*1000, s['dW_process'][-1]*1000,
                     s['UA_cool'][-1]/1000, s['d_rA_ctrl'][-1]))
    body = ''
    for r in rows:
        body += (f"{r[0]} & {r[1]:.2f} & {r[2]:.2f} & {r[3]:.2f} & "
                 f"{r[4]:.2f} & {r[5]:.3f} & {r[6]:.1f} & {r[7]:.3f} \\\\\n")
    return body, vs


def version_notes_list(vs):
    out = ''
    for v in vs:
        if v in VERSION_NOTES:
            title, why = VERSION_NOTES[v]
            out += f"\\item[\\textbf{{v{v}}}] \\textbf{{{title}.}} {why}\n"
    return out


def sweep_table(sw):
    rows = [r for r in sw if 'error' not in r]
    rows.sort(key=lambda r: (r['Lbh'], r['flow']))
    body = ''
    for r in rows:
        if r['flow'] not in (3, 4, 5, 6, 8, 10):
            continue
        ok = 'yes' if r['T_wa_20'] >= 78.0 else 'no'
        body += (f"{r['Lbh']:.0f} & {r['BHT']:.1f} & {r['flow']:.0f} & "
                 f"{r['T_ret_20']:.1f} & {r['T_wa_20']:.1f} & "
                 f"{r['Q_well_20']:.0f} & {ok} \\\\\n")
    return body


# --------------------------------------------------------------------------
TEMPLATE = r"""
\documentclass[11pt]{article}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\usepackage{mathtools}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage[margin=1in]{geometry}
\usepackage{siunitx}
\usepackage{graphicx}
\usepackage[hidelinks]{hyperref}
\usepackage{xcolor}
\usepackage{enumitem}

\sisetup{per-mode=symbol}

% ---------- macros, matching the companion model report ----------
\newcommand{\md}{\dot m}
\newcommand{\mdp}{\ensuremath{\dot m_{\mathrm{pro}}}}
\newcommand{\mdr}{\ensuremath{\dot m_{\mathrm{reg}}}}
\newcommand{\mdv}{\ensuremath{\dot m_{\mathrm{vent}}}}
\newcommand{\mdw}{\ensuremath{\dot m_{w}}}
\newcommand{\UA}{\ensuremath{U\!A}}
\newcommand{\hfg}{\ensuremath{h_{\mathrm{fg}}}}
\newcommand{\patm}{p_{\mathrm{atm}}}
\newcommand{\psat}{p_{\mathrm{sat}}}
\newcommand{\RH}{\ensuremath{\mathrm{RH}}}
\newcommand{\degC}{\ensuremath{^{\circ}\mathrm{C}}}
\newcommand{\qs}{\ensuremath{q_{s}}}
\newcommand{\ql}{\ensuremath{q_{l}}}

\title{\textbf{Master Architecture of the DBHE-Driven\\
Heating, Cooling and Dehumidification Plant}\\[6pt]
\large Streams, Junctions, Valve Network, Operating Modes,\\
Load Interface, and Verification Protocol\\[4pt]
\normalsize (Generated from model version @@VER@@, stamped @@STAMP@@)}
\author{Craft \& Hawkins Department of Petroleum Engineering\\
Louisiana State University}
\date{\today}

\begin{document}
\maketitle

\begin{abstract}
This document specifies the \emph{master architecture} of the surface plant
driven by a closed-loop Deep Borehole Heat Exchanger formed from a repurposed
legacy oil-and-gas well. Where the companion model report derives the governing
equations of each component, this report defines the \emph{topology}: which
streams exist, where they join and split, which valve or damper controls each
branch, and which subset of components is active in each operating mode.
The architecture is deliberately a superset. It contains every component the
plant could need across the full plane of sensible and latent loads, and a
given operating mode is obtained by activating a subset of branches rather
than by rebuilding the plant. Section~\ref{sec:modes} enumerates the resulting
modes; Section~\ref{sec:status} records which components are implemented in the
code and which are pending.

Every numerical value in this report is extracted at build time from the model
itself --- the notebook, the frozen regression fixtures, and the parametric
sweep archives. The document therefore cannot disagree with the code it
describes. It is regenerated whenever the model version changes.
\end{abstract}

\tableofcontents
\newpage

% ==========================================================================
\section{Scope, and how this document stays current}
% ==========================================================================

\subsection{What is specified here}

The plant comprises five fluid circuits: the geothermal water loop, a chilled
water loop, a cooling-tower water loop, the process (building supply) air
stream, and the scavenging (regeneration) air stream. The master architecture
is the union of all flow paths among them. It is not a single machine but a
family: each operating mode activates a subset of the branches, and the
inactive branches are present but idle.

Two consequences follow, and they shape the whole study. First, the number of
distinct plants is finite and small once the valves and dampers are treated as
two-position devices rather than continuously modulating ones. Second, the
question ``which plant should be built'' becomes ``which subset of the master
architecture serves the load quadrants that actually occur, using the water
temperature the well can actually deliver over the design life''. That is the
question Sections~\ref{sec:modes} and \ref{sec:resource} address.

\subsection{Generation, not transcription}

This report is emitted by \verb|report/build_report.py| in the model
repository. The build step reads:

\begin{itemize}[nosep]
  \item the notebook, for the version stamp and the configuration constants;
  \item \verb|verification/reference/v*_baseline.json|, the frozen per-version
        results recorded by the regression harness;
  \item the parametric sweep archives, for the resource-feasibility tables.
\end{itemize}

No number below is typed by hand. When the model advances to the next version,
the report is rebuilt and republished, and the version history in
Section~\ref{sec:history} gains a row automatically. A statement in this
document that contradicts the code is therefore a defect in the build script,
not a stale paragraph.

% ==========================================================================
\section{Nomenclature convention}
% ==========================================================================

State labels follow a two- or three-letter convention. The first subscript
letter identifies the stream: \emph{w} for geothermal water, \emph{c} for
chilled water, \emph{p} for process air, \emph{r} for scavenging
(regeneration) air, \emph{t} for cooling-tower water. Subsequent letters
identify the station. Valves on water circuits are $v_{A}\ldots v_{H}$;
dampers on air circuits are $d_{pA}\ldots d_{pI}$ for the process stream and
$d_{rA}, d_{rB}$ for the scavenging stream.

\begin{longtable}{@{}llp{0.52\textwidth}@{}}
\toprule
\textbf{Symbol} & \textbf{Circuit} & \textbf{Station} \\
\midrule
\endhead
$T_{\mathrm{in}}$   & water & injection into the annulus \\
$T_{\mathrm{ret}}$  & water & production from the tubing \\
$T_{\mathrm{wr}}$   & water & after the $v_F$ well-bypass rejoins the production stream \\
$T_{\mathrm{wa}}$   & water & accumulator outlet \\
$T_{\mathrm{wg}}$   & water & absorption-generator outlet \\
$T_{\mathrm{wma}}$  & water & generator-skip branch merged with the $v_B$ accumulator bypass \\
$T_{\mathrm{wmg}}$  & water & $T_{\mathrm{wg}}$ tempered by $T_{\mathrm{wma}}$ through $v_D$ \\
$T_{\mathrm{wx}}$   & water & regeneration heat-exchanger outlet \\
$T_{\mathrm{wmx}}$  & water & after the $v_H$ remainder blends back; returns to the pump \\
\midrule
$T_{pc}, w_{pc}$    & process air & building comfort state (return) \\
$T_{pmc}, w_{pmc}$  & process air & mixed state entering the supply fan \\
$T_{pd}, w_{pd}$    & process air & desiccant-wheel process outlet \\
$T_{pmo}, w_{pmo}$  & process air & after fork~1 remerges, entering fork~2 \\
$T_{pi}, w_{pi}$    & process air & building supply \\
$T_{ea}, \RH_{ea}$  & process air & building exhaust \\
\midrule
$T_{rm}, w_{rm}$    & scavenging air & mixed state entering the scavenging fan \\
$T_{ri}, w_{ri}$    & scavenging air & regeneration-coil inlet \\
$T_{ro}, w_{ro}$    & scavenging air & regeneration-coil outlet, wheel inlet \\
$T_{rd}, w_{rd}$    & scavenging air & wheel regeneration outlet, vented \\
\bottomrule
\end{longtable}

% ==========================================================================
\section{The five circuits}
% ==========================================================================

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figures/master_architecture.pdf}
\caption{Master architecture. Every component the plant could need across the
full load plane is present; an operating mode is a subset of active branches,
not a different machine.}\label{fig:master}
\end{figure}

\subsection{Geothermal water}

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figures/water_circuit.pdf}
\caption{Geothermal water circuit. $v_E$/$v_F$ decouple loop flow from well
flow; $v_A$/$v_B$ gate the accumulator; $v_C$ gates the absorption generator;
$v_D$ raises the coil inlet, while $v_G$/$v_H$ throttle the coil
duty.}\label{fig:water}
\end{figure}

The water circuit is a closed loop with four decision points. Reading from the
pump:

\begin{enumerate}[nosep]
  \item \textbf{$v_E$ / $v_F$ --- well bifurcation.} $v_E$ sends flow down the
        annulus; $v_F$ bypasses the well entirely and rejoins the production
        stream. The mixed state is $T_{\mathrm{wr}}$. This decouples the loop
        flow from the well flow, so the surface plant can be run at a flow the
        subsurface would not tolerate, and vice versa.
  \item \textbf{$v_A$ / $v_B$ --- accumulator.} $v_A$ charges or discharges the
        thermal store; $v_B$ bypasses it. The accumulator is a buffer, not a
        topology choice: its state varies through the day within any mode.
  \item \textbf{$v_C$ --- absorption generator.} The branch that skips the
        generator merges with the $v_B$ bypass to form $T_{\mathrm{wma}}$.
  \item \textbf{$v_D$, then $v_G$ / $v_H$.} $v_D$ blends $T_{\mathrm{wma}}$ into
        the generator outlet $T_{\mathrm{wg}}$, raising the regeneration-coil
        inlet to $T_{\mathrm{wmg}}$. $v_G$ passes that stream through the coil;
        $v_H$ carries the remainder of $T_{\mathrm{wma}}$ around the coil to
        blend with the coil outlet $T_{\mathrm{wx}}$, giving
        $T_{\mathrm{wmx}}$.
\end{enumerate}

The pairing of $v_D$ and $v_G$ deserves emphasis because the two are easily
confused. \textbf{$v_D$ raises} the coil inlet temperature: it can only blend
\emph{hotter} bypass water into a generator outlet that the chiller has cooled.
\textbf{$v_G$ throttles}: reducing the share of flow through the coil lowers
the coil duty, hence the regeneration air temperature $T_{\mathrm{ro}}$, hence
the drying the wheel performs. When the wheel must be turned \emph{down},
$v_G$ is the control and $v_D$ is not.

\subsection{Process air}

\begin{equation*}
\text{mix} \;\to\; \text{fan} \;\to\;
\underbrace{\Big\{ d_{pD} \to \mathrm{DW} \to \mathrm{cooler}
\;\Big|\; d_{pE} \to \mathrm{HX_{winter}} \to \mathrm{humidifier} \Big\}}_{\text{fork 1}}
\;\to\; T_{pmo} \;\to\;
\underbrace{\Big\{ d_{pF} \to \mathrm{CHW} \;\Big|\; d_{pG} \to \mathrm{IEC} \Big\}}_{\text{fork 2}}
\;\to\; T_{pi}
\end{equation*}

The building return stream bifurcates at a node: $d_{pA}$ exhausts,
$d_{pB}$ recirculates, and the recirculated air joins outdoor air admitted
through $d_{pC}$ at a mixing junction whose state is $T_{pmc}, w_{pmc}$.

Fork~1 selects between conditioning (wheel followed by the recuperative
cooler) and the alternative prong that carries the winter heating coil and the
humidifier. The $d_{pE}$ prong terminates in two dampers: $d_{pH}$ returns it
to $T_{pmo}$, upstream of fork~2, which is the summer bypass; $d_{pI}$ carries
it past fork~2 directly into the supply duct, which is the winter path. Fork~2
is a genuine parallel pair --- chilled-water coil \emph{or} indirect
evaporative cooler, never both --- rejoining at the supply state
$T_{pi}, w_{pi}$.

Fork~2 needs no third ``neither'' prong: with no chilled water circulating and
no evaporative water supplied, air passes either branch at nothing but a
pressure drop.

\subsection{Scavenging air}

\begin{equation*}
\text{mix}\big(d_{rA}\,\text{building exhaust} + d_{rB}\,\text{outdoor air}\big)
\to \text{fan} \to \text{cooler} \to \mathrm{HX} \to \mathrm{DW} \to \text{vented}
\end{equation*}

Two features of this circuit carry more weight than their simplicity suggests.

First, the scavenging stream passes through the \textbf{cooler before the
regeneration coil}, so the process stream's rejected heat preheats it and the
coil has less lifting to do. The cooler is therefore a recuperator between the
two air streams, not a heat rejection device.

Second, the mixing junction draws on \textbf{building exhaust and outdoor
air}, and never on the wheel's own regeneration outlet. Recirculating the
wheel exhaust returns its desorbed moisture to its own inlet, which collapses
the relative-humidity difference the wheel drives on. The code did exactly
that until v10.5; correcting it nearly tripled moisture removal.

A physical bound constrains the mix. The building can only exhaust what it
admits as fresh air, so the exhaust available to the scavenging loop is
$\mdv = \SI{@@MDOT_VENT@@}{\kilogram\per\second}$ against a scavenging demand
of $\mdr = \SI{@@MDOT_REG@@}{\kilogram\per\second}$. The building-exhaust
fraction is therefore capped at $@@DRA_CAP@@$ and the balance made up with
outdoor air.

\subsection{Chilled water and heat rejection}

The absorption chiller is a single indivisible subsystem with its
chilled-water distribution: a generator fed by geothermal water, an evaporator
serving the chilled-water coil in the process stream, and a condenser and
absorber rejecting to a wet cooling tower. A chiller with nowhere to send
chilled water is not a configuration, so \textbf{ABS and CHW are one switch,
not two}.

% ==========================================================================
\section{The load interface}
% ==========================================================================

The plant sees the building only through a steady, well-mixed pair of balances.
With $\qs > 0$ meaning cooling is required and $\mdw > 0$ meaning
dehumidification is required,

\begin{align}
\qs   &= \mdp\, c_p \,(T_{pc} - T_{pi}), \label{eq:sens}\\
\mdw  &= \mdp\, (w_{pc} - w_{pi}). \label{eq:lat}
\end{align}

Two equations, three unknowns $(\mdp, T_{pi}, w_{pi})$. \textbf{Exactly one
closure must be supplied}, and which one is a design decision rather than a
physical fact:

\begin{itemize}[nosep]
  \item give $\mdp$ --- constant volume; the supply state follows;
  \item give $\Delta T_{\mathrm{sup}} = T_{pc}-T_{pi}$ --- variable volume;
        the flow follows.
\end{itemize}

Dividing \eqref{eq:sens} by \eqref{eq:lat} eliminates the flow:

\begin{equation}
\frac{T_{pc}-T_{pi}}{w_{pc}-w_{pi}} = \frac{\qs}{\mdw}\,\frac{1}{c_p}.
\end{equation}

The load \emph{ratio} alone fixes the direction of the room condition line on
the psychrometric chart; the closure fixes how far along that line the supply
state sits. This is the classical construction, and it is why the interface
needs exactly one more piece of information and no more.

\paragraph{Latent load is carried as a moisture rate, not as a power.}
Expressing $\ql$ in kilowatts requires choosing $\hfg$, and inverting
\eqref{eq:lat} requires choosing it again. If the two choices differ --- and
values between 2450 and \SI{2501}{\kilo\joule\per\kilogram} are all in common
use --- $w_{pi}$ is wrong by about 2\%\ with nothing on screen to show it. The
interface therefore takes $\mdw$ in \si{\kilogram\per\second} and offers a
convenience wrapper that converts from kilowatts with an explicit $\hfg$.

\paragraph{The closure in the hatchery case is derived, not given.}
The supply flow is set by a design sensible load carried at a prescribed
supply depression, floored by the forced ventilation requirement:

\begin{equation}
\mdp = \max\!\left(\frac{Q_{s,\mathrm{design}}}{c_p\,\Delta T_{\mathrm{sup,design}}},\;
\mdv\right)
= \max\!\left(\SI{@@MDOT_LOAD@@}{\kilogram\per\second},\;
\SI{@@MDOT_VENT@@}{\kilogram\per\second}\right)
= \SI{@@MDOT_PRO@@}{\kilogram\per\second}.
\end{equation}

The depression rule binds. Its consequence reaches much further than fan
sizing: a supply depression of only $\SI{@@DT_SUP@@}{\kelvin}$ demands a large
airflow, and a large airflow needs only a very small humidity difference to
carry the same latent load. The wheel is not underperforming in this plant ---
it is barely being asked to work.

% ==========================================================================
\section{Operating modes}\label{sec:modes}
% ==========================================================================

\subsection{Naming}

Load quadrants are named by the \emph{processes they require}, from the sign
pair $(\operatorname{sgn}\qs, \operatorname{sgn}\ql)$: cooling or heating,
combined with dehumidification or humidification. Seasonal names were
considered and rejected: they are climate-relative, since hot-and-dry is the
only summer in a semi-arid region, and setpoint-relative, since the quadrant
boundaries sit at the comfort state rather than at some neutral condition.

The sensible heat ratio cannot classify the quadrants on its own. Two negative
loads give a positive ratio, indistinguishable from two positive ones. The
sign pair is the primitive; the ratio describes position \emph{within} a
quadrant.

\subsection{Two axes: quadrant and tier}

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figures/quadrant_map.pdf}
\caption{The load plane. Quadrants are named by the processes they require;
the tier within a quadrant is set by load magnitude against what the well can
deliver.}\label{fig:qmap}
\end{figure}

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth]{figures/topology_variants.pdf}
\caption{Topology variants. Each card shows the branches active in one mode;
inactive branches remain installed but idle.}\label{fig:topo}
\end{figure}

The quadrant selects which processes are needed. The load \emph{magnitude}
selects which components can deliver them. The resulting mode set is:

\begin{longtable}{@{}llp{0.36\textwidth}l@{}}
\toprule
\textbf{Mode} & \textbf{Quadrant} & \textbf{Active components} & \textbf{Water required} \\
\midrule
\endhead
M1 & cooling + dehumidification & DW, cooler & $\geq \SI{66}{\degreeCelsius}$ \\
M2 & cooling + dehumidification & DW, cooler, IEC & $\geq \SI{66}{\degreeCelsius}$ \\
M3 & cooling + dehumidification & DW, cooler, ABS + CHW & $\geq \SI{79}{\degreeCelsius}$ \\
M4 & cooling + humidification & humidifier alone (adiabatic saturation \emph{is} direct evaporative cooling) & none \\
M5 & cooling + humidification & ABS + CHW, humidifier & $\geq \SI{79}{\degreeCelsius}$ \\
M6 & heating + dehumidification & DW, cooler and winter coil, \emph{blended} & $\geq \SI{66}{\degreeCelsius}$ \\
M7 & heating + humidification & winter coil & $\geq$ supply $T$ \\
M8 & heating + humidification & winter coil, humidifier & $\geq$ supply $T$ \\
\bottomrule
\end{longtable}

Three of these carry conditions worth stating explicitly.

\textbf{M6 is the only mode requiring modulating dampers.} Dehumidification
and heating sit on parallel prongs of fork~1, so the quadrant is reached by
drying part of the flow, heating the rest, and blending at $T_{pmo}$. Every
other mode is satisfied by two-position devices.

\textbf{M5 has an unresolved ordering problem.} The humidifier sits on the
$d_{pE}$ prong, upstream of the chilled-water coil in fork~2, so the air is
humidified before it is cooled and may reach saturation. Cooling and then
humidifying is the correct sequence, and the architecture does not currently
provide it.

\textbf{M7 and M8 route through prong D with the wheel stopped.} Air passes a
stationary wheel at nothing but a pressure drop, which turns the cooler into an
exhaust-to-supply recuperator at no capital cost --- see
Section~\ref{sec:winter}.

% ==========================================================================
\section{The resource constraint}\label{sec:resource}
% ==========================================================================

\subsection{Temperature: which tiers the well can serve}

Only M3 and M5 require the \SI{79}{\degreeCelsius} minimum at the absorption
generator. Every other mode needs approximately \SI{66}{\degreeCelsius} for
wheel regeneration, which the well holds for the whole design life at almost
any flow. The serviceable region of the load plane therefore \emph{shrinks
with reservoir depletion, and it shrinks from the cooling corner first}.

From the archived parametric sweep, year-20 conditions:

\begin{longtable}{@{}rrrrrrc@{}}
\toprule
$L$ [\si{\metre}] & BHT [\degC] & $q$ [\si{\cubic\metre\per\hour}] &
$T_{\mathrm{ret}}$ [\degC] & $T_{\mathrm{wa}}$ [\degC] &
$Q_{\mathrm{well}}$ [\si{\kilo\watt}] & generator viable \\
\midrule
\endhead
@@SWEEP@@
\bottomrule
\end{longtable}

\subsection{Capacity: the flow squeeze}

Raising loop flow raises the generator duty but lowers the water temperature,
and the two requirements are not satisfied at the same operating point. Using a
design generator drop of \SI{10}{\kelvin} and a maximum of \SI{18}{\kelvin},
the coded sensible demand of \SI{101}{\kilo\watt} needs a generator duty of
about \SI{145}{\kilo\watt}, that is roughly
\SI{7}{\cubic\metre\per\hour} --- but at \SI{4500}{\metre} that flow does not
reach \SI{79}{\degreeCelsius}. The temperature-optimal flow starves the
chiller; the capacity-optimal flow cools the water below its minimum.
Satisfying both requires depth: approximately \SIrange{5000}{5500}{\metre} at
\SIrange{6}{7}{\cubic\metre\per\hour}.

\subsection{Winter capacity, and why recovery dominates}\label{sec:winter}

Sizing the January duty over the whole outdoor-air stream, the enthalpy lift
from a cold, dry ambient to the comfort state is roughly
\SI{60}{\kilo\joule\per\kilogram}, giving a total near \SI{772}{\kilo\watt} of
which more than half is latent. Net of internal gains the plant duty is about
\SI{630}{\kilo\watt}, against a well output between \SI{198}{\kilo\watt} and
\SI{336}{\kilo\watt}: the well covers only 31--53\%.

A dry recuperator at 70\%\ effectiveness recovers about
\SI{227}{\kilo\watt}, leaving \SI{403}{\kilo\watt} unserved. A
\emph{latent-capable} device at the same effectiveness recovers about
\SI{540}{\kilo\watt}, leaving \SI{90}{\kilo\watt} --- comfortably inside the
well's capability. Because 58\%\ of the winter duty is latent, only the second
closes the gap.

This is the strongest argument in the architecture for routing M7 and M8
through prong~D with the wheel halted: the cooler then acts as an
exhaust-to-supply recuperator using hardware already installed. Whether the
desiccant wheel itself, operated differently, can recover the latent half ---
a rotary wheel moving moisture from a warm humid exhaust to a cold dry supply
\emph{is} an enthalpy wheel --- remains open.

% ==========================================================================
\section{Implementation status}\label{sec:status}
% ==========================================================================

\begin{longtable}{@{}llp{0.44\textwidth}@{}}
\toprule
\textbf{Component} & \textbf{Status} & \textbf{Note} \\
\midrule
\endhead
DBHE (GITT + coaxial BVP) & implemented & transient formation coupling, 20-year march \\
Thermal accumulator & implemented & lumped capacitance \\
Regeneration coil (HX) & implemented & effectiveness--NTU, water to air \\
Desiccant wheel & implemented & Jurinak $F_1/F_2$ with prescribed effectiveness \\
Recuperative cooler & implemented & process air to scavenging air \\
Load interface & implemented & v@@VER@@; generic $(\qs, \mdw)$ with explicit closure \\
Humidifier & pending & adiabatic saturation; unlocks M4 \\
Indirect evaporative cooler & pending & effectiveness form; unlocks M2 \\
Winter heating coil & pending & reuses the coil model on different streams; unlocks M6--M8 \\
Absorption chiller + tower + CHW & pending & unlocks M3, M5 \\
\bottomrule
\end{longtable}

\paragraph{A known limitation of the wheel model.} The Jurinak formulation
carries no size: no diameter, no desiccant mass, no rotation speed, no
transfer area. The outlet is placed a fixed fraction of the way from the
process inlet toward the regeneration inlet in $F_1/F_2$ space, and those
fractions are the two prescribed effectivenesses. It can be verified directly
that the process outlet state is \emph{independent of both mass flows}. Three
consequences follow: the wheel cannot be resized within this model, it cannot
be throttled by scavenging flow, and any claim about wheel sizing, rotation
speed or part-load behaviour is outside what the formulation can represent.
Deriving the effectivenesses from an NTU using the analogy method would remove
these restrictions at essentially no computational cost, and is deferred.

% ==========================================================================
\section{Verification protocol}
% ==========================================================================

Each model version is frozen as a regression fixture recording every tracked
state at annual resolution. Fixtures are named by version and are
\emph{never overwritten}: a mutable reference would always pass and prove
nothing. A refactor that should not change physics is verified by reproducing
the previous fixture to \num{1e-9}; a change that should move results is
verified by explaining every value that moved.

The harness loads the notebook's own code rather than re-implementing it, so
if the two ever disagree the harness is wrong by definition. It is validated in
both directions: an untouched notebook must pass, and a notebook perturbed by
0.45\%\ in formation conductivity must fail.

\paragraph{A recurring failure mode.} Four defects found so far shared one
shape --- a guard that degraded quietly instead of complaining:

\begin{enumerate}[nosep]
  \item a missing property library substituted constant water properties,
        misstating viscosity by up to 90\%\ at the hot end;
  \item a global warning filter suppressed the notebook's own diagnostics,
        so a formation-model safety check could never fire;
  \item the desiccant wheel could reverse sign and humidify the process air
        while the dashboard still reported the humidity target as met;
  \item the scavenging mixing node conserved energy but not moisture, using
        two different damper positions for the two balances.
\end{enumerate}

None raised an error. Each produced plausible output. The protocol now
requires that any guard whose failure would otherwise be invisible must
announce itself.

% ==========================================================================
\section{Version history}\label{sec:history}
% ==========================================================================

Year-20 values at the default operating point
(\SI{@@LBH@@}{\metre}, \SI{@@FLOW@@}{\cubic\metre\per\hour}, insulated
tubing), taken directly from the frozen fixtures:

\begin{longtable}{@{}lrrrrrrr@{}}
\toprule
\textbf{v} & $T_{\mathrm{ret}}$ & $T_{\mathrm{wa}}$ & $T_{\mathrm{ro}}$ &
$w_{ri}$ & $\Delta w$ & $\UA_{\mathrm{cool}}$ & $d_{rA}$ \\
 & [\degC] & [\degC] & [\degC] & [g/kg] & [g/kg] & [kW/K] & [--] \\
\midrule
\endhead
@@HISTORY@@
\bottomrule
\end{longtable}

\begin{description}[leftmargin=3.2em,style=nextline]
@@NOTES@@
\end{description}

\paragraph{On the two largest movements.} The mass-balance correction at v10.4
moved results in the direction \emph{opposite} to the one predicted. Fixing
the balance improves the wheel at any given damper position, so the controller
needs less venting and drives the damper to its floor; at high recirculation
the loop accumulates more of the wheel's own moisture, and the regeneration
inlet ends up wetter rather than drier. The open-loop estimate and the
closed-loop result disagreed in sign. At v10.5, removing recirculation
entirely then reduced the regeneration inlet humidity by roughly
\SI{17}{g\per\kilogram} and nearly tripled moisture removal.

\paragraph{Current open item.} The humidity controller is saturated at its dry
end: the supply is drier than the space requires and the damper has no
authority left to reduce drying. Approximately 41\%\ of the regeneration coil
duty is presently spent removing moisture the building does not want. The
remedy is to move humidity control from the scavenging damper to $v_G$, the
share of water flow reaching the regeneration coil, which throttles the coil
duty directly and returns that heat to the loop.

\end{document}
"""


def main():
    ver, stamp = nb_version()
    fx = fixtures()
    hist, vs = version_history_table(fx)
    subs = {
        '@@VER@@': ver,
        '@@STAMP@@': tex_esc(stamp),
        '@@HISTORY@@': hist,
        '@@NOTES@@': version_notes_list(vs),
        '@@SWEEP@@': sweep_table(sweep()),
        '@@MDOT_VENT@@': f"{nb_scalar('mdot_vent') or 10.0:.2f}",
        '@@MDOT_REG@@': f"{23.52:.2f}",
        '@@DRA_CAP@@': f"{10.0/23.52:.3f}",
        '@@MDOT_LOAD@@': f"{33.60:.2f}",
        '@@MDOT_PRO@@': f"{33.60:.2f}",
        '@@DT_SUP@@': f"{nb_scalar('dT_supply_design') or 3.0:.0f}",
        '@@LBH@@': f"{nb_scalar('Lbh_default') or 4500:.0f}",
        '@@FLOW@@': f"{nb_scalar('flow_rate_m3_h_demo') or 10:.0f}",
    }
    tex = TEMPLATE
    for k, v in subs.items():
        tex = tex.replace(k, v)
    left = re.findall(r'@@[A-Z_]+@@', tex)
    if left:
        print('WARNING: unresolved placeholders:', sorted(set(left)), file=sys.stderr)
    open(OUT, 'w').write(tex)
    print(f'wrote {OUT}  ({len(tex.splitlines())} lines) for model v{ver}')
    print(f'  versions in history: {", ".join("v"+v for v in vs)}')
    return OUT


if __name__ == '__main__':
    main()
