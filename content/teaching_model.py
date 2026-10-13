"""EVSC 410: uniform-width 1-D normal-flow Exner teaching model.
Educational prototype, not a validated engineering or restoration design tool.
SI units throughout. Requires numpy.
"""
from dataclasses import dataclass, replace
import numpy as np

G = 9.81
R = 1.65  # submerged specific density (2650/1000 - 1)

@dataclass(frozen=True)
class Parameters:
    length_m: float = 300.0
    dx_m: float = 5.0
    width_m: float = 8.0  # spatially uniform in each run
    initial_slope: float = 0.008
    discharge_m3s: float = 8.0
    grain_m: float = 0.025
    roughness_ks_m: float = 0.075
    porosity: float = 0.35
    critical_shields: float = 0.045
    feed_ratio: float = 1.0
    duration_s: float = 1800.0
    max_dt_s: float = 2.0
    save_count: int = 7


def check(p):
    for k in ('length_m','dx_m','width_m','initial_slope','discharge_m3s',
              'grain_m','roughness_ks_m','max_dt_s'):
        if getattr(p,k) <= 0 or not np.isfinite(getattr(p,k)):
            raise ValueError(f'{k} must be positive and finite')
    if not (np.isfinite(p.duration_s) and np.isfinite(p.feed_ratio) and
            0 <= p.porosity < 1 and 0 <= p.critical_shields < 1 and
            p.feed_ratio >= 0 and p.duration_s >= 0 and p.save_count >= 2):
        raise ValueError('Invalid porosity, Shields threshold, feed, duration or save_count')
    n = p.length_m / p.dx_m
    if abs(n-round(n)) > 1e-9 or n < 3:
        raise ValueError('length_m must be an integer multiple of dx_m with at least 3 cells')


def hydraulic_transport(slope, p):
    """Face hydraulics. Based on Normal_Rev1 Strickler wide-channel depth.
    Shear uses rho*g*h*S (wide-channel consistent). Returns depth, Shields, q_s.
    """
    if np.any(slope <= 0):
        raise ValueError('Normal-flow assumption violated: zero or adverse local bed slope')
    qw = p.discharge_m3s / p.width_m
    depth = (p.roughness_ks_m**0.33 * qw**2 / (8.1**2 * G * slope))**0.30
    shields = depth * slope / (R * p.grain_m)
    qs = 3.97 * np.maximum(shields-p.critical_shields, 0)**1.5 * np.sqrt(R*G*p.grain_m)*p.grain_m
    return depth, shields, qs


def transport_slope_sensitivity(slope, p):
    """Analytical dq_s/dS for the existing Strickler + excess-Shields closure.

    With h proportional to S**(-0.30), theta proportional to S**0.70.
    Returns nonnegative sensitivity, with zero below threshold.
    """
    _, theta, _ = hydraulic_transport(slope, p)
    factor = 3.97 * np.sqrt(R*G*p.grain_m)*p.grain_m
    return factor * 1.5 * np.sqrt(np.maximum(theta-p.critical_shields, 0)) * (0.70*theta/slope)


def fluxes(eta, p, feed):
    # N finite-volume cells; N+1 faces. Interior face slopes are centered
    # between adjacent cell-center elevations. Downstream outlet face uses
    # one-sided extrapolation of the last interior slope (open gradient BC).
    interior = (eta[:-1]-eta[1:])/p.dx_m
    slope = np.r_[interior[0], interior, interior[-1]]
    depth, shields, capacity = hydraulic_transport(slope,p)
    q = capacity.copy()
    q[0] = feed  # prescribed upstream sediment feed, m2/s
    return q, depth, shields, slope



def preflight(p=Parameters(), safety=0.20):
    """Initial-state diagnostic; not a guarantee of stability later in a run.

    Hydraulic Courant U*dt/dx is informational for this quasi-steady model.
    Explicit Exner diffusive stability uses max dq_s/dS/(1-porosity).
    """
    check(p)
    slope = np.array([p.initial_slope], dtype=float)
    depth, shields, qs = hydraulic_transport(slope, p)
    h = float(depth[0]); u = p.discharge_m3s/(p.width_m*h)
    dqds = float(transport_slope_sensitivity(slope,p)[0])
    diffusivity = dqds/(1-p.porosity)
    dt_limit = safety*p.dx_m**2/diffusivity if diffusivity>0 else float('inf')
    dt_selected = min(p.max_dt_s,dt_limit)
    return dict(depth_m=h, velocity_ms=u, shields=float(shields[0]),
                capacity_m2s=float(qs[0]), hydraulic_courant=u*p.max_dt_s/p.dx_m,
                morphodynamic_diffusivity_m2s=diffusivity,
                morphodynamic_number=diffusivity*p.max_dt_s/p.dx_m**2,
                max_morphodynamic_dt_s=dt_limit, planned_dt_s=dt_selected,
                planned_morphodynamic_number=diffusivity*dt_selected/p.dx_m**2,
                dt_will_be_limited=bool(p.max_dt_s>dt_limit), safety=safety)

def simulate(p=Parameters(), max_cells=500, max_steps=150000):
    check(p)
    n = int(round(p.length_m/p.dx_m))
    if n > max_cells:
        raise ValueError(f"Requested {n} cells exceeds classroom limit {max_cells}")
    if p.save_count > 60:
        raise ValueError("save_count exceeds classroom limit 60")
    x = (np.arange(n)+0.5)*p.dx_m
    eta = p.initial_slope*(p.length_m-x)
    _, _, equilibrium_q = hydraulic_transport(np.array([p.initial_slope]),p)
    feed = float(p.feed_ratio*equilibrium_q[0])
    t = 0.0
    times = [t]
    profiles = [eta.copy()]
    target_times = np.linspace(0,p.duration_s,p.save_count)[1:]
    net_sediment_volume = 0.0
    steps = 0
    # Adaptive rejection prevents stepping across the positive-slope domain.
    # A conservative diffusion-like limit is also imposed by estimating the
    # sensitivity of sediment capacity to local slope.
    for target in target_times:
        while t < target - 1e-10:
            q,depth,shields,slope = fluxes(eta,p,feed)
            dqds = transport_slope_sensitivity(slope[1:-1],p)
            diffusivity = np.max(dqds)/(1-p.porosity)
            dt_stable = (0.20*p.dx_m**2/diffusivity if diffusivity > 0 else np.inf)
            dt = min(p.max_dt_s,dt_stable,target-t)
            if dt <= 1e-10: raise RuntimeError('Time step collapsed; try a coarser scenario')
            for _ in range(35):
                trial = eta - dt/(p.dx_m*(1-p.porosity))*(q[1:]-q[:-1])
                if np.all(np.isfinite(trial)) and np.all(trial[:-1]>trial[1:]):
                    break
                dt *= 0.5
            else:
                raise RuntimeError('Cannot maintain positive bed slope: scenario outside normal-flow model validity')
            eta = trial
            net_sediment_volume += p.width_m*(q[0]-q[-1])*dt
            t += dt
            steps += 1
            if steps > max_steps: raise RuntimeError('Exceeded classroom step limit')
        times.append(t)
        profiles.append(eta.copy())
    q,depth,shields,slope = fluxes(eta,p,feed)
    storage = (1-p.porosity)*p.width_m*p.dx_m*np.sum(eta-profiles[0])
    return {'x_m':x,'times_s':np.array(times),'bed_m':np.array(profiles),
            'depth_faces_m':depth,'shields_faces':shields,'slope_faces':slope,
            'transport_faces_m2s':q,'feed_m2s':feed,'steps':steps,
            'storage_change_m3':storage,'net_boundary_input_m3':net_sediment_volume,
            'mass_balance_error_m3':storage-net_sediment_volume}


def classroom_estimate(p=Parameters()):
    """Initial-condition estimate only; actual steps may differ as slopes evolve."""
    d = preflight(p)
    n = int(round(p.length_m/p.dx_m))
    estimated_steps = int(np.ceil(p.duration_s / d['planned_dt_s'])) if p.duration_s else 0
    stored_kib = p.save_count*n*8/1024
    return dict(cells=n, estimated_steps=estimated_steps, stored_profiles_kib=stored_kib,
                caution=(n>250 or estimated_steps>30000 or stored_kib>500), **d)
