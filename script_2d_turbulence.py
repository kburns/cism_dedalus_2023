
import numpy as np
np.seterr(over="raise")
import dedalus.public as d3
import logging
logger = logging.getLogger(__name__)


# Domain parameters
L = 2 * np.pi
N = 4096
mesh = None
dtype = np.float64
dealias = 3 / 2

# Forcing parameters
epsilon = 1     # Energy injection rate
kf = 50         # Forcing wavenumber
kfw = 2         # Forcing bandwidth
seed = None     # Random seed

# Problem parameters
L_diss = L / N  # Dissipation scale
L_fric = L      # Friction scale

# Solver parameters
timestepper = d3.RK443
stop_sim_time = 5

# Analysis parameters
snapshots_dt = 0.1
scalars_dt = 0.01

# Timestepping parameters
dx = L / N                          # Grid spacing
U = epsilon**(1/3) * L_fric**(1/3)  # Friction velocity
safety = 0.5                        # CFL safety factor
max_dt = safety * dx / U            # Timestep

# Derived parameters
eta = epsilon * kf**2  # Enstrophy injection rate
nu = L_diss**2 * eta**(1/3)                 # Viscosity
alpha = epsilon**(1/3) * L_fric**(-2/3)     # Friction

# Domain
coords = d3.CartesianCoordinates('x', 'y')
dist = d3.Distributor(coords, mesh=mesh, dtype=dtype)
xbasis = d3.RealFourier(coords[0], N, bounds=(0, L), dealias=dealias)
ybasis = d3.RealFourier(coords[1], N, bounds=(0, L), dealias=dealias)

# Fields
psi = dist.Field(name='psi', bases=(xbasis, ybasis))
c = dist.Field(name='c')

# Substitutions
u = -d3.skew(d3.grad(psi))  # velocity vector: [dy(psi), -dx(psi)]
w = -d3.lap(psi)            # vorticity: dx(uy) - dy(ux)
e = (u@u) / 2               # energy density
z = (w*w) / 2               # enstrophy density

# Forcing field and derived parameters
Fw = dist.Field(name='Fw', bases=(xbasis, ybasis))
kx = xbasis.wavenumbers[dist.local_modes(xbasis)]
ky = ybasis.wavenumbers[dist.local_modes(ybasis)]
dkx = dky = 2 * np.pi / L

# Forcing function
rand = np.random.RandomState(seed)

def draw_gaussian_random_field():
    """Create Gaussian random field concentrating on a ring in Fourier space with unit variance."""
    k = (kx**2 + ky**2)**0.5
    # 1D power spectrum: normalized Gaussian, no mean
    P1 = np.exp(-(k-kf)**2/2/kfw**2) / np.sqrt(kfw**2 * np.pi / 2) * (k != 0)
    # 2D power spectrum: divide by polar Jacobian
    P2 = P1 / 2 / np.pi / (k + (k==0))
    # 2D coefficient poewr spectrum: divide by mode power
    Pc = P2 / 2**((kx == 0).astype(float) + (ky == 0).astype(float) - 2)
    # Forcing amplitude, including division between sine and cosine
    f_amp = (Pc / 2 * dkx * dkx)**0.5
    # Forcing with random phase
    f = f_amp * rand.randn(*k.shape)
    return f

def set_vorticity_forcing(timestep):
    """Set vorticity forcing field from scaled Gaussian random field."""
    # Set forcing to normalized Gaussian random field
    Fw['c'] = draw_gaussian_random_field()
    # Rescale by forcing rate, including factor for 1/2 in kinetic energy
    Fw['c'] *= (2 * eta / timestep)**0.5

# Problem
problem = d3.IVP([psi, c], namespace=locals())
problem.add_equation("dt(w) - nu*lap(w) + alpha*w + c = -u@grad(w) + Fw")
problem.add_equation("integ(psi) = 0");

# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=snapshots_dt, max_writes=10, mode='overwrite')
snapshots.add_task(psi, name='psi')
snapshots.add_task(w, name='vorticity')

scalars = solver.evaluator.add_file_handler('scalars', sim_dt=scalars_dt, mode='overwrite')
ave = d3.Average
scalars.add_task(ave(e), name='E')
scalars.add_task(ave(z), name='Z')
scalars.add_task(ave(-alpha * 2 * e), name='E friction')
scalars.add_task(ave(-alpha * 2 * z), name='Z friction')
scalars.add_task(ave(nu * u @ d3.lap(u)), name='E viscosity')
scalars.add_task(ave(nu * w * d3.lap(w)), name='Z viscosity')

# CFL
CFL = d3.CFL(solver, initial_dt=max_dt, cadence=10, safety=safety, max_change=1.5, min_change=0.5, max_dt=max_dt, threshold=0.05)
CFL.add_velocity(u)

# Main loop
try:
    logger.info('Starting loop')
    while solver.proceed:
        timestep = CFL.compute_timestep()
        set_vorticity_forcing(timestep)
        solver.step(timestep)
        if (solver.iteration-1) % 10 == 0:
            logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()
