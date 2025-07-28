```python
"""
EIDOS: Emergent Information Dynamics of Observational Spacetime

A parameter-free framework for explaining cosmological anomalies through
staged quantum decoherence.
"""

__version__ = "1.0.0"
__author__ = "EIDOS Collaboration"
__email__ = "contact@eidos-cosmology.org"

# Core imports
from .core import EIDOSModel, entropy_field, power_spectrum_modification
from .analysis import (
    analyze_cmb,
    analyze_lss,
    compute_chi_squared,
    parameter_forecast
)
from .visualization import (
    plot_power_spectrum,
    plot_cmb_spectra,
    plot_parameter_constraints,
    make_triangle_plot
)
from .utils import (
    download_data,
    load_planck_data,
    load_boss_data,
    cosmology_calculator
)

# Convenience imports
from .core import EIDOSModel as Model

__all__ = [
    # Core
    "EIDOSModel",
    "Model",
    "entropy_field",
    "power_spectrum_modification",
    
    # Analysis
    "analyze_cmb",
    "analyze_lss",
    "compute_chi_squared",
    "parameter_forecast",
    
    # Visualization
    "plot_power_spectrum",
    "plot_cmb_spectra",
    "plot_parameter_constraints",
    "make_triangle_plot",
    
    # Utils
    "download_data",
    "load_planck_data",
    "load_boss_data",
    "cosmology_calculator",
]