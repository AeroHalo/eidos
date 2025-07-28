"""
Visualization tools for EIDOS framework.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns
from typing import Optional, List, Tuple, Dict
import corner
from .core import EIDOSModel


def setup_plotting():
    """Set up matplotlib parameters for publication-quality plots."""
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 11
    rcParams['axes.labelsize'] = 12
    rcParams['axes.titlesize'] = 14
    rcParams['xtick.labelsize'] = 10
    rcParams['ytick.labelsize'] = 10
    rcParams['legend.fontsize'] = 10
    rcParams['figure.figsize'] = (8, 6)
    rcParams['figure.dpi'] = 100
    rcParams['savefig.dpi'] = 300
    rcParams['savefig.bbox'] = 'tight'
    rcParams['axes.grid'] = True
    rcParams['grid.alpha'] = 0.3
    

def plot_power_spectrum(model: Optional[EIDOSModel] = None,
                       k_range: Tuple[float, float] = (1e-4, 1),
                       show_components: bool = True,
                       save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot power spectrum modifications.
    
    Parameters
    ----------
    model : EIDOSModel, optional
        Model instance
    k_range : tuple
        Range of k values to plot
    show_components : bool
        Whether to show individual components
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.Figure
        Figure object
    """
    if model is None:
        model = EIDOSModel()
        
    setup_plotting()
    
    # Create k array
    k = np.logspace(np.log10(k_range[0]), np.log10(k_range[1]), 1000)
    
    # Get modifications
    S_mixed = model.mixed_entropy(k)
    P_mod = S_mixed**2
    
    # Create figure
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
    
    # Top panel: Entropy field
    ax1.semilogx(k, S_mixed, 'k-', lw=2, label='Mixed S(k)')
    
    if show_components:
        S_g = model.entropy_field(k, model.k_t1)
        S_m = model.entropy_field(k, model.k_t2)
        ax1.semilogx(k, S_g, 'b--', lw=1.5, alpha=0.7, label='Gravity')
        ax1.semilogx(k, S_m, 'r--', lw=1.5, alpha=0.7, label='Matter')
    
    # Mark transitions
    ax1.axvline(model.k_t1, color='b', linestyle=':', alpha=0.5)
    ax1.axvline(model.k_t2, color='r', linestyle=':', alpha=0.5)
    
    ax1.set_ylabel('Entropy Field S(k)')
    ax1.set_ylim(-0.05, 1.05)
    ax1.legend(loc='upper left')
    ax1.set_title('EIDOS Entropy Field Evolution')
    
    # Bottom panel: Power spectrum modification
    ax2.semilogx(k, P_mod, 'k-', lw=2)
    ax2.fill_between(k, P_mod, 1, alpha=0.3, color='gray')
    
    # Add observational scales
    ax2.axvspan(1e-4, 1e-3, alpha=0.2, color='orange', label='CMB')
    ax2.axvspan(1e-3, 1e-2, alpha=0.2, color='green', label='LSS')
    ax2.axvspan(1e-2, 1e-1, alpha=0.2, color='blue', label='Galaxy')
    
    ax2.axhline(1, color='gray', linestyle='--', alpha=0.5)
    ax2.set_xlabel('k [h/Mpc]')
    ax2.set_ylabel('P(k) / P$_{\\Lambda\\rm CDM}$(k)')
    ax2.set_ylim(0, 1.1)
    ax2.legend(loc='lower right')
    ax2.set_title('Power Spectrum Modification')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        
    return fig


def plot_cmb_spectra(model: Optional[EIDOSModel] = None,
                    ell_max: int = 2500,
                    show_data: bool = True,
                    save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot CMB angular power spectra.
    
    Parameters
    ----------
    model : EIDOSModel, optional
        Model instance
    ell_max : int
        Maximum multipole
    show_data : bool
        Whether to show mock data points
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.Figure
        Figure object
    """
    if model is None:
        model = EIDOSModel()
        
    setup_plotting()
    
    ell = np.arange(2, ell_max + 1)
    mods = model.cmb_cl_modifications(ell)
    
    # Create mock ΛCDM spectrum (simplified)
    D_ell_lcdm = 5800 * (ell/220)**(-0.4) * np.exp(-ell/1800)
    D_ell_lcdm[ell < 30] = 1200 * np.exp(-((ell[ell < 30]-10)/20)**2)
    
    # Apply modifications
    D_ell_TT = D_ell_lcdm * mods['TT']
    
    # Create figure with subplots
    fig = plt.figure(figsize=(12, 10))
    
    # TT spectrum
    ax1 = plt.subplot(2, 2, 1)
    ax1.loglog(ell, D_ell_lcdm, 'gray', lw=2, alpha=0.5, label='ΛCDM')
    ax1.loglog(ell, D_ell_TT, 'b-', lw=2, label='EIDOS')
    
    if show_data:
        # Add mock error bars
        ell_data = np.logspace(np.log10(2), np.log10(2000), 30)
        ell_data = ell_data.astype(int)
        D_data = D_ell_TT[ell_data - 2]
        errors = D_data * 0.05 + 50
        ax1.errorbar(ell_data, D_data, yerr=errors, fmt='ko', 
                    markersize=4, alpha=0.5, label='Mock data')
    
    ax1.set_xlabel('$\\ell$')
    ax1.set_ylabel('$\\ell(\\ell+1)C_\\ell^{TT}/2\\pi$ [$\\mu$K$^2$]')
    ax1.legend()
    ax1.set_title('Temperature Power Spectrum')
    ax1.set_xlim(2, ell_max)
    
    # Low-ell zoom
    ax2 = plt.subplot(2, 2, 2)
    ell_low = ell[ell < 50]
    ax2.plot(ell_low, D_ell_lcdm[ell < 50], 'gray', lw=2, alpha=0.5)
    ax2.plot(ell_low, D_ell_TT[ell < 50], 'b-', lw=2)
    
    # Highlight quadrupole
    ax2.scatter([2], [D_ell_TT[0]], color='red', s=100, zorder=5)
    ax2.annotate(f'{(1-mods["TT"][0])*100:.0f}% suppression', 
                xy=(2, D_ell_TT[0]), xytext=(5, 800),
                arrowprops=dict(arrowstyle='->', color='red'))
    
    ax2.set_xlabel('$\\ell$')
    ax2.set_ylabel('$\\ell(\\ell+1)C_\\ell^{TT}/2\\pi$ [$\\mu$K$^2$]')
    ax2.set_title('Low-$\\ell$ Suppression')
    ax2.set_xlim(2, 50)
    
    # BB spectrum
    ax3 = plt.subplot(2, 2, 3)
    # Mock BB spectrum
    D_ell_BB_lcdm = 0.05 * np.exp(-((ell - 80)/50)**2)
    D_ell_BB = D_ell_BB_lcdm * mods['BB']
    
    ax3.semilogy(ell, D_ell_BB_lcdm, 'gray', lw=2, alpha=0.5, label='ΛCDM')
    ax3.semilogy(ell, D_ell_BB, 'r-', lw=2, label='EIDOS')
    
    # Highlight enhancement region
    mask = (ell > 80) & (ell < 120)
    ax3.fill_between(ell[mask], D_ell_BB_lcdm[mask], D_ell_BB[mask],
                    alpha=0.3, color='red')
    
    ax3.set_xlabel('$\\ell$')
    ax3.set_ylabel('$\\ell(\\ell+1)C_\\ell^{BB}/2\\pi$ [$\\mu$K$^2$]')
    ax3.legend()
    ax3.set_title('B-mode Enhancement')
    ax3.set_xlim(2, 300)
    
    # TE phase
    ax4 = plt.subplot(2, 2, 4)
    ax4.plot(ell, mods['TE_phase'], 'g-', lw=2)
    ax4.axhline(0, color='gray', linestyle='--', alpha=0.5)
    
    ax4.set_xlabel('$\\ell$')
    ax4.set_ylabel('TE Phase Shift')
    ax4.set_title('TE Correlation Phase')
    ax4.set_xlim(2, 1000)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path)
        
    return fig


def plot_parameter_constraints(model: Optional[EIDOSModel] = None,
                             include_planck: bool = True,
                             include_sh0es: bool = True,
                             save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot parameter constraints in H0-S8 plane.
    
    Parameters
    ----------
    model : EIDOSModel, optional
        Model instance
    include_planck : bool
        Show Planck contours
    include_sh0es : bool
        Show SH0ES measurement
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.Figure
        Figure object
    """
    if model is None:
        model = EIDOSModel()
        
    setup_plotting()
    
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # EIDOS prediction
    params = model.predict_parameters()
    ax.scatter(params['H0'], params['S8'], color='red', s=200, 
              marker='*', zorder=5, label='EIDOS')
    
    # Error ellipse for EIDOS
    from matplotlib.patches import Ellipse
    ellipse = Ellipse((params['H0'], params['S8']), 
                     width=2*params['H0_err'], 
                     height=2*params['S8_err'],
                     angle=0, facecolor='red', alpha=0.3)
    ax.add_patch(ellipse)
    
    if include_planck:
        # Planck constraints (simplified)
        H0_planck = 67.4
        S8_planck = 0.834
        ax.scatter(H0_planck, S8_planck, color='blue', s=100, 
                  marker='o', label='Planck 2018')
        
        # Planck contours
        theta = np.linspace(0, 2*np.pi, 100)
        for sigma in [1, 2]:
            H0_contour = H0_planck + sigma * 0.5 * np.cos(theta)
            S8_contour = S8_planck + sigma * 0.02 * np.sin(theta)
            ax.plot(H0_contour, S8_contour, 'b-', alpha=0.5-0.2*sigma)
    
    if include_sh0es:
        # SH0ES
        H0_shoes = 73.2
        ax.axvline(H0_shoes, color='green', linestyle='--', alpha=0.5)
        ax.fill_betweenx([0.7, 0.9], H0_shoes - 1.3, H0_shoes + 1.3,
                        alpha=0.2, color='green', label='SH0ES')
    
    # Weak lensing
    S8_lensing = 0.759
    ax.axhline(S8_lensing, color='orange', linestyle='--', alpha=0.5)
    ax.fill_between([64, 76], S8_lensing - 0.025, S8_lensing + 0.025,
                   alpha=0.2, color='orange', label='Weak Lensing')
    
    ax.set_xlabel('$H_0$ [km/s/Mpc]')
    ax.set_ylabel('$S_8 = \\sigma_8\\sqrt{\\Omega_m/0.3}$')
    ax.set_xlim(64, 76)
    ax.set_ylim(0.7, 0.9)
    ax.legend(loc='upper right')
    ax.set_title('EIDOS Resolves Cosmological Tensions')
    
    if save_path:
        plt.savefig(save_path)
        
    return fig


def make_triangle_plot(samples: np.ndarray,
                      labels: List[str],
                      truths: Optional[List[float]] = None,
                      save_path: Optional[str] = None) -> plt.Figure:
    """
    Make corner plot for parameter constraints.
    
    Parameters
    ----------
    samples : array_like
        MCMC samples (n_samples, n_params)
    labels : list
        Parameter labels
    truths : list, optional
        True parameter values
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.Figure
        Figure object
    """
    setup_plotting()
    
    fig = corner.corner(samples, labels=labels, truths=truths,
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True, title_kwargs={"fontsize": 12})
    
    if save_path:
        plt.savefig(save_path)
        
    return fig


def plot_future_forecasts(experiments: List[str] = ['euclid', 'cmb-s4', 'desi'],
                         model: Optional[EIDOSModel] = None,
                         save_path: Optional[str] = None) -> plt.Figure:
    """
    Plot detection forecasts for future experiments.
    
    Parameters
    ----------
    experiments : list
        List of experiments to include
    model : EIDOSModel, optional
        Model instance
    save_path : str, optional
        Path to save figure
        
    Returns
    -------
    fig : matplotlib.Figure
        Figure object
    """
    if model is None:
        model = EIDOSModel()
        
    setup_plotting()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    exp_names = []
    significances = []
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    
    for i, exp in enumerate(experiments):
        forecast = model.forecast_detection(exp)
        exp_names.append(exp.upper())
        significances.append(forecast.get('detection_sigma', 0))
    
    # Bar plot
    bars = ax.bar(exp_names, significances, color=colors[:len(experiments)])
    
    # Add significance levels
    ax.axhline(3, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(5, color='gray', linestyle='--', alpha=0.5)
    ax.text(len(experiments)-0.5, 3.1, '3σ', fontsize=10)
    ax.text(len(experiments)-0.5, 5.1, '5σ', fontsize=10)
    
    # Add values on bars
    for bar, sig in zip(bars, significances):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.1,
               f'{sig:.1f}σ', ha='center', va='bottom')
    
    ax.set_ylabel('Detection Significance')
    ax.set_title('EIDOS Detection Forecasts for Future Experiments')
    ax.set_ylim(0, max(significances) * 1.2)
    
    if save_path:
        plt.savefig(save_path)
        
    return fig