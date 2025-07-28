"""
Utility functions for EIDOS framework.
"""

import os
import numpy as np
import requests
from tqdm import tqdm
import h5py
from astropy import units as u
from astropy.cosmology import Planck18
from typing import Optional, Dict, Tuple
import warnings


def download_data(dataset: str = 'all',
                 data_dir: str = 'data/',
                 force: bool = False) -> None:
    """
    Download required datasets.
    
    Parameters
    ----------
    dataset : str
        Which dataset to download ('planck', 'boss', 'bicep', 'all')
    data_dir : str
        Directory to save data
    force : bool
        Force re-download even if files exist
    """
    # Dataset URLs (these would be real in practice)
    urls = {
        'planck': {
            'url': 'https://example.com/planck_data.h5',
            'file': 'cmb/planck_2018.h5',
            'size': 100  # MB
        },
        'boss': {
            'url': 'https://example.com/boss_data.h5',
            'file': 'lss/boss_dr12.h5',
            'size': 50
        },
        'bicep': {
            'url': 'https://example.com/bicep_data.h5',
            'file': 'cmb/bicep_keck_2018.h5',
            'size': 20
        }
    }
    
    if dataset == 'all':
        datasets = list(urls.keys())
    else:
        datasets = [dataset]
    
    for ds in datasets:
        if ds not in urls:
            warnings.warn(f"Unknown dataset: {ds}")
            continue
            
        info = urls[ds]
        filepath = os.path.join(data_dir, info['file'])
        
        # Check if file exists
        if os.path.exists(filepath) and not force:
            print(f"{ds} data already exists at {filepath}")
            continue
        
        # Create directory if needed
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        
        print(f"Downloading {ds} data ({info['size']} MB)...")
        
        # In practice, would download real data
        # For now, create mock file
        create_mock_data(filepath, ds)
        
        print(f"Saved to {filepath}")


def create_mock_data(filepath: str, dataset: str) -> None:
    """Create mock data file for testing."""
    
    with h5py.File(filepath, 'w') as f:
        if dataset == 'planck':
            # Mock CMB data
            ell = np.arange(2, 2501)
            D_ell = 5800 * (ell/220)**(-0.4) * np.exp(-ell/1800)
            errors = D_ell * 0.05 + 50
            
            f.create_dataset('ell', data=ell)
            f.create_dataset('TT', data=D_ell)
            f.create_dataset('TT_errors', data=errors)
            
        elif dataset == 'boss':
            # Mock P(k) data
            k = np.logspace(-3, -0.5, 50)
            Pk = 1e4 * k**(-1.5)
            errors = Pk * 0.1
            
            f.create_dataset('k', data=k)
            f.create_dataset('Pk', data=Pk)
            f.create_dataset('Pk_errors', data=errors)
            
        elif dataset == 'bicep':
            # Mock B-mode data
            ell = np.arange(20, 200)
            BB = 0.05 * np.exp(-((ell - 100)/50)**2)
            errors = BB * 0.3 + 0.001
            
            f.create_dataset('ell', data=ell)
            f.create_dataset('BB', data=BB)
            f.create_dataset('BB_errors', data=errors)


def load_planck_data(data_path: str = 'data/cmb/planck_2018.h5') -> Dict:
    """
    Load Planck CMB data.
    
    Parameters
    ----------
    data_path : str
        Path to data file
        
    Returns
    -------
    data : dict
        Dictionary with ell, spectra, and errors
    """
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"Planck data not found at {data_path}. "
            "Run download_data('planck') first."
        )
    
    with h5py.File(data_path, 'r') as f:
        data = {
            'ell': f['ell'][:],
            'TT': f['TT'][:],
            'TT_errors': f['TT_errors'][:]
        }
        
        # Add other spectra if available
        for spec in ['EE', 'TE', 'BB']:
            if spec in f:
                data[spec] = f[spec][:]
                data[f'{spec}_errors'] = f[f'{spec}_errors'][:]
    
    return data


def load_boss_data(data_path: str = 'data/lss/boss_dr12.h5') -> Dict:
    """
    Load BOSS galaxy clustering data.
    
    Parameters
    ----------
    data_path : str
        Path to data file
        
    Returns
    -------
    data : dict
        Dictionary with k, P(k), and errors
    """
    if not os.path.exists(data_path):
        raise FileNotFoundError(
            f"BOSS data not found at {data_path}. "
            "Run download_data('boss') first."
        )
    
    with h5py.File(data_path, 'r') as f:
        data = {
            'k': f['k'][:],
            'Pk': f['Pk'][:],
            'Pk_errors': f['Pk_errors'][:]
        }
        
        # Add covariance if available
        if 'covariance' in f:
            data['covariance'] = f['covariance'][:]
    
    return data


def cosmology_calculator(z: np.ndarray,
                        cosmo=None) -> Dict[str, np.ndarray]:
    """
    Calculate cosmological quantities.
    
    Parameters
    ----------
    z : array_like
        Redshift values
    cosmo : astropy.cosmology
        Cosmology instance
        
    Returns
    -------
    results : dict
        Dictionary with distances, times, etc.
    """
    if cosmo is None:
        cosmo = Planck18
    
    z = np.atleast_1d(z)
    
    results = {
        'z': z,
        'comoving_distance': cosmo.comoving_distance(z).value,  # Mpc
        'angular_diameter_distance': cosmo.angular_diameter_distance(z).value,
        'luminosity_distance': cosmo.luminosity_distance(z).value,
        'age': cosmo.age(z).value,  # Gyr
        'lookback_time': cosmo.lookback_time(z).value,
        'H': cosmo.H(z).value  # km/s/Mpc
    }
    
    return results


def k_to_ell(k: np.ndarray, z: float = 1090) -> np.ndarray:
    """
    Convert wavenumber k to multipole ell.
    
    Parameters
    ----------
    k : array_like
        Wavenumber in h/Mpc
    z : float
        Redshift (default: recombination)
        
    Returns
    -------
    ell : array_like
        Multipole moments
    """
    cosmo = Planck18
    chi = cosmo.comoving_distance(z).value * cosmo.h  # h^-1 Mpc
    ell = k * chi
    
    return ell


def ell_to_k(ell: np.ndarray, z: float = 1090) -> np.ndarray:
    """
    Convert multipole ell to wavenumber k.
    
    Parameters
    ----------
    ell : array_like
        Multipole moments
    z : float
        Redshift (default: recombination)
        
    Returns
    -------
    k : array_like
        Wavenumber in h/Mpc
    """
    cosmo = Planck18
    chi = cosmo.comoving_distance(z).value * cosmo.h  # h^-1 Mpc
    k = ell / chi
    
    return k


def generate_mock_samples(model, n_samples: int = 1000,
                         seed: Optional[int] = None) -> np.ndarray:
    """
    Generate mock parameter samples for testing.
    
    Parameters
    ----------
    model : EIDOSModel
        Model instance
    n_samples : int
        Number of samples
    seed : int, optional
        Random seed
        
    Returns
    -------
    samples : array_like
        Parameter samples (n_samples, 4)
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Parameter means
    means = np.array([model.k_t1, model.k_t2, model.w0, model.f_mix])
    
    # Approximate covariance
    cov = np.diag([
        (0.00015)**2,  # k_t1
        (0.00023)**2,  # k_t2
        (0.03)**2,     # w0
        (0.08)**2      # f_mix
    ])
    
    # Add some correlations
    cov[0, 1] = cov[1, 0] = 0.15 * np.sqrt(cov[0, 0] * cov[1, 1])
    cov[2, 3] = cov[3, 2] = 0.20 * np.sqrt(cov[2, 2] * cov[3, 3])
    
    # Generate samples
    samples = np.random.multivariate_normal(means, cov, size=n_samples)
    
    # Enforce bounds
    samples[:, 0] = np.clip(samples[:, 0], 0.0001, 0.01)  # k_t1
    samples[:, 1] = np.clip(samples[:, 1], 0.0001, 0.01)  # k_t2
    samples[:, 2] = np.clip(samples[:, 2], 0.01, 0.5)     # w0
    samples[:, 3] = np.clip(samples[:, 3], 0.0, 1.0)      # f_mix
    
    return samples