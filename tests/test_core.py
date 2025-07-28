"""
Unit tests for EIDOS core functionality.
"""

import pytest
import numpy as np
from eidos.core import EIDOSModel, entropy_field, power_spectrum_modification


class TestEIDOSModel:
    """Test the main EIDOS model class."""
    
    def test_initialization(self):
        """Test model initialization."""
        model = EIDOSModel()
        
        # Check parameters are set
        assert hasattr(model, 'k_t1')
        assert hasattr(model, 'k_t2')
        assert hasattr(model, 'f_mix')
        assert hasattr(model, 'w0')
        
        # Check parameter values
        assert 0.0005 < model.k_t1 < 0.002
        assert 0.001 < model.k_t2 < 0.005
        assert 0.5 < model.f_mix < 0.9
        assert 0.1 < model.w0 < 0.3
        
    def test_entropy_field(self):
        """Test entropy field computation."""
        model = EIDOSModel()
        
        # Test at transition
        k = np.array([model.k_t1])
        S = model.entropy_field(k, model.k_t1)
        assert np.abs(S[0] - 0.5) < 0.1  # Should be ~0.5 at transition
        
        # Test limits
        k_low = np.array([1e-5])
        k_high = np.array([1.0])
        
        S_low = model.entropy_field(k_low, model.k_t1)
        S_high = model.entropy_field(k_high, model.k_t1)
        
        assert S_low[0] < 0.1  # Should be ~0 at low k
        assert S_high[0] > 0.9  # Should be ~1 at high k
        
    def test_mixed_entropy(self):
        """Test mixed entropy field."""
        model = EIDOSModel()
        
        k = np.logspace(-4, 0, 100)
        S = model.mixed_entropy(k)
        
        # Check bounds
        assert np.all(S >= 0)
        assert np.all(S <= 1)
        
        # Check monotonicity
        assert np.all(np.diff(S) >= 0)
        
    def test_power_spectrum_modification(self):
        """Test power spectrum modification."""
        model = EIDOSModel()
        
        k = np.logspace(-4, 0, 100)
        P_mod = model.power_spectrum_modification(k)
        
        # Check bounds
        assert np.all(P_mod >= 0)
        assert np.all(P_mod <= 1)
        
        # Check specific scales
        k_001 = 0.001
        idx = np.argmin(np.abs(k - k_001))
        suppression = 1 - P_mod[idx]
        
        assert 0.05 < suppression < 0.15  # Should be ~8% suppression
        
    def test_parameter_predictions(self):
        """Test cosmological parameter predictions."""
        model = EIDOSModel()
        params = model.predict_parameters()
        
        # Check H0
        assert 'H0' in params
        assert 68 < params['H0'] < 72
        
        # Check S8
        assert 'S8' in params
        assert 0.75 < params['S8'] < 0.85
        
        # Check A_L
        assert 'A_L' in params
        assert 1.0 < params['A_L'] < 1.3
        
    def test_chi_squared_improvement(self):
        """Test chi-squared calculation."""
        model = EIDOSModel()
        delta_chi2, sigma = model.chi_squared_improvement()
        
        assert delta_chi2 > 20  # Should be significant
        assert sigma > 3  # Should be >3σ
        
    def test_cmb_modifications(self):
        """Test CMB spectrum modifications."""
        model = EIDOSModel()
        
        ell = np.arange(2, 2500)
        mods = model.cmb_cl_modifications(ell)
        
        # Check quadrupole suppression
        assert mods['TT'][0] < 0.8  # Should be suppressed
        
        # Check B-mode enhancement
        bb_100 = mods['BB'][(ell > 95) & (ell < 105)].mean()
        assert bb_100 > 1.05  # Should be enhanced
        
    def test_forecast_detection(self):
        """Test future experiment forecasts."""
        model = EIDOSModel()
        
        # Test Euclid
        euclid = model.forecast_detection('euclid')
        assert euclid['detection_sigma'] > 5
        
        # Test CMB-S4
        cmbs4 = model.forecast_detection('cmb-s4')
        assert cmbs4['detection_sigma'] > 4


class TestStandaloneFunctions:
    """Test standalone functions."""
    
    def test_entropy_field_function(self):
        """Test standalone entropy field."""
        k = np.array([0.001, 0.01, 0.1])
        k_t = 0.001
        w = 0.15
        
        S = entropy_field(k, k_t, w)
        
        assert len(S) == len(k)
        assert np.all(S >= 0)
        assert np.all(S <= 1)
        
    def test_power_spectrum_modification_function(self):
        """Test standalone power spectrum modification."""
        k = np.logspace(-4, 0, 10)
        P_mod = power_spectrum_modification(k)
        
        assert len(P_mod) == len(k)
        assert np.all(P_mod >= 0)
        assert np.all(P_mod <= 1)


if __name__ == "__main__":
    pytest.main([__file__])