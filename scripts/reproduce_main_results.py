#!/usr/bin/env python
"""
Reproduce all main results from the EIDOS paper.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import eidos
from eidos.analysis import analyze_cmb, analyze_lss, parameter_forecast
from eidos.visualization import (
    setup_plotting, plot_power_spectrum, plot_cmb_spectra,
    plot_parameter_constraints, plot_future_forecasts
)


def main():
    """Run full EIDOS analysis pipeline."""
    
    print("="*60)
    print("EIDOS: Reproducing Main Results")
    print("="*60)
    
    # Create output directory
    os.makedirs('figures', exist_ok=True)
    os.makedirs('results', exist_ok=True)
    
    # Initialize model
    model = eidos.EIDOSModel()
    
    # 1. Report theoretical parameters
    print("\n1. Theoretical Parameters (no fitting!):")
    print(f"   k_t1 = {model.k_t1:.5f} h/Mpc  (gravitational)")
    print(f"   k_t2 = {model.k_t2:.5f} h/Mpc  (matter)")
    print(f"   f_mix = {model.f_mix:.2f}      (mixing fraction)")
    print(f"   w0 = {model.w0:.2f}           (base width)")
    
    # 2. CMB analysis
    print("\n2. CMB Analysis:")
    try:
        cmb_results = analyze_cmb(model=model)
        
        print(f"   Quadrupole suppression: {cmb_results['quadrupole_suppression']:.1%}")
        print(f"   Octopole suppression: {cmb_results['octopole_suppression']:.1%}")
        print(f"   TE phase amplitude: {cmb_results['te_phase_amplitude']:.3f}")
        print(f"   B-mode enhancement (ℓ~100): {cmb_results['bmode_enhancement']:.1%}")
        print(f"   CMB Δχ² = {cmb_results['total_chi2_improvement']:.1f}")
    except Exception as e:
        print(f"   Warning: CMB analysis failed - {e}")
        cmb_results = {}
    
    # 3. Large-scale structure
    print("\n3. Large-Scale Structure:")
    try:
        lss_results = analyze_lss(model=model)
        
        print(f"   P(k) suppression at k=0.01: {lss_results['pk_suppression_001']:.1%}")
        print(f"   BAO shift: {lss_results['bao_shift']:.1%}")
        print(f"   S8 modification: {lss_results['s8_modification']:.3f}")
        print(f"   LSS Δχ² = {lss_results['chi2_improvement']:.1f}")
    except Exception as e:
        print(f"   Warning: LSS analysis failed - {e}")
        lss_results = {}
    
    # 4. Parameter predictions
    print("\n4. Cosmological Parameters:")
    params = model.predict_parameters()
    
    print(f"   H₀ = {params['H0']:.1f} ± {params['H0_err']:.1f} km/s/Mpc")
    print(f"   S₈ = {params['S8']:.3f} ± {params['S8_err']:.3f}")
    print(f"   A_L = {params['A_L']:.3f} ± {params['A_L_err']:.3f}")
    print(f"   r_s ratio = {params['r_s_ratio']:.3f}")
    
    # 5. Statistical significance
    print("\n5. Statistical Significance:")
    delta_chi2, sigma = model.chi_squared_improvement()
    
    print(f"   Total Δχ² = {delta_chi2:.1f}")
    print(f"   Significance = {sigma:.1f}σ")
    print(f"   p-value = {2 * (1 - eidos.analysis.stats.norm.cdf(sigma)):.2e}")
    print(f"   (With 0 free parameters!)")
    
    # 6. Generate figures
    print("\n6. Generating Figures...")
    
    # Figure 1: Power spectrum modifications
    fig1 = plot_power_spectrum(model, show_components=True)
    fig1.savefig('figures/figure1_power_spectrum.pdf', dpi=300)
    print("   Created: figures/figure1_power_spectrum.pdf")
    
    # Figure 2: CMB angular power spectra
    fig2 = plot_cmb_spectra(model, show_data=True)
    fig2.savefig('figures/figure2_cmb_spectra.pdf', dpi=300)
    print("   Created: figures/figure2_cmb_spectra.pdf")
    
    # Figure 3: Parameter constraints
    fig3 = plot_parameter_constraints(model)
    fig3.savefig('figures/figure3_parameters.pdf', dpi=300)
    print("   Created: figures/figure3_parameters.pdf")
    
    # Figure 4: Future forecasts
    fig4 = plot_future_forecasts(['euclid', 'cmb-s4', 'desi'], model)
    fig4.savefig('figures/figure4_forecasts.pdf', dpi=300)
    print("   Created: figures/figure4_forecasts.pdf")
    
    # 7. Future predictions
    print("\n7. Predictions for Future Experiments:")
    
    experiments = ['euclid', 'cmb-s4', 'desi']
    forecasts = parameter_forecast(model, experiments)
    
    print("\n   Detection Significances:")
    for _, row in forecasts.iterrows():
        if row['experiment'] != 'combined':
            print(f"   {row['experiment'].upper()}: {row['detection_sigma']:.1f}σ")
    
    # Save results
    print("\n8. Saving Results...")
    
    results = {
        'parameters': {
            'k_t1': model.k_t1,
            'k_t2': model.k_t2,
            'f_mix': model.f_mix,
            'w0': model.w0
        },
        'predictions': params,
        'significance': {
            'delta_chi2': delta_chi2,
            'sigma': sigma
        },
        'cmb': cmb_results,
        'lss': lss_results
    }
    
    # Save as text file
    with open('results/main_results.txt', 'w') as f:
        f.write("EIDOS Main Results\n")
        f.write("=" * 50 + "\n\n")
        
        for category, values in results.items():
            f.write(f"{category.upper()}:\n")
            if isinstance(values, dict):
                for key, val in values.items():
                    if isinstance(val, (int, float)):
                        f.write(f"  {key}: {val:.5g}\n")
                    else:
                        f.write(f"  {key}: {val}\n")
            f.write("\n")
    
    print("   Saved: results/main_results.txt")
    
    print("\n" + "="*60)
    print("Analysis Complete!")
    print("="*60)
    
    # Show plots
    plt.show()


if __name__ == "__main__":
    main()