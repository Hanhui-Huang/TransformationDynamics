import numpy as np
from scipy.optimize import root_scalar


def solve_Bethe(f_A, f_B, s, p_end, alpha_array):
    """
    Solves the recursive Bethe approximation (Macosko-Miller) for the A-group 
    extinction probability (Q_A) across an array of degradation states (alpha).
    """
    
    # 1. Resolve the effective reaction probabilities based on stoichiometry
    # Mimics the logic of p_net(), p_A(), and p_B() from the legacy Igor script
    if s > 0.5:
        p_net = p_end * 2.0 * (1.0 - s)
    else:
        p_net = p_end * 2.0 * s
        
    p_A = 0.0 if s == 0 else p_net / (2.0 * s)
    p_B = 0.0 if s == 1 else p_net / (2.0 * (1.0 - s))
    
    # 2. Prepare the output tensor for Q_A
    Q_A_array = np.zeros_like(alpha_array, dtype=float)
    
    # 3. Define the recursive topological expectation equation
    def bethe_eq(x, alpha):
        # x represents the trial extinction probability Q_A
        term1 = 1.0 - p_A
        term2 = (1.0 - alpha) * (x**(f_A - 1.0)) + alpha * (x**(f_A / 2.0 - 1.0))
        term3 = (1.0 - p_B) + p_B * term2
        return term1 + p_A * (term3**(f_B - 1.0)) - x

    # 4. Numerically solve the roots for each decomposition state
    for i, alpha in enumerate(alpha_array):
        # Edge case: Perfect, complete network
        if p_A == 1.0 and p_B == 1.0:
            Q_A_array[i] = 0.0
            continue
            
        # The Macosko-Miller equation always has a trivial root at x = 1.
        # We must check if a lower physical root exists (f(1 - eps) < 0).
        eps = 1e-9
        if bethe_eq(1.0 - eps, alpha) < 0:
            # A physical gel root exists. Isolate it using Brent's method.
            sol = root_scalar(bethe_eq, args=(alpha,), bracket=[0.0, 1.0 - eps], method='brentq')
            Q_A_array[i] = sol.root
        else:
            # No lower root exists. The network is entirely soluble (sol fraction = 1).
            Q_A_array[i] = 1.0 
            
    return Q_A_array






# --- Execution Test ---
if __name__ == "__main__":
    # Generate 101 values from 0 to 1, exactly matching the Igor simulation scale
    alpha_test = np.linspace(0, 1, 101)
    
    # Calculate the extinction probabilities for a standard trial
    # (e.g., 8-arm and 2-arm crosslinkers, s=0.5, p_end=0.9)
    Q_A_results = solve_Bethe(f_A=8, f_B=2, s=0.5, p_end=0.9, alpha_array=alpha_test)
    
    # Display the first 5 results to verify the numerical engine is functioning
    print("Alpha values (first 5):", alpha_test[:5])
    print("\n")
    print("Q_A results (first 5):", Q_A_results[:5])