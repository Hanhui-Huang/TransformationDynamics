import numpy as np
import pandas as pd
from scipy.optimize import root_scalar, curve_fit
from scipy.special import comb

class TopologicalTransformation:
    def __init__(self, f_A, f_B, s, M_A, M_B, c_total, Q_eq_0):
        """
        Fixed parameters
        """
        self.f_A = float(f_A)
        self.f_B = float(f_B)
        self.s = float(s)
        self.M_A = float(M_A)
        self.M_B = float(M_B)
        self.c_total = float(c_total)
        self.Q_eq_0 = float(Q_eq_0)

    def _solve_extinction(self, alpha_array, p_end):
        """
        Numerical solver for the A-group extinction probability (Q_A).
        """

        # 1. Calculate p_net, p_A, p_B based on s, p_end
        if self.s > 0.5:
            p_net = p_end * 2.0 * (1.0 - self.s)
        else:
            p_net = p_end * 2.0 * self.s
            
        p_A = 0.0 if self.s == 0 else p_net / (2.0 * self.s)
        p_B = 0.0 if self.s == 1 else p_net / (2.0 * (1.0 - self.s))
        
        # 2. Initialization of Q_A_array
        Q_A_array = np.zeros_like(alpha_array, dtype=float)
        
        # 3. Define the Bethe equation
        def bethe_eq(x, alpha):
            term1 = 1.0 - p_A
            term2 = (1.0 - alpha) * (x**(self.f_A - 1.0)) + alpha * (x**(self.f_A / 2.0 - 1.0))
            term3 = (1.0 - p_B) + p_B * term2
            return term1 + p_A * (term3**(self.f_B - 1.0)) - x

        for i, alpha in enumerate(alpha_array):
            if p_A == 1.0 and p_B == 1.0:
                Q_A_array[i] = 0.0
                continue
            
            # Check if there is a root less than 1
            eps = 1e-9
            if bethe_eq(1.0 - eps, alpha) < 0:
                sol = root_scalar(bethe_eq, args=(alpha,), bracket=[0.0, 1.0 - eps], method='brentq')
                Q_A_array[i] = sol.root
            else:
                Q_A_array[i] = 1.0 
                
        return Q_A_array, p_B

    def _calc_active_fractions(self, f, Q_array):
        """
        given Q_A/Q_B, calculate the sum of X_A_f / X_A_f*(f/2) [node density/ strand density ]
        """

        # If the molecule has fewer than 3 functional groups, it cannot form a network node
        if f < 3:
            return np.zeros_like(Q_array), np.zeros_like(Q_array)
            
        X_active = np.zeros_like(Q_array)
        X_strands = np.zeros_like(Q_array)
        
        # A network node requires at least 3 connected paths (i >= 3)
        for i in range(3, int(f) + 1):
            # Binomial probability of having exactly 'i' connected arms
            prob = comb(f, i) * ((1.0 - Q_array)**i) * (Q_array**(f - i))
            X_active += prob
            X_strands += (i / 2.0) * prob
            
        return X_active, X_strands

    def _pure_physics(self, alpha_array, p_end):

        # 1. Solve the extinction probabilities
        Q_A, p_B = self._solve_extinction(alpha_array, p_end)
        Q_B = 1.0 - p_B + p_B * ((1.0 - alpha_array) * (Q_A**(self.f_A - 1.0)) + alpha_array * (Q_A**(self.f_A / 2.0 - 1.0)))
        
        # 2. Compute absolute number densities (N_A, N_B, N_C)
        mass_den = (self.M_A * self.f_B * self.s) + (self.M_B * self.f_A * (1.0 - self.s))

        N_A = self.c_total * self.f_B * self.s * (1.0 - alpha_array) / mass_den
        N_B = self.c_total * self.f_A * (1.0 - self.s) / mass_den
        N_C = self.c_total * self.f_B * self.s * 2.0 * alpha_array / mass_den
        
        # 3. Soluble fraction calculations (c_out)
        X_A0 = Q_A**(self.f_A)
        X_B0 = Q_B**(self.f_B)
        X_C0 = Q_A**(self.f_A / 2.0)
        c_out = (N_A * self.M_A * X_A0) + (N_B * self.M_B * X_B0) + (N_C * (self.M_A / 2.0) * X_C0)
        
        # 4. Extract active nodes and strands via combinatorics
        X_A_active, X_A_strands = self._calc_active_fractions(self.f_A, Q_A)
        X_B_active, X_B_strands = self._calc_active_fractions(self.f_B, Q_B)
        X_C_active, X_C_strands = self._calc_active_fractions(self.f_A / 2.0, Q_A)
        
        # 5. Calculate number densities of the network topology
        mu = (X_A_active * N_A) + (X_B_active * N_B) + (X_C_active * N_C)
        nu = (X_A_strands * N_A) + (X_B_strands * N_B) + (X_C_strands * N_C)
        xi = nu - mu
        
        # 6. Swelling Ratio
        c_gel = self.c_total - c_out
        c_gel_0 = self.c_total - c_out[0]
        xi_0 = xi[0]
        
        # Default the entire array to physical infinity (the sol state)
        Q = np.full_like(alpha_array, np.inf) 
        
        # Prevent division by zero if the initial state is somehow ungelled
        if c_gel_0 > 0 and xi_0 > 0:
            # Mask to isolate only the data points where the macroscopic network exists
            gel_mask = xi > 0
            
            # Normalize concentration and cycle rank
            norm_c_out = c_gel[gel_mask] / c_gel_0
            norm_xi = xi[gel_mask] / xi_0
            

            flory_nu = 0.588
            exp1 = (3 * flory_nu) / (3 * flory_nu - 1)
            norm_xi_c_out = norm_xi / (norm_c_out**exp1)
            
            exp2 = 1 / ((6 * flory_nu + 1) / (9 * flory_nu - 3))
            
            # Apply the swelling ratio only to the gelled portion
            Q[gel_mask] = self.Q_eq_0 / (norm_xi_c_out**exp2)
            
        return c_out, mu, nu, xi, Q
    
    def simulate(self, p_end, k_alpha, num_points=101):
        """
        USE CASE 1: The Forward Model. Generates the fully stamped DataFrame.
        """
        
        alpha_array = np.linspace(0, 1, num_points)
        
        # Protect against division by zero if k_alpha is 0
        if k_alpha > 0:
            t_array = np.zeros_like(alpha_array)
            
            # At alpha = 1.0 (fully degraded), time is infinity
            t_array[alpha_array == 1.0] = np.inf

            # Calculate for alpha < 1.0
            mask = alpha_array < 1.0
            t_array[mask] = -np.log(1.0 - alpha_array[mask]) / k_alpha
            
        else:
            t_array = np.zeros_like(alpha_array)
            
        c_out, mu, nu, xi, Q = self._pure_physics(alpha_array, p_end)
        
        df = pd.DataFrame({
            'alpha': alpha_array,
            't': t_array,
            'c_out': c_out,
            'mu': mu,
            'nu': nu,
            'xi': xi,
            'Q': Q
        })
        
        df['p_end'] = p_end
        df['k_alpha'] = k_alpha
        df['s'] = self.s
        df['f_A'] = self.f_A
        df['f_B'] = self.f_B
        
        return df

    def fit_experiment(self, t_exp, Q_exp, p_end_guess, k_alpha_guess):
        def objective_function(t, p_end_trial, k_alpha_trial):
            alpha_trial = 1.0 - np.exp(-k_alpha_trial * t)
            _, _, _, _, Q_simulated = self._pure_physics(alpha_trial, p_end_trial)
            return Q_simulated

        best_params, covariance = curve_fit(
            objective_function, 
            t_exp, 
            Q_exp, 
            p0=[p_end_guess, k_alpha_guess],
            bounds=([0.0, 0.0], [1.0, np.inf])
        )
        return best_params[0], best_params[1]

# --- Execution Block (Use Case 1) ---
if __name__ == "__main__":
    
    # 1. Bolt the fixed parameters to the chassis, now including Q_eq_0.
    gel_system = TopologicalTransformation(
        f_A=8.0, 
        f_B=2.0, 
        s=0.7, 
        M_A=40000.0, 
        M_B=10000.0, 
        c_total=90.0,
        Q_eq_0=3.204 # Provide your initial experimental swelling ratio here
    )
    
    print("Running theoretical simulation with swelling integration...")
    results_df = gel_system.simulate(p_end=0.665, k_alpha=0.828677349, num_points=101)
    
    print("\n--- Simulation Results (First 20 Rows) ---")
    print(results_df.head(101).to_string())