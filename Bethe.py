import numpy as np
import pandas as pd
from scipy.optimize import root_scalar, curve_fit
from scipy.special import comb

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import warnings

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
        """
        USE CASE 2: The Inverse Model. 
        Assuming the experimental data flawlessly begins at t=0.
        """
        t_exp = np.asarray(t_exp)
        Q_exp = np.asarray(Q_exp)
        
        def objective_function(t_array, p_end_trial, k_alpha_trial):
            # Reverse-engineer alpha strictly from the provided experimental time
            alpha_trial = 1 - np.exp(-k_alpha_trial * t_array)
            
            _, _, _, _, Q_sim = self._pure_physics(alpha_trial, p_end_trial)
            
            return Q_sim

        best_params, _ = curve_fit(
            objective_function, 
            t_exp, 
            Q_exp, 
            p0=[p_end_guess, k_alpha_guess],
            bounds=([0.0, 0.0], [1.0, 5])
        )
        
        return best_params[0], best_params[1]


"""

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



if __name__ == "__main__":
    
    try:
        # 1. Ingest the raw data FIRST to establish physical reality
        # (Assuming your CSV has headers. If it does not, add header=None to read_csv)
        df_exp = pd.read_csv("experimental_data.csv")
        
        # Define how many physically invalid points at the end of the experiment you wish to murder.
        points_to_drop = 2
        
        # Amputate the arrays. 
        # The syntax [:-4] instructs Python to extract everything EXCEPT the last 4 items.
        if points_to_drop > 0:
            t_raw = df_exp.iloc[:-points_to_drop, 0].values
            Q_exp = df_exp.iloc[:-points_to_drop, 1].values
        else:
            t_raw = df_exp.iloc[:, 0].values
            Q_exp = df_exp.iloc[:, 1].values

        # Apply the mandatory manuscript scaling standard
        t_exp = t_raw / 1e5
        
        # Dynamically extract Q_eq_0 from the very first data point
        extracted_Q0 = Q_exp[0]
        
        print(f"Data ingested. Extracted baseline Q_eq_0: {extracted_Q0:.3f}")
        
        # 2. Initialize the apparatus using the real initial condition
        gel_system = TopologicalTransformation(
            f_A=8.0, 
            f_B=2.0, 
            s=0.7, 
            M_A=40000.0, 
            M_B=10000.0, 
            c_total=90.0,
            Q_eq_0=extracted_Q0  
        )
        
        print("Fitting, Awaiting convergence...")
        
        # 3. Fit the model
        p_end_opt, k_alpha_opt = gel_system.fit_experiment(
            t_exp=t_exp, 
            Q_exp=Q_exp, 
            p_end_guess=0.7, 
            k_alpha_guess=0.8
        )
        
        print(f"\n--- Convergence Achieved ---")
        print(f"Optimal p_end:   {p_end_opt:.5f}")
        print(f"Optimal k_alpha: {k_alpha_opt:.3e}")

        # --- NEW EXPORT SEQUENCE ---
        print("\nGenerating high-resolution theoretical curve...")
        
        # 4. Feed the optimized parameters back into the forward model
        # Forcing the engine to calculate 1001 points across the degradation gradient
        fitted_df = gel_system.simulate(
            p_end=p_end_opt, 
            k_alpha=k_alpha_opt, 
            num_points=1001
        )
        
        # 5. Export the flat matrix to a CSV file
        output_filename = "fitted_theoretical_curve.csv"
        
        # The index=False argument is strictly mandatory. 
        # Without it, Pandas acts like a bureaucratic accountant and permanently 
        # appends a completely useless column of row numbers (0 to 1000) to your data.
        fitted_df.to_csv(output_filename, index=False)
        
        print(f"Data successfully exiled to {output_filename}.")

        
    except FileNotFoundError:
        print("Error: The experimental_data.csv file is missing from this directory.")
    except Exception as e:
        print(f"The optimizer encountered a fatal reality check. Reason: {e}")


"""


class GelThermodynamicsGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Macosko-Miller Topological Dashboard")
        self.root.geometry("450x700")
        
        # 强制静音 Pandas 的未来弃用警告，免得污染后台
        warnings.simplefilter(action='ignore', category=FutureWarning)
        
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        self.tab_sim = ttk.Frame(self.notebook)
        self.tab_fit = ttk.Frame(self.notebook)
        
        self.notebook.add(self.tab_sim, text="Theoretical Simulation")
        self.notebook.add(self.tab_fit, text="Experimental Fitting")
        
        # 全局参数（被两个标签页共享，以防你输入两次不同的分子量）
        self.shared_vars = {
            'f_A': tk.StringVar(value="8.0"),
            'f_B': tk.StringVar(value="2.0"),
            's': tk.StringVar(value="0.7"),
            'M_A': tk.StringVar(value="40000.0"),
            'M_B': tk.StringVar(value="10000.0"),
            'c_total': tk.StringVar(value="90.0")
        }
        
        self._build_sim_tab()
        self._build_fit_tab()

    def _create_input_row(self, parent, label_text, text_var, row):
        """批量生产输入框的无情流水线"""
        ttk.Label(parent, text=label_text).grid(row=row, column=0, padx=5, pady=5, sticky='e')
        ttk.Entry(parent, textvariable=text_var, width=18).grid(row=row, column=1, padx=5, pady=5, sticky='w')

    def _build_shared_inputs(self, parent_frame):
        frame = ttk.LabelFrame(parent_frame, text="Global Structural Constants")
        frame.pack(fill='x', padx=10, pady=5)
        row = 0
        for key, var in self.shared_vars.items():
            self._create_input_row(frame, f"{key}:", var, row)
            row += 1
        return frame

    def _build_sim_tab(self):
        self._build_shared_inputs(self.tab_sim)
        
        frame_sim = ttk.LabelFrame(self.tab_sim, text="Simulation Parameters")
        frame_sim.pack(fill='x', padx=10, pady=5)
        
        self.sim_Q0 = tk.StringVar(value="3.204")
        self.sim_p = tk.StringVar(value="0.665")
        self.sim_k = tk.StringVar(value="0.828")
        self.sim_pts = tk.StringVar(value="1001")
        
        self._create_input_row(frame_sim, "Q_eq_0 (Initial Swelling):", self.sim_Q0, 0)
        self._create_input_row(frame_sim, "p_end:", self.sim_p, 1)
        self._create_input_row(frame_sim, "k_alpha:", self.sim_k, 2)
        self._create_input_row(frame_sim, "Number of Points:", self.sim_pts, 3)
        
        ttk.Button(self.tab_sim, text="Run Simulation", command=self.run_simulation).pack(pady=15)

    def _build_fit_tab(self):
        self._build_shared_inputs(self.tab_fit)
        
        frame_file = ttk.LabelFrame(self.tab_fit, text="Experimental Data Source")
        frame_file.pack(fill='x', padx=10, pady=5)
        
        self.csv_path = tk.StringVar()
        ttk.Entry(frame_file, textvariable=self.csv_path, state='readonly').pack(side='left', fill='x', expand=True, padx=5, pady=5)
        ttk.Button(frame_file, text="Browse CSV...", command=self._browse_file).pack(side='right', padx=5, pady=5)
        
        frame_fit = ttk.LabelFrame(self.tab_fit, text="Optimization Parameters")
        frame_fit.pack(fill='x', padx=10, pady=5)
        
        self.fit_p_guess = tk.StringVar(value="0.7")
        self.fit_k_guess = tk.StringVar(value="0.8")
        self.fit_pts = tk.StringVar(value="1001")
        self.fit_drop = tk.StringVar(value="2")
        
        self._create_input_row(frame_fit, "p_end_guess:", self.fit_p_guess, 0)
        self._create_input_row(frame_fit, "k_alpha_guess:", self.fit_k_guess, 1)
        self._create_input_row(frame_fit, "Export Resolution:", self.fit_pts, 2)
        self._create_input_row(frame_fit, "Points to Drop (Terminal):", self.fit_drop, 3)
        
        ttk.Button(self.tab_fit, text="Amputate & Fit", command=self.run_fitting).pack(pady=15)

    def _browse_file(self):
        filename = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if filename:
            self.csv_path.set(filename)

    def _get_shared_floats(self):
        return {k: float(v.get()) for k, v in self.shared_vars.items()}

    def run_simulation(self):
        try:
            params = self._get_shared_floats()
            p_end = float(self.sim_p.get())
            k_alpha = float(self.sim_k.get())
            Q_eq_0 = float(self.sim_Q0.get())
            pts = int(self.sim_pts.get())
            
            # 实例化你的物理引擎
            gel = TopologicalTransformation(
                f_A=params['f_A'], f_B=params['f_B'], s=params['s'], 
                M_A=params['M_A'], M_B=params['M_B'], c_total=params['c_total'], Q_eq_0=Q_eq_0
            )
            
            df = gel.simulate(p_end=p_end, k_alpha=k_alpha, num_points=pts)
            
            # 严格按照需求生成包含关键参数的文件名
            filename = f"sim_s{params['s']}_p{p_end}_k{k_alpha:.3f}_Q{Q_eq_0}.csv"
            df.to_csv(filename, index=False)
            
            messagebox.showinfo("Task Complete", f"Theoretical curve generated.\nExiled to: {filename}")
            
        except ValueError:
            messagebox.showerror("Parsing Error", "Stop typing non-numeric characters into the text fields.")
        except Exception as e:
            messagebox.showerror("Catastrophic Failure", str(e))

    def run_fitting(self):
        if not self.csv_path.get():
            messagebox.showerror("Logic Error", "The algorithm cannot optimize a void. Select a CSV file.")
            return
            
        try:
            params = self._get_shared_floats()
            p_guess = float(self.fit_p_guess.get())
            k_guess = float(self.fit_k_guess.get())
            pts = int(self.fit_pts.get())
            drop = int(self.fit_drop.get())
            
            import pandas as pd
            df_exp = pd.read_csv(self.csv_path.get())
            
            # 数字断头台：切除扩散区的幽灵数据[cite: 5]
            if drop > 0:
                t_raw = df_exp.iloc[:-drop, 0].values
                Q_exp = df_exp.iloc[:-drop, 1].values
            else:
                t_raw = df_exp.iloc[:, 0].values
                Q_exp = df_exp.iloc[:, 1].values
                
            # 强制执行你实验室的 1e5 时间缩放标准[cite: 5]
            t_exp = t_raw / 1e5
            Q_eq_0_real = Q_exp[0]
            
            gel = TopologicalTransformation(
                f_A=params['f_A'], f_B=params['f_B'], s=params['s'], 
                M_A=params['M_A'], M_B=params['M_B'], c_total=params['c_total'], Q_eq_0=Q_eq_0_real
            )
            
            # 警告：由于没有多线程，你的窗口将会假死。接受现实。
            print("SciPy optimizer engaged. UI will freeze. Do not panic.")
            
            p_opt, k_opt = gel.fit_experiment(t_exp, Q_exp, p_guess, k_guess)
            
            # 生成高分辨率拟合曲线
            df_fit = gel.simulate(p_end=p_opt, k_alpha=k_opt, num_points=pts)
            
            # 将优化结果烙印在文件名上
            filename = f"fit_s{params['s']}_p{p_opt:.4f}_k{k_opt:.3f}.csv"
            df_fit.to_csv(filename, index=False)
            
            msg = (f"Convergence achieved.\n\n"
                   f"Extracted Q_eq_0: {Q_eq_0_real:.3f}\n"
                   f"Optimized p_end: {p_opt:.5f}\n"
                   f"Optimized k_alpha: {k_opt:.3e}\n\n"
                   f"High-res curve saved to: {filename}")
            messagebox.showinfo("Optimization Complete", msg)
            
        except Exception as e:
            messagebox.showerror("Optimizer Choked", f"The physics engine rejected your reality. Reason:\n{e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = GelThermodynamicsGUI(root)
    root.mainloop()