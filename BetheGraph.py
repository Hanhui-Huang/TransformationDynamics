import numpy as np
import pandas as pd
from scipy.optimize import root_scalar, curve_fit
from scipy.special import comb

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import warnings

import matplotlib.pyplot as plt

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
        """
        This function return c_out, mu, nu, xi, Q arrays based on the input alpha_array and p_end.
        """

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

    def one_parameter_fit_experiment(self, t_exp, Q_exp, p_end_guess, k_alpha_fix):
        """
        USE CASE 3: The Inverse Model. But only fit for p_end_guess.
        The kinetic rate is frozen in the outer scope.
        """
        t_exp = np.asarray(t_exp)
        Q_exp = np.asarray(Q_exp)
        
        # Remove k_alpha_fix from the signature. SciPy is no longer confused.
        def objective_function(t_array, p_end_trial):
            # The function naturally inherits k_alpha_fix from the surrounding environment.
            alpha_trial = 1 - np.exp(-k_alpha_fix * t_array)
            
            _, _, _, _, Q_sim = self._pure_physics(alpha_trial, p_end_trial)
            
            return Q_sim

        best_pend, _ = curve_fit(
            objective_function, 
            t_exp, 
            Q_exp, 
            p0=p_end_guess,
            bounds=(0.0, 1.0)
        )
        
        # Extract the solitary float from the numpy array to prevent downstream type errors.
        return best_pend[0]



class GelThermodynamicsGUI:
    def __init__(self, root):
        """
        The constructor for the graphical shell. Initializes the main window, 
        suppresses annoying warnings, and builds the tabbed architecture.
        """
        self.root = root
        self.root.title("Macosko-Miller Topological Dashboard")
        
        # Fixed window dimensions to prevent the layout from collapsing on itself.
        self.root.geometry("500x780")
        
        # Brutally silence Pandas' future deprecation warnings. 
        # Without this, the terminal will be polluted with irrelevant upgrade nags.
        warnings.simplefilter(action='ignore', category=FutureWarning)
        
        # Construct the primary tabbed notebook architecture using ttk.
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Initialize the three distinct functional frames (Workshops).
        self.tab_sim = ttk.Frame(self.notebook)
        self.tab_fit = ttk.Frame(self.notebook)
        self.tab_plot = ttk.Frame(self.notebook) 
        
        # Register the frames as tabs within the notebook.
        self.notebook.add(self.tab_sim, text="Theoretical Simulation")
        self.notebook.add(self.tab_fit, text="Experimental Fitting")
        self.notebook.add(self.tab_plot, text="Data Visualization")
        
        # Global structural constants stored as Tkinter StringVars.
        # Storing them globally ensures that if you change the molecular weight
        # in the Simulation tab, it automatically updates in the Fitting tab.
        # This prevents catastrophic user-induced inconsistencies.
        self.shared_vars = {
            'f_A': tk.StringVar(value="8.0"),
            'f_B': tk.StringVar(value="2.0"),
            's': tk.StringVar(value="0.7"),
            'M_A': tk.StringVar(value="40000.0"),
            'M_B': tk.StringVar(value="10000.0"),
            'c_total': tk.StringVar(value="90.0")
        }
        
        # Rosters to temporarily hold the file paths of data you intend to plot.
        self.exp_roster = []
        self.theo_roster = []
        
        # Execute the construction methods for each individual tab.
        self._build_sim_tab()
        self._build_fit_tab()
        self._build_plot_tab()

    def _create_input_row(self, parent, label_text, text_var, row):
        """
        A ruthless assembly line method to generate standardized input fields.
        Saves roughly 40 lines of repetitive Tkinter grid boilerplate.
        """
        ttk.Label(parent, text=label_text).grid(row=row, column=0, padx=5, pady=5, sticky='e')
        ttk.Entry(parent, textvariable=text_var, width=18).grid(row=row, column=1, padx=5, pady=5, sticky='w')

    def _build_shared_inputs(self, parent_frame):
        """
        Bolts the global structural constants panel to the top of the provided frame.
        """
        frame = ttk.LabelFrame(parent_frame, text="Global Structural Constants")
        frame.pack(fill='x', padx=10, pady=5)
        
        # Iterate through the shared dictionary to build the entry fields dynamically.
        row = 0
        for key, var in self.shared_vars.items():
            self._create_input_row(frame, f"{key}:", var, row)
            row += 1
        return frame

    def _build_sim_tab(self):
        """
        Constructs Tab 1: The Forward Model Simulation interface.
        """
        self._build_shared_inputs(self.tab_sim)
        
        # Frame specific to simulation kinetics.
        frame_sim = ttk.LabelFrame(self.tab_sim, text="Simulation Parameters")
        frame_sim.pack(fill='x', padx=10, pady=5)
        
        # Local variables specifically governing the theoretical timeline.
        self.sim_Q0 = tk.StringVar(value="3.204")
        self.sim_p = tk.StringVar(value="0.665")
        self.sim_k = tk.StringVar(value="0.828")
        self.sim_pts = tk.StringVar(value="1001")
        
        self._create_input_row(frame_sim, "Q_eq_0 (Initial Swelling):", self.sim_Q0, 0)
        self._create_input_row(frame_sim, "p_end (Terminal Conversion):", self.sim_p, 1)
        self._create_input_row(frame_sim, "k_alpha (Degradation Rate):", self.sim_k, 2)
        self._create_input_row(frame_sim, "Export Resolution (Rows):", self.sim_pts, 3)
        
        # The execution button linking to the run_simulation method.
        ttk.Button(self.tab_sim, text="Run Forward Simulation", command=self.run_simulation).pack(pady=15)


    def _build_fit_tab(self):
        """
        Constructs Tab 2: The Experimental Fitting and Optimization interface.
        """
        self._build_shared_inputs(self.tab_fit)
        
        frame_file = ttk.LabelFrame(self.tab_fit, text="Experimental Data Source")
        frame_file.pack(fill='x', padx=10, pady=5)
        
        self.csv_path = tk.StringVar()
        ttk.Entry(frame_file, textvariable=self.csv_path, state='readonly').pack(side='left', fill='x', expand=True, padx=5, pady=5)
        ttk.Button(frame_file, text="Browse CSV...", command=self._browse_file).pack(side='right', padx=5, pady=5)
        
        frame_fit = ttk.LabelFrame(self.tab_fit, text="Optimization Parameters & Trimming")
        frame_fit.pack(fill='x', padx=10, pady=5)
        
        self.fit_p_guess = tk.StringVar(value="0.7")
        self.fit_k_guess = tk.StringVar(value="0.8")
        self.fit_pts = tk.StringVar(value="1001")
        self.fit_drop = tk.StringVar(value="2")
        
        self._create_input_row(frame_fit, "p_end_guess (Initial Anchor):", self.fit_p_guess, 0)
        self._create_input_row(frame_fit, "k_alpha_guess (Or Fixed Value):", self.fit_k_guess, 1)
        self._create_input_row(frame_fit, "Export Resolution (Rows):", self.fit_pts, 2)
        self._create_input_row(frame_fit, "Terminal Points to Drop (Guillotine):", self.fit_drop, 3)
        
        # Dual execution controls
        frame_buttons = ttk.Frame(self.tab_fit)
        frame_buttons.pack(fill='x', padx=10, pady=15)
        
        ttk.Button(frame_buttons, text="Amputate & Fit (Both Parameters)", command=self.run_fitting).pack(side='left', expand=True, fill='x', padx=5)
        ttk.Button(frame_buttons, text="Amputate & Fit (Fix k_alpha, Fit p_end)", command=self.run_one_parameter_fitting).pack(side='right', expand=True, fill='x', padx=5)


    """
    def _build_fit_tab(self):

        self._build_shared_inputs(self.tab_fit)
        
        # File selection area.
        frame_file = ttk.LabelFrame(self.tab_fit, text="Experimental Data Source")
        frame_file.pack(fill='x', padx=10, pady=5)
        
        self.csv_path = tk.StringVar()
        # Read-only entry to display the selected path; prevents manual typo corruption.
        ttk.Entry(frame_file, textvariable=self.csv_path, state='readonly').pack(side='left', fill='x', expand=True, padx=5, pady=5)
        ttk.Button(frame_file, text="Browse CSV...", command=self._browse_file).pack(side='right', padx=5, pady=5)
        
        # Optimizer parameter bounds and data amputation settings.
        frame_fit = ttk.LabelFrame(self.tab_fit, text="Optimization Parameters & Trimming")
        frame_fit.pack(fill='x', padx=10, pady=5)
        
        self.fit_p_guess = tk.StringVar(value="0.7")
        self.fit_k_guess = tk.StringVar(value="0.8")
        self.fit_pts = tk.StringVar(value="1001")
        self.fit_drop = tk.StringVar(value="2")
        
        self._create_input_row(frame_fit, "p_end_guess (Initial Anchor):", self.fit_p_guess, 0)
        self._create_input_row(frame_fit, "k_alpha_guess (Initial Anchor):", self.fit_k_guess, 1)
        self._create_input_row(frame_fit, "Export Resolution (Rows):", self.fit_pts, 2)
        self._create_input_row(frame_fit, "Terminal Points to Drop (Guillotine):", self.fit_drop, 3)
        
        ttk.Button(self.tab_fit, text="Amputate Data & Engage Optimizer", command=self.run_fitting).pack(pady=15)
    """


    def _build_plot_tab(self):
        """
        Constructs Tab 3: The Data Visualization and Rendering module.
        Acts as a staging area for multiple files before throwing them into Matplotlib.
        """
        # Listbox for raw experimental scatter points.
        frame_exp = ttk.LabelFrame(self.tab_plot, text="Raw Experimental Data (Scatter)")
        frame_exp.pack(fill='x', padx=10, pady=5)
        
        self.list_exp = tk.Listbox(frame_exp, height=4)
        self.list_exp.pack(fill='x', padx=5, pady=5)
        ttk.Button(frame_exp, text="Load Experimental CSV...", command=self._load_exp_file).pack(pady=2)

        # Listbox for theoretical/fitted solid curves.
        frame_theo = ttk.LabelFrame(self.tab_plot, text="Theoretical Data (Lines)")
        frame_theo.pack(fill='x', padx=10, pady=5)
        
        self.list_theo = tk.Listbox(frame_theo, height=4)
        self.list_theo.pack(fill='x', padx=5, pady=5)
        ttk.Button(frame_theo, text="Load Theoretical CSV...", command=self._load_theo_file).pack(pady=2)

        # Execution controls for the renderer.
        frame_ctrl = ttk.Frame(self.tab_plot)
        frame_ctrl.pack(fill='x', padx=10, pady=15)
        
        ttk.Button(frame_ctrl, text="Clear Staged Roster", command=self._clear_plot_roster).pack(side='left', expand=True, fill='x', padx=5)
        ttk.Button(frame_ctrl, text="Render Architecture", command=self.render_plot).pack(side='right', expand=True, fill='x', padx=5)

    def _browse_file(self):
        """Opens a system dialog to retrieve a single CSV path for the fitting tab."""
        filename = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if filename:
            self.csv_path.set(filename)

    def _load_exp_file(self):
        """Appends a selected experimental CSV to the scatter plot roster."""
        filename = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if filename:
            self.exp_roster.append(filename)
            # Display only the base filename in the UI, not the massive absolute path.
            self.list_exp.insert(tk.END, os.path.basename(filename))

    def _load_theo_file(self):
        """Appends a selected theoretical CSV to the line plot roster."""
        filename = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if filename:
            self.theo_roster.append(filename)
            self.list_theo.insert(tk.END, os.path.basename(filename))

    def _clear_plot_roster(self):
        """Purges all memory of loaded files from both the backend lists and UI listboxes."""
        self.exp_roster.clear()
        self.theo_roster.clear()
        self.list_exp.delete(0, tk.END)
        self.list_theo.delete(0, tk.END)

    def _get_shared_floats(self):
        """
        Mercilessly casts fragile Tkinter string inputs into cold, hard floats.
        Returns a dictionary for immediate injection into the physics engine.
        """
        return {k: float(v.get()) for k, v in self.shared_vars.items()}

    def run_simulation(self):
        """
        Executes the pure forward-thermodynamics pipeline.
        """
        try:
            # 1. Gather all required parameters from the GUI fields.
            params = self._get_shared_floats()
            p_end = float(self.sim_p.get())
            k_alpha = float(self.sim_k.get())
            Q_eq_0 = float(self.sim_Q0.get())
            pts = int(self.sim_pts.get())
            
            # 2. Instantiate the physical apparatus.
            gel = TopologicalTransformation(
                f_A=params['f_A'], f_B=params['f_B'], s=params['s'], 
                M_A=params['M_A'], M_B=params['M_B'], c_total=params['c_total'], Q_eq_0=Q_eq_0
            )
            
            # 3. Command the engine to generate the degradation matrix.
            df = gel.simulate(p_end=p_end, k_alpha=k_alpha, num_points=pts)
            
            # 4. Forge the highly specific output filename to prevent data amnesia.
            filename = f"sim_s{params['s']}_p{p_end}_k{k_alpha:.3f}_Q{Q_eq_0}.csv"
            df.to_csv(filename, index=False) # index=False prevents Pandas from printing useless row numbers.
            
            # Notify the user of survival.
            messagebox.showinfo("Task Complete", f"Theoretical curve generated.\nExiled to: {filename}")
            
        except ValueError:
            # Catch the inevitable moment you accidentally type a letter instead of a number.
            messagebox.showerror("Parsing Error", "Stop typing non-numeric characters into the text fields.")
        except Exception as e:
            messagebox.showerror("Catastrophic Failure", str(e))

    def run_fitting(self):
        """
        Executes the SciPy Non-Linear Least Squares optimizer against reality.
        """
        if not self.csv_path.get():
            messagebox.showerror("Logic Error", "The algorithm cannot optimize a void. Select a CSV file.")
            return
            
        try:
            # 1. Harvest parameters and initial guesses.
            params = self._get_shared_floats()
            p_guess = float(self.fit_p_guess.get())
            k_guess = float(self.fit_k_guess.get())
            pts = int(self.fit_pts.get())
            drop = int(self.fit_drop.get())
            
            # 2. Ingest the raw experimental reality.
            df_exp = pd.read_csv(self.csv_path.get())
            
            # 3. The Digital Guillotine: Amputate non-physical diffusion tails.
            # Negative slicing ([:-drop]) isolates the solid elastic regime.
            if drop > 0:
                t_raw = df_exp.iloc[:-drop, 0].values
                Q_exp = df_exp.iloc[:-drop, 1].values
            else:
                t_raw = df_exp.iloc[:, 0].values
                Q_exp = df_exp.iloc[:, 1].values
                
            # 4. Enforce the arbitrary 1e5 time-scaling standard demanded by your manuscript.
            t_exp = t_raw / 1e5
            
            # 5. Extract absolute reality: Equilibrium swelling exactly at t=0.
            Q_eq_0_real = Q_exp[0]
            
            # 6. Build the physics engine around that baseline equilibrium.
            gel = TopologicalTransformation(
                f_A=params['f_A'], f_B=params['f_B'], s=params['s'], 
                M_A=params['M_A'], M_B=params['M_B'], c_total=params['c_total'], Q_eq_0=Q_eq_0_real
            )
            
            # Warning: UI thread will completely freeze during SciPy optimization. 
            # Multi-threading this is an over-engineered waste of time. Accept reality.
            print("SciPy optimizer engaged. UI will freeze. Do not panic.")
            
            # 7. Unleash the optimizer.
            p_opt, k_opt = gel.fit_experiment(t_exp, Q_exp, p_guess, k_guess)
            
            # 8. Feed the victorious parameters back into the forward model for a high-res curve.
            df_fit = gel.simulate(p_end=p_opt, k_alpha=k_opt, num_points=pts)
            df_fit['drop'] = drop

            # 9. Brand the filename with the newly discovered kinetics.
            filename = f"fit_s{params['s']}_p{p_opt:.4f}_k{k_opt:.3f}_drop{drop}.csv"
            df_fit.to_csv(filename, index=False)
            
            msg = (f"Convergence achieved.\n\n"
                   f"Extracted Q_eq_0: {Q_eq_0_real:.3f}\n"
                   f"Optimized p_end: {p_opt:.5f}\n"
                   f"Optimized k_alpha: {k_opt:.3e}\n\n"
                   f"Dropped points: {drop}\n\n"
                   f"High-res curve saved to: {filename}")
            messagebox.showinfo("Optimization Complete", msg)
            
        except Exception as e:
            messagebox.showerror("Optimizer Choked", f"The physics engine rejected your reality. Reason:\n{e}")


    def run_one_parameter_fitting(self):
        """
        Executes a restricted optimization, holding kinetics constant while forcing 
        the algorithm to solve strictly for the terminal topology.
        """
        if not self.csv_path.get():
            messagebox.showerror("Logic Error", "The algorithm cannot optimize a void. Select a CSV file.")
            return
            
        try:
            params = self._get_shared_floats()
            p_guess = float(self.fit_p_guess.get())
            # In this context, the guess field is treated as absolute law.
            k_fixed = float(self.fit_k_guess.get()) 
            pts = int(self.fit_pts.get())
            drop = int(self.fit_drop.get())
            
            df_exp = pd.read_csv(self.csv_path.get())
            
            if drop > 0:
                t_raw = df_exp.iloc[:-drop, 0].values
                Q_exp = df_exp.iloc[:-drop, 1].values
            else:
                t_raw = df_exp.iloc[:, 0].values
                Q_exp = df_exp.iloc[:, 1].values
                
            t_exp = t_raw / 1e5
            Q_eq_0_real = Q_exp[0]
            
            gel = TopologicalTransformation(
                f_A=params['f_A'], f_B=params['f_B'], s=params['s'], 
                M_A=params['M_A'], M_B=params['M_B'], c_total=params['c_total'], Q_eq_0=Q_eq_0_real
            )
            
            print("SciPy optimizer engaged in restricted 1D mode. UI will freeze. Do not panic.")
            
            # Unleash the restricted optimizer
            p_opt = gel.one_parameter_fit_experiment(t_exp, Q_exp, p_guess, k_fixed)
            
            # Feed the victorious parameter back into the forward model
            df_fit = gel.simulate(p_end=p_opt, k_alpha=k_fixed, num_points=pts)
            df_fit['drop'] = drop
            # Brand the filename differently so you know it was a restricted fit
            filename = f"fit1D_s{params['s']}_p{p_opt:.4f}_k{k_fixed:.3f}_drop{drop}.csv"
            df_fit.to_csv(filename, index=False)
            
            msg = (f"Partial convergence achieved.\n\n"
                   f"Extracted Q_eq_0: {Q_eq_0_real:.3f}\n"
                   f"Fixed k_alpha: {k_fixed:.3e}\n"
                   f"Optimized p_end: {p_opt:.5f}\n\n"
                   f"Dropped points: {drop}\n\n"
                   f"High-res curve saved to: {filename}")
            messagebox.showinfo("Optimization Complete", msg)
            
        except Exception as e:
            messagebox.showerror("Optimizer Choked", f"The physics engine rejected your restricted reality. Reason:\n{e}")


    def render_plot(self):
        """
        The visual projector. Consumes the staged rosters and commands 
        Matplotlib to render the overlay architectures.
        """
        if not self.exp_roster and not self.theo_roster:
            messagebox.showerror("Void Error", "There is absolutely nothing to plot. Load some files.")
            return

        # Initialize the figure canvas and axis with high DPI for acceptable publication quality.
        fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
        
        # 1. Project raw experimental data as hollow scatter points.
        for exp_file in self.exp_roster:
            try:
                df = pd.read_csv(exp_file)
                
                # We strictly enforce the 1e5 scale dynamically before plotting.
                t_scaled = df.iloc[:, 0] / 1e5 
                Q_raw = df.iloc[:, 1]
                
                # Strip the path and extension to create a clean legend label.
                file_label = os.path.basename(exp_file).replace('.csv', '')
                ax.scatter(t_scaled, Q_raw, marker='o', s=40, alpha=0.6, label=f"Exp: {file_label}")
            except Exception as e:
                print(f"Failed to render experimental file {exp_file}: {e}")

        # 2. Project the perfect theoretical degradation matrices as solid lines.
        for theo_file in self.theo_roster:
            try:
                df = pd.read_csv(theo_file)
                
                # Extract embedded parameters directly from the top row of the DataFrame.
                # This bypasses the need to parse strings out of filenames.
                s_val = df['s'].iloc[0]
                p_val = df['p_end'].iloc[0]
                k_val = df['k_alpha'].iloc[0]
                
                # Filter out mathematically infinite sol-state swelling arrays.
                # If np.inf touches the axis, Matplotlib's scaling algorithms will immediately self-destruct.
                drop_val = int(df['drop'].iloc[0])
                df_valid = df[df['Q'] != np.inf]
                
                legend_tag = f"Model: s={s_val}, p={p_val:.3f}, k={k_val:.2e}, drop={drop_val}"
                ax.plot(df_valid['t'], df_valid['Q'], linewidth=2.5, label=legend_tag)
                
            except Exception as e:
                print(f"Failed to render theoretical file {theo_file}: {e}")

        # 3. Apply superficial aesthetic decorations and labeling.
        ax.set_xlabel(r"Degradation Time ($10^5$ s)", fontsize=12, fontweight='bold')
        ax.set_ylabel("Swelling Ratio ($Q$)", fontsize=12, fontweight='bold')
        ax.set_title("Gel Degradation Kinetics", fontsize=14)
        
        ax.set_ylim(0, 50)

        # Enable gridlines, because you will inevitably attempt to estimate half-lives by eye.
        ax.grid(True, linestyle='--', alpha=0.5)
        
        # Force the legend into the optimal empty space to prevent data occlusion.
        ax.legend(loc='best', framealpha=0.9, edgecolor='black')
        
        # tight_layout() prevents axis labels from being decapitated by the window border.
        plt.tight_layout()
        
        # Command the OS to render the graphical window and halt execution until it is closed.
        plt.show()

# --- Engine Ignition Sequence ---
if __name__ == "__main__":
    # Initialize the base Tkinter interpreter shell.
    root = tk.Tk()
    # Attach our complex architecture to the shell.
    app = GelThermodynamicsGUI(root)
    # Engage the infinite event-monitoring loop.
    root.mainloop()