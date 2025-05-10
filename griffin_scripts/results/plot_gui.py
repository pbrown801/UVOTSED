import tkinter as tk
from tkinter import ttk, simpledialog, messagebox
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# Try importing mplcursors for interactive tooltips
try:
    import mplcursors
except ImportError:
    mplcursors = None

# Load model data
df_model = pd.read_csv("summary.csv")
model_params = ['Spectrum', 'lb', 'rv', 'av16', 'av31', 'redshift', 'dmb', 'f11']
filter_names = ['U', 'B', 'V', 'UVW2', 'UVM2', 'UVW1']
color_indices = [f"{f1}-{f2}" for f1 in filter_names for f2 in filter_names if f1 != f2]

df_model[model_params] = df_model[model_params].fillna('Any')

# Load observational data
obs_path = "../../snbrown_output_20240426.csv"
try:
    df_obs = pd.read_csv(obs_path, skipinitialspace=True)
    obs_filters = {
        'UVW2': ('W2_MBmax', 'W2_eMBmax'),
        'UVM2': ('M2_MBmax', 'M2_eMBmax'),
        'UVW1': ('W1_MBmax', 'W1_eMBmax'),
        'U': ('U_MBmax', 'U_eMBmax'),
        'B': ('B_Mmax', 'B_eMmax'),
        'V': ('V_MBmax', 'V_eBMmax')
    }
except Exception as e:
    df_obs = None
    obs_filters = {}
    print(f"Failed to load observational data: {e}")

class ColorPlotGUI:
    def __init__(self, master):
        self.master = master
        master.title("UVOT Color Plot GUI")

        self.param_vars = {}
        self.combos = {}

        # Axis selectors
        tk.Label(master, text="X-axis Color:").grid(row=0, column=0)
        self.x_choice = ttk.Combobox(master, values=color_indices)
        self.x_choice.grid(row=0, column=1)

        tk.Label(master, text="Y-axis Color:").grid(row=1, column=0)
        self.y_choice = ttk.Combobox(master, values=color_indices)
        self.y_choice.grid(row=1, column=1)

        # Parameter filters
        for i, param in enumerate(model_params):
            tk.Label(master, text=f"{param}:").grid(row=i+2, column=0)
            unique_vals = sorted(df_model[param].dropna().astype(str).unique())
            values = ['Any'] + unique_vals
            combo = ttk.Combobox(master, values=values)
            combo.set('Any')
            combo.grid(row=i+2, column=1)
            self.param_vars[param] = combo

        # Buttons and checkbox
        self.plot_button = tk.Button(master, text="Plot", command=self.make_plot)
        self.plot_button.grid(row=len(model_params)+2, column=0, columnspan=2)

        self.save_button = tk.Button(master, text="Save Plot", command=self.save_plot)
        self.save_button.grid(row=len(model_params)+3, column=0, columnspan=2)

        self.reset_button = tk.Button(master, text="Reset Filters", command=self.reset_filters)
        self.reset_button.grid(row=len(model_params)+4, column=0, columnspan=2)

        self.obs_var = tk.BooleanVar()
        self.obs_checkbox = tk.Checkbutton(master, text="Show Observational Data", variable=self.obs_var)
        self.obs_checkbox.grid(row=len(model_params)+5, column=0, columnspan=2)

        # Plot canvas
        self.figure = plt.Figure(figsize=(6, 5), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master)
        self.canvas.get_tk_widget().grid(row=0, column=2, rowspan=30)

    def get_filtered_data(self):
        filtered_df = df_model.copy()
        for param, combo in self.param_vars.items():
            val = combo.get()
            if val != 'Any':
                filtered_df = filtered_df[filtered_df[param].astype(str) == val]
        return filtered_df

    def make_plot(self):
        x_color = self.x_choice.get()
        y_color = self.y_choice.get()
        if not x_color or not y_color:
            messagebox.showerror("Input Error", "Please select both X and Y color indices.")
            return

        self.ax.clear()
        filtered_df = self.get_filtered_data()

        try:
            x1, x2 = x_color.split("-")
            y1, y2 = y_color.split("-")

            x_vals = filtered_df[x1] - filtered_df[x2]
            y_vals = filtered_df[y1] - filtered_df[y2]

            self.ax.scatter(x_vals, y_vals, label="Model", alpha=0.7)
        except KeyError as e:
            messagebox.showerror("Data Error", f"Missing filter: {e}")
            return

        if self.obs_var.get() and df_obs is not None:
            try:
                x_obs = df_obs[obs_filters[x1][0]] - df_obs[obs_filters[x2][0]]
                y_obs = df_obs[obs_filters[y1][0]] - df_obs[obs_filters[y2][0]]
                xerr = np.sqrt(df_obs[obs_filters[x1][1]]**2 + df_obs[obs_filters[x2][1]]**2)
                yerr = np.sqrt(df_obs[obs_filters[y1][1]]**2 + df_obs[obs_filters[y2][1]]**2)

                points = self.ax.errorbar(x_obs, y_obs, xerr=xerr, yerr=yerr,
                                          fmt='o', label="Observations", color='red', alpha=0.7)

                if mplcursors:
                    cursor = mplcursors.cursor(self.ax, hover=True)
                    @cursor.connect("add")
                    def on_add(sel):
                        i = sel.index
                        if i is not None and i < len(df_obs):
                            sel.annotation.set_text(df_obs.iloc[i]['snname2'])

            except KeyError as e:
                messagebox.showerror("Observation Error", f"Missing observation filter: {e}")
                return

        self.ax.set_xlabel(x_color)
        self.ax.set_ylabel(y_color)
        self.ax.set_title("Color-Color Plot")
        self.ax.legend()
        self.canvas.draw()

    def save_plot(self):
        filename = simpledialog.askstring("Save Plot", "Enter filename (e.g., plot.png):")
        if filename:
            try:
                self.figure.savefig(filename)
                messagebox.showinfo("Success", f"Plot saved as {filename}")
            except Exception as e:
                messagebox.showerror("Save Error", f"Failed to save plot:\n{e}")

    def reset_filters(self):
        self.x_choice.set("")
        self.y_choice.set("")
        for combo in self.param_vars.values():
            combo.set("Any")
        self.obs_var.set(False)

if __name__ == "__main__":
    root = tk.Tk()
    app = ColorPlotGUI(root)
    root.mainloop()
