import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class MRACSimulation:
    def __init__(self):
        self.t = 0
        self.dt = 0.1
        self.actual_trajectory = np.array([0.0, 0.0])
        self.reference_trajectory = lambda t: np.array([np.cos(0.1 * t), np.sin(0.1 * t)])

    def run_mrac(self, error):
        """MRAC control algorithm (proportional control for simplicity)"""
        k = 1.0  # Gain (can be dynamic)
        return k * error

    def update_trajectory(self):
        """Updates the trajectory based on the MRAC control algorithm"""
        self.t += self.dt
        ref = self.reference_trajectory(self.t)
        error = ref - self.actual_trajectory
        control_input = self.run_mrac(error)
        self.actual_trajectory += control_input * self.dt
        return self.actual_trajectory, ref

class MRACGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MRAC Trajectory Simulation")
        self.geometry("600x500")

        # MRAC simulation object
        self.mrac = MRACSimulation()

        # Labels and Buttons
        self.start_button = ttk.Button(self, text="Start Simulation", command=self.start_simulation)
        self.start_button.pack(pady=10)

        self.stop_button = ttk.Button(self, text="Stop Simulation", command=self.stop_simulation)
        self.stop_button.pack(pady=10)

        # Setup matplotlib figure for trajectory visualization
        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Placeholder for the line plot (initial empty plot)
        self.actual_line, = self.ax.plot([], [], 'b-', label='Actual Trajectory')
        self.reference_line, = self.ax.plot([], [], 'r--', label='Reference Trajectory')
        self.ax.legend()
        self.ax.set_xlim(-2, 2)
        self.ax.set_ylim(-2, 2)

        self.simulation_running = False

    def start_simulation(self):
        """Start the simulation"""
        self.simulation_running = True
        self.update_simulation()

    def stop_simulation(self):
        """Stop the simulation"""
        self.simulation_running = False

    def update_simulation(self):
        """Update the UI with the new trajectory from MRAC"""
        if not self.simulation_running:
            return

        actual_traj, ref_traj = self.mrac.update_trajectory()

        # Update the plot with new data
        self.actual_line.set_data(np.append(self.actual_line.get_xdata(), actual_traj[0]),
                                  np.append(self.actual_line.get_ydata(), actual_traj[1]))
        self.reference_line.set_data(np.append(self.reference_line.get_xdata(), ref_traj[0]),
                                     np.append(self.reference_line.get_ydata(), ref_traj[1]))

        self.ax.set_xlim(-2, 2)
        self.ax.set_ylim(-2, 2)
        self.canvas.draw()

        # Continue the simulation
        self.after(100, self.update_simulation)

# Run the GUI application
if __name__ == '__main__':
    app = MRACGUI()
    app.mainloop()
