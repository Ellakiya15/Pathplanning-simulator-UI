import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.interpolate import interp1d


class MRACSimulation:
    def __init__(self):
        self.t = 0
        self.dt = 0.1
        self.actual_trajectory = np.array([0.0, 0.0])
        self.reference_trajectory = []  # List to store user-drawn trajectory points
        self.ref_index = 0  # Index to track the reference point in the trajectory

    def set_reference_trajectory(self, trajectory):
        """Set the user-drawn reference trajectory"""
        self.reference_trajectory = trajectory
        self.ref_index = 0

    def run_mrac(self, error):
        """MRAC control algorithm (proportional control for simplicity)"""
        k = 1.0  # Gain (can be dynamic)
        return k * error

    def update_trajectory(self):
        """Updates the trajectory based on the MRAC control algorithm"""
        self.t += self.dt
        if len(self.reference_trajectory) > 0 and self.ref_index < len(self.reference_trajectory):
            ref = self.reference_trajectory[self.ref_index]
            self.ref_index += 1
        else:
            ref = np.array([0, 0])  # Default if trajectory is not available
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

        # Labels, Canvas, and Buttons
        self.create_widgets()

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
        self.drawn_points = []  # List to store the drawn reference trajectory points

    def create_widgets(self):
        """Create UI components"""
        # Drawing canvas
        self.drawing_canvas = tk.Canvas(self, bg="white", height=200)
        self.drawing_canvas.pack(pady=5)
        self.drawing_canvas.bind("<B1-Motion>", self.draw_reference_trajectory)

        # Start/Stop buttons
        self.start_button = ttk.Button(self, text="Start Simulation", command=self.start_simulation)
        self.start_button.pack(pady=10)

        self.stop_button = ttk.Button(self, text="Stop Simulation", command=self.stop_simulation)
        self.stop_button.pack(pady=10)

        # Clear drawing button
        self.clear_button = ttk.Button(self, text="Clear Trajectory", command=self.clear_trajectory)
        self.clear_button.pack(pady=5)

    def draw_reference_trajectory(self, event):
        """Capture the user's drawing on the canvas as the reference trajectory"""
        # Get the (x, y) coordinates of the drawn point
        x, y = event.x, event.y
        self.drawing_canvas.create_oval(x - 2, y - 2, x + 2, y + 2, fill='red')

        # Normalize the coordinates (canvas size is 600x200, adjust to match matplotlib plot)
        norm_x = (x / 300) - 1  # Normalize x between -1 and 1
        norm_y = 1 - (y / 100)  # Normalize y between -1 and 1 (invert y-axis for plotting)

        # Append the point to the drawn trajectory
        self.drawn_points.append(np.array([norm_x, norm_y]))

    def resample_trajectory(self, points, num_points=100):
        """Resample the drawn trajectory to have evenly spaced points"""
        points = np.array(points)
        distances = np.sqrt(np.sum(np.diff(points, axis=0) ** 2, axis=1))
        cumulative_distances = np.insert(np.cumsum(distances), 0, 0)

        total_distance = cumulative_distances[-1]
        interpolated_distance = np.linspace(0, total_distance, num_points)

        # Interpolate x and y separately
        interp_x = interp1d(cumulative_distances, points[:, 0], kind='linear')
        interp_y = interp1d(cumulative_distances, points[:, 1], kind='linear')

        resampled_points = np.column_stack((interp_x(interpolated_distance), interp_y(interpolated_distance)))
        return resampled_points

    def start_simulation(self):
        """Start the simulation with the user-drawn reference trajectory"""
        if self.drawn_points:
            # Resample the user-drawn trajectory to have evenly spaced points
            resampled_trajectory = self.resample_trajectory(self.drawn_points)
            self.mrac.set_reference_trajectory(resampled_trajectory)

            self.simulation_running = True
            self.update_simulation()

    def stop_simulation(self):
        """Stop the simulation"""
        self.simulation_running = False

    def clear_trajectory(self):
        """Clear the drawing canvas and reset the trajectory"""
        self.drawing_canvas.delete("all")
        self.drawn_points = []  # Clear stored points

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
