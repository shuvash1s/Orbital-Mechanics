from astropy import units as u
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.plotting.static import StaticOrbitPlotter
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as spc
from matplotlib.animation import FuncAnimation
import tkinter as tk

def simulate():
        isp = float(isp_entry.get()) * u.s
        r1 = float(r1_entry.get()) * u.km
        r2 = float(r2_entry.get()) * u.km

        g0 = 9.81 * u.m / u.s**2
        ve = isp * g0

        mass_earth = 5.97219e24 * u.kg
        G = spc.gravitational_constant * u.m**3 / u.kg / u.s**2 * u.Unit("")

        # Velocities at circular orbits
        vc1 = np.sqrt(G * mass_earth / r1)
        vc2 = np.sqrt(G * mass_earth / r2)
        a = (r1 + r2) / 2 # Semi-major axis
        ecc = (r2 - r1) / (r1 + r2)

        # Velocities at burns
        vb1 = np.sqrt(G * mass_earth * (2 / r1 - 1 / a))
        vb2 = np.sqrt(G * mass_earth * (2 / r2 - 1 / a))

        delta_v1 = (vb1 - vc1).to(u.m/u.s)
        delta_v2 = (vc2 - vb2).to(u.m/u.s)
        delta_v = delta_v1 + delta_v2

        mass_ratio1 = np.exp(delta_v1 / ve)
        mass_ratio2 = np.exp(delta_v2 / ve)
        mass_ratio_total = np.exp(delta_v / ve)

        print(f"Mass ratio for Burn 1: {mass_ratio1:.3f}")
        print(f"Mass ratio for Burn 2: {mass_ratio2:.3f}")
        print(f"Total mass ratio: {mass_ratio_total:.3f} \n")

        print(f"delta_V1 (Burn 1): {delta_v1:.2f}")
        print(f"delta_V2 (Burn 2): {delta_v2:.2f}")
        print(f"Total delta_V: {delta_v:.2f}")

        epoch = Time("2025-01-01", scale="tdb")

        leo = Orbit.circular(Earth, alt=(r1 - Earth.R).to(u.km), epoch=epoch)
        geo = Orbit.circular(Earth, alt=(r2 - Earth.R).to(u.km), epoch=epoch)
        hohmann = Orbit.from_classical(Earth, a=a, ecc=ecc, inc=0*u.deg, raan=0*u.deg, argp=0*u.deg, nu=0*u.deg, epoch=epoch)

        num_frames = 200
        times = hohmann.epoch + np.linspace(0, hohmann.period.to_value(u.s)/2, num_frames) * u.s
        positions = np.array([hohmann.propagate(t - hohmann.epoch).r.to(u.km).value for t in times])

        fig, ax = plt.subplots(figsize=(10, 10))
        plotter = StaticOrbitPlotter(ax)
        plotter.plot(leo, label="Initial Orbit", color="green")
        plotter.plot(geo, label="Target Orbit", color="blue")
        plotter.plot(hohmann, label="Hohmann Transfer", color="red")

        # Burn annotations
        ax.text(r1.value, 0, f"delta_V1 = {delta_v1.value:.2f} m/s", ha='right', va='bottom', color='red', fontsize=10,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
        ax.text(-r2.value, 0, f"delta_V2 = {delta_v2.value:.2f} m/s", ha='left', va='top', color='red', fontsize=10,
                bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

        # Rocket marker
        rocket_dot, = ax.plot([], [], 'ro', markersize=8, label="Rocket")

        ax.set_title("Hohmann Transfer", pad=20)
        ax.legend(loc="upper left")
        plt.grid()

        path_line, = ax.plot([], [], '', linewidth=1, label="Traveled Path", color='black')

        def update(frame):
                x, y, _ = positions[frame]
                rocket_dot.set_data(x, y)
                path_line.set_data(positions[:frame+1, 0], positions[:frame+1, 1])
                return rocket_dot, path_line

        anim = FuncAnimation(fig, update, frames=num_frames, interval=30, blit=True)
        plt.show()

root = tk.Tk()
root.title("Hohmann Transfer Inputs")

tk.Label(root, text="Specific Impulse (s):").pack()
isp_entry = tk.Entry(root)
isp_entry.insert(0, "0")  # default
isp_entry.pack()

tk.Label(root, text="Initial Orbit Radius (km):").pack()
r1_entry = tk.Entry(root)
r1_entry.insert(0, "0")
r1_entry.pack()

tk.Label(root, text="Final Orbit Radius (km):").pack()
r2_entry = tk.Entry(root)
r2_entry.insert(0, "0")
r2_entry.pack()

tk.Button(root, text="Run Simulation", command=simulate).pack(pady=10)
root.geometry("300x300")

root.mainloop()