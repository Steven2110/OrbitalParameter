import numpy as np
import matplotlib.pyplot as plt

# Given parameters for the GPS object
GM = 3.986004418e5  # Gravitational parameter in km^3/s^2
a = 26559.94123  # Semi-major axis in km
e = 0.00202  # Eccentricity
T = 43077.61  # Orbital period in seconds
i = np.radians(54.71892)  # Inclination in radians
Ω = np.radians(83.15839)  # Right ascension of ascending node in radians
ω = np.radians(27.00210)  # Argument of periapsis in radians
M0 = np.radians(-165.3855)  # Mean anomaly at epoch in radians

# Earth's radius in km
earth_radius = 6371.0  

# Time array for one complete orbit
time_steps = np.linspace(0, T, num=int(T))

# Function to solve Kepler's equation for E
def solve_kepler(M, e):
    E = M  # Initial guess
    for _ in range(100):  # Iterative solving
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    return E

# Calculate positions over time
x_vals, y_vals, z_vals = [], [], []
for t in time_steps:
    M = M0 + np.sqrt(GM / a**3) * t  # Mean anomaly at time t
    E = solve_kepler(M, e)  # Solve for eccentric anomaly
    # True anomaly
    ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    r = a * (1 - e * np.cos(E))  # Distance from central body

    # Calculate position in orbital plane
    x_orb = r * np.cos(ν)
    y_orb = r * np.sin(ν)
    z_orb = 0

    # Rotate to the 3D inertial frame
    x = (np.cos(Ω) * np.cos(ω + ν) - np.sin(Ω) * np.sin(ω + ν) * np.cos(i)) * x_orb
    x += (-np.cos(Ω) * np.sin(ω + ν) - np.sin(Ω) * np.cos(ω + ν) * np.cos(i)) * y_orb
    
    y = (np.sin(Ω) * np.cos(ω + ν) + np.cos(Ω) * np.sin(ω + ν) * np.cos(i)) * x_orb
    y += (-np.sin(Ω) * np.sin(ω + ν) + np.cos(Ω) * np.cos(ω + ν) * np.cos(i)) * y_orb
    
    z = (np.sin(i) * np.sin(ω + ν)) * x_orb + (np.sin(i) * np.cos(ω + ν)) * y_orb

    x_vals.append(x)
    y_vals.append(y)
    z_vals.append(z)

# Plot the orbit and Earth in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot the GPS object's orbit
ax.plot(x_vals, y_vals, z_vals, label='GPS Orbit', color='red')

# Create a sphere to represent Earth
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x_sphere = earth_radius * np.outer(np.cos(u), np.sin(v))
y_sphere = earth_radius * np.outer(np.sin(u), np.sin(v))
z_sphere = earth_radius * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the Earth sphere
ax.plot_surface(x_sphere, y_sphere, z_sphere, color='blue', alpha=0.3)

# Set labels and title
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('Orbit of GPS Object around Earth')
ax.legend()

# Adjust the aspect ratio to make the Earth look like a real sphere
max_radius = a  # Semi-major axis is the furthest extent of the orbit
ax.set_xlim([-max_radius, max_radius])
ax.set_ylim([-max_radius, max_radius])
ax.set_zlim([-max_radius, max_radius])

# Set equal scaling
ax.set_box_aspect([1, 1, 1])  # Equal scaling for all axes (aspect ratio = 1)

plt.show()
