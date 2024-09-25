# calculation.py
import numpy as np

# Calculation of position and velocity
def calculate_position_velocity(E, a, e, Omega, omega, i, GM):
  # Calculate the value of 𝛏 and 𝛈
  ksi = a * (np.cos(E) - e)
  eta = a * np.sqrt(1-e**2) * np.sin(E)

  # Calculate all sin & cos of some variables
  sin_omega = np.sin(omega)
  cos_omega = np.cos(omega)
  sin_Omega = np.sin(Omega)
  cos_Omega = np.cos(Omega)
  sin_i = np.sin(i)
  cos_i = np.cos(i)

  # Calculate coordinates position (x, y, z)
  x = ksi * (cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i) + eta * (-sin_omega * cos_Omega - cos_omega * sin_Omega * cos_i)
  y = ksi * (cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i) + eta * (-sin_omega * sin_Omega + cos_omega * cos_Omega * cos_i)
  z = ksi * sin_omega * sin_i + eta * cos_omega * sin_i
  
  r_vec_inertial = np.array([x, y, z])

  # Calculate the p & r
  p = a * (1 - e**2)
  r = a * (1 - e * np.cos(E)) # or r = np.linalg.norm(r_vec_inertial)

  # Calculate u
  sin_u = z / (r * sin_i)
  cos_u = x / r * cos_Omega + y / r * sin_Omega

  # Calculate the true anomaly ν
  sin_nu = sin_u * cos_omega - cos_u * sin_omega
  cos_nu = cos_u * cos_omega + sin_u * sin_omega

  # Calculate Vr and Vn
  Vr = np.sqrt(GM/p) * e * sin_nu
  Vn = np.sqrt(GM/p) * (1 + e * cos_nu)

  # Calculate velocity
  v_x = x/r * Vr + (-sin_u * cos_Omega - cos_u * sin_Omega * cos_i) * Vn
  v_y = y/r * Vr + (-sin_u * sin_Omega + cos_u * cos_Omega * cos_i) * Vn
  v_z = z/r * Vr + cos_u * sin_i * Vn
  
  v_vec_inertial = np.array([v_x, v_y, v_z])

  return r_vec_inertial, v_vec_inertial


# Newton's Method
def solve_kepler_newton(M, e, tol=1e-15, max_iter=1000):
  E = M  # Initial guess
  iteration = 0
  for _ in range(max_iter):
    f_E = E - e * np.sin(E) - M  # f(E)
    f_prime_E = 1 - e * np.cos(E)  # f'(E)
    
    # Update E using the iterative formula
    E_next = E - f_E / f_prime_E
    iteration += 1
    # Check the convergence condition
    if abs(E - E_next) < tol:
      return E_next, iteration  # Accept solution if convergence criterion is met
    
    E = E_next  # Update E for the next iteration
  
  return E, iteration  # Return the final value if max_iter is reached

# Danby's Method
def solve_kepler_danby(M, e, tol=1e-10, max_iter=1000):
  E = M + 0.85 * e  # Initial guess according to the formula
  iteration = 0
  
  for _ in range(max_iter):
    f_E = M + e * np.sin(E) - E  # (M + e * sin(E_n) - E_n)
    f_prime_E = E - 2 * (M + e * np.sin(E)) + M + e * np.sin(M + e * np.sin(E))  # Denominator part
    
    # Update E using the iterative formula
    E_next = E - (f_E**2) / f_prime_E
    iteration += 1
    # Check for convergence
    if abs(E - E_next) < tol:
      return E_next, iteration  # Return when the convergence condition is met
    
    E = E_next  # Update E for the next iteration
  
  return E, iteration  # Return the final value if max_iter is reached

def solve_kepler_bisection(M, e, tol=1e-15, max_iter=1000):
  # Since we are using lower and upper bound based on the limit of E value,
  # we need to convert our M if it's negative to the positive radians.
  # If we don't convert it to positive radians we need to change the initial
  # bounds to be [-2π rad;2π rad].
  if M < 0:
    M = M + np.radians(360)

  # Set initial bounds
  E_h = 0  # Lower bound
  E_k = 2 * np.pi  # Upper bound

  # Ensure f(E_h) * f(E_k) < 0
  f_h = E_h - e * np.sin(E_h) - M
  f_k = E_k - e * np.sin(E_k) - M

  if f_h * f_k >= 0:
    raise ValueError("Bisection method requires that f(E_h) and f(E_k) have opposite signs")

  iteration = 0

  while abs(E_k - E_h) > tol and iteration < max_iter:
    iteration += 1
    # Find the midpoint
    E_cp = (E_h + E_k) / 2
    f_cp = E_cp - e * np.sin(E_cp) - M

    # Check if the midpoint is a root or if we should adjust the interval
    if f_cp == 0 or abs(E_k - E_h) / 2 < tol:
      return E_cp, iteration

    # Update the interval based on the sign
    if f_h * f_cp < 0:
      E_k = E_cp
      f_k = f_cp
    else:
      E_h = E_cp
      f_h = f_cp

  # Return the midpoint as the best estimate
  return (E_h + E_k) / 2, iteration

def solve_kepler_muller(M, e, tol=1e-15, max_iter=1000):
  def f(E):
    return E - e * np.sin(E) - M

  def divided_differences(E0, E1):
    f_E0_E1 = (f(E1) - f(E0)) / (E1 - E0)
    return f_E0_E1


  # Initial guesses for Muller's method
  E_k = 2 * np.pi #Ek
  E_k1 = M #Ek-1
  E_k2 = 0 #Ek-2

  iteration = 0

  while iteration < max_iter:
    iteration += 1
    # Divided differences
    f_Ek_Ek1 = divided_differences(E_k, E_k1) # f[Ek;Ek-1]
    f_Ek1_Ek2 = divided_differences(E_k1, E_k2) # f[Ek-1;Ek-2]
    f_Ek_Ek2 = divided_differences(E_k, E_k2) # f[Ek;Ek-2]

    f_Ek_Ek1_Ek2 = (f_Ek1_Ek2 - f_Ek_Ek1) / (E_k2 - E_k) # f[Ek;Ek-1;Ek-2]

    omega = f_Ek_Ek1 + f_Ek_Ek2 - f_Ek1_Ek2

    D = omega ** 2 - 4 * f(E_k) * f_Ek_Ek1_Ek2 # Discriminant

    if D < 0:
      print(f"Discriminant became negative at iteration {iteration}. Stopping the method.")
      return np.nan, iteration

    if abs(omega + np.sqrt(D)) > abs(omega - np.sqrt(D)):
      denominator = omega + np.sqrt(D)
    else:
      denominator = omega - np.sqrt(D)

    E_next = E_k - (2 * f(E_k)) / denominator

    if abs(E_next - E_k) < tol:
      return E_next, iteration

    E_k2, E_k1, E_k = E_k1, E_k, E_next

  return E_k, iteration

def solve_kepler_simple_iteration(M, e, tol=1e-15, max_iter=1000):
    E = M  # Initial guess
    iteration = 0
    
    while iteration < max_iter:
      iteration += 1
      E_new = M + e * np.sin(E)
      if abs(E_new - E) < tol:
          return E_new, iteration
      E = E_new
    
    return E, iteration

# Export calculation methods as dictionary
methods = {
  "Newton": solve_kepler_newton,
  "Danby": solve_kepler_danby,
  "Bisection": solve_kepler_bisection,
  "Muller": solve_kepler_muller,
  "Simple Iteration": solve_kepler_simple_iteration,
}
