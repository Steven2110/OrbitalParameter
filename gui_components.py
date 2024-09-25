# gui_components.py
import tkinter as tk
from tkinter import ttk
import numpy as np
from tkinter import messagebox
from calculation import methods, calculate_position_velocity

results = {}

# Function to create the input frame with input fields and the calculate button
def create_input_frame(parent, calculate_callback):
    input_frame = ttk.Frame(parent, padding="10", relief=tk.GROOVE)
    input_frame.grid_propagate(False)  # Prevent the frame from resizing
    input_frame.update_idletasks()  # Update the frame to get its current size
    input_frame.configure(width=1400, height=350)  # Adjust these dimensions as needed

    # Input fields label
    input_label = ttk.Label(input_frame, text="Input Parameters", font=("Helvetica", 24, "bold"))
    input_label.grid(row=0, column=0, columnspan=3, pady=(0, 10))

    # Define unit variables
    global inclination_unit, omega_unit, Omega_unit, M0_unit
    inclination_unit = tk.StringVar()
    omega_unit = tk.StringVar()
    Omega_unit = tk.StringVar()
    M0_unit = tk.StringVar()

    # Create input fields with default values
    global eccentricity_entry, inclination_entry, semi_major_axis_entry, asc_node_entry, arg_periapsis_entry, mean_anomaly_entry, period_entry, grav_param_entry
    eccentricity_entry = create_input_field(input_frame, "Eccentricity (e):", 1, tk.StringVar(), "0.00202")
    inclination_entry = create_input_field(input_frame, "Inclination (i):", 2, inclination_unit, "54.71892", with_unit=True)
    semi_major_axis_entry = create_input_field(input_frame, "Semi-major axis (a) [km]:", 3, tk.StringVar(), "26559.94123")
    asc_node_entry = create_input_field(input_frame, "Right Ascension (Ω):", 4, Omega_unit, "83.15839", with_unit=True)
    arg_periapsis_entry = create_input_field(input_frame, "Argument of Periapsis (ω):", 5, omega_unit, "27.00210", with_unit=True)
    mean_anomaly_entry = create_input_field(input_frame, "Mean Anomaly (M0):", 6, M0_unit, "-165.3855", with_unit=True)
    period_entry = create_input_field(input_frame, "Orbital Period (T) [sec]:", 7, tk.StringVar(), "43077.61")
    grav_param_entry = create_input_field(input_frame, "Gravitational Parameter (GM) [km^3/s^2]:", 8, tk.StringVar(), "398600.4418")

    # Create a button to trigger calculations
    calculate_button = ttk.Button(input_frame, text="Calculate", command=calculate_callback)
    calculate_button.grid(row=9, columnspan=3, pady=10)

    return input_frame

# Helper function to create an individual input field with units
def create_input_field(parent, label_text, row, unit_var, default_value, with_unit = False, default_unit="degrees"):
    label = ttk.Label(parent, text=label_text)
    label.grid(row=row, column=0, sticky=tk.W, padx=5)
    
    entry = ttk.Entry(parent, width=20)
    entry.grid(row=row, column=1, padx=5)
    entry.insert(0, default_value)
    
    if with_unit:
        unit_options = ["degrees", "radians"]
        unit_menu = ttk.OptionMenu(parent, unit_var, default_unit, *unit_options)
        unit_menu.grid(row=row, column=2, padx=5)
    
    return entry

# Function to handle calculations
def calculate_orbital_elements(results_frame):
    try:
        # Retrieve values from input fields
        e = float(eccentricity_entry.get())
        a = float(semi_major_axis_entry.get())
        T = float(period_entry.get())
        GM = float(grav_param_entry.get())

        inclination = float(inclination_entry.get())
        omega = float(arg_periapsis_entry.get())
        Omega = float(asc_node_entry.get())
        M0 = float(mean_anomaly_entry.get())

        # Convert degrees to radians if necessary
        if inclination_unit.get() == "degrees":
            inclination = np.radians(inclination)
        if omega_unit.get() == "degrees":
            omega = np.radians(omega)
        if Omega_unit.get() == "degrees":
            Omega = np.radians(Omega)
        if M0_unit.get() == "degrees":
            M0 = np.radians(M0)

        # Perform calculations for each method
        for method_name, method in methods.items():
            E, iteration = method(M0, e)
            r_inertial, v_inertial = calculate_position_velocity(E, a, e, Omega, omega, inclination, GM)

            # Store the results in the results dictionary
            results[method_name] = {
                "Eccentric Anomaly (E)": E,
                "Iterations": iteration,
                "Position": r_inertial.tolist(),
                "Velocity": v_inertial.tolist(),
            }

        # Display the results in the GUI
        display_results(results_frame)
        
    except ValueError:
        messagebox.showerror("Input Error", "Please enter valid numbers for all parameters.")

# Create the results frame
def create_results_frame(parent):
    results_frame = ttk.Frame(parent, padding="10", relief=tk.GROOVE)
    return results_frame

# Display calculation results
def display_results(results_frame):
    # Clear any existing widgets in results_frame
    for widget in results_frame.winfo_children():
        widget.destroy()

    row = 0
    column = 0
    
    # Display the first three methods on the left
    for method_name in ["Newton", "Danby", "Bisection"]:
        add_method_display(results_frame, method_name, results[method_name], row, column)
        row += 6

    # Display the remaining two methods on the right
    row = 0
    column = 2  # Start second column for right side
    for method_name in ["Muller", "Simple Iteration"]:
        add_method_display(results_frame, method_name, results[method_name], row, column)
        row += 6
        
    row += 5
    ttk.Separator(results_frame, orient='horizontal').grid(row=row, column=column, sticky="ew", pady=(10, 5), columnspan=2)

def add_method_display(frame, method_name, result, row, column):
    ttk.Label(frame, text=f"{method_name} Method", font=("Helvetica", 18, "bold")).grid(row=row, column=column, sticky=tk.W, columnspan=2, padx=10, pady=(10, 0))
    row += 1
    
    ttk.Label(frame, text=f"Eccentric Anomaly (E) (radians):").grid(row=row, column=column, sticky=tk.W, padx=10)
    ttk.Label(frame, text=f"{result['Eccentric Anomaly (E)']}").grid(row=row, column=column + 1, sticky=tk.W, padx=10)
    row += 1
    
    ttk.Label(frame, text=f"Iterations:").grid(row=row, column=column, sticky=tk.W, padx=10)
    ttk.Label(frame, text=f"{result['Iterations']}").grid(row=row, column=column + 1, sticky=tk.W, padx=10)
    row += 1
    
    ttk.Label(frame, text=f"Coordinates (km):").grid(row=row, column=column, sticky=tk.W, padx=10)
    ttk.Label(frame, text=f"{result['Position']}").grid(row=row, column=column + 1, sticky=tk.W, padx=10)
    row += 1
    
    ttk.Label(frame, text=f"Velocity (km/s):").grid(row=row, column=column, sticky=tk.W, padx=10)
    ttk.Label(frame, text=f"{result['Velocity']}").grid(row=row, column=column + 1, sticky=tk.W, padx=10)
    row += 1
    
    # Add a divider
    ttk.Separator(frame, orient='horizontal').grid(row=row, column=column, sticky="ew", pady=(10, 5), columnspan=2)
