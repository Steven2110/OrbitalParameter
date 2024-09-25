# main.py
import tkinter as tk
from tkinter import ttk
import numpy as np
from tkinter import messagebox
from calculation import methods, calculate_position_velocity
from gui_components import create_input_frame, create_results_frame, calculate_orbital_elements

# Create main window
root = tk.Tk()
root.title("Orbital Elements Calculator")
root.geometry("1500x800")

# Create a main frame
main_frame = ttk.Frame(root, padding="10")
main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create input frame
input_frame = create_input_frame(main_frame, lambda: calculate_orbital_elements(results_frame))
input_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Create a frame for displaying results
results_frame = create_results_frame(main_frame)
results_frame.grid(row=1, column=0, sticky=(tk.W, tk.E))

# Start the Tkinter main loop
root.mainloop()
