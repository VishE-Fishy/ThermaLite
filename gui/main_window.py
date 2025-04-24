from PyQt6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
                               QLabel, QLineEdit, QPushButton, QGroupBox, QFormLayout,
                               QMessageBox, QDoubleSpinBox, QComboBox, QStatusBar)
from PyQt6.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np
import os
import sys

# Add parent directory to Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from thermal.simulator import ThermalSimulator

class ValidatedDoubleSpinBox(QDoubleSpinBox):
    def __init__(self, min_val, max_val, default_val, decimals=3, step=0.1):
        super().__init__()
        self.setDecimals(decimals)
        self.setMinimum(min_val)
        self.setMaximum(max_val)
        self.setValue(default_val)
        self.setSingleStep(step)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ThermaLite - On-Orbit Thermal-Cycle Simulator")
        self.setMinimumSize(1200, 800)
        
        # Initialize material presets first
        self._setup_material_presets()
        
        # Create central widget and main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        
        # Create input panel
        input_panel = self._create_input_panel()
        main_layout.addWidget(input_panel)
        
        # Create plot panel
        plot_panel = self._create_plot_panel()
        main_layout.addWidget(plot_panel)

        # Add status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        
    def _setup_material_presets(self):
        self.material_presets = {
            "Select Material...": {"α": 0.0, "ε": 0.0},
            "Black Paint": {"α": 0.98, "ε": 0.89},
            "White Paint": {"α": 0.25, "ε": 0.85},
            "Polished Aluminum": {"α": 0.15, "ε": 0.05},
            "Solar Cell": {"α": 0.92, "ε": 0.85},
            "Gold Coating": {"α": 0.30, "ε": 0.03},
            "OSR (Optical Solar Reflector)": {"α": 0.10, "ε": 0.80}
        }
        
    def _create_input_panel(self):
        panel = QGroupBox("Input Parameters")
        layout = QVBoxLayout()
        
        # Geometry inputs
        geometry_group = QGroupBox("Satellite Geometry")
        geometry_layout = QFormLayout()
        
        self.length_input = ValidatedDoubleSpinBox(0.1, 10.0, 1.0)
        self.width_input = ValidatedDoubleSpinBox(0.1, 10.0, 1.0)
        self.height_input = ValidatedDoubleSpinBox(0.1, 10.0, 1.0)
        
        geometry_layout.addRow("Length (m):", self.length_input)
        geometry_layout.addRow("Width (m):", self.width_input)
        geometry_layout.addRow("Height (m):", self.height_input)
        geometry_group.setLayout(geometry_layout)
        
        # Material properties
        material_group = QGroupBox("Material Properties")
        material_layout = QFormLayout()
        
        # Add material preset selector
        self.material_preset = QComboBox()
        self.material_preset.addItems(self.material_presets.keys())
        self.material_preset.currentTextChanged.connect(self._update_material_properties)
        material_layout.addRow("Preset Materials:", self.material_preset)
        
        self.absorptivity_input = ValidatedDoubleSpinBox(0.0, 1.0, 0.3)
        self.emissivity_input = ValidatedDoubleSpinBox(0.0, 1.0, 0.8)
        material_layout.addRow("Solar Absorptivity (α):", self.absorptivity_input)
        material_layout.addRow("Thermal Emissivity (ε):", self.emissivity_input)
        material_group.setLayout(material_layout)
        
        # Orbit parameters
        orbit_group = QGroupBox("Orbit Parameters")
        orbit_layout = QFormLayout()
        
        self.altitude_input = ValidatedDoubleSpinBox(200, 36000, 500)
        self.inclination_input = ValidatedDoubleSpinBox(0, 180, 51.6)  # ISS inclination as default
        self.eccentricity_input = ValidatedDoubleSpinBox(0, 0.9, 0.0, 4)
        
        orbit_layout.addRow("Altitude (km):", self.altitude_input)
        orbit_layout.addRow("Inclination (deg):", self.inclination_input)
        orbit_layout.addRow("Eccentricity:", self.eccentricity_input)
        orbit_group.setLayout(orbit_layout)
        
        # Add run button with professional styling
        run_button = QPushButton("Run Simulation")
        run_button.setStyleSheet("""
            QPushButton {
                background-color: #2196F3;
                color: white;
                padding: 8px;
                font-weight: bold;
                border-radius: 4px;
            }
            QPushButton:hover {
                background-color: #1976D2;
            }
        """)
        run_button.clicked.connect(self._run_simulation)
        
        # Add all groups to panel
        layout.addWidget(geometry_group)
        layout.addWidget(material_group)
        layout.addWidget(orbit_group)
        layout.addWidget(run_button)
        layout.addStretch()
        
        panel.setLayout(layout)
        return panel

    def _update_material_properties(self, material_name):
        if material_name in self.material_presets:
            props = self.material_presets[material_name]
            self.absorptivity_input.setValue(props["α"])
            self.emissivity_input.setValue(props["ε"])
            
    def _create_plot_panel(self):
        panel = QGroupBox("Results")
        layout = QVBoxLayout()
        
        # Create matplotlib figure with improved styling
        plt.style.use('bmh')  # Using a built-in style instead of seaborn
        self.figure = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)
        
        panel.setLayout(layout)
        return panel
    
    def _run_simulation(self):
        """Run the thermal simulation with current parameters"""
        try:
            # Get geometry parameters
            geometry = {
                'length': self.length_input.value(),
                'width': self.width_input.value(),
                'height': self.height_input.value()
            }
            
            # Get material properties
            material_props = {
                'absorptivity': self.absorptivity_input.value(),
                'emissivity': self.emissivity_input.value()
            }
            
            # Get orbit parameters
            orbit_params = {
                'altitude': self.altitude_input.value(),
                'inclination': self.inclination_input.value(),
                'eccentricity': self.eccentricity_input.value()
            }
            
            # Create and run simulator
            simulator = ThermalSimulator(geometry, material_props, orbit_params)
            time, temp = simulator.run_simulation()
            
            # Plot results
            self.plot_results(time, temp)
            self.status_bar.showMessage("Simulation completed successfully")
            
        except ValueError as e:
            # Display the specific validation error message
            self.status_bar.showMessage(f"Input Error: {str(e)}")
        except Exception as e:
            # Handle any other unexpected errors
            self.status_bar.showMessage(f"Simulation failed: {str(e)}")

    def plot_results(self, time, temp):
        # Plot results with enhanced styling
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        
        # Plot temperature curve
        ax.plot(time / 60, temp - 273.15, 
               color='#2196F3', linewidth=2, label='Temperature')
        
        # Add eclipse indication
        eclipse_duration = 0.3 * time[-1]
        ax.axvspan(0, eclipse_duration / 60, color='#9E9E9E', 
                  alpha=0.3, label='Eclipse')
        
        # Calculate key metrics
        max_temp = np.max(temp) - 273.15
        min_temp = np.min(temp) - 273.15
        delta_t = max_temp - min_temp
        
        # Add annotations
        ax.annotate(f'Max: {max_temp:.1f}°C', 
                   xy=(time[-1]/60, max_temp),
                   xytext=(10, 10), textcoords='offset points')
        ax.annotate(f'Min: {min_temp:.1f}°C',
                   xy=(time[0]/60, min_temp),
                   xytext=(10, -10), textcoords='offset points')
        
        # Enhance plot styling
        ax.set_xlabel('Time (minutes)', fontsize=10)
        ax.set_ylabel('Temperature (°C)', fontsize=10)
        ax.set_title('Surface Temperature vs. Time\n'
                    f'ΔT = {delta_t:.1f}°C', fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.7)
        ax.legend(loc='upper right')
        
        # Adjust layout
        self.figure.tight_layout()
        self.canvas.draw() 
