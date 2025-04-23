# ThermaLite

# ThermaLite - On-Orbit Thermal-Cycle Simulator

ThermaLite is a professional-grade thermal simulation tool for spacecraft designers. It provides comprehensive thermal analysis of satellite surfaces during orbital operations, including sophisticated modeling of sun/eclipse transitions and environmental heat loads.

## Advanced Features

### Thermal Modeling
- Complete multi-source thermal analysis:
  - Direct solar radiation
  - Earth albedo effects
  - Earth IR radiation
  - Deep space radiation
  - View factor calculations
- Sophisticated eclipse modeling with beta angle calculations
- Penumbra transition effects
- Accurate orbital mechanics using astronomical constants
- Dynamic thermal mass estimation

### Material Properties
- Built-in presets for common spacecraft materials:
  - Black Paint
  - White Paint
  - Polished Aluminum
  - Solar Cells
  - Gold Coating
  - Optical Solar Reflector (OSR)
- Customizable absorptivity and emissivity values

### Orbital Parameters
- Comprehensive orbit configuration:
  - Altitude range: 200km - 36,000km
  - Full inclination range: 0° - 180°
  - Eccentricity support: 0 - 0.9
  - Beta angle calculations
  - View factor computations
- ISS orbit parameters available as defaults

### Analysis & Visualization
- Real-time temperature plotting
- Key metrics display:
  - Maximum/minimum temperatures
  - Temperature delta (ΔT)
  - Eclipse duration
- Professional-grade plotting:
  - Eclipse period visualization
  - Temperature annotations
  - Grid customization
  - Clear axis labeling
- Status updates and progress tracking

## Installation

1. Ensure you have Python 3.8 or newer installed
2. Clone this repository
3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

1. Launch the application:
```bash
python main.py
```

2. Input Parameters:
   - Geometry:
     - Length (0.1m - 10m)
     - Width (0.1m - 10m)
     - Height (0.1m - 10m)
   
   - Material Properties:
     - Select from preset materials OR
     - Custom properties:
       - Solar Absorptivity (α): 0.0 - 1.0
       - Thermal Emissivity (ε): 0.0 - 1.0
   
   - Orbit Parameters:
     - Altitude: 200km - 36,000km
     - Inclination: 0° - 180°
     - Eccentricity: 0.0 - 0.9

3. Click "Run Simulation" to perform the analysis

## Technical Details

### Physical Models
- Solar radiation: 1361 W/m²
- Earth albedo: 0.3 (average)
- Earth IR: 237 W/m²
- Stefan-Boltzmann constant: 5.67e-8 W/m²·K⁴

### Numerical Methods
- RK45 integration method
- Adaptive step size
- Enhanced numerical stability:
  - Relative tolerance: 1e-6
  - Absolute tolerance: 1e-6

### Thermal Calculations
- Dynamic view factor computation
- Sophisticated eclipse modeling
- Multi-source heat flux integration
- Volume-based mass estimation
- Material-specific heat capacity

## Limitations

- Assumes simplified box geometry
- Beta angle calculation assumes sun in ecliptic plane
- No internal heat generation modeling
- Simplified thermal conduction between surfaces
- Single-node thermal model

## Contributing

Contributions are welcome! Areas for potential enhancement:
- Multi-node thermal modeling
- Internal heat source support
- Complex geometry handling
- Thermal contact modeling
- Additional material presets

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

- Astronomical constants from Astropy
- Orbital mechanics calculations using SkyField
- Scientific computing powered by NumPy and SciPy
- Visualization using Matplotlib
- GUI framework using PyQt6 
