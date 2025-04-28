import numpy as np
from scipy.integrate import solve_ivp
from skyfield.api import load, wgs84
from astropy import units as u
from astropy.constants import sigma_sb, G, M_earth

class ThermalSimulator:
    def __init__(self, geometry, material_props, orbit_params):
        """
        Initialize the thermal simulator with satellite parameters
        
        Parameters:
        -----------
        geometry : dict
            Contains 'length', 'width', 'height' in meters
        material_props : dict
            Contains 'absorptivity', 'emissivity'
        orbit_params : dict
            Contains 'altitude' (km), 'inclination' (deg), 'eccentricity'
        """
        # Validate inputs
        self._validate_inputs(geometry, material_props, orbit_params)
        
        self.geometry = geometry
        self.material_props = material_props
        self.orbit_params = orbit_params
        
        # Physical constants
        self.solar_constant = 1361.0  # W/m²
        self.earth_radius = 6371000.0  # meters
        self.stefan_boltzmann = sigma_sb.value
        self.earth_albedo = 0.3  # Earth's average albedo
        self.earth_IR = 237.0  # W/m² (Earth's average IR emission)
        
        # Calculate derived parameters
        self.surface_area = self._calculate_surface_area()
        self.mass = self._estimate_mass()
        self.specific_heat = 900.0  # J/kg·K (aluminum)
        self.orbit_period = self._calculate_orbit_period()
        self.beta_angle = self._calculate_beta_angle()

    def _validate_inputs(self, geometry, material_props, orbit_params):
        """Validate all input parameters"""
        # Validate geometry
        for dim in ['length', 'width', 'height']:
            if dim not in geometry:
                raise ValueError(f"Missing {dim} in geometry parameters")
            if not isinstance(geometry[dim], (int, float)):
                raise ValueError(f"{dim} must be a number")
            if geometry[dim] <= 0:
                raise ValueError(f"{dim} must be positive")
            if geometry[dim] > 10.0:
                raise ValueError(f"{dim} must be less than 10 meters")

        # Validate material properties
        for prop in ['absorptivity', 'emissivity']:
            if prop not in material_props:
                raise ValueError(f"Missing {prop} in material properties")
            if not isinstance(material_props[prop], (int, float)):
                raise ValueError(f"{prop} must be a number")
            if not 0 <= material_props[prop] <= 1:
                raise ValueError(f"{prop} must be between 0 and 1")

        # Validate orbit parameters
        if 'altitude' not in orbit_params:
            raise ValueError("Missing altitude in orbit parameters")
        if not isinstance(orbit_params['altitude'], (int, float)):
            raise ValueError("Altitude must be a number")
        if not 200 <= orbit_params['altitude'] <= 36000:
            raise ValueError("Altitude must be between 200 and 36000 km")

        if 'inclination' not in orbit_params:
            raise ValueError("Missing inclination in orbit parameters")
        if not isinstance(orbit_params['inclination'], (int, float)):
            raise ValueError("Inclination must be a number")
        if not 0 <= orbit_params['inclination'] <= 180:
            raise ValueError("Inclination must be between 0 and 180 degrees")

        if 'eccentricity' not in orbit_params:
            raise ValueError("Missing eccentricity in orbit parameters")
        if not isinstance(orbit_params['eccentricity'], (int, float)):
            raise ValueError("Eccentricity must be a number")
        if not 0 <= orbit_params['eccentricity'] < 0.9:
            raise ValueError("Eccentricity must be between 0 and 0.9")

    def _calculate_surface_area(self):
        """Calculate the total surface area of the satellite"""
        length = self.geometry['length']
        width = self.geometry['width']
        height = self.geometry['height']
        
        # Calculate area of each face
        front_back = 2 * length * width
        top_bottom = 2 * length * height
        sides = 2 * width * height
        
        return front_back + top_bottom + sides

    def _estimate_mass(self):
        """Estimate satellite mass based on volume and typical density"""
        length = self.geometry['length']
        width = self.geometry['width']
        height = self.geometry['height']
        
        # Calculate volume in m³
        volume = length * width * height
        
        # Assume typical CubeSat density of 1000 kg/m³
        density = 1000  # kg/m³
        return volume * density

    def _calculate_orbit_period(self):
        """Calculate orbital period using Kepler's Third Law"""
        # Convert altitude to meters and add Earth radius
        orbit_radius = (self.orbit_params['altitude'] * 1000) + self.earth_radius
        
        # Calculate period using T = 2π√(r³/μ)
        # where μ = GM (standard gravitational parameter)
        mu = G.value * M_earth.value
        period = 2 * np.pi * np.sqrt(orbit_radius**3 / mu)
        return period

    def _calculate_beta_angle(self):
        """Calculate beta angle (angle between orbit plane and sun vector)"""
        # Convert inclination to radians
        inclination = np.radians(self.orbit_params['inclination'])
        
        # Assume worst-case beta angle for thermal analysis
        # This is typically the maximum beta angle possible for the given inclination
        beta = np.abs(inclination)
        return beta

    def run_simulation(self, duration=None):
        """
        Run the thermal simulation
        
        Parameters:
        -----------
        duration : float, optional
            Simulation duration in seconds. Defaults to one orbit period.
            
        Returns:
        --------
        time : ndarray
            Time points in seconds
        temperatures : ndarray
            Temperature at each time point in Kelvin
        """
        if duration is None:
            duration = self.orbit_period

        # Initial temperature (assume room temperature)
        T0 = 293.15  # 20°C in Kelvin

        # Define the thermal differential equation
        def dT_dt(t, T):
            # Solar heat input (consider eclipse)
            if t % self.orbit_period < 0.7 * self.orbit_period:  # In sunlight
                Q_solar = self.solar_constant * self.material_props['absorptivity'] * (self.surface_area / 6)
                Q_albedo = self.solar_constant * self.earth_albedo * self.material_props['absorptivity'] * (self.surface_area / 6)
            else:  # In eclipse
                Q_solar = 0
                Q_albedo = 0

            # Earth IR heat input (constant)
            Q_earth = self.earth_IR * self.material_props['absorptivity'] * (self.surface_area / 6)

            # Radiative heat loss
            Q_rad = self.material_props['emissivity'] * self.stefan_boltzmann * self.surface_area * T**4

            # Net heat flow
            Q_net = Q_solar + Q_albedo + Q_earth - Q_rad

            # Temperature change rate
            dT = Q_net / (self.mass * self.specific_heat)
            return dT

        # Solve the differential equation
        t_eval = np.linspace(0, duration, 1000)
        solution = solve_ivp(dT_dt, (0, duration), [T0], t_eval=t_eval, method='RK45')

        return solution.t, solution.y[0]

    # ... rest of the existing methods stay the same ... 
