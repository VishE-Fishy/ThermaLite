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
        
    def _calculate_surface_area(self):
        """Calculate total surface area of the satellite"""
        l, w, h = (self.geometry['length'], self.geometry['width'], 
                   self.geometry['height'])
        return 2 * (l*w + l*h + w*h)
    
    def _estimate_mass(self):
        """Estimate satellite mass based on volume and typical density"""
        l, w, h = (self.geometry['length'], self.geometry['width'], 
                   self.geometry['height'])
        volume = l * w * h
        density = 100.0  # kg/m³ (typical satellite density)
        return volume * density
    
    def _calculate_orbit_period(self):
        """Calculate orbital period using Kepler's Third Law"""
        alt = self.orbit_params['altitude'] * 1000  # convert to meters
        semi_major_axis = self.earth_radius + alt
        period = 2 * np.pi * np.sqrt(semi_major_axis**3 / (G.value * M_earth.value))
        return period
    
    def _calculate_beta_angle(self):
        """Calculate beta angle (angle between orbit plane and sun vector)"""
        inclination = np.radians(self.orbit_params['inclination'])
        # Simplified beta angle calculation (assumes sun in ecliptic plane)
        return np.arcsin(np.sin(inclination) * np.sin(np.radians(23.5)))
    
    def _calculate_view_factors(self, orbit_position):
        """Calculate view factors to Earth and Sun based on orbit position"""
        alt = self.orbit_params['altitude'] * 1000
        earth_angle = np.arcsin(self.earth_radius / (self.earth_radius + alt))
        
        # Simplified view factors
        F_earth = (1 - np.cos(earth_angle)) / 2
        F_space = 1 - F_earth
        
        return F_earth, F_space
    
    def _calculate_eclipse_factor(self, t, period):
        """Calculate eclipse factor (0 = full shadow, 1 = full sun)"""
        # More sophisticated eclipse model using beta angle
        eclipse_duration = period * (1 - np.abs(np.sin(self.beta_angle)))
        
        # Smooth transition for penumbra
        phase = (t % period) / period * 2 * np.pi
        if phase < np.pi:
            return np.clip(np.cos(phase * 2), 0, 1)
        return 1.0
    
    def _thermal_derivative(self, t, T, solar_flux):
        """Define the thermal differential equation with all heat sources"""
        # Get view factors
        F_earth, F_space = self._calculate_view_factors(t / self.orbit_period * 2 * np.pi)
        
        # Solar heat input (direct)
        eclipse_factor = self._calculate_eclipse_factor(t, self.orbit_period)
        q_solar = solar_flux * self.material_props['absorptivity'] * eclipse_factor
        
        # Earth albedo heat input
        q_albedo = (self.solar_constant * self.earth_albedo * F_earth * 
                   self.material_props['absorptivity'] * eclipse_factor)
        
        # Earth IR heat input
        q_earth_ir = self.earth_IR * F_earth * self.material_props['emissivity']
        
        # Heat radiation to space and Earth
        q_out = (self.material_props['emissivity'] * self.stefan_boltzmann * 
                self.surface_area * T**4)
        
        # Net heat flow
        q_net = q_solar + q_albedo + q_earth_ir - q_out
        
        # Temperature change rate
        dT_dt = q_net / (self.mass * self.specific_heat)
        
        return dT_dt
    
    def run_simulation(self, duration=None, time_steps=1000):
        """
        Run the thermal simulation for one orbit
        
        Parameters:
        -----------
        duration : float, optional
            Simulation duration in seconds
        time_steps : int, optional
            Number of time steps for simulation
            
        Returns:
        --------
        tuple : (time_points, temperatures)
        """
        if duration is None:
            duration = self.orbit_period
            
        # Create time points
        t = np.linspace(0, duration, time_steps)
        
        # Calculate solar flux considering orbit geometry
        solar_flux = np.ones(time_steps) * self.solar_constant
        
        # Initial temperature guess (based on equilibrium with Earth IR)
        T0 = (self.earth_IR / (self.stefan_boltzmann * 
              self.material_props['emissivity']))**0.25
        
        # Solve the thermal differential equation
        solution = solve_ivp(
            lambda t, y: self._thermal_derivative(t, y, 
                np.interp(t, t, solar_flux)),
            [0, duration],
            [T0],
            t_eval=t,
            method='RK45',
            rtol=1e-6,
            atol=1e-6
        )
        
        return solution.t, solution.y[0] 
