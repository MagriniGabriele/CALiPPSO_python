import numpy as np
import utils

# This class defines a sphere as an object with its radius and coordinates
# TODO: maybe implement a neighbour list for each sphere?
class Sphere:

    def __init__(self, radius, x = None, y=None):
        self.radius = radius
        self.x = x
        self.y = y
        self.diameter = 2*radius
        self.coordinates = np.array([x,y])
    
    def get_x(self):
        return self.coordinates[0]
    
    def get_y(self):
        return self.coordinates[1]
    
    def get_radius(self):
        return self.radius
    
    def get_diameter(self):
        return self.diameter
    
    def get_coordinates(self):
        return self.coordinates