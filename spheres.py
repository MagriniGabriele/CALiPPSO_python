import numpy as np
import utils

# This class defines a sphere as an object with its radius and coordinates
# TODO: maybe implement a neighbour list for each sphere?
class Sphere:

    def __init__(self, radius, x = None, y=None, z = None):
        self.radius = radius
        self.x = x
        self.y = y
        self.z = z
        self.diameter = 2*radius
        list_coord = [x,y,z]
        # Remove the None values from the list - in this case, to remove 
        # the non-used dimensions out of the 3 possible
        list_coord = utils.remove_none(list_coord)
        self.coordinates = np.array(list_coord)
    
    def get_x(self):
        return self.x
    
    def get_y(self):
        return self.y
    
    def get_z(self):
        return self.z
    
    def get_radius(self):
        return self.radius
    
    def get_diameter(self):
        return self.diameter
    
    def get_coordinates(self):
        return self.coordinates