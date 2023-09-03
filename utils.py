import numpy as np
import math
import matplotlib.pyplot as plt

# This function computes the volume of a sphere in d dimensions - works from 2d upwards
def compute_volume(radius, dimensions):
  return (math.pi**(dimensions/2))*(radius**dimensions)/math.gamma(dimensions/2 + 1)

def compute_particles_volume(dimensions, particles):
  volumes_sum = 0
  for i in range(len(particles)):
    volumes_sum += compute_volume(radius = (particles[i].get_radius()), dimensions=dimensions)
  return volumes_sum

# Computes the density of the system(phi) to check if it is jammed
def compute_density(dimensions,particles, volume):
  return compute_particles_volume(dimensions, particles)/volume

# This computes the distance with the nearest image convention (as in the original paper)
def compute_distance(A, B, box_size):
    dx = np.abs(A.get_x() - B.get_x())
    dy = np.abs(A.get_y() - B.get_y())

    # Apply nearest image convention
    if dx > 0.5 * box_size:
        dx = box_size - dx
    if dy > 0.5 * box_size:
        dy = box_size - dy

    #distance = np.sqrt(dx**2 + dy**2)
    return [dx, dy]

# This computes the distance without the nearest image convention
def compute_easy(A,B):
    return np.linalg.norm(A.get_coordinates()-B.get_coordinates())

# This computes the distance on each axis - necessary for the constraints
def compute_alt(A,B):
     return[A.x-B.x, A.y-B.y]


# This function computes the cutoff l(phi) for a given gamma and bounds - as in the julia code
def bounds_and_cutoff(sqr_gamma, radius, l0, d, thresholds=(5e-4, 1e-5), sbound=0.01):
    th1, th2 = thresholds

    if (sqr_gamma > 1 + th1):
        s_bound = (0.5 * l0 - max(1.5, sqr_gamma) * radius) / (d ** 0.5)
        return l0, s_bound
    elif (th1 >= sqr_gamma - 1 > th2):
        return min(4 * radius, l0), 0.1 * radius
    else:
        return min(2.7 * radius, l0), sbound * radius
    

# Simple function to check if two spheres overlap
def check_overlap(new_sphere, existing_spheres, radius):
    for sphere in existing_spheres:
        if np.linalg.norm(new_sphere - sphere) < 2 * radius:
            return True
    return False

# This function generates a random packing of spheres in a container of size "container_size"
# NB: this works through rejection sampling, so it might take a while to generate a packing
def generate_packing(num_spheres, radius, container_size):
    spheres = []
    for _ in range(num_spheres):
        new_sphere = np.random.rand(2) * (container_size - 2 * radius) + radius
        while check_overlap(new_sphere, spheres, radius):
            new_sphere = np.random.rand(2) * (container_size - 2 * radius) + radius
        spheres.append(new_sphere)
    return spheres

# This function visualizes a packing of spheres - for 2d
def visualize_packing(sphere_list, container_size):
    fig, ax = plt.subplots()
    ax.set_xlim(0, container_size)
    ax.set_ylim(0, container_size)
    ax.set_aspect('equal')
    
    for sphere in sphere_list:
        circle = plt.Circle(sphere.coordinates, sphere.radius, color='b', fill=False)
        ax.add_artist(circle)
    
    plt.show()

# This function removes the None values from a list - necessary in sphere object initialization
def remove_none(list):
    return [x for x in list if x is not None]