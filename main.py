import utils
import my_callipso
import my_callipso_remake
import spheres
import argparse

def main():

    parser = argparse.ArgumentParser(description="CALiPPSO for python!")

    # Add arguments with default values
    parser.add_argument('--container_size', default=1, help='Size of the container')
    parser.add_argument('--radius', default=0.1, help='Radius of the spheres')
    parser.add_argument('--num_spheres', default=10, help='Number of spheres to generate and jam')
    parser.add_argument('--plot', default=True, help='Whether to plot the initial and final configurations')

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Generate a random initial configuration of particles.
    container_size = args.container_size
    num_spheres =  args.num_spheres
    radius = args.radius
    if(args.plot == 1):
        plot = True
    else:
        plot = False

    # Printing info about the parameters
    #print("Container size:", container_size)
    print("Number of spheres:", num_spheres)
    print("Radius:",radius)
    print("Plot:", plot)
    print("Container size:", container_size)

    # Generate a random initial configuration of particles
    packing = utils.generate_packing(num_spheres, radius, container_size)

    # Create a list of spheres as objects
    sphere_list = []
    for i in range(len(packing)):
        sphere_list.append(spheres.Sphere(radius, packing[i][0], packing[i][1]))

    # Visualize the initial configuration
    if plot:
        utils.visualize_packing(sphere_list, container_size)

    # Producing a jammed configuration
    #new_spheres, contact_vect, force_magnitude = my_callipso.my_calippso(sphere_list, [0.0000000001, 0.001], 1, 2)
    new_spheres= my_callipso_remake.my_calippso(sphere_list, [0.0000000001, 0.001], 1, 2)

    # Now collect the coordinates of the spheres
    packing = []
    for i in range(len(new_spheres)):   
        packing.append(new_spheres[i].get_coordinates())
    for i in range(len(sphere_list)):
            print("Sphere", i, "coordinates:", sphere_list[i].get_coordinates())

    # Visualize the final configuration
    if plot:
        utils.visualize_packing(sphere_list,container_size)

    # Print the contact vector and the force magnitude
    print("Contact vector:", contact_vect)
    print("Force magnitude:", force_magnitude)

if __name__ == "__main__":
    main()