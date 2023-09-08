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
    parser.add_argument('--debug', default=False, help='To print info necessary for debugging')
    parser.add_argument('--max_iter', default=1000, help='Maximum number of iterations for the algorithm')

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Generate a random initial configuration of particles.
    container_size = int(args.container_size)
    num_spheres =  int(args.num_spheres)
    radius = float(args.radius)
    debug = args.debug
    max_iter = int(args.max_iter)
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

    # Check the initial density - if it is too low, warn the user
    initial_density = utils.compute_density(2, sphere_list, container_size)
    if(initial_density < 0.5):
        print(f"Warning!!! Initial density({initial_density}) is too low, convergence may be slow and not assured")
        print("\n")
        print("Want to continue? (y/n)")
        answer = input()
        if(answer == "n"):
            exit()

    # Producing a jammed configuration
    new_spheres, contact_vect, force_magnitude= my_callipso_remake.my_calippso(sphere_list, container_size, 2, debug=debug, max_iter=max_iter)

    # Now collect the coordinates of the spheres
    packing = []
    for i in range(len(new_spheres)):   
        packing.append(new_spheres[i].get_coordinates())
    for i in range(len(sphere_list)):
            print("Final Sphere", i, "coordinates:", sphere_list[i].get_coordinates())

    # Visualize the final configuration
    if plot:
        utils.visualize_packing(sphere_list,container_size)

    # Print the contact vector and the force magnitude
    print("\n")
    print("Contact vector:", contact_vect)
    print("\n")
    print("Force magnitude:", force_magnitude)

if __name__ == "__main__":
    main()