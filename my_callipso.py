import utils
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import math

# THIS FILE IS DEPRECATED - USE MY_CALLIPSO_REMAKE.PY INSTEAD
#
# This function implements the CALiPPSO algorithm
#
# Parameters:
#   -spheres: a list of Sphere objects
#   -tolerances: a list of two values, the first one is the 
#               tolerance for the maximum overlap between spheres, 
#               the second one is the tolerance for the density - not used here, we adopt the furobi default ones
#   -volume: the volume of the container
#   -dim: the dimension of the system
#
def my_calippso(spheres, tolerances, volume,dim):
    
    # Compute initial density, phi and cutoff l(phi)
    dens = utils.compute_density(dim, spheres, volume)
    print("Initial Density:", dens)
    # Initial cutoff as in the original julia code
    cutoff = 4*spheres[0].get_radius()
    print("Initial cutoff:", cutoff)
    
    # Configuring dimension and number of spheres
    num_spheres = len(spheres)
    d = dim
    iter = 0
    converged = False
    # Initial sbound as in the original julia code
    sbound = 0.01
    while(converged == False):
        # Defining opt model, fresh at each iteration
        model = gp.Model("MyLPModel")

        # Set Info visible about unbounded rays - OPTIONAL FOR DEBUGGING
        model.setParam("InfUnbdInfo", 1)

        # Setting the barrier method as in the original paper is suggested
        model.setParam('Method', 2)

        # If sbound is too small, set it to 0.01 to avoid numerical issues and instabilities
        print("Sbound:", sbound)
        sbound = 0.01
        # Create decision variables S[i, j] for i in range(d) and j in range(n) (a matrix of size d x n)
    
        S = model.addVars((d), (num_spheres), name="S", lb=-sbound, ub=sbound)
        # Create a variable for gamma - here, all the spheres will share the same gamma as in the paper and julia code
        gamma = model.addVar(name="Gamma", lb=1)

        iter += 1
        # Construct the neighbours lists
        neighbours_lists = []
        for i in range(num_spheres):
            neighbours_i = []
            for j in range(num_spheres):
                if i != j:
                    # Compute the distance between the two spheres to check if they are neighbours
                    dist_vect = utils.compute_distance(spheres[i], spheres[j], 1)
                    dist = np.linalg.norm(dist_vect)
                    if dist <= cutoff:
                        neighbours_i.append(j)
            neighbours_lists.append(neighbours_i)

        # Print the neighbours lists and add the constraints
        print("Neighbours lists:", neighbours_lists)
        for i in range(num_spheres):
                    for j in ((neighbours_lists[i])):
                        if i<j:
                            #centre_dist =  utils.compute_distance(spheres[i], spheres[j],1)
                            centre_dist =  utils.compute_alt(spheres[i], spheres[j])
                            
                            #print("Centre distance:", centre_dist)
                            diameter_1 = spheres[i].get_diameter()
                            diameter_2 = spheres[j].get_diameter()
                            radius_diff = (diameter_1 + diameter_2)/2
                            
                            # Imposing the possible positions of S as constraints
                            var_diff = [S[0,i]-S[0,j], S[1,i]-S[1,j]]
                            
                            #Define constraints
                            constraint_expr = -2* np.dot(centre_dist ,var_diff) + gamma*(radius_diff)**2 - np.linalg.norm(centre_dist)**2
                            # Finally, add the constraint to the model
                            model.addConstr(constraint_expr <= 0, name=f"{i,j}")

                            # Now we add the closed boundaries constraints ( As explained at the end of page 3 in the original paper )
                            
                            for k in range(d):
                                dim_constraint =  S[k,i] + spheres[i].get_coordinates()[k]
                                model.addConstr(dim_constraint <= 1, name=f"Constraint_{i}_sbound")
                                model.addConstr(dim_constraint >= 0, name=f"Constraint_{i}_sbound")

        # Set the objective to maximize gamma
        model.setObjective(gamma, GRB.MAXIMIZE)

        # Optimize the model
        model.optimize()

        print("Model STATUS ",model.status)

        # Check if the model is infeasible
        if model.status == GRB.INFEASIBLE:
            print("The model is infeasible")
            return spheres
        if model.status == GRB.UNBOUNDED:
            print("The model is unbounded")
            unbounded_ray = model.unbdRay
            print("Unbounded ray:", unbounded_ray)
            return spheres
        
        # Get the optimal gamma and S values
        optimal_gamma = gamma.X
        optimal_S = {(i, j): S[i, j].X for i in range(d) for j in range(num_spheres)}

        # Print the results
        print("Optimal gamma:", optimal_gamma)
        """
        print("Optimal values for S:")
        for i in range(d):
            for j in range(num_spheres):
                print(f"S[{i}, {j}]: {optimal_S[(i, j)]}")

        for i, constr in enumerate(model.getConstrs()):
            print(f"Constraint {i + 1}: {constr.ConstrName}")
        """
        #print(model.getConstrs())
            #print(f"Dual for {constr.ConstrName}: {constr.Pi}")

        # Update particles' positions
        for i in range(num_spheres):
            for j in range(d):
                spheres[i].coordinates[j] += optimal_S[(j, i)]

        # Update radius
        for i in range(num_spheres):
            spheres[i].diameter = spheres[i].diameter * math.sqrt(optimal_gamma)
            spheres[i].radius = spheres[i].diameter/2


        # Updating the cutoff, density and new bounds
        dens = utils.compute_density(dim, spheres, volume)
        print("New density:", dens)
        cutoff, sbound = utils.bounds_and_cutoff(math.sqrt(optimal_gamma), spheres[0].get_radius(), cutoff, dim)
        sbound = abs(sbound)
       
        print("New cutoff:", cutoff)
        print("New sbound:", sbound)

        # Finally, we store contact indices
        non_zero_dual_indices = []
        for i, constr in enumerate(model.getConstrs()):
            dual_value = constr.Pi
            if dual_value != 0.0:
                non_zero_dual_indices.append(i)
                print(f"Index {i} has a non-zero dual value: {i}")
        # And save the values of the duals
        duals = [model.getConstrs()[i].Pi for i in non_zero_dual_indices]
        print("Duals:", duals)

        max_Si = max(optimal_S.values())
        if (math.sqrt(optimal_gamma) - 1 <= 0.000000000001) and (max_Si <= 0.00001):
            # As suggested in the original julia code, we dispose of the opt model and re-create it afresh for the next iteration
            converged = True
        else:
            model.dispose()
    
    # Defining the jammed configuration
    jammed_config = spheres

    # And finally, we compute the contact vectors and force magnitudes
    
    contact_vects = []
    force_magnitudes = []
    if len(non_zero_dual_indices) == 0:
        print("No contacts found")
        return jammed_config, contact_vects, force_magnitudes
    else:
        print("Contacts found")
        for couple in non_zero_dual_indices:
            i,j = model.getConstrs()[couple].ConstrName.split("_")[0].split(",")
            i = int(i[1])
            j = int(j[1])
            #print("i:", i)
            #print("j:", j)
            coord_diff = utils.compute_alt(spheres[i], spheres[j])
            radius_diff = (spheres[i].get_diameter() + spheres[j].get_diameter())/2
            contact_vects.append(np.array(coord_diff)/radius_diff)
            index_non_zero = non_zero_dual_indices.index(couple)
            force_magnitudes.append(duals[index_non_zero]/radius_diff)
    
    return jammed_config, contact_vects, force_magnitudes