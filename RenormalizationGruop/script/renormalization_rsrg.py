import numpy as np
import time 
import renormalization_utils as ut
import matplotlib.pyplot as plt



def rsrg():
    
    # Hilbert space dimension
    dimension = ut.dimension
 
    # Size of the system = # of particles
    N_particles = ut.N_particles   
    
    # Bins for the strenght (lambda) parameter
    binlambda = ut.binlambda
    
    # Precision for the evaluation of the results
    precision = ut.precision
    
    # Max value of parameter lambda for the variation of the strenght parameter
    p_lambda = ut.p_lambda
        
    # Discretization of the strenght parameter
    dlambda = ut.dlambda
    
    # Initialization of the Pauli matrices
    sigma_x, sigma_z = ut.pauli()
    
    # Interaction part of the Isign model
    interaction_hamiltonian_x = ut.x_hamiltonian(sigma_x, N_particles)
    
    # Non-Interaction part of the Isign model
    n_interaction_hamiltonian_z = ut.z_hamiltonian(sigma_z, N_particles)
    
    # Array for store the partial time of the computation    
    time_spent = []
    
    # Identity matrix of dimension d**(N-1) x d**(N-1)    
    identity_Nm1xNm1 = np.eye(dimension**(N_particles-1))
    
    # Identity matrix of dimension d**(N) x d**(N)
    identity_NxN = np.eye(dimension**N_particles)
        
    # Storing the starting time
    start = time.time()
    
    # Initializing the arrays that will be used to store the results of the computation for the plotting part.
    step_lambda_0 = []
    level_lambda_0 = []
    step_lambda_max = []
    level_lambda_max = []
    steps_convergence = []
    lambda_convergence = []
    gs_convergence = []
    
    for i in range(binlambda+1):
        
        # Taking the partial of the computation for index 'i'
        partial = time.time()
        
        # Initialization of the full Hamiltonian of the Isign model
        full_hamiltonian = dlambda*i*n_interaction_hamiltonian_z + interaction_hamiltonian_x
        
        # ========================================================== #
        # ======================DOUBLING PART======================= #
        # ========================================================== #
        
        # Left-interaction term
        left_interaction_part = np.kron(identity_Nm1xNm1, sigma_x)

        # Right interaction term
        right_interaction_part = np.kron(sigma_x, identity_Nm1xNm1)
        
        # Storing the number of step to compute the renormalization
        step = 0
                
        # Storing the ground state for the convergence
        levels = []
        
        while True:
            
            # Defining the current size of the system
            N_particles = ut.N_particles * (2**(step+1))
            
            # Defining the doubled hamiltonian as 1 system + 1 system + the interaction between the two
            # The Hamiltonian of the system doubled.
            doubled_full_hamiltonian = np.kron(full_hamiltonian, identity_NxN) + np.kron(identity_NxN, full_hamiltonian) +  np.kron(left_interaction_part, right_interaction_part)
            
            # Computation of eigenvalues and eigevectors
            eigvals, eigvects = np.linalg.eigh(doubled_full_hamiltonian)
            
            # Storing the ground state
            levels.append(eigvals[0]/N_particles)
            
            if i == 0:
                step_lambda_0.append(step+1)
                level_lambda_0.append(levels[step])
            elif i == binlambda:
                step_lambda_max.append(step+1)
                level_lambda_max.append(levels[step])
            
            # Checking the convergence
            if step>0:
                
                # If the ground state are similar, then convergence
                if abs(levels[step] - levels[step-1]) < precision:
                    
                    lambda_convergence.append(dlambda*i)
                    gs_convergence.append(levels[step])
                    steps_convergence.append(step+1)
                    break
                            
            # If there isn't convergence, the system will be doubled
            
            # Defining the projector
            projector = np.zeros((dimension**(2*ut.N_particles), dimension**ut.N_particles))
            
            for j in range(len(projector)):
                for k in range(len(projector[0])):
                    projector[j,k] = eigvects[j,k]
                    
            # Projecting the doubled hamiltonian on the full hamiltonian
            full_hamiltonian = ut.projection(projector, doubled_full_hamiltonian)
            
            # Projecting the left-interaction part 
            left_interaction_part = ut.projection(projector, np.kron(identity_NxN, left_interaction_part))
            
            # Projecting the right-interaction part
            right_interaction_part = ut.projection(projector, np.kron(right_interaction_part, identity_NxN))
            
            # Increase the step
            step += 1
        
        # Adding the partial time to the array
        time_spent.append(time.time()-partial) 
        
        # Printing the progress bar
        ut.progress_bar(i+1, binlambda+1, ut.expected_time(time_spent, i+1, binlambda+1))
        
    # Storing the time at the end of the computation.
    end = time.time()
    
    # Printing the time spent in the computation.
    print(f"Time spent: {int((end-start)/60)} minutes and {int((end - start) - int((end-start)/60)*60)} seconds")
    
    # Plotting the results of this computation
    fig = plt.figure()
    fig.set_figheight(9)
    fig.set_figwidth(15)
    
    ax1 = plt.subplot2grid(shape=(4,12), loc=(0,0), colspan=4, rowspan=2)
    ax1.plot(step_lambda_0, level_lambda_0)
    ax1.set_title("Step of convergence vs GS of convergence for $\lambda=0$")
    ax1.set_xlabel("Step")
    ax1.set_ylabel("GS/$N_{particles}$")
    
    ax2 = plt.subplot2grid((4,12), (0,4), colspan=4, rowspan=2)
    ax2.plot(step_lambda_max, level_lambda_max)
    ax2.set_title(f"Step of convergence vs GS of convergence for $\lambda={p_lambda}$")
    ax2.set_xlabel("Step")
    ax2.set_ylabel("GS/$N_{particles}$")
    
    ax3 = plt.subplot2grid((4,12), (0,8), colspan=4, rowspan=2)
    ax3.plot(lambda_convergence, gs_convergence)
    ax3.set_title("$\lambda$ of covergence vs GS of convergence")
    ax3.set_xlabel("$\lambda$")
    ax3.set_ylabel("GS/$N_{particles}$")
    
    ax4 = plt.subplot2grid((4,12), (2,2), colspan=4, rowspan=2)
    ax4.plot(lambda_convergence, steps_convergence)
    ax4.set_title("$\lambda$ of covergence vs Step of convergence")
    ax4.set_xlabel("$\lambda$")
    ax4.set_ylabel("Step")
    
    ax5 = plt.subplot2grid((4,12), (2,6), colspan=4, rowspan=2)
    ax5.plot(steps_convergence, gs_convergence)
    ax5.set_title("Step of convergence vs GS of convergence")
    ax5.set_xlabel("Step")
    ax5.set_ylabel("GS/$N_{particles}$")
    
    plt.subplots_adjust(hspace=0.5, wspace=15)
    
    plt.tight_layout()
    
    plt.show()

if __name__ == '__main__':
    rsrg()