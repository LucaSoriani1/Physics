import renormalization_utils as ut
import numpy as np
import matplotlib.pyplot as plt
import time

def infinite_dmrg():
    
    # Hilbert space dimension
    dimension = ut.dimension
    
    # Size of the system = # of particles
    N_particles = ut.N_particles
    
    # Bins for the strenght (lambda) parameter
    binlambda = ut.binlambda
    
    # Max value of parameter lambda for the variation of the strenght parameter
    p_lambda = ut.p_lambda
   
    # Precision for the evaluation of the results
    precision = ut.precision
    
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
        full_hamiltonian = dlambda*i*n_interaction_hamiltonian_z + interaction_hamiltonian_x #4x4
        
        
        # Just initializing the left and right part of the full Hamiltonian.
        left_part_of_full_hamiltonian = full_hamiltonian #4x4
        right_part_of_full_hamiltonian = full_hamiltonian #4x4
        
        # Creation of the left-part of the interaction
        # LEFT_(2^N x 2^N) = 1_(2^(N-1) x 2^(N-1)) x S_x_(2x2)
        left_interaction_part = np.kron(identity_Nm1xNm1, sigma_x) #4x4

        # Creation of the right-part of the interaction
        # RIGHT_(2^N x 2^N) ) = S_x_(2x2) x 1_(2^(N-1) x 2^(N-1))
        right_interaction_part = np.kron(sigma_x, identity_Nm1xNm1)
        
        # Storing the number of step to compute the renormalization
        step = 0
                
        # Storing the ground state for the convergence
        levels = []
               
        while True:
            
            # Defining the current size of the system
            N_particles = 2*ut.N_particles + 2*(step+1)
            
            #===============================================================!
            #======================NON-INTERACTION PART=====================!
            #===============================================================!
            
            # It is stored in order to compute the projection
            # PROJECT1_(2^(N+1)  x  2^(N+1)) = H_tot_(2^N x 2^N)  x  1_(2x2)
            non_interaction_to_project_1 = np.kron(left_part_of_full_hamiltonian, np.eye(2)) #8x8
            
            # H1_(2^(2N+2) x 2^(2N+2)) = (H_tot_(2^N x 2^N)  x  1_(2 x 2)) x (1_(2x2) x 1_(2^N x 2^N))
            non_interaction_1 = np.kron(
                np.kron(left_part_of_full_hamiltonian, np.eye(2)), #8x8
                np.kron(np.eye(2), identity_NxN) #8x8
                ) #64x64
            
            # It is stored in order to compute the projection
            # PROJECT2_(2^(N+1)  x  2^(N+1)) = dl * i * (1_(2^N x 2^N) x Sigma_z)                      
            non_interaction_to_project_2 = dlambda * i * np.kron(identity_NxN, sigma_z) #8x8
            
            # H2_(2^(2N+2)  x  2^(2N+2)) = (1_(2^N x 2^N)  x  Sigma_z)  x (1_(2 x 2) x 1_(2^N x 2^N))
            non_interaction_2 = np.kron(
                dlambda*i*np.kron(identity_NxN, sigma_z), #8x8
                np.kron(np.eye(2), identity_NxN) #8x8
                ) #64x64

            # It is stored in order to compute the projection
            # PROJECT3_(2^(N+1)  x  2^(N+1)) = dl * i * (Sigma_z x 1_(2^N x 2^N))           
            non_interaction_to_project_3 = dlambda*i*np.kron(sigma_z, identity_NxN) #8x8
            
            #H3_(2^(2N+2)  x  2^(2N+2)) = (1_(2 x 2) x 1_(2^N x 2^N))  x  (Sigma_z  x 1_(2^N x 2^N))
            non_interaction_3 = np.kron(
                np.kron(np.eye(2), identity_NxN), #8x8
                dlambda * i * np.kron(sigma_z, identity_NxN) #8x8
                ) #64x64

            # It is stored in order to compute the projection
            # PROJECT4_(2^(N+1)  x  2^(N+1)) = 1_(2 x 2) x H_tot_(2^N x 2^N)          
            non_interaction_to_project_4 = np.kron(np.eye(2), right_part_of_full_hamiltonian) #8x8
            
            # H4_(2^(2N+2) x 2^(2N+2)) = (1_(2 x 2) x 1_(2^N x 2^N)) x (1_(2 x 2) x H_tot_(2^N x 2^N))            
            non_interaction_4 = np.kron(
                np.kron(np.eye(2), identity_NxN), #8x8
                np.kron(np.eye(2), right_part_of_full_hamiltonian) #8x8
                ) #64x64
            
            # =============================================================== #
            # ========================INTERACTION PART======================= #
            # =============================================================== #
            
            # It is stored in order to compute the projection
            # PROJECT12_(2^(N+1)  x  2^(N+1)) = LEFT_(2^N x 2^N) x Sigma_x 
            interaction_to_project_12 = np.kron(left_interaction_part, sigma_x) #8x8
            
           # H_int_12_(2^(2N+2) x 2^(2N+2)) = (LEFT_(2^N x 2^N) x Sigma_x_(2x2)) x (1_(2 x 1) x 1_(2^N x 2^N))
            interaction_12 = np.kron(
                                    np.kron(left_interaction_part, sigma_x), #8x8
                                    np.kron(np.eye(2), identity_NxN) #8x8
                                    ) #64
    
            # H_int_23_(2^(2N+2) x 2^(2N+2)) = (1_(2^N x 2^N) x Sigma_x) x (Sigma_x x 1_(2^N x 2^N))
            interaction_23 = np.kron(
                np.kron(identity_NxN, sigma_x), #8x8
                np.kron(sigma_x, identity_NxN) #8x8
            ) #64x64
            
            # It is stored in order to compute the projection
            # PROJECT34_(2^(N+1)  x  2^(N+1)) = Sigma_x x RIGHT_(2^N x 2^N)
            interaction_to_project_34 = np.kron(sigma_x, right_interaction_part) #8x8
            
            # H_int_34_(2^(2N+2) x 2^(2N+2)) = (1_(2 x 2) x 1_(2^N x 2^N)) x (Sigma_x x RIGHT_(2^N x 2^N))
            interaction_34 = np.kron(
                np.kron(np.eye(2), identity_NxN), #8x8
                np.kron(sigma_x, right_interaction_part) #8x8
            ) #64x64
            
            
            # Adding all the terms of the Hamiltonian together.
            new_full_hamiltonian = non_interaction_1 + non_interaction_2 + non_interaction_3 + non_interaction_4 + interaction_12 + interaction_23 + interaction_34
            
            # ============================================================== #
            # ===============DIAGONALIZATION OF H_(2N+2)==================== #
            # ============================================================== #
            
            # Finding the eigenvalues and eigenvectors of the new_full_hamiltonian matrix.
            eigvals, eigvects = np.linalg.eigh(new_full_hamiltonian)
            
            # Calculating the ground state energy of the system as the ground state dived by the number of particles.
            levels.append(eigvals[0]/N_particles)
            
            # Appending the step and level to the list.
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
                    
                    # Appending the steps, levels and the lambdas to the list.
                    lambda_convergence.append(dlambda*i)
                    gs_convergence.append(levels[step])
                    steps_convergence.append(step+1)
                    break
            # ----------------------Density matrix for the GS---------------- #
            # Defining the ket as a 2N+2 x 1 matrix
            ket = np.zeros((dimension**(2*ut.N_particles+2), 1))
            
            # Defining the bra as 1 x 2N+2 matrix
            bra = np.zeros((1, dimension**(2*ut.N_particles+2)))
            
            for j in range(dimension**(2*ut.N_particles+2)):
                bra[0, j] = eigvects[j, 0]
                ket[j,0] = eigvects[j,0]
            
            # Tensor product between the bra and the ket in order to get the density matrix
            density_matrix = np.kron(bra, ket)
            
            #---------------------REDUCED DENSITY MATRIX-------------------!
            #--------------------------------L-----------------------------!
            
            #Computing the reduced density matrix for the L-subsystem
            reduced_left_density_matrix = ut.reduced(density_matrix, 1, 8)
            
            # Diagonalization of the reduced density matrix for L-sub-system
            eigvals, eigvects = np.linalg.eigh(reduced_left_density_matrix)
            
            # Setting the small eigenvalue to zero in order to have more precision
            for j in range(len(eigvals)):
                if abs(eigvals[j]) < 1e-8:
                    eigvals[j] = +0
                    
            # Building the projector in the descending order
            projector = np.zeros((len(reduced_left_density_matrix), dimension**ut.N_particles))
    
            for j in range(len(projector)):
                for k in range(len(projector[0])):
                    projector[j, k] = eigvects[j, dimension**(ut.N_particles+1)-k-1]
                    
            # Defining the L-part of the total Hamiltonian
            left_part_of_full_hamiltonian = non_interaction_to_project_1 + non_interaction_to_project_2 + interaction_to_project_12
            
            # Projecting the L-part of the total Hamiltonian
            left_part_of_full_hamiltonian = ut.projection(projector, left_part_of_full_hamiltonian)
            
            # LEFT_(2^(N+1) x 2^(N+1)) = 1_(2^N x 2^N) x Sigma_x
            left_interaction_part = np.kron(identity_NxN, sigma_x)
            
            # Projecting the LEFT matrix
            left_interaction_part = ut.projection(projector, left_interaction_part)
            
            #--------------------------------R-----------------------------!
            
            # Computing the reduced density matrix for the R-subsystem
            reduced_right_density_matrix = ut.reduced(density_matrix, 2, 8)
            
            # Diagonalization of the reduced density matrix for L-sub-system
            eigvals, eigvects = np.linalg.eigh(reduced_right_density_matrix)
            
            # Setting the small eigenvalue to zero in order to have more precision
            for j in range(len(eigvals)):
                if abs(eigvals[j]) < 1e-8:
                    eigvals[j] = +0
            
            # Building the projector in the descending order                    
            for j in range(len(projector)):
                for k in range(len(projector[0])):
                    projector[j,k] = eigvects[j, dimension**(ut.N_particles+1)-k-1]
            # Defining the R-part of the total Hamiltonian
            right_part_of_full_hamiltonian = non_interaction_to_project_3 + non_interaction_to_project_4 + interaction_to_project_34
            
            # Projecting the R-part of the total Hamiltonian
            right_part_of_full_hamiltonian = ut.projection(projector, right_part_of_full_hamiltonian)
            
            # RIGHT_(2^(N+1) x 2^(N+1)) = Sigma_x x 1_(2^N x 2^N) 
            right_interaction_part = np.kron(sigma_x, identity_NxN)
            
            # Projecting the Right matrix
            right_interaction_part = ut.projection(projector, right_interaction_part)
            
            # Increasing the step
            step +=1
            
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
            
if __name__ == "__main__":
    infinite_dmrg()