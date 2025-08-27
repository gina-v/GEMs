import numpy as np
from scipy.sparse import lil_matrix

    # define a function called iMAT that returns a tissueModel and has the input arguments: model, expressionRxns, threshold_lb, threshold_ub, tol=None, core=None, logfile='MILPlog', runtime=60, epsilon=1
def iMAT(model, expressionRxns, threshold_lb, threshold_ub, tol=None, core=None, logfile='MILPlog', runtime=60, epsilon=1):
    # A docstring to describe the function when help(iMAT) is called.
    """
    Returns a context-specific tissue model using the iMAT algorithm.
    """
    # the result of the function is a variable called tissueModel.

    # CTB: tissueModel is not defined, and you don't want to return immediately
    # anyway. Note to self @CTB: look at the matlab code to grok why this is here.
    #CTB return tissueModel # ** what is tissueModel? It's not an object and not assigned to a class...

# check if the model has a field called C or E. If so, return a warning using hasattr (check if it has attribute)
# ** what is C and E? ** 
    if hasattr(model, 'C') or hasattr(model, 'E'):
        print("Warning: iMAT does not handle the additional constraints and variables defined in the model structure (fields 'C' and 'E'). It will only use the stoichiometry provided.")

# check if the function has the attribute. Python can check by name whereas matlab checks by the attribute's position. 
    if epsilon is None:
        epsilon = 1

    if runtime is None:
        runtime = 60

    if logfile is None:
        logfile = 'MILPlog'

    if core is None:
        core = []

    if tol is None:
        tol = 1e-8

        # create two numpy arrays of the positions of the reactions that are above or below a set threshold. 
        # numpy.where(condition, [x, y, ]/)
        # np.where is for multdimensional arrays. End in [0] to return the first element of the tuple with the indices. 
        RHindex = np.where(expressionRxns >= threshold_ub)[0]
        RLindex = np.where(expressionRxns <= threshold_lb)[0]

        # Ensure the predefined list of core reactions is always included in the model even if they have low expression. 
        # Three steps: 
        # 1. loop through a list of core reactions that that we want to keep
        # 2. Add them to the model manually if they aren't already in the RHindex
        # 3. Check if any core reaction isn't in the model 
    if core: 
        for rxn in core:
            if rxn in model.rxns:
                rxn_location = model.rxns.index(rxn)
                if rxn_location not in RHindex:
                    RHindex.append(rloc) # ** why would rloc be in RHindex? They're both lists of positions? **
                elif rxn_location not in model.rxns:
                    print(f"Warning: core reaction {rxn} not found in the model reactions.")

    # Purpose: Store S matrix in easy to use variable, extract reaction upper bounds and lower bounds. 
    #S = model.S # model is an object of the cobra.Model class. model.s is an attribute of the model object that stores the SciPy sparse stoichiometric matrix.
    # CTB: I _think_, looking at https://cobrapy.readthedocs.io/en/latest/getting_started.html#Reactions, that Python's model.reactions is the equivalent of model.S.
    S = model.reactions

    # concept: list comprehension
    # rxn.lower_bound is an attribute of the cobra.Reaction class that stores the lower bound of a reaction and vice versa. 
    lb = [rxn.lower_bound for rxn in model.reactions] # list of lower bounds for each reaction made by looping through the rxns in model.reactions list
    ub = [rxn.upper_bound for rxn in model.reactions] # list of upper bounds for each reaction  made by looping through the rxns in model.reactions list
    # Convert lb and ub to numpy arrays as a faster alternative: 
        # lb = np.array([rxn.lower_bound for rxn in model.reactions])
        # ub = np.array([rxn.upper_bound for rxn in model.reactions])

    # CTB note: n_rows and n_cols need to exist BEFORE they are used to
    # create the lil_matrix A.
    # Purpose: Calculate the number of rows and columns in the constraint matrix A.
            # A contains all equations and inequalties that:
                # Enforce mass balance in S matrix
                # Enforce high expression reactions to carry flux 
                # Enforce low expression reactions to carry no flux

    # CTB: need to define n_mets. Presumably number of metabolites?
    # CTB: need to define n_rxns. Presumably number of reactions?
    n_mets = len(model.metabolites)
    n_rxns = len(model.reactions)

    n_rows = n_mets + 2 * len(RHindex) + 2 * len(RLindex)
    # n_mets: the row for mass balance equations in the S matrix (number of metabolites)
    # 2 rows of constraints per each highly expressed reaction (enforce activation)
    # 2 rows of constraints per each lowly expressed reaction (enforce suppression)
    n_cols = n_rxns + 2 * len(RHindex) + len(RLindex)
    # n_rxns: one column for each reaction flux 
    # 2 columns of binary variables (like 0 and 1) to enforce activation of highly expressed reactions
    # 1 column of binary variables to enforce suppression of lowly expressed reactions

    # Purpose: Get size of the S matrix
    # Concept: Tuple unpacking - a shortcut to unpack a tuple to multiple variables at once 
    # S.shape = (number of rows, number of columns)
    # rows = metabolites, columns = reactions
    # Concept: Use a "List of Lists", LIL, format to make a sparse matrix because it is efficient for row-by-row assignment
    # ** Look into the difference between a sparse and dense matrix. Why are most of the values zero? What are the values?
    # CTB: a sparse matrix is one where most of the elements are 0, and such
    # matrices can be stored more efficiently than a dense matrix.
    # CTB: the values are (probably) interactions between
    # CTB: n_rows and n_cols are the number of rows and columns here.
    A = lil_matrix((n_rows, n_cols)) # ** why do you have to specify the names of the rows and columns? Why can't it just be A = lil_matrix()??
    n_RH = len(RHindex)
    n_RL = len(RLindex)
    # Copy S into the top-left corner of A
    # Concept: Sparce matrix slicing
    # CTB: here we need to figure out 'S' better :)
    A[:n_mets, :n_rxns] = S 
    # take the first n_mets rows and first n_rxns columns of A and fill this block with the values from S. 

# More info on A: 
    # The matrix A is used in the constraint: A * x = b
    # where A is the constraint coefficient matrix, x is a vector of variables, and b is a vector of bounds. 
    # ** make a visual diagram of what the A matrix looks like. 
# How to assign a value to a cell in a matrix: A[row, col] = value

    # create constraints for highly expressed reactions
    for i, rxn in enumerate(RHindex): # counter, i, tells you how many times you've looped through or position. rxn is the value at that position
        A[n_mets + i, rxn] = 1 # For each reaction in RHindex, go one row below n_mets each time and place a 1 in the correct column.
        A[n_mets + i, n_rxns + i] = lb[rxn] - epsilon # In the row for the minimum flux constraint of reaction rxn, set the coefficient for the binary variable y_i equal to the reaction’s lower bound minus epsilon.
        A[n_mets + n_RH + i, rxn] = 1
        A[n_mets + n_RH + i, n_rxns + len(RHindex) + n_RL + i] = ub[rxn] + epsilon

    for i, rxn in enumerate(RLindex):
        A[n_mets + i + 2 * n_RH, rxn] = 1
        A[n_mets + i + 2 * n_RH, n_rxns + i + n_RH] = lb[rxn]
        A[n_mets + i + 2 * n_RH + n_RL, rxn] = 1
        A[n_mets + i + 2 * n_RH + n_RL, n_rxns + i + n_RH] = ub[rxn]
    

    # Creating csense
    # Concept: constraint sense vector
    csense1 = ['E'] * n_mets
    csense2 = ['G'] * n_RH
    csense3 = ['L'] * n_RH
    csense4 = ['G'] * n_RL
    csense5 = ['L'] * n_RL

    csense = csense1 + csense2 + csense3 + csense4 + csense5

    # Creating lb and ub
    lb_y = np.zeros(2*n_RH + n_RL)
    ub_y = np.ones(2*n_RH + n_RL)
    lb = np.concatenate([lb, lb_y])
    ub = np.concatenate([ub, ub_y])

    # Creating c
    n_rxns = S.shape[1]
    c_v = np.zeros(n_rxns)
    c_y = np.ones(2*n_RH + n_RL)

    c = np.concatenate([c_v, c_y])

    # Creating b
    n_mets = S.shape[0]

    b_s = np.zeros(n_mets)
    lb_rh = lb[RHindex]
    ub_rh = ub[RHindex]
    lb_rl = lb[RLindex]
    ub_rl = ub[RLindex]

    b = np.concatenate([b_s, lb_rh, ub_rh, lb_rl, ub_rl])

    # Creating vartype - 
    MILPproblem = {
        "A": A,               # scipy.sparse OK
        "b": b,               # numpy 1D array
        "c": c,               # numpy 1D array
        "lb": lb,             # numpy 1D array
        "ub": ub,             # numpy 1D array
        "csense": csense,     # list of 'E','G','L'
        "vartype": vartype,   # list of 'C','B'
        "osense": -1,         # maximize
        "x0": None,
}


        # meeting with Titus 08.07.25 next meeting 08.21.25
        # code stability test - test that when you run this multiple times you get the same result
        # make sample files to test this code and check if the outputs are what's expected, do this a fedw times with fake files
        # what is a solver wrapper 
    

    # put on github
    # look at cobrapy that has an algorithim or just the documentation
    # imat paper 

    # get the syntax working 
