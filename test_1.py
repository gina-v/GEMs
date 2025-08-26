import iMAT
import numpy as np

import cobra
from cobra.io import load_model

model = load_model("textbook")

# expressionRxns should be a numpy array holding the (normalized?) expression
# values for a bunch of genes. Let's make a fake one!
expressionRxns = np.array([0, .1, .2, .3, .4, .5])

iMAT.iMAT(model, expressionRxns, 0, 1)
