
# coding: utf-8

# # 1.0 Single Qubit Operations
# 
# In this notebook we'll explore some properties of single qubit operations
# $\newcommand{\ket}[1]{\left|{#1}\right\rangle}$
# $\newcommand{\bra}[1]{\left\langle{#1}\right|}$

# In[2]:


from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit, BasicAer, execute
from qiskit.tools.visualization import plot_histogram

import matplotlib.pyplot as plt
import numpy as np


# ### 1.1 Pauli X
# - Initialise a single quantum bit and classical register
# - confirm the properties of the X operation `x(q[0])`, what happens when you put $\ket{0}$ through, and when you put $\ket{1}$ through?

# In[4]:


# Initliaise quantum and classical register


# apply the X operation



# apply the measurement



# Set the initial state (either |0> or |1>)


# Load backend QasmSimulator and run the job


# get the counts (how many events in each bin) and print them


# plot histogram


# ### 1.2 Hadamard
# Lets look at the properties of the Hadamard gate `h(q[0])`.  What happens when you put $\ket{0}$ in, then $\ket{1}$ in?

# In[5]:


# Initliaise quantum and classical register


# Apply the H operation


# Apply the measurement


# Set the initial state (either |0> or |1>)


# Load backend QasmSimulator and run the job


# get the counts (how many events in each bin)


# plot


# Eeek! The quanutm statistics are the same.  Lets try putting in different states to see whats going on:
# $$\ket{+}=(\ket{0}+\ket{1})/\sqrt{2}$$
# $$\ket{-}=(\ket{0}-\ket{1})/\sqrt{2}$$

# In[6]:


# Set the initial state (either |+> or |->)


# Load backend QasmSimulator and run the job


# get the counts (how many events in each bin)


# plot


# So you can figure out what the gate is doing by inputting different states and making measurments.  Well done, youve just done quantum tomography!
