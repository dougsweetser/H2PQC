
# coding: utf-8

# # 2.1 Entangling Two Remote Spins
# 
# In this notebook we look at a protocol to entangle two remote spins.

# In[1]:


from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit, Aer, execute, BasicAer
from qiskit.tools.visualization import plot_histogram

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi


# Initialize the register (I'm going to give you this because the indexing can be a bit annoying)

# In[12]:


# Initialize Alice registers
q_alice_spin = QuantumRegister(1, 'q_alice_spin')
q_alice_photon = QuantumRegister(1, 'q_alice_photon')
c_alice_spin = ClassicalRegister(1, 'c_alice_spin')
c_alice_photon = ClassicalRegister(1, 'c_alice_photon')

# Build Alices Quantum Circuit
qc_alice = QuantumCircuit(q_alice_spin, q_alice_photon, c_alice_spin, c_alice_photon)

# Initialize Bob registers
q_bob_spin = QuantumRegister(1, 'q_bob_spin')
q_bob_photon = QuantumRegister(1, 'q_bob_photon')
c_bob_spin = ClassicalRegister(1, 'c_bob_spin')
c_bob_photon = ClassicalRegister(1, 'c_bob_photon')

# Build Bobs Quantum Circuit
qc_bob = QuantumCircuit(q_bob_spin, q_bob_photon, c_bob_spin, c_bob_photon)


# Perform the operations:

# In[13]:


# Alices Opertaions, HINT: you can call each qubit via e.g. q_alice_spin[0] in Alices and Bobs seperate quantum circuits
qc_alice.h(q_alice_spin[0])
qc_alice.cx(q_alice_spin[0], q_alice_photon[0])
qc_alice.draw()
# Bobs Operations



# join alice and bobs circuits together (HINT: you can use plain old '+' for this)  

# CX between photons

# Hadmard Alices photon


# add barrier and Draw just to check


# In[14]:


qc_bob.h(q_bob_spin[0])
qc_bob.cx(q_bob_spin[0], q_bob_photon[0])
qc_bob.draw()


# In[19]:


qc_alice_bob = qc_alice + qc_bob


# In[20]:


qc_alice_bob.cx(q_alice_photon[0], q_bob_photon[0])
qc_alice_bob.draw()


# Make measurements:

# In[7]:


# make sure the quantum registers match up with the correct classical registers

#draw just to check


# Lets look at the counts:

# In[21]:


# Load backend QasmSimulator and run the job

# select the number of shots (repeats) of the experiment, and run the job

# get the counts (how many events in each bin) + print

# plot
qc_alice_bob.measure(q, c)

# Load backend QasmSimulator and run the job
backend = BasicAer.get_backend('qasm_simulator')
job = execute(qc_alice_bob, backend, shots=1024)
result = job.result()

# get the counts (how many events in each bin)
counts = result.get_counts(qc_alice_bob)
print(counts)

# plot
plot_histogram(counts)


# Note. Basis goes {Bobs Photon, Bobs Spin, Alices Photon, Alices Spin}
# 
# **Questions:**
# - What is the state of the spins when Alice measures 0, and Bob measures 0?
# - What about 01, 10, 11 ?
