
# coding: utf-8

# $\newcommand{\ket}[1]{\left|{#1}\right\rangle}$
# $\newcommand{\bra}[1]{\left\langle{#1}\right|}$
# # 0.1 H2PQC
# 
# The following describes the basic functionality of Qiskit and how to:
# 
# (1) Initialize qubits and classical registers
# 
# (2) Perform operations
# 
# (3) Make measurements
# 
# In addition we cover different simulation backends and plotting.  Lots of this is also in the IBM Qiskit Terra tutorial.

# ### Code Import

# In[1]:


from qiskit import QuantumRegister, ClassicalRegister
from qiskit import QuantumCircuit, Aer, execute, BasicAer
from qiskit.tools.visualization import plot_histogram

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi


# ### Initialize quanutm and classical registers

# In[2]:


# quantum circuit with 2 qubits initliased in the |0> state
q = QuantumRegister(2)

# classical register for readout with 2 bits of representation
c = ClassicalRegister(2)


# ### Build the quantum circuit and apply some operations
# 
# A list of all the availible single and multi-qubit operations can be found here: https://qiskit.org/documentation/summary_of_quantum_operations.html.  For the moment lets just use the ones we've learnt about so far, let's use Y, Z and this mystery H.

# In[3]:


# build the circuit
qc = QuantumCircuit(q, c)

# apply Y to qubit 0
qc.y(q[0])

# apply Z to qubit 1
qc.z(q[1])

# apply Hadamard to qubit 0
qc.h(q[0])

# apply Hadamard to qubit 1
qc.h(q[1])


# ### Visualisation
# 
# To check we applied the correct operations we can draw it!

# In[4]:


qc.draw()


# ### Backend
# 
# Once we have the circuit we need to choose our backend.  This will determine how we run the circuit.  You can either do this numerically (so you output the state vector or unitary exactly), simulate a quantum computer (so you output fake, noisey, data) or run it on a real quanutm computer.  For our purpose we will moslty be generating fake data, but lets check out the other two:

# In[5]:


BasicAer.backends()


# #### Numerical Simulation (statevector)

# In[6]:


# load the backend
backend = BasicAer.get_backend('statevector_simulator')

# create a program to run
job = execute(qc, backend)

# get the results
result = job.result()

# get the statevector and print it
outputstate = result.get_statevector()
print(outputstate)


# Lets plot the probabilities $|a_{ij}|^2$ of each element in the basis:

# In[7]:


plt.bar(range(len(outputstate)), np.abs(outputstate)**2)
plt.xticks(range(len(outputstate)), ['|00>','|01>', '|10>', '|11>' ])
plt.ylabel('probabitiy')
plt.xlabel('state')
plt.show()


# #### Numerical Simulation (unitary)

# In[8]:


# load the backend
backend = BasicAer.get_backend('unitary_simulator')

# create a program to run
job = execute(qc, backend)

# get the results
result = job.result()

# Print the unitary as an array
print(result.get_unitary())


# #### Draw Samples
# 
# We're mostly going to be interested in the output from a quanutm computer where we don't find out what the full state is, rather, we get a single output with probaility goverend by the respective probability amplitude.  Forutnately Qiskit let's us do that!

# In[9]:


# we first need to measure the qubits
qc.barrier(q)
qc.measure(q, c)

# lets check this
qc.draw()


# In[11]:


# load the QASM backend
backend = BasicAer.get_backend('qasm_simulator')

# select the number of shots (repeats) of the experiment, and run the job
job = execute(qc, backend, shots=1000)
result = job.result()

# get the counts (how many events in each bin)
counts = result.get_counts(qc)

print(counts)
# plot
plot_histogram(counts)


# Yay! Your first quanutm data!

# ### Different input states
# 
# You can also change the input state, **note** that for two qubits the input array is in binary ordering:
# 
# $\{\ket{00}, \ket{01}, \ket{10}, \ket{11}\}$
# 
# with the indexing $\ket{q_{1} q_0}$ (qubit 1, qubit 0)
# 

# In[12]:


# Build Circuit
q = QuantumRegister(2)
c = ClassicalRegister(2)
qc = QuantumCircuit(q, c)

# apply X to qubit 0
qc.x(q[0])

# apply Y to qubit 1
qc.y(q[1])

# Measure the qubit
qc.measure(q, c)

# Set the initial state, two qubit basis is [00, 01, 10, 11]
opts = {"initial_statevector": np.array([0,1,0,0])}

# Load backend QasmSimulator and run the job
backend = BasicAer.get_backend('qasm_simulator')

# select the number of shots (repeats) of the experiment, and run the job
job = execute(qc, backend, shots=1024, backend_options=opts)
result = job.result()

# get the counts (how many events in each bin)
counts = result.get_counts(qc)
print(counts)

# plot
plot_histogram(counts)

