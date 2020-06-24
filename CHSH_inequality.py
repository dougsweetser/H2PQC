
# coding: utf-8

# # Quaternion Series QM Proof of CHSH inequalites

# Doug Sweetser <sweetser@alum.mit.edu>

# Apparently there is work by Joy Christian claiming that quaternions can be used to do Bell's inequality incorrectly. That is unfortunate. Let's just _do the right thing_.
# 
# Quaternions are a normed division algebra. A quaternion series is a vector space of quaternions with two additional pieces of information, integers for rows and columns. Quaternion series are neither normed nor a division algebra. Two quaternion series can be orthogonal, showing the norm is not preserved. There is exactly one additive identity, zero for each of the n states. There are $2^n$ multiplicative inverses for a quaternion series with each state being either zero or one. Quaternion series can be thought of as a semi-group with inverses.
# 
# The next section is cut and pasted from the notebook e91_quantum_key_distribution_protocol.ipynb and has the relations that my tools must necessarily recreate.

# ## *CHSH inequality*

# The entangled wave function used:
# 
# $$\lvert\psi_s\rangle =
#   \frac{1}{\sqrt{2}}(\lvert0\rangle_A\otimes\lvert1\rangle_B - \lvert1\rangle_A\otimes\lvert0\rangle_B) =
#   \frac{1}{\sqrt{2}}(\lvert01\rangle - \lvert10\rangle),$$
# 
# In the framework of classical physics, it is impossible to create a correlation inherent in the singlet state $\lvert\psi_s\rangle$.
# Indeed, let us measure the observables $X$, $Z$ for qubit *A* and observables $W = \frac{1}{\sqrt{2}} (X + Z)$, $V = \frac{1}{\sqrt{2}} (-X + Z)$ for qubit *B*.
# Performing joint measurements of these observables, the following expectation values can be obtained:
# \begin{eqnarray*}
#  \langle X \otimes W \rangle_{\psi_s} &= -\frac{1}{\sqrt{2}}, \quad 
#  \langle X \otimes V \rangle_{\psi_s} &= \frac{1}{\sqrt{2}}, \qquad\qquad (2) \\
#  \langle Z \otimes W \rangle_{\psi_s} &= -\frac{1}{\sqrt{2}}, \quad
#  \langle Z \otimes V \rangle_{\psi_s} &= -\frac{1}{\sqrt{2}}.
# \end{eqnarray*}
# 
# $\textbf{Exercise:}$ Given the singlet state described in the previous section, show that 
# \begin{eqnarray*}
# \langle X \otimes W \rangle_{\psi_s} &= -\frac{1}{\sqrt{2}}
# \end{eqnarray*}
# 
# Now we can costruct the *Clauser-Horne-Shimony-Holt (CHSH) correlation value*:
# 
# $$C =
# \langle X\otimes W \rangle - \langle X \otimes V \rangle + \langle Z \otimes W \rangle + \langle Z \otimes V \rangle =
# -2 \sqrt{2}. \qquad\qquad (3)$$
# 
# The [local hidden variable theory](https://en.wikipedia.org/wiki/Local_hidden_variable_theory) which was developed in an attempt to explain the quantum correlations with a classical theory gives that $\lvert C \rvert \leqslant 2$.
# But [Bell's theorem](https://en.wikipedia.org/wiki/Bell's_theorem) states that "no physical theory of local hidden variables can ever reproduce all of the predictions of quantum mechanics."
# Thus, the violation of the [CHSH inequality](https://en.wikipedia.org/wiki/Bell's_theorem#Bell_inequalities_are_violated_by_quantum_mechanical_predictions) (i.e. $C = -2 \sqrt{2}$ for the singlet state), which is a generalized form of Bell's inequality, can serve as an *indicator of quantum entanglement*.
# This fact finds its application in the E91 protocol.

# ## A complex-valued analysis of the CHSH correlation value

# I have written a library, Q_tools, to do the work with quaternions. There are two important classes:
# 1. QH, for quaternions the normed divison algebra
# 2. QHStates, for quaternion series, the semi-group with inverses
# 
# In the work of this section, although the tools can take quaternion values, the third and fourth slots of the quaternions are always equal to zero making the quaternions formally identical to complex numbers. This should show that the libraries do as they promise.
# 
# Load the needed resources.

# In[1]:


get_ipython().run_cell_magic('capture', '', '%matplotlib inline\nimport numpy as np\nimport sympy as sp\nimport matplotlib.pyplot as plt\nimport math\n\n# To get equations the look like, well, equations, use the following.\nfrom sympy.interactive import printing\nprinting.init_printing(use_latex=True)\nfrom IPython.display import display\n\n# Tools for manipulating quaternions.\nimport Q_tools as qt\nfrom IPython.core.display import display, HTML, Math, Latex\ndisplay(HTML("<style>.container { width:100% !important; }</style>"))')


# The assignment does not dictate what basis should be used for the spin states |0> and |1>. Of course there is one exceptionally popular one:

# In[2]:


q_0, q_1, q_i, q_j, q_k = qt.QH().q_0(), qt.QH().q_1(), qt.QH().q_i(), qt.QH().q_j(), qt.QH().q_k()

u = qt.QHStates([q_1, q_0])
d = qt.QHStates([q_0, q_1])

u.print_state("|u>")
d.print_state("|d>")


# To demonstrate these tools are robust, another basis will also be used that uses imaginary values.

# In[3]:


q_2 = qt.QHStates( [ qt.QH([sp.sqrt(1/2), 0, 0, 0]) ] )
q_2i = qt.QHStates([qt.QH([0, sp.sqrt(1/2), 0, 0])])

i = q_2.product(u).add(q_2i.product(d)).ket()
o = q_2.product(u).dif(q_2i.product(d)).ket()

i.print_state("|i>")
o.print_state("|o>")


# A technical note about the Q_tools... Everything is a quaternion or a series of quaternions. That makes everything look, well, the same. To tell what has happened to create a particular quaternion, I created the idea of a "qtype" which is the string that follows that has a rough record of how a particular quaternion was created. It is a little bit fun and useless when the qtype gets particularly long. The printing of qtype can be suppressed.

# Now define the four operators to be used as QHStates.

# In[4]:


# The numbers
q_1r2 = qt.QH().q_1(sp.sqrt(1/2))
q_1r2_n = qt.QH().q_1(- sp.sqrt(1/2))
q_12 = qt.QH().q_1(1/2)
q_12_n = qt.QH().q_1(-1/2)

# The operators
σ_x = qt.QHStates([q_0, q_1r2, q_1r2, q_0], qs_type="op")
σ_z = qt.QHStates([q_1r2, q_0, q_0, q_1r2_n], qs_type="op")
W = qt.QHStates([q_12, q_12, q_12, q_12_n], qs_type="op")
V = qt.QHStates([q_12, q_12_n, q_12_n, q_12_n], qs_type="op")

# Print out
σ_x.print_state("σ_x")
σ_z.print_state("σ_z")
W.print_state("W")
V.print_state("V")

σ_x.norm_squared().print_state("norm of σ_x")
W.norm_squared().print_state("norm of W")


# All the needed players are here: a way to represent the spin 2 states and the operators. The additional function one needs is called "bracket" which works like so:

# In[5]:


qt.QHStates().bracket(u, qt.QHStates().identity(2, operator=True), u).print_state("<u|I|u>")
qt.QHStates().bracket(u, σ_x, u).print_state("<u|σ_x|u>")


# Operators mix things around! 
# 
# Here is the first task:
# 
# Given a wave function:
# 
# $$\lvert\psi_s\rangle =
#   \frac{1}{\sqrt{2}}(\lvert0\rangle_A\otimes\lvert1\rangle_B - \lvert1\rangle_A\otimes\lvert0\rangle_B) =
#   \frac{1}{\sqrt{2}}(\lvert01\rangle - \lvert10\rangle),$$
#   
# show:
# 
# \begin{eqnarray*}
# \langle X \otimes W \rangle_{\psi_s} &= -\frac{1}{\sqrt{2}}
# \end{eqnarray*}

# There are eight brackets that have to be calculated:
# 1. $<0 X 0><1 W 1>$
# 1. $<1 X 1><0 W 0>$
# 1. $-<0 X 1><1 W 0>$
# 1. $-<1 X 0><0 W 1>$
# 
# These will have to be repeated a bunch of times too. Time to write a function for the task.

# In[63]:


bracket = qt.QHStates().bracket

def bracket_bracket(ket_0, ket_1, op_1, op_2, verbose=False):
    """Form the inner product for a superposition of states."""
    
    b_010 = bracket(ket_0, op_1, ket_0)
    b_111 = bracket(ket_1, op_1, ket_1)
    b_011 = bracket(ket_0, op_1, ket_1)
    b_110 = bracket(ket_1, op_1, ket_0)

    b_121 = bracket(ket_1, op_2, ket_1)
    b_020 = bracket(ket_0, op_2, ket_0)
    b_120 = bracket(ket_1, op_2, ket_0)
    b_021 = bracket(ket_0, op_2, ket_1)
    
    b_0011 = b_010.product(b_121)
    b_1100 = b_111.product(b_020)
    b_0110 = b_011.product(b_120)
    b_1001 = b_110.product(b_021)
    
    bb = b_0011.add(b_1100).dif(b_0110).dif(b_1001)
    
    if verbose:
        b_010.print_state("b_010", quiet=True)
        b_011.print_state("b_011", quiet=True)
        b_110.print_state("b_110", quiet=True)
        b_111.print_state("b_111", quiet=True)
        b_121.print_state("b_121", quiet=True)
        b_120.print_state("b_120", quiet=True)
        b_021.print_state("b_021", quiet=True)
        b_020.print_state("b_020", quiet=True)
        
        b_0011.print_state("b_0011", quiet=True)
        b_1100.print_state("b_1100", quiet=True)
        b_0110.print_state("b_0110", quiet=True)
        b_1001.print_state("b_1001", quiet=True)
    
    return bb


# In[64]:


XxW = bracket_bracket(d, u, σ_x, W)
XxW.print_state("XxW")


# Bingo, bingo. This was not a trivial calculation given the long qtype.
# 
# Demonstrate the expression is true even if the basis is switched to |i> and |o> which use imaginary values.

# In[65]:


bracket_bracket(i, o, σ_x, W).print_state("XxW", quiet=True)


# Now calculate the other 3: XxV, ZxW, and ZxV.

# In[66]:


XxV = bracket_bracket(i, o, σ_x, V)
XxV.print_state("XxV", quiet=True)

ZxW = bracket_bracket(i, o, σ_z, W)
ZxW.print_state("ZxV", quiet=True)

ZxV = bracket_bracket(i, o, σ_z, V)
ZxV.print_state("ZxV", quiet=True)


# Add 'em all up but XxV which gets subtracted.

# In[67]:


CHSH = XxW.add(ZxW).add(ZxV).dif(XxV)
CHSH.print_state("CHSH", quiet=True)
-2 * sp.sqrt(2.0)


# This is the correct answer, a good thing.

# ## Redo with quaternion-valued states

# Take the state vectors |i> and |o> and add in a non-zero $j$ and $k$. The normalization has to be tweeked so that the states remain ortho-normal.

# In[87]:


n = sp.sqrt(6)
i3 = qt.QHStates( [ qt.QH([sp.sqrt(1/2), 0, 0, 0]), qt.QH([0, 1/n, 1/n, 1/n])])
o3 = qt.QHStates( [ qt.QH([sp.sqrt(1/2), 0, 0, 0]), qt.QH([0, -1/n, -1/n, -1/n])])

i3.norm_squared().print_state("i3 norm")
o3.norm_squared().print_state("o3 norm")
o3.bra().Euclidean_product(i3).print_state("<o3|i3>")


# There is a bit of rounding error, but I will count this as a quaternion-valued ortho-normal spin 2 states. Calculate the XxV as before.

# In[88]:


XxV = bracket_bracket(i3, o3, σ_x, W)
XxV.print_state("XxV", quiet=True)


# There better not be anything special about the direction 1, 1, 1. Show a different direction also works.

# In[92]:


n = sp.sqrt(28)
i123 = qt.QHStates( [ qt.QH([sp.sqrt(1/2), 0, 0, 0]), qt.QH([0, 1/n, 2/n, 3/n])])
o123 = qt.QHStates( [ qt.QH([sp.sqrt(1/2), 0, 0, 0]), qt.QH([0, -1/n, -2/n, -3/n])])

i123.norm_squared().print_state("i123 norm")
o123.norm_squared().print_state("o123 norm")
o123.bra().Euclidean_product(i123).print_state("<o123|i123>")

XxV = bracket_bracket(i123, o123, σ_x, W)
XxV.print_state("XxV", quiet=True)


# The rounding error is worse, but that is trivial. 

# ## Why this *had* to work

# Complex numbers are a subgroup of quaternions. Nature and mathematics are logically consistent, so it follows that any expression written using complex numbers can be rewritten using quaternions that have a pair of zeros. The majority of this notebook showed the trivial double-zero quaternion process worked.
# 
# If one interprets the imaginary numbers of quaternions as a spatial thing, then one can argue on physical grounds that space is homogeneous, so it should not matter what direction one points in: $i$, $j$, $k$, or any combination of those. There are details that have to be done right, and it is easy to mess up those details. 
