# H2PQC
Jupyter Lab Tutorials for IAP2019 Class 'How to Program a Quantum Computer'.  The easiest way to access these tutorials is via the [binder image](https://mybinder.org/v2/gh/jacquescarolan/H2PQC.git/master) but if you want to install them on your own machine you can:
- Install (conda)[https://conda.io/docs/index.html]
- Create a conda environment for Qiskit and install packages with the `environment.yml` file


```
cd H2PQC
conda env create -f environment.yml
```
- initiate the environment `source activate Qiskitenv`

Will see if I can get any of my own quaternion based tools to work within
this context.

1. Class QH: does all the work needed for the normed division algebra of quaternions.
2. Class QHStates: This is a non-normed semi-group for quaternion-valued wave
   functions. Basically it is an ordered array of quaternion values. I am
   trying to think of it as a little bit more than just a vector space since
   inverses and products exist. Two quaternion series can be orthogonal, so
   that is why the norm is not perserved. For multiplication, there are $2^n$
   possible inverses because each element might be zero or non-zero having an
   identity of 0 or 1.
3. Class QHbits: Multiple QHStates. This class will deal with how to form
   products of multiple quaternion series. The idea is based off the fact that
   |a><b| for quaternion series is a tensor product, having the size
   dim(a)*dim(b). So I can imagine 

    <a|Op_a|a><b|Op_b|b>...

    yeilding a quaternion of the form: (p, 0, 0, 0). Everyone points out the **p**
    as the odds of an observer seeing an event. I like to point out the three zeros
   that show the spatial location is precisely where the observer has to be.

[Scott Aaronson in blog Why are amplitudes complex?](https://www.scottaaronson.com/blog/?p=4021)
writes:

    Namely: quantum mechanics over the quaternions is a flaming garbage fire, 
    which would’ve been rejected at an extremely early stage of God and the 
    angels’ deliberations about how to construct our universe.

He justifies the position by showing that superluminal information transfer
is possible with quaternions as players. The way he did the calcuation, that is
true so is not woth anyone's time. I hope to develop a different approach that
does not have this fatal flaw.
