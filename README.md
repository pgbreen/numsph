Numsph
============================

Calculates and returns all spherical harmonics up to a given degree (including all orders <= degree), faster than Scipy for large arrays.

Example:

```python
import numpy as np
import numsph

az = np.random.uniform(0, 2.*np.pi, 10000)
pol = np.random.uniform(0, np.pi, 10000)

#calculate shperical harmonics up to degree 4
sphdic = numsph.sph(4,az,pol)
print(sphdic[(2,4)])

[ 0.29883719+0.02571323j  0.36657912+0.11878028j -0.38797512+0.00997551j
 ...  0.19399634+0.38149024j -0.18006638-0.00831573j
 -0.01988602+0.08946939j]


# With derivative=True, two additional dictionaries are returned: derivatives with respect to azimuthal angle and polar, i.e.
sphdic,dsphaz,dsphpol = numsph.sph(4,az,pol,derivative=True)

# return associated Legendre polynomials and its derivative evaluted at x
lpdic = numsph.alp(4,x)
lpdic,dlpdic = numsph.alp(4,x,derivative=True)


```
 - Returns dictionary containing spherical harmonics up to a given degree, key tuple (order,degree)
 - Aslo calculates the derivatives of the spherical harmonics with respect to both angles, e.g. sph(l,az,pol,derivative=True)
 - Better preformance than scipy for large arrays (by about a factor 3 for a 100000 element array using l=4)
 - Associated Legendre polynomials with optional derivatives, same key i.e. (order,degree)


Notation
---------
In order to hopefully reduce some of the confusion caused by different conventions for spherical coordinates (nicely covered on [Wolfram Mathword](http://mathworld.wolfram.com/SphericalCoordinates.html)), in the code the angles are referred to as azimuthal and polar rather than adopting the physics or mathematics convention, for degree l (or li in loops) is used and for order m (or mi in loops).


Gegenbauer ultraspherical polynomials
--------------------------------------


```python
#  Gegenbauer ultraspherical polynomials up n 
gegdic = numsph.gegenbauer(n, alpha, x)
gegdic, dgegdic = numsph.gegenbauer(n, alpha, x, derivative=True)
```



