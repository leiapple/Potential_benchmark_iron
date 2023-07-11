#!/usr/bin/env python
# ------original comments-----
# Fit equations of state (EOS) to energy vs volume curves and optionally plot the results
# Davide Ceresoli <dceresoli@gmail.com>, 03/22/2011
# Inspired by ase.utils.eosase2
# TODO: Keane
# ------modified (degenerated) by leizhang@rug.nl

import sys, numpy, math
from scipy.optimize import leastsq

# Murnaghan equation of state
def eos_murnaghan(params, vol):
    'From Phys. Rev. B 28, 5480 (1983)'
    E0, B0, Bp, V0 = params 
    E = E0 + B0/Bp * vol * ((V0/vol)**Bp/(Bp-1.0)+1.0) - V0*B0/(Bp-1.0)
    return E

# Birch-Murnaghan equation of state
def eos_birch_murnaghan(params, vol):
    'From Phys. Rev. B 70, 224107'
    E0, B0, Bp, V0 = params 
    eta = (vol/V0)**(1.0/3.0)
    E = E0 + 9.0*B0*V0/16.0 * (eta**2-1.0)**2 * (6.0 + Bp*(eta**2-1.0) - 4.0*eta**2)
    return E

# Birch equation of state
def eos_birch(params, vol):
    '''
    From Intermetallic compounds: Principles and Practice, Vol. I: Princples
    Chapter 9 pages 195-210 by M. Mehl. B. Klein, D. Papaconstantopoulos
    '''
    E0, B0, Bp, V0 = params 
    E = (E0 + 9.0/8.0*B0*V0*((V0/vol)**(2.0/3.0) - 1.0)**2
            + 9.0/16.0*B0*V0*(Bp-4.)*((V0/vol)**(2.0/3.0) - 1.0)**3)
    return E

# Vinet equation of state
def eos_vinet(params, vol):
    'From Phys. Rev. B 70, 224107'
    E0, B0, Bp, V0 = params 
    eta = (vol/V0)**(1.0/3.0)
    E = E0 + 2.0*B0*V0/(Bp-1.0)**2 * \
        (2.0 - (5.0 + 3.0*Bp*(eta-1.0) - 3.0*eta)*numpy.exp(-3.0*(Bp-1.0)*(eta-1.0)/2.0))
    return E


# Customized input with default and accepted values
def myinput(prompt, default, accepted):
    while True:
        res = input(prompt + " [default=%s]: " % (default))
        if res == '': res = default
        if res in accepted:
            break
        else:
            print("accepted values:", accepted)
    return res


print("Welcome to eos-fit.py")
print()

fname = 'volume.dat'
f = open(fname, 'rt')

# read data from file
print()
print("Data read from file:")
vol = []
ene = []
while True:
    line = f.readline().strip()
    if line == '': break
    if line[0] == '#' or line[0] == '!': continue
    v, e = [float(x) for x in line.split()[:2]]
    vol.append(v)
    ene.append(e)
    print(v, e)
print()
f.close()

# transform to numpy arrays
vol = numpy.array(vol)
ene = numpy.array(ene)

print()

# fit a parabola to the data and get inital guess for equilibirum volume
# and bulk modulus
a, b, c = numpy.polyfit(vol, ene, 2)
V0 = -b/(2*a)
E0 = a*V0**2 + b*V0 + c
B0 = 2*a*V0
Bp = 4.0

# initial guesses in the same order used in the Murnaghan function
x0 = [E0, B0, Bp, V0]

def print_params(label, params):
    E0, B0, Bp, V0 = params
    print(label, ": E0 = %f eV" % (E0))
    print(label, ": B0 = %f GPa" % (B0*160.21765))
    print(label, ": Bp = %f" % (Bp))
    print(label, ": V0 = %f angstrom^3" % (V0))
    print(label, ": a0 = %f angstrom" % ((V0*2)**(1/3)))
    print()

# write to Birch_Murnaghan EOS the results.txt file: ../results.txt
def write_para(label, params):
    E0, B0, Bp, V0 = params
    with open('./data/results.txt','a') as fw1:
        fw1.write("============================================ \n")
        fw1.write("%s \n" % label)
        fw1.write("E0 = %f eV\n" % (E0))
        fw1.write("B0 = %f GPa\n" %  (B0*160.21765))
        fw1.write("Bp = %f\n" %  (Bp))
        fw1.write("V0 = %f angstrom^3\n" %  (V0))
        fw1.write("a0 = %f angstrom\n" % ((V0*2)**(1/3)))
    fw1.close()

    
# fit the equations of state
target = lambda params, y, x: y - eos_murnaghan(params, x)
murn, ier = leastsq(target, x0, args=(ene,vol))
print_params("Murnaghan", murn)

target = lambda params, y, x: y - eos_birch_murnaghan(params, x)
birch_murn, ier = leastsq(target, x0, args=(ene,vol))
print_params("Birch-Murnaghan", birch_murn)
write_para("Birch-Murnaghan", birch_murn)

target = lambda params, y, x: y - eos_birch(params, x)
birch, ier = leastsq(target, x0, args=(ene,vol))
print_params("Birch", birch)

target = lambda params, y, x: y - eos_vinet(params, x)
vinet, ier = leastsq(target, x0, args=(ene,vol))
print_params("Vinet", vinet)
