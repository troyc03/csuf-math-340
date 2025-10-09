# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 09:15:07 2025

@author: WINDOWS
"""

import numpy as np
import matplotlib.pyplot as plt

# Question 1

# Recursive method
def my_bisection(f, a, b, tol):
    if np.sign(f(a)) == np.sign(f(b)):
        raise Exception("The scalars a and b do not bound a root.")
    
    # Get midpoint
    m = (a+b)/2
    
    if np.abs(f(m)) < tol:
        return m
    elif np.sign(f(a)) == np.sign(f(m)):
        return my_bisection(f, m, b, tol) # Make recursive call with a = m
    elif np.sign(f(b)) == np.sign(f(m)):
        return my_bisection(f, a, m, tol) # Make recursive call with b = m
    

f = lambda x: 3*x**3+x**2-x-5
r1 = my_bisection(f, -1, 4, 0.1)
print(f"r1 = {r1}")
r2 = my_bisection(f, -1, 4, 0.01)
print(f"r2 = {r2}")
    
print(f"f(r1) = {f(r1)}")
print(f"f(r2) = {f(r2)}")

print('-----------------------------------')

# Iterative method
f = lambda x: (np.cos(x))**2-x+6
a = 6
b = 7
etol = 0.01
k = 0

# Initial check to ensure a root exists within the interval [a, b]
if f(a) * f(b) >= 0:
    print("Error: The function does not change sign within the given interval [a, b].")
else:
    while abs(b - a) / 2 > etol:  # The loop continues as long as the interval width is greater than etol
        x = (a + b) / 2
        
        # Check if x is the root or if the root lies in [a, x] or [x, b]
        if f(x) == 0:
            break # Found the exact root
        elif f(a) * f(x) < 0: # Root is in the left subinterval [a, x]
            b = x
        else: # Root is in the right subinterval [x, b]
            a = x
        k += 1

    print(f"Approximate root: {x}")
    print(f"Number of iterations: {k}")
    print(f"Final interval: [{a}, {b}]")
    
print('-----------------------------------')

# Question 4
def f(x):
    return (1-(3/4*x))**(1/3)
def f_prime(x):
    u = (1-(3/4*x))
    return (1/4*x**2)*u**(-2/3)

def newtons_method(x0, n_iter=50):
    xs = [x0]
    x= x0
    for i in range(n_iter):
        x = 4*x*(1-x)
        xs.append(x)
    return np.array(xs)

x0 = 0.6
n_iter = 50
xs = newtons_method(x0, n_iter)

for i, val in enumerate(xs):
    print(f"Iteration {i}: {i}")
    
plt.figure(figsize=(8,5))
plt.plot(range(n_iter+1), xs, marker='o', label=f"x0={x0}")
plt.axhline(0.75, color='red', linestyle='--', label="Root x=0.75")
plt.xlabel("Iteration")
plt.ylabel("x_n")
plt.title("Newton Iterations for f(x) = (1 - 3/(4x))^(1/3)")
plt.legend()
plt.grid(True)
plt.show()    

print('-----------------------------------')

def secant_method(f, x0, x1, tol=1e-6, max_iter=100):

    
    for i in range(max_iter):
        fx0 = f(x0)
        fx1 = f(x1)

        if abs(fx1) < tol:  # Check for convergence
            return x1

        if abs(fx1 - fx0) < 1e-9:  # Avoid division by zero if f(x1) approx f(x0)
            print("Warning: Division by zero risk. Try different initial guesses.")
            return None

        # Secant method formula
        x_new = x1 - fx1 * (x1 - x0) / (fx1 - fx0)

        # Update x0 and x1 for the next iteration
        x0 = x1
        x1 = x_new

    print(f"Secant method did not converge within {max_iter} iterations.")
    return None

# Example usage:
def my_function(x):
    return np.exp(x)+np.sin(x)-4

# Initial guesses
x_initial_0 = 1.0
x_initial_1 = 2.0

root = secant_method(my_function, x_initial_0, x_initial_1)

if root is not None:
    print(f"The approximate root is: {root}")
    print(f"f({root}) = {my_function(root)}")
        
    

