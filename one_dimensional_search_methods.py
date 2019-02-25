from __future__ import division
import numpy as np
import math

def f1(x):
    return np.log(x * x - 1.6 * x + 3)
def f2(x):
    return (x - 1) * math.exp(-2 * x * x + 4 * x)
def f1prime(x):
    return (2 * x - 1.6) / (x * x - 1.6 * x + 3)
def f2prime(x):
    return math.exp(-2 * x * x + 4 * x) + math.exp(-2 * x * x + 4 * x) * (-4 * x * x + 8 * x - 4)
def f1primeprime(x):
    return (2 * (x * x - 1.6 * x + 3) - (2 * x - 1.6) * (2 * x - 1.6)) / (x * x - 1.6 * x + 3) * (x * x - 1.6 * x + 3)
def f2primeprime(x):
    return math.exp(-2 * x * x + 4 * x) * (-4 * x + 4) + math.exp(-2 * x * x + 4 * x) * (-4 * x + 4) * (-4 * x * x + 8 * x - 4) + math.exp(-2 * x * x + 4 * x) * (-8 * x + 8)

def bisection_method(f, fprime, a, b, N):
    while N > 0:
        midpoint = (a + b) / 2
        if fprime(midpoint) > 0:
            b = midpoint
        elif fprime(midpoint) < 0:
            a = midpoint
        else:
            break
        N = N - 1
    print 'Minimizer is', (a+b)/2

bisection_method(f1,f1prime,0,3,2)
bisection_method(f1,f1prime,0,3,4)
bisection_method(f1,f1prime,0,3,8)
bisection_method(f1,f1prime,0,3,16)

bisection_method(f2,f2prime,0,3,2)
bisection_method(f2,f2prime,0,3,4)
bisection_method(f2,f2prime,0,3,8)
bisection_method(f2,f2prime,0,3,16)

bisection_method(f2,f2prime,-1,1,2)
bisection_method(f2,f2prime,-1,1,4)
bisection_method(f2,f2prime,-1,1,8)
bisection_method(f2,f2prime,-1,1,16)

def newtons_method(f, fprime, fprimeprime, a, b, n):
    x0 = (a + b)/2
    while n > 0:
        x0 = x0 - fprime(x0) / fprimeprime(x0)
        n = n - 1
    print 'Minimizer is', x0

newtons_method(f1,f1prime,f1primeprime,0,3,2)
newtons_method(f1,f1prime,f1primeprime,0,3,4)
newtons_method(f1,f1prime,f1primeprime,0,3,8)
newtons_method(f1,f1prime,f1primeprime,0,3,16)

newtons_method(f2,f2prime,f2primeprime,0,3,2)
newtons_method(f2,f2prime,f2primeprime,0,3,4)
newtons_method(f2,f2prime,f2primeprime,0,3,8)
newtons_method(f2,f2prime,f2primeprime,0,3,16)

newtons_method(f2,f2prime,f2primeprime,-1,1,2)
newtons_method(f2,f2prime,f2primeprime,-1,1,4)
newtons_method(f2,f2prime,f2primeprime,-1,1,8)
newtons_method(f2,f2prime,f2primeprime,-1,1,16)

def fibonacci(n):
    if n < 1:
        return 0
    elif n == 1:
        return 1
    else:
        return fibonacci(n - 1) + fibonacci(n - 2)

def fibonacci_search(f, initiala, initialb, n):
    a1 = initiala
    b1 = initialb
    while n > 0:
        x0 = fibonacci(n) / fibonacci(n + 1)
        a2 = a1 + (1 - x0) * (b1 - a1)
        b2 = a1 + x0 * (b1 - a1)
        if f(a1) < f(b1):
            b1 = b2
        else:
            a1 = a2
        n = n - 1
    print 'Minimizer is', (a1 + b1) / 2

fibonacci_search(f1,0,3,2)
fibonacci_search(f1,0,3,4)
fibonacci_search(f1,0,3,8)
fibonacci_search(f1,0,3,16)

fibonacci_search(f2,0,3,2)
fibonacci_search(f2,0,3,4)
fibonacci_search(f2,0,3,8)
fibonacci_search(f2,0,3,16)

fibonacci_search(f2,-1,1,2)
fibonacci_search(f2,-1,1,4)
fibonacci_search(f2,-1,1,8)
fibonacci_search(f2,-1,1,16)

def golden_section_search(f,initiala,initialb,n):
    a1 = initiala
    b1 = initialb
    while n > 0:
        a2 = a1 + (3 - np.sqrt(5))/2 * (b1 - a1)
        b2 = a1 + (-1 + np.sqrt(5))/2 * (b1 - a1)
        if f(a1) < f(b1):
            b1 = b2
        else:
            a1 = a2
        n = n - 1
    print 'Minimizer is', (a1 + b1) / 2

golden_section_search(f1,0,3,2)
golden_section_search(f1,0,3,4)
golden_section_search(f1,0,3,8)
golden_section_search(f1,0,3,16)

golden_section_search(f2,0,3,2)
golden_section_search(f2,0,3,4)
golden_section_search(f2,0,3,8)
golden_section_search(f2,0,3,16)

golden_section_search(f2,-1,1,2)
golden_section_search(f2,-1,1,4)
golden_section_search(f2,-1,1,8)
golden_section_search(f2,-1,1,16)
