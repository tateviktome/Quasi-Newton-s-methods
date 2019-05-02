from __future__ import division
import autograd.numpy as np
from autograd import grad
from scipy.optimize import line_search

def printMessage(algorithmName, solution, iterationCount, alphasUsed, hks):
    print algorithmName + ' solution: ', solution, '  IterationCount: ', iterationCount, '  Alphas Used: ', alphasUsed, '  matrix approximations: ', hks;

def objcFunction(x1,x2):
    # return 2*x1**2 + x2**2 + 2*x1*x2+ x1-x2;
    return np.exp(x1-1) + np.exp(1-x2) + (x1-x2)**2;
    # return 2.5*x1**2 - 3*x1*x2 + x2**2-x2+5;

def objcFunction1(params):
    # return 2*params[0]**2 + params[1]**2 + 2*params[0]*params[1]+ params[0]-params[1];
    return np.exp(params[0] - 1) + np.exp(1 - params[1]) + (params[0] - params[1]) ** 2;
    # return 2.5 * params[0] ** 2 - 3 * params[0] * params[1] + params[1] ** 2 - params[1] + 5;

def identityMatrix():
    return np.array([[1,0], [0,1]]);

def makeGradient(func, values):
    grad_func_x = grad(func, 0);
    grad_func_y = grad(func, 1);
    return np.array([grad_func_x(values[0], values[1]), grad_func_y(values[0], values[1])]);

def makeHessian(func, values):
    hessian_gradFunc_x_x = grad(grad(func, 0), 0);
    hessian_gradFunc_x_y = grad(grad(func, 0), 1);
    hessian_gradFunc_y_x = grad(grad(func, 1), 0);
    hessian_gradFunc_y_y = grad(grad(func, 1), 1);

    return np.array([[hessian_gradFunc_x_x(values[0], values[1]), hessian_gradFunc_x_y(values[0], values[1])],
                          [hessian_gradFunc_y_x(values[0], values[1]), hessian_gradFunc_y_y(values[0], values[1])]]);

def bfgs_matrix_approximation(hk, xSub, gradSub):
    rk = 1 / np.dot(np.transpose(xSub), gradSub);
    return (np.dot((identityMatrix() - rk * np.dot(xSub[np.newaxis].T, gradSub[np.newaxis])), np.dot(hk, (identityMatrix() - rk * np.dot(gradSub[np.newaxis].T, xSub[np.newaxis])))) + rk * np.dot(xSub[np.newaxis].T, xSub[np.newaxis]));

# BFGS
def bfgs_method(x0, iterationCount):
    count = 1;
    grad_objcFunc = makeGradient(objcFunction, x0);
    h0 = identityMatrix()
    xk = x0
    hk = h0

    #for reporting
    alphas = [];
    hks = [];
    while count <= iterationCount:
         searchDir = np.dot(-hk, grad_objcFunc);
         alpha = line_search(objcFunction1, grad(objcFunction1), xk, searchDir)[0]
         alphas.append(alpha);
         xKPlus1 = xk + alpha * searchDir;
         gradSub = makeGradient(objcFunction, xKPlus1) - grad_objcFunc;
         xSub = xKPlus1 - xk;
         hk = bfgs_matrix_approximation(hk, xSub, gradSub);
         hks.append(hk);
         grad_objcFunc = makeGradient(objcFunction, xKPlus1);
         xk = xKPlus1;
         count += 1;
    print printMessage('BFGS', xk, iterationCount, alphas, hks);

bfgs_method([0.0, 0.0], 4);

def dfp_matrix_approximation(hk, newSearchDir, gradSub):
    return hk + (np.dot(newSearchDir[np.newaxis].T, newSearchDir[np.newaxis]) / np.dot(np.transpose(newSearchDir), gradSub)) - (np.dot((np.dot(hk, gradSub))[np.newaxis].T, (np.dot(hk, gradSub))[np.newaxis]) / np.dot(np.transpose(gradSub), np.dot(hk, gradSub)));

def dfp_method(x0, iterationCount):
    count = 1;
    h0 = identityMatrix();
    hk = h0;
    xk = x0;

    #For reporting
    alphas = [];
    hks = [];
    grad_objcFunc = makeGradient(objcFunction, xk);
    while count <= iterationCount:
        searchDir = np.dot(-hk, grad_objcFunc);
        alpha = line_search(objcFunction1, grad(objcFunction1), xk, searchDir)[0]
        alphas.append(alpha);
        xkPlus1 = xk + alpha * searchDir;
        newSearchDir = alpha * searchDir;
        gradSub = makeGradient(objcFunction, xkPlus1) - grad_objcFunc;
        hk = dfp_matrix_approximation(hk, newSearchDir, gradSub);
        hks.append(hk);
        grad_objcFunc = makeGradient(objcFunction, xkPlus1);
        xk = xkPlus1;
        count += 1;
    print printMessage('DFP', xk, iterationCount, alphas, hks);

dfp_method([0.0, 0.0], 4);

def newton_method(x0, iterationCount):
    count = 1;
    xk = x0;
    while count <= iterationCount:
        xkPlus1 = xk - np.dot(np.linalg.inv(makeHessian(objcFunction, xk)), makeGradient(objcFunction, xk));
        xk = xkPlus1;
        count += 1;
    print printMessage('Newton', xk, iterationCount, [], []);

newton_method([0.0, 0.0], 4);

newton_method([0.0, 0.0], 30);