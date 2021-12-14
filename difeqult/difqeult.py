#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math

#<dimensions and data>
x_low = -5
x_high = 5
x_num = 20

x_data = np.linspace (x_low, x_high, x_num)
#</dimensions and data>

#<functions>
def eq1 (x_data, f_data):
    return np.exp (x_data)

def eq1_ref (x_data):
    return np.exp (x_data) - 1

def eq2 (x_data, f_data):
    return np.sin (x_data) + x_data * x_data

def eq2_ref (x_data):
    return -np.cos (x_data) + x_data**3 / 3 + 1

def eq3 (x_data, f_data):
    return np.exp (-np.sin (x_data)) - f_data * np.cos (x_data)

def eq3_ref (x_data):
    return x_data * np.exp (-np.sin (x_data))

def eq3_dif (x_data, f_data):
    return -np.cos (x_data)
#</functions>

#<equation>
# df/dx = F(x)
class difequation:
    def __init__ (self, F, x0, f0, reference_F, name, dif_F = None):
        self.F = F
        self.reference_F = reference_F
        self.f0 = f0
        self.x0 = x0
        self.name = name
        self.dif_F = dif_F
    
    def data_reference_F (self, x_data):
        return self.reference_F (x_data)
    
    def solve_euler (self, x_data):
        assert (isinstance (x_data, np.ndarray))
        assert (x_data[0] < self.x0  and x_data[x_data.size - 1] > self.x0)
        f_data = np.empty (x_data.size)
        i = 0
        while x_data[i] < self.x0:
            i += 1
        for i_up in range (i, x_data.size):
            if i_up == i and (x_data[i_up] - self.x0) > (x_data[i_up] - x_data[i_up - 1]) / 5:
                f_data[i_up] = self.f0 + (x_data[i_up] - self.x0) * self.F(self.x0, self.f0)
            elif i_up == i:
                f_data[i_up] = self.f0
            else:
                f_data[i_up] = f_data[i_up - 1] + (x_data[i_up] - x_data[i_up - 1]) * self.F(x_data[i_up - 1], f_data[i_up - 1]) 
        for i_down in range (i - 1, -1, -1):
            if i_down == i - 1 and (-x_data[i_down] + self.x0) > (x_data[i_down + 1] - x_data[i_down]) / 5:
                f_data[i_down] = self.f0 + (x_data[i_down] - self.x0) * self.F(self.x0, self.f0)
            elif i_down == i:
                f_data[i_down] = self.f0
            else:
                f_data[i_down] = f_data[i_down + 1] + (x_data[i_down] - x_data[i_down + 1]) * self.F(x_data[i_down + 1], f_data[i_down + 1])
        return f_data
    
    def solve_euler_recalc (self, x_data):
        assert (isinstance (x_data, np.ndarray))
        assert (x_data[0] < self.x0 and x_data[x_data.size - 1] > self.x0)
        f_data = np.empty (x_data.size)
        i = 0
        while x_data[i] < self.x0:
            i += 1       
        for i_up in range (i, x_data.size):
            if i_up == i and (x_data[i_up] - self.x0) > (x_data[i_up] - x_data[i_up - 1]) / 5:
                f_data[i_up] = self.f0 + (x_data[i_up] - self.x0) * self.F(self.x0, self.f0)
                f_data[i_up] = self.f0 + (x_data[i_up] - self.x0) * (self.F(self.x0, self.f0) + self.F(x_data[i_up], f_data[i_up])) / 2
            elif i_up == i:
                f_data[i_up] = self.f0
            else:
                f_data[i_up] = f_data[i_up - 1] + (x_data[i_up] - x_data[i_up - 1]) * self.F(x_data[i_up - 1], f_data[i_up - 1])
                f_data[i_up] = f_data[i_up - 1] + (x_data[i_up] - x_data[i_up - 1]) * (self.F(x_data[i_up - 1], f_data[i_up - 1]) + self.F(x_data[i_up], f_data[i_up])) / 2
        for i_down in range (i - 1, -1, -1):
            if i_down == i - 1 and (-x_data[i_down] + self.x0) > (x_data[i_down + 1] - x_data[i_down]) / 5:
                f_data[i_down] = self.f0 + (x_data[i_down] - self.x0) * self.F(self.x0, self.f0)
                f_data[i_down] = f_data[i_down + 1] + (x_data[i_down] - x_data[i_down + 1]) * (self.F(self.x0, self.f0) + self.F(x_data[i_down], f_data[i_down])) / 2
            elif i_down == i:
                f_data[i_down] = self.f0
            else:
                f_data[i_down] = f_data[i_down + 1] + (x_data[i_down] - x_data[i_down + 1]) * self.F(x_data[i_down + 1], f_data[i_down + 1])
                f_data[i_down] = f_data[i_down + 1] + (x_data[i_down] - x_data[i_down + 1]) * (self.F(x_data[i_down + 1], f_data[i_down + 1]) + self.F(x_data[i_down], f_data[i_down])) / 2
        return f_data
    
    def resolve_nonlinear (self, x, h, y0):
        #newton if available, secants otherwise
        y_1 = y0
        y = y_1 + 10 * h
        if self.dif_F != None:
            while abs(y - y_1) > h * h:
                y_1 = y
                y = y_1 - (y_1 - h * self.F (x, y_1) - y0) / (1 - h * self.dif_F (x,y_1))
        else:
            y_2 = y - h
            while abs(y - y_1) > h * h:
                y_1 = y
                y = y_1 - (y_1 - h * self.F (x, y_1) - y0) * (y_1 - y_2) / ((y_1 - h * self.F (x, y_1) - y0) - (y_2 - h * self.F (x, y_2) - y0))
                y_2 = y_1
        return y
    
    def solve_euler_nonlinear (self, x_data):
        assert (isinstance (x_data, np.ndarray))
        assert (x_data[0] < self.x0 and x_data[x_data.size - 1] > self.x0)
        f_data = np.empty (x_data.size)
        i = 0
        while x_data[i] < self.x0:
            i += 1       
        for i_up in range (i, x_data.size):
            if i_up == i and (x_data[i_up] - self.x0) > (x_data[i_up] - x_data[i_up - 1]) / 5:
                f_data[i_up] = self.resolve_nonlinear (x_data[i_up], x_data[i_up] - self.x0, self.f0)
            elif i_up == i:
                f_data[i_up] = self.f0
            else:
                f_data[i_up] = self.resolve_nonlinear (x_data[i_up], x_data[i_up] - x_data[i_up - 1], f_data[i_up - 1])
        for i_down in range (i - 1, -1, -1):
            if i_down == i - 1 and (-x_data[i_down] + self.x0) > (x_data[i_down + 1] - x_data[i_down]) / 5:
                f_data[i_down] = self.resolve_nonlinear (x_data[i_down], x_data[i_down] - self.x0, self.f0)
            elif i_down == i:
                f_data[i_down] = self.f0
            else:    
                f_data[i_down] = self.resolve_nonlinear (x_data[i_down], x_data[i_down] - x_data[i_down + 1], f_data[i_down + 1])
        return f_data
    
    def solve_runge_cutt_4 (self, x_data):
        assert (isinstance (x_data, np.ndarray))
        assert (x_data[0] < self.x0 and x_data[x_data.size - 1] > self.x0)
        f_data = np.empty (x_data.size)
        i = 0
        while x_data[i] < self.x0:
            i += 1       
        for i_up in range (i, x_data.size):
            if i_up == i and (x_data[i_up] - self.x0) > (x_data[i_up] - x_data[i_up - 1]) / 5:
                h = (x_data[i_up] - self.x0)
                k1 = self.F (self.x0, self.f0)
                k2 = self.F (self.x0 + h / 2, self.f0 + k1 * h / 2)
                k3 = self.F (self.x0 + h / 2, self.f0 + k2 * h / 2)
                k4 = self.F (self.x0 + h, self.f0 + k3 * h)
                f_data[i_up] = self.f0 + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
            elif i_up == i:
                f_data[i_up] = self.f0
            else:
                h = (x_data[i_up] - x_data[i_up - 1])
                k1 = self.F (x_data[i_up - 1], f_data[i_up - 1])
                k2 = self.F (x_data[i_up - 1] + h / 2, f_data[i_up - 1] + k1 * h / 2)
                k3 = self.F (x_data[i_up - 1] + h / 2, f_data[i_up - 1] + k2 * h / 2)
                k4 = self.F (x_data[i_up - 1] + h, f_data[i_up - 1] + k3 * h)
                f_data[i_up] = f_data[i_up - 1] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
        for i_down in range (i - 1, -1, -1):
            if i_down == i - 1 and (-x_data[i_down] + self.x0) > (x_data[i_down + 1] - x_data[i_down]) / 5:
                h = (x_data[i_down] - self.x0)
                k1 = self.F (self.x0, self.f0)
                k2 = self.F (self.x0 + h / 2, self.f0 + k1 * h / 2)
                k3 = self.F (self.x0 + h / 2, self.f0 + k2 * h / 2)
                k4 = self.F (self.x0 + h, self.f0 + k3 * h)
                f_data[i_down] = self.f0 + (k1 + 2 + k2 + 2 * k3 + k4) * h / 6
            elif i_down == i:
                f_data[i_down] = self.f0
            else:
                h = (x_data[i_down] - x_data[i_down + 1])
                k1 = self.F (x_data[i_down + 1], f_data[i_down + 1])
                k2 = self.F (x_data[i_down + 1] + h / 2, f_data[i_down + 1] + k1 * h / 2)
                k3 = self.F (x_data[i_down + 1] + h / 2, f_data[i_down + 1] + k2 * h / 2)
                k4 = self.F (x_data[i_down + 1] + h, f_data[i_down + 1] + k3 * h)
                f_data[i_down] = f_data[i_down + 1] + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6
        return f_data
#</equation>

#<plotting>
#</plotting>

#<main>
if __name__ == "__main__":
    eqvs = []
    eqvs.append (difequation (eq1, 0, 0, eq1_ref, "df/dx = exp(x); f(0) = 0"))
    eqvs.append (difequation (eq2, 0, 0, eq2_ref, "df/dx = sin(x) + x^2; f(0) = 0"))
    eqvs.append (difequation (eq3, 0, 0, eq3_ref, "df/dx = sin(x) + x^2; f(0) = 0", eq3_dif))

    subplot = math.ceil (math.sqrt (float (len (eqvs)))) * 10 + math.ceil (float (len (eqvs)) - math.sqrt (float (len (eqvs)))) * 100
    for i in range (0, len (eqvs)):
        plt.subplot (subplot + i + 1)
        plt.grid (True)
        plt.title (eqvs[i].name)
        plt.plot (x_data, eqvs[i].data_reference_F (x_data), "b-", label = eqvs[i].name + " reference")
        plt.plot (x_data, eqvs[i].solve_euler (x_data), "g--", label = eqvs[i].name + " euler")
        plt.plot (x_data, eqvs[i].solve_euler_recalc (x_data), "c--", label = eqvs[i].name + " euler recalc")
        plt.plot (x_data, eqvs[i].solve_euler_nonlinear (x_data), "m--", label = eqvs[i].name + " euler nonlinear")
        plt.plot (x_data, eqvs[i].solve_runge_cutt_4 (x_data), "r--", label = eqvs[i].name + " runge cutt")
        plt.plot (eqvs[i].x0, eqvs[i].f0, "ro", label = eqvs[i].name + " (x0,f0)")
        plt.legend ()
    plt.show ()
#</main>