#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math
import sys

#<dimensions and data>
x_low = 0.5
x_high = 2
x_num = 100

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

def eq4 (x_data, f_data):
    return np.exp (-f_data)

def eq4_ref (x_data):
    return np.log (x_data)

def eq4_dif (x_data, f_data):
    return -np.exp (-f_data)
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

    def mistake_graph (self, num_start, num_end):
        mistake_graphs_list = []
        mist_data_1 = np.empty (num_end - num_start)
        mist_data_2 = np.empty (num_end - num_start)
        mist_data_3 = np.empty (num_end - num_start)
        mist_data_4 = np.empty (num_end - num_start)
        for num in range (num_start, num_end):
            temp_x_data = np.linspace (x_low, x_high, num)
            temp_f_data = self.solve_euler (temp_x_data)
            temp_ref_data = self.reference_F (temp_x_data)
            max = 0
            for i in range (0, temp_x_data.size):
                if abs (temp_f_data[i] - temp_ref_data[i]) > max and temp_x_data[i] > self.x0:
                    max = abs (temp_f_data[i] - temp_ref_data[i])
            mist_data_1[num - num_start] = max
        mistake_graphs_list.append (mist_data_1)
        for num in range (num_start, num_end):
            temp_x_data = np.linspace (x_low, x_high, num)
            temp_f_data = self.solve_euler_nonlinear (temp_x_data)
            temp_ref_data = self.reference_F (temp_x_data)
            max = 0
            for i in range (0, temp_x_data.size):
                if abs (temp_f_data[i] - temp_ref_data[i]) > max and temp_x_data[i] > self.x0:
                    max = abs (temp_f_data[i] - temp_ref_data[i])
            mist_data_2[num - num_start] = max
        mistake_graphs_list.append (mist_data_2)
        for num in range (num_start, num_end):
            temp_x_data = np.linspace (x_low, x_high, num)
            temp_f_data = self.solve_euler_recalc (temp_x_data)
            temp_ref_data = self.reference_F (temp_x_data)
            max = 0
            for i in range (0, temp_x_data.size):
                if abs (temp_f_data[i] - temp_ref_data[i]) > max and temp_x_data[i] > self.x0:
                    max = abs (temp_f_data[i] - temp_ref_data[i])
            mist_data_3[num - num_start] = max
        mistake_graphs_list.append (mist_data_3)
        for num in range (num_start, num_end):
            temp_x_data = np.linspace (x_low, x_high, num)
            temp_f_data = self.solve_runge_cutt_4 (temp_x_data)
            temp_ref_data = self.reference_F (temp_x_data)
            max = 0
            for i in range (0, temp_x_data.size):
                if abs (temp_f_data[i] - temp_ref_data[i]) > max and temp_x_data[i] > self.x0:
                    max = abs (temp_f_data[i] - temp_ref_data[i])
            mist_data_4[num - num_start] = max
        mistake_graphs_list.append (mist_data_4)
        return mistake_graphs_list
            

#</equation>

#<plotting>
#</plotting>

#<main>
if __name__ == "__main__":
    if len (sys.argv) != 2:
        print ("Usage: " + sys.argv[0] + "<solution, mistake>")
        sys.exit (1)
    if not (sys.argv[1] == "solution" or sys.argv[1] == "mistake"):
        print ("Usage: " + sys.argv[0] + "<solution, mistake>")
        sys.exit (1)

    if sys.argv[1] == "solution":
        eqvs = []
        eqvs.append (difequation (eq1, 1, 1, eq1_ref, "df/dx = exp(x); f(0) = 0"))
        eqvs.append (difequation (eq2, 1, 1, eq2_ref, "df/dx = sin(x) + x^2; f(0) = 0"))
        eqvs.append (difequation (eq3, 1, 1, eq3_ref, "df/dx = sin(x) + x^2; f(0) = 0", eq3_dif))
        eqvs.append (difequation (eq4, 1, 0, eq4_ref, "df/dx = exp (-f); f(1) = 0", eq4_dif))

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
    if sys.argv[1] == "mistake":
        num_max = 100
        num_min = 20
        eq_m = difequation (eq4, 1, 0, eq4_ref, "df/dx = exp (-f); f(1) = 0", eq4_dif)
        mistake_graphs = eq_m.mistake_graph (num_min, num_max)
        mistake_n_data = np.arange (num_min, num_max)
        mistake_x_data =np.empty (mistake_n_data.size)
        for i in range (0, mistake_n_data.size):
            mistake_x_data[i] = (x_high - x_low) / mistake_n_data[i]
        plt.grid (True)
        plt.xscale ('log')
        plt.yscale ('log')
        plt.plot (mistake_x_data, mistake_graphs[0], "g", label = "euler")
        plt.plot (mistake_x_data, mistake_graphs[1], "m", label = "euler nonlin")
        plt.plot (mistake_x_data, mistake_graphs[2], "c", label = "euler recalc")
        plt.plot (mistake_x_data, mistake_graphs[3], "r", label = "runge cutt")
        k_euler = (math.log (mistake_graphs[0][num_max - num_min - 1]) - math.log (mistake_graphs[0][0])) / (math.log (mistake_x_data[num_max - num_min - 1]) - math.log (mistake_x_data[0]))
        k_nonlin = (math.log (mistake_graphs[1][num_max - num_min - 1]) - math.log (mistake_graphs[1][0])) / (math.log (mistake_x_data[num_max - num_min - 1]) - math.log (mistake_x_data[0]))
        k_recalc = (math.log (mistake_graphs[2][num_max - num_min - 1]) - math.log (mistake_graphs[2][0])) / (math.log (mistake_x_data[num_max - num_min - 1]) - math.log (mistake_x_data[0]))
        k_runge = (math.log (mistake_graphs[3][num_max - num_min - 1]) - math.log (mistake_graphs[3][0])) / (math.log (mistake_x_data[num_max - num_min - 1]) - math.log (mistake_x_data[0]))
        plt.text (mistake_x_data[num_max - num_min - 1], mistake_graphs[3][0], 'k euler = {:.6f}'.format(k_euler) + '\nk nonlin = {:.6f}'.format(k_nonlin) + '\nk recalc = {:.6f}'.format(k_recalc) + '\nk runge = {:.6f}'.format(k_runge))
        plt.legend ()
        plt.show ()
#</main>