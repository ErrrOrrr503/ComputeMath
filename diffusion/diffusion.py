#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import subprocess

#<dimensions and data>
x_low = 0
x_high = math.pi / 8
x_num = 20

eq_x0 = x_low
eq_x1 = x_high

#<1 kind>
eq_y0 = 1
eq_y1 = 0
#</1 kind>

#<2 kind>
eq_dy0 = -math.sqrt (2) * math.sin (math.pi / 8) * math.exp (-math.pi / 8) + 0.5
eq_dy1 = (-(1 + (1 / math.sqrt (2)) * math.sin (math.pi / 8) * math.exp (-math.pi / 8)) * 2 * math.sqrt (2) * math.exp (math.pi / 4)) + math.exp (math.pi / 8) * (math.sin (math.pi / 8) + math.cos (math.pi / 8)) / 2
#</2 kind>

mistake_num_max = 100
mistake_num_min = 20

x_data = np.linspace (x_low, x_high, x_num)
#</dimensions and data>