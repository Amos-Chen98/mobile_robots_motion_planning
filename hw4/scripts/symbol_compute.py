'''
Author: Yicheng Chen - BUAA
Date: 2022-08-01 16:37:37
LastEditTime: 2022-08-01 18:20:20
FilePath: /hw4/scripts/symbol_compute.py
'''

# !/usr/bin/python3
# -*- coding: UTF-8 -*-
import sympy

T = sympy.Symbol('T')

px0 = sympy.Symbol('px0')
py0 = sympy.Symbol('py0')
pz0 = sympy.Symbol('pz0')

vx0 = sympy.Symbol('vx0')
vy0 = sympy.Symbol('vy0')
vz0 = sympy.Symbol('vz0')

pxf = sympy.Symbol('pxf')
pyf = sympy.Symbol('pyf')
pzf = sympy.Symbol('pzf')

vxf = 0.0
vyf = 0.0
vzf = 0.0

delta_p_x = pxf - vx0 * T - px0
delta_p_y = pyf - vy0 * T - py0
delta_p_z = pzf - vz0 * T - pz0
delta_v_x = vxf - vx0
delta_v_y = vyf - vy0
delta_v_z = vzf - vz0

alpha1 = -(12/T**3)*delta_p_x + (6/T**2)*delta_v_x
alpha2 = -(12/T**3)*delta_p_y + (6/T**2)*delta_v_y
alpha3 = -(12/T**3)*delta_p_z + (6/T**2)*delta_v_z
beta1 = (6/T**2)*delta_p_x - (2/T)*delta_v_x
beta2 = (6/T**2)*delta_p_y - (2/T)*delta_v_y
beta3 = (6/T**2)*delta_p_z - (2/T)*delta_v_z

J = T + \
    ((1/3)*alpha1**2*T**3 + alpha1*beta1*T**2 + beta1**2*T) +\
    ((1/3)*alpha2**2*T**3 + alpha2*beta2*T**2 + beta2**2*T) +\
    ((1/3)*alpha3**2*T**3 + alpha3*beta3*T**2 + beta3**2*T)

J_expand = sympy.expand(J)
J_expr = sympy.collect(J_expand, T)
print(J_expr)
print('-------------')

J_diff = sympy.diff(J, T)
# J_diff = J_diff * T**4
J_diff_expand = sympy.expand(J_diff)
J_diff_expr = sympy.collect(J_diff_expand, T)
print(J_diff_expr)
