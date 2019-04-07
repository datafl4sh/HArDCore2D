#!/bin/bash
#

# Edge and cell degrees
k=1
l=1

# Boundary conditions (0=Dirichlet, 1=Neumann)
bc=0

# Test case
tcsol=1
tcdiff=2

# Solver
solver_type="bicgstab"
#solver_type="ma41"

# Meshes (without .typ2)
#mesh[1]="cart10x10"
#mesh[2]="cart20x20"
#mesh[3]="cart40x40"
#mesh[4]="cart80x80"

#mesh[1]="mesh1_2"
#mesh[2]="mesh1_3"
#mesh[3]="mesh1_4"
#mesh[4]="mesh1_5"

mesh[1]="hexa1_1"
mesh[2]="hexa1_2"
mesh[3]="hexa1_3"
mesh[4]="hexa1_4"
#mesh[5]="hexa1_5"

#mesh[1]="mesh3_2"
#mesh[2]="mesh3_3"
#mesh[3]="mesh3_4"
#mesh[4]="mesh3_5"
