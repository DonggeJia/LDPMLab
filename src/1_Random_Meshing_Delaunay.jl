cd(@__DIR__)
using LinearAlgebra
using DelimitedFiles
using Statistics
using CSV 
using Tables
using TetGen
using PlotlyJS
mutable struct RawTetGenIO{T}
    pointlist::Matrix
end
# read values of the Geometrical parameters

Positions, Size = pointsresult, pointsdiameterfinal
#
## Delaunay
#
input =TetGen.RawTetGenIO{Cdouble}(pointlist=Positions')
triangulation = TetGen.tetrahedralize(input,"Q")
TRI=triangulation.tetrahedronlist'


#
## 
#
n_nodes=size(Positions,1)                                                                # 6 degrees for each nodes 
gdl_n = 6;                                                                               # 6 degrees for each nodes 
gdl = gdl_n*n_nodes;                                                                     # Total Number of DOFs


#
## check volume is 0
#

n_tet = size(TRI,1);
V=ones(n_tet)
for i in 1:n_tet;
    A = Positions[TRI[i,1],:];                                                                                                      # 1st Node of the Connection
    B111 = Positions[TRI[i,2],:];                                                                                                      # 2nd Node of the Connection
    C = Positions[TRI[i,3],:];                                                                                                      # 2nd Node of the Connection
    D = Positions[TRI[i,4],:];    
    b_1 = A-B111; # Coordinates of P1 in the Reference System for the Current Point
    c_1 = A-C; # Coordinates of P2 in the Reference System for the Current Point
    d_1 = A-D; # Coordinates of P3 in the Reference System for the Current Point
    V[i] = abs(dot(b_1,cross(c_1,d_1)))/6; # mass of the Tetrahedron
end

V=sort(V)
V_total=sum(V)






