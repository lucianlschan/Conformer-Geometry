# Conformer-Geometry

"""
SMARTS patterns for experimental torsion angle preferences
taken from RDkit torsionPreference-v2 file

torsion-angle potential form:
V = V1*(1 + s1*cos(1x)) + V2*(1 + s2*cos(2x)) + V3*(1 + s3*cos(1x))
     + V4*(1 + s4*cos(1x)) + V5*(1 + s5*cos(1x)) + V6*(1 + s6*cos(1x))

format: [SMARTS, s1, V1, s2, V2, s3, V3, s4, V4, s5, V5, s6, V6]

Total: 365 Smarts Pattern 
"""
