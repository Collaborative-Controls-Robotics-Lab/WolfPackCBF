This document briefly describes each of the files:

WP_EvaderCBF
This is the current file for single-run simulations. Has decentralized CBFs as well as evader CBF

WP_Data, CoverageData, PP_Data
These files are for collecting 100 run data for the paper. Simulate 100 random runs with pursuer positions uniformly sampled in the half
ellipse.

WPTrajectories, CoverageTrajectories, PPTrajectories
Generate trajectories for the paper

WolfPACK_CentralCBF
Old code for the centralized CBF. No evader CBF

WolfPack_TypeII
My attempt at reproducing Type II results

circle, CircleIntersection, CircleMinDist, ComputeArc:
All have to do with Apollonius circles, getting function handles for circles, find intersection points of two circles, 
finding the closest point to a circle, and the arclength along a circle

closestEllipse, ellipse, ellipseData, lineEllipse
All have to do with the CRS. ellipse gives the function handle of the ellipse, ellipseData gives the values, A,B,C,D,E,F of the ellipse 
in terms of Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0, closestEllipse finds the point on the ellipse closest to a given point, and lineEllipse
finds the point on the ellipse intersecting the line closest to the agent

dxstardt, dxstardxe, dxstardxi
Computes partial derivatives of x^i*, which are relatively nasty

getpatches
Generates a patch for the pursuers and evader in terms of position and heading. Just makes simulation look pretty

quadobj, quadconstr, quadhess
Generates quadratic functions for the QCQP optimization. quadobj is the objective function, quadconstr are the constraint functions, and 
quadhess are the hessians of the constraint functions

