# Truss3DClass
PHP Class for analysing 3 dimensional pin truss by using FEM


 # Truss3DClass
 
 ## Class Inheritance
 * [FEMSolver]
 *  -> [Truss3D]
 * [FEMSolver] is base class and includes several matrix operation for the standard FEM solutions.
 * [Truss3D] is a class for the FEM solution process and include data structure of 3 dimensional pin jointed truss.
 * Solution process run for the loaded truss to analyse deformations, reactions and element forces.
 * Multiple load cases will be solved simultaneously.
 * Html result tables are generated during the process.
 * Model of 3d truss can be generated within class by assigning values to variables.
 * Or model can be created by loading CSV file.
 * Model csv file is very simple comma separated text file in-wich the properties of FEM element, boundary conditions and loads are written.

# Technical Reference

## Assumptions

1 Structure behave as a linear system.
2 Stress and strain inside of the nembers are small enough to be in the range of elastic portion.
3 Displacements of joints/nodes are small enough sothat secondary effects will be negalected.
4 Members are large enough to prevent bucklings.

## FEM model

Each element/member connected to 2 joints/nodes.
Joint has 3 degrees of freedom, ux, uy and uz.
where
ux is translation along x axis.
uy is translation along x axis.
uz is translation along x axis.

Element move linearly respectively to the end nodes.
Pin-jointed truss can be idealized as a Pintruss3D.

## Applicable field
It can solve structural mechanic problems such as bridges, transmission tower, trustal elevated tower for supporting storage tank.
Many material such as steel, alluminium and wood can be assigned with the appropriate physical properties parameter.

## Contact
 * This class is free for the educational use as long as maintain this header together with this class.
 * Author: Win Aung Cho
 * Contact winaungcho@gmail.com
 * version 1.0
 * Date: 30-9-2020
 *
