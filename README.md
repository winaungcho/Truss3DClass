# Truss3DClass
PHP Class for analysing 3 dimensional pin truss by using FEM


 * Truss3DClass
 *
 * Class Inheritance
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
 * 
 * This class is free for the educational use as long as maintain this header together with this class.
 * Author: Win Aung Cho
 * Contact winaungcho@gmail.com
 * version 1.0
 * Date: 30-9-2020
 *
