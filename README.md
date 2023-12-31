# Simulating the 2D Schrodinger Equation Using the Crank-Nicolson scheme

<p>
This repository contains the tools necessary to simulate the wave equation of a
particle in a box. This will be your guide for utilizing this code.
</p>

## Dependencies

<p>
For this code to run, you will need the following programs:

<ul>
  <li>GNU Make <href>https://www.gnu.org/software/make/#download</href></li>
  <li>GNU Compiler Collection <href>https://gcc.gnu.org/</href></li>
  <li>Armadillo 12.2 or greater <href>https://arma.sourceforge.net/download.html</href>
  <ul>
    <li>CMake <href>https://cmake.org/download/</href></li>
    <li>OpenBLAS <href>https://www.openblas.net/</href></li>
    <li>ARPACK <href>https://www.arpack.org/download</href></li>
    <li>SuperLU <href>https://portal.nersc.gov/project/sparse/superlu/</href></li>
  </ul>
  </li>
  <li>Python 3 or greater <href>https://www.python.org/downloads/</href></li>
  <li>NumPy <href>https://numpy.org/install/</href></li>
  <li>Matplotlib <href>https://matplotlib.org/3.8.2/index.html</href></li>
  <li>PyArma <href>https://pyarma.sourceforge.io/</href></li>
</ul>
</p>

## Building the code

<p>
To compile and link the code, simply type in </br>
  <code>make build</code></br>
in your terminal. If this worked, congratulations! You can skip to the next section, or if you are really impatient.
</br> If you get an issue with for example Armadillo when running this command, I reccommend finding it's directory, and manually link it. You can do this by adding </br>
<code>-I/directory/to/your/library/include</code>  or  <code>-L/directory/to/your/library/lib</code></br>
when linking depending on whether the library has a /include folder, og /lib folder. I recommend building and linking in one, as this removes the issue caused by having multiple versions of a library installed. You can do this by </br>
<code>g++ -o main.exe main.cpp -std=c++11 -lsuperlu -lblas -larmadillo</code></br> and changing out the library that is acting up with one of the commands given above.
</p>

## Running the code
<p>
 If you want to run a default run, meaning a run without any special settings, run the command <code>make run</code>.
  This command will run the simulation for 1, 2 and 3 slits at default settings, and then run plotter.py to plot different aspects of the waves. </br>

  To run the simulation with different settings, run the command </br>
  <code>./main.exe [h] [dt] [T] [x_c] [sigma_x] [p_x] [y_c] [sigma_y] [p_y] [v_0] [N] [y/n]</code>

  ### Parameters:
 
| Parameter name |                        Description                             |         Range        | Default value |
|:---------------|:---------------------------------------------------------------|:--------------------:|--------------:|
|    h           | The steplength, sets the spatial resolution                    | (0,1)                | 0.005         |
|    dt          | The timestep length                                            | (0, infinity)        | 2.5e-5        |
|    T           | Total time of simulation                                       | (0, infinity)        | 0.002         |
|    x_c         | Centering of the wave packet in the x-axis                     | [0,1]                | 0.25          |
|    sigma_x     | Spread of the wave packet in the x-axis                        | [0,1]                | 0.05          |
|    p_x         | momentum in the x-axis                                         | none                 | 200           |
|    y_c         | centering of the wave packet in the y-axis                     | [0,1]                | 0.5           |
|    sigma_y     | spread of the wave packet in the y-axis                        | [0,1]                | 0.2           |
|    p_y         | momentum in the y-axis                                         | none                 | 0             |
|    v_0         | the value of the potential in the slit wall                    | none                 | 1e10          |
|    N           | number of slits                                                | [-1, 2*h/slit width] | 2             |
|    y/n         | Do you want to further modify the parameters of the slit wall? | y / n                | n             |
--------------------------------------------------------------------------------------------------------------------------

  if you pick <code>y</code> on the last parameter, you will be met with </br><code>Enter slit width, slit center, wall width and wall center: </code></br>
  The next step is trivial: you fill in the respective values. All of them have a range = (0,1). The respective default values are 0.05, 0.5, 0.02 and 0.5.

  Picking a the number of slits to be -1 will remove the wall.
  </p>

</p>



## I will provide a description of the code later.
