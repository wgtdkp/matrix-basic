		matrix exercises to study numeric analysis

This is my numeric analysis exercises. It includes some basical 
elements about playing with matrix.

To make the code straightforward, i did not make these optimizations 
to reduce its memory. So that, you can easily understand what the 
code is doing with any numeric analysis textbook.
	
It maybe not very efficient when the input matrix is too big. But 
i think it's quick enough to give the answer of your numeric analysis 
exercise or you can use it to check if your hand work is correct.
	
HOW TO INSTALL IT
	0. clone it to your linux machine.
	1. navigate to the root of the project and type command:
		$:make
	2. type command:
		$:./mym 
       to make sure it is installed successfully.

HOW TO USE IT
	0. type command:
		$:./mym
	   or 
		%:./mym help
	   to show the usage.
    1. it supports (main element) gauss elimination and (main element) 
	   triangular decomposition.
	2. as a study tool, it supports to give out the calculation steps 
	   with a optional argument "step", like:
		%:./mym gauss .\input\m.txt step
	   will give out the A matrix and B matrix at each elimiaation step.

WHAT IS THE INPUT FILE
	The input matrix file is a simple representation of the linear equations.
	0. it starts with a integer n stands for the dimensionality.
	1. the next n line(each line includes n floating) is the coefficient matrix.
 	2. then a line includes n floating is the constant term.


