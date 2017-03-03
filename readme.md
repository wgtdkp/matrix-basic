## matrix exercises to study numeric analysis ##

This is my numeric analysis exercises. It includes some basical 
elements about playing with matrix.

To make the code straightforward, i did not make these optimizations 
to reduce its memory usage. So, you can easily understand what the 
code is doing with any numeric analysis textbook.
	
It maybe not very efficient when the input matrix is too big. But 
i think it's quick enough to give the answer of your numeric analysis 
exercise or you can use it to verify your hand work.

Notes that when the numbers come to be very big, these is some 
calculation error because of the machine representation of floating number.

**HOW TO INSTALL IT**
 * clone it to your linux machine.
	
	```$ git clone https://github.com/wgtdkp/matrix-basic.git```
 * navigate to the root of the project and type command: 
	
	 ```$ make```
 * type command: 		
	
	```$ ./mym ```

    to make sure it is installed successfully.

**HOW TO USE IT**
 * type command:
    
    ```$ ./mym```
	
	or
	
    ```$ ./mym help```

    to show the usage.
 * it supports (main element) gauss elimination and (main element) 
   triangular decomposition.
 * as a study tool, it supports to give out the calculation steps 
	 with an optional argument "step", like:

	```$ ./mym gauss ./input/m.txt step```
	
	 will give out the A matrix and B matrix at each elimiaation step.

**WHAT IS THE INPUT FILE**<br />
The input matrix file is a simple representation of the linear equations.

 * it starts with a integer n stands for the dimensionality.
 * the next n line(each line includes n floating) is the coefficient matrix.
 * then a line includes n floating is the constant term.
