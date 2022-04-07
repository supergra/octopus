# octopus
Exploration of solvability of various sets of equations

## The question

Is a system with a given "connectivity" between equations and variables
"likely" to be solvable?

## Motivation

In aircraft applications, you want to answer a series of interrelated questions based on the same set of equations and variables. For example, consider the range equation:

``` range = V/g * (1/SFC) * L/D * ln(Winit/Wfinal) ```

We can imagine being interested in a series of questions, depending on our goal:
- [Envelope] : What is the maximum range, given V, SFC, L/D and weight fraction?
- [Operations] : Given V, SFC, L/D, how much fuel do we need to achieve a range?
- [Design] : What weight fraction do we need to achieve to get a required range at a specific velocity V?

The answers are found by symbolically rearranging this equation. For more complex problems, there can be many interlinked equations, and the math can get quite tedious.

This code was an experiment to see if I could encode common flight or design equations, and then automatically solve them.

In the end, it was more trouble than it was worth, but it was an interesting experiment.

## Notes


- It might not be solvable due to nonlinearities and data-dependence
- It might not be solvable even if linear, due to redundant equations
- Even if not clearly solvable, there may be infinite solutions for
   "special cases".

e.g. consider a connectivity chart of:
    [1 1 1]
C = [1 1 0]
    [0 1 1]

The ones indicate non-zero dependence, but could be arbitrary sensitivities.
e.g.
x+   y +z = 0
x+0.5y    = 0
  0.5y +z = 0

[1   1   1]
[1  0.5  0]
[0  0.5  1]

--> This is unsolvable, but an *arbitrary* set of non-zero values is
    almost always solvable.


## References

Unit systems:
https://pint.readthedocs.io/en/latest/faq.html#you-mention-other-similar-python-libraries-can-you-point-me-to-those

Uncertainties:
https://pythonhosted.org/uncertainties/numpy_guide.html

Symbolic:
https://docs.sympy.org/latest/modules/numeric-computation.html
