These are the source codes used in

Cressman JR, Ullah G, Ziburkus J, Schiff SJ, Barreto E (2009) The
influence of sodium and potassium dynamics on excitability, seizures,
and the stability of persistent states: I. Single neuron dynamics. J
Comput Neurosci 26:159-70

The code for the full model is written in C and uses routines from
Numerical Recipes in C.

The code for the reduced model runs on XPP.


To run, compile with:

`gcc -o neuron_model neuron_model.c nrutil.c -lm`

and then run with:

`./neuron_model.c`

Then, plot with:

`python plotter.py`

or 

`python3 plotter.py`
