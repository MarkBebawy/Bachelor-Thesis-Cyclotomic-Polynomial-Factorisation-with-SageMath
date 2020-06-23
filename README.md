# Bachelor Thesis Mathematics & Computer Science: Higher-order Fourier analysis and Locally Decodable Codes
## Cyclotomic Polynomial Factorisation with SageMath
Read the thesis: [URL]

For my bachelor's thesis in mathematics and computer science, I studied the applications of locally decodable codes to additive combinatorics and I analysed a recent refutation of a conjecture from higher-order Fourier analysis concerning dual functions and their approximations in terms of polynomial phase functions. This code constitutes my contributions to this refutation, as cyclotomic polynomial factorisation is tightly linked to the search for stronger and more general counterexamples. My results motivated the conjecture that there are counterexamples for all prime numbers strictly larger than 5.

### Files content
#### Cyclotomic Polynomial Factorisation-v4.ipynb
- This notebook contains all code used for factoring cyclotomic polynomials, running large experiments and analysing the strength of the results. This file also contains a systematic approach of how I analysed the results and what I conjectured/concluded from these results.

#### Cyclotomic Polynomial Factorisation-v4-CODE_ONLY.ipynb
- This notebook contains the same code as in the notebook above, but without the elaborate discussion about the results. Use this file to run your own cylotomic polynomial factorisation experiments.

#### cycloFact.py
- Same code as above, but without example and just a python file, for the people that don't prefer to use jupyter notebooks.

#### 'images/'
- overviewThesis.png: the image used above

#### 'csvFiles/'
- results of 30+ hours of cyclotomic polynomial factorisations (see the notebooks above for more details).

### Usage
The usage is explained in the comments of the Jupyter notebooks. There is one experiment class and one strength analysis class.
The header comments of both classes contain an elaborate explanation of the way to initialise and use these classes.
Furthermore, both notebooks contain usage examples (the complete notebook with the results contains the most examples).

### Abstract of thesis
In this thesis, we study locally decodable codes (LDCs), a special type of error-correcting codes, and we study their
applications to cryptography and to additive combinatorics.
We present how LDCs can be translated into private information retrieval schemes
and we study the known bounds on the length of the codeword of an LDC as a function
of the message length.
Matching vector codes are currently the best-known LDCs and these constructions
have been used in a recent paper by Jop Briët and Farrokh Labib
to refute a conjecture from additive combinatorics. Roughly, this conjecture
says that 'dual functions' from $\FF_p^n$ to $\CC$ can
be approximated well in terms of 'polynomial phase functions' in the
$L^\infty$-norm. As such, we study an application of computer science to additive
combinatorics, whereas typically, applications go in the opposite direction.
The conjecture has strong connections to a stochastic refinement
of a celebrated theorem of Szemerédi (in a finite field setting), where the common
differences of the arithmetic progressions are restricted to be in a random set.
If the conjecture were true, a lower bound (due to Altman) on the minimal sample
probability of this random set would be optimal.
%
The proof of Bri\"et and Labib relies on constructions
of matching vector codes and on the factorisation
of cyclotomic polynomials over finite fields. We give an overview of all
components needed in
their proof and we study the limitations of their proof method.
Furthermore, we build upon their counterexamples by
searching for more general ones and for stronger ones.
We created a computer program that searches for counterexamples, and we have
found these for the first 1000 primes $p$, except for $p \in \set{2, 3, 5}$.
We prove that these three cases indeed give no counterexamples.
However, as the strength of the counterexamples tends to increase as $p$ gets larger,
we conjecture that for all primes $p \geq 7$ there are counterexamples.
In particular, we have found that stronger counterexamples
emerge from factoring cyclotomic polynomials over smaller finite fields $\FF_q$,
especially when the order of $q$ in $\FF_p^\ast$ equals $\frac{p-1}{2}$.

### Schematic overview of thesis
![Thesis Overview](/images/overviewThesis.png)
The code in this repository concerns itself around the cyclotomic polynomial node and its link to matching vector codes.
