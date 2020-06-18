# Bachelor's Thesis: Higher-order Fourier analysis and Locally Decodable Codes
## Cyclotomic Polynomial Factorisation with SageMath
Read the thesis: [URL]

For my bachelor's thesis, I studied the applications of locally decodable codes to additive combinatorics and I analysed a recent refutation of a conjecture from higher-order Fourier analysis concerning dual functions and their approximations in terms of polynomial phase functions. This code constitutes my contributions to this refutation, as cyclotomic polynomial factorisation is tightly linked to the search for stronger and more general counterexamples. My results motivated the conjecture that there are counterexamples for all prime numbers strictly larger than 5.

### Abstract of thesis
[Abstract]

### Files content
#### Cyclotomic Polynomial Factorisation-v4.ipynb
- This notebook contains all code used for factoring cyclotomic polynomials, running large experiments and analysing the strength of the results. This file also contains a systematic approach of how I analysed the results and what I conjectured/concluded from these results.

#### Cyclotomic Polynomial Factorisation-v4-CODE_ONLY.ipynb
- This notebook contains the same code as in the notebook above, but without the elaborate discussion about the results. Use this file to run your own cylotomic polynomial factorisation experiments.

### cycloFact.py
- Same code as above, but without example and just a python file, for the people that don't prefer to use jupyter notebooks.

#### 'images/'
- overviewThesis.png: the image used above

#### 'csvFiles/'
- results of 30+ hours of cyclotomic polynomial factorisations (see the notebooks above for more details).

### Usage
The usage is explained in the comments of the Jupyter notebooks. There is one experiment class and one strength analysis class.
The header comments of both classes contain an elaborate explanation of the way to initialise and use these classes.
Furthermore, both notebooks contain usage examples (the complete notebook with the results contains the most examples).

### Schematic overview of thesis
![Thesis Overview](/images/overviewThesis.png)
The code in this repository concerns itself around the cyclotomic polynomial node and its link to matching vector codes.
