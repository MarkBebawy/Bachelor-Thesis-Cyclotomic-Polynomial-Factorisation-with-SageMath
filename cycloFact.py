from sage.rings.polynomial.cyclotomic import cyclotomic_coeffs
import matplotlib.pyplot as plt
import csv
import os

# -------------------------- CYCLOTOMIC FACTORISATION CLASS --------------------------
class CycloFactExperiment:
    """This class allows the execution and analysis of numerical experiments concerning the factorisation
    of p-th cyclotomic polynomials (p prime) over finite fields (see docstring of __init__ function
    for details).
    Dependency:
        - SageMath (https://www.sagemath.org/)
    Usage:
        - Create a folder called 'csvFiles' in the same directory as this file.
        - Initialise a cycloFactExperiment object:
            cfexp = CycloFactExperiment(fileName, read, numPrimes, primeList, numFields, fieldList,
            showTrace) see docstring of __init__ for the meaning of these parameters.
        - Call cfexp.runExperiment()
            if read=False, this will execute the calculations and create two files in the
            folder csvFiles/
                - fileNameSPARSEST.csv, this is a list of the counterexamples with sparsest
                    irreducible factors for each prime
                - fileNameALL.csv, this is a list of all counterexamples
            if read=True, this will read the data from the two
            files in csvFiles/ (the files must exist)

            The rows in each file consist of (p, qMin, minSupp, order, bool) where the {p}-th
            cyclotomic polynomial is factored over FF_{qMin} into irreducible factors, one of
            which was a sparse factor, namely sparsity {minSupp}. Furthermore, {order} is
            ord_{p}({qMin}) and {bool} is True if and only if this is a valid counterexample
            (i.e., if {minSupp} <= ord_{p}({qMin})).
        - cfexp.plotSparsest() and cfexp.plotAll() shows a graph of {qMin} plotted against {p} for the
          sparsest (respectively all) counterexamples.
        - cfexp.plotStrength(ylim, strengthFunc, ylabel)
          plots the strength of all counterexamples against the primes p. The definition of strength is
          determined by strengthFunc, which should be a function of the order, the sparsity and
          the prime p. Default is strengthFunc=(lambda order, spars, p: order - spars). The parameter 
          ylabel should be the same formula in words, to be printed as the ylabel
          (default '$ord_p(q)$ - sparsity'). ylim is optional (default None) and can be used to set
          the yrange of the plot.
        - cfexp.searchHighOrder(read, fileName, primeList) tests for every prime p in primeList
          whether the first prime q encountered with ord_p(q) = (p - 1)/2 entails a counterexample.
          It returns a list of counterexamples and writes that to csvFiles/fileNameOrderTest.csv (or reads
          counterexamples from that file and returns the list if read==True).
    """
    def __filterMersenne(self, primeList):
        """This function removes all Mersenne primes from a list of prime numbers (primeList).
        (This actually removes all numbers of the form 2^t - 1 from a list of integers.)"""
        return [p for p in primeList if not log(p + 1, 2).n() in ZZ]

    
    def __init__(self, fileName='ram', read=True, numPrimes=100, primeList=None, numFields=10,
                 fieldList=None, showTrace=True):
        """
            Initialise a cycloFact object to conduct numerical experiments. For every prime p in primeList,
            this class allows to factor the p-th cyclotomic polynomial into irreducible factors
            over finite fields, with order q specified by fieldList, and store the q that gives the sparsest
            factor with sparsity at most the multiplicative order of q modulo p. A table is generated that
            stores the 'optimal' q for every p, with the corresponding sparsity, as well as a table that
            stores all factorisations with sparsity <= order. For p Mersenne, the behaviour
            is understood well, hence these primes are filtered out from the list of primes.

            Parameters:
            - fileName: name of the csv-file in which data will be stored or from which data will be read.
                        NOTE: the files are stored in and read from the folder 'csvFiles/'. Make sure this
                        folder exists as a subfolder in the folder this file is in. The extension '.csv' is
                        added automatically. Thus, if fileName='ram', data will be stored in (or read from)
                        csvFiles/ramALL.csv and csvFiles/ramSPARSEST.csv.
                        When writing, the file does not have to exist already,
                        but the folder 'csvFiles/' does have to exist. When reading, the files must exist.
            - read: boolen indicating whether to read from or to write to the file specified by fileName.
            - numPrimes: if no primeList is given, the first 'numPrimes' primes will be considered for p.
            - primeList: list of primes to consider for p (if this is not None, this will ignore numPrimes).
            - numFields: if no fieldList is given, the first 'numFields' primes will be considered for q.
            - fieldList: list of primes to consider for q (if this is not None, this will ignore numFields).
            - showTrace: prints which steps are being done during experiments when set to True.
        """
        # Iterable list of prime numbers
        self.primes = Primes()
        self.read = read
        self.fileName = f"csvFiles/{fileName}"
        
        # Filter out Mersenne primes, as these are well-understood cases
        self.primeList = [self.primes.unrank(i) for i in range(numPrimes)] if primeList is None else primeList
        self.primeList = self.__filterMersenne(self.primeList)
        self.fieldList = [self.primes.unrank(i) for i in range(numFields)] if fieldList is None else fieldList
        self.showTrace = showTrace
    
        self.bestCounterExmpls = list()
        self.allCounterExmpls = list()

    
    def multOrder(self, p, r):
        """Calculate ord_p(r), (i.e. the multiplicative order of r in (Z/pZ)^*) or return p + 1
        when this is not well-defined."""
        Zp = Integers(p)
        order = p + 1 if Zp(r) == Zp(0) else Zp(r).multiplicative_order()
        return order


    def support(self, pol):
        """This function calculates the support of a polynomial (pol).
        That is, it counts the number of non-zero coefficients of pol.
        This function assumes pol is a polynomial in one variable."""
        return len(pol.coefficients())
    
    
    def write_csv(self, fileName, counterExmpls):
        """Store list of four-tuples (counterExmpls) in a csv file
        called 'csvFiles/fileName.csv'."""
        with open(fileName + '.csv', 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            # Add titles
            titles = ['p-th Cyclo pol', '#F_q (pol factored over F_q)', 'sparsity',
                      'ord_p(q)', 'sparsity <= ord_p(q)?']
            writer.writerow(titles)
            for (p, qMin, minSup, order) in counterExmpls:
                writer.writerow([p, qMin, minSup, order, minSup <= order and order < p])


    def read_csv(self, fileName):
        """Loads counterexamples from csv file called 'fileName.csv' and returns
        result as list of tuples."""
        sparsity = list()

        with open(fileName + '.csv', newline='') as csvfile:
            reader = list(csv.reader(csvfile))

            # Skip title row
            for row in reader[1:]:
                sparsity.append((int(row[0]), int(row[1]), int(row[2]), int(row[3])))

        return sparsity
    
    
    def factorCyclo(self, pCyclo, qFactor):
        """Factorise the {pCyclo}-th cyclotomic polynomial
        over the finite field with {qFactor} elements.

        Return value (three tuple):
        - qFactor: same as input, number of elements in the field
        - factorsSupp: list of tuples (fac, supp), where fac is an irreducible factor
            of the pCyclo-th cyclotomic polynomial and supp is the support size of fac.
        - minSupp: the support size of the sparsest irreducible factor.
        """
        # Only prime numbers are allowed
        assert(pCyclo in self.primes)
        assert(qFactor in self.primes)
        
        # Create field and polynomial ring
        field = GF(qFactor)
        polRing = PolynomialRing(field, 'x')
        x = polRing.gen()
        
        # Create and factorise cyclotomic polynomial. The list cycloFactor consists of
        # tuples (pol, pow) where pol is an irreducible polynomial and pow is its multiplicity
        # (thus \prod_{cycloFactor} pol^pow equals the {pCyclo}-th cyclotomic polynomial)
        cycloFactor = list(polRing(cyclotomic_coeffs(pCyclo)).factor())
        
        # Small program verification step
        def verifyAndRaise(pol, power):
            """From theory it is known that we should have distinct irreducible factors,
            thus all powers should be 1. We assert this and return the factor."""
            if not power == 1:
                print(f"ERROR: Factorising {pCyclo}-th pol over {field} gives {cycloFactor}")
            return pow(pol, power)
        factors = [verifyAndRaise(pol, power) for (pol, power) in cycloFactor]
        factorsSupp = [(fac, self.support(fac)) for fac in factors]
        minSupp = min([supp for (pol, supp) in factorsSupp])

        return (qFactor, factorsSupp, minSupp)
    
    
    def sparsestField(self, p, facFields):
        """This function factors the {p}-th cyclotomic polynomial
        over different finite fields (fields with number of elements in {facFields})
        and returns the field and the factorisation that has the most sparse
        irreducible factor. All valid counterexamples found
        are appended to self.allCounterExmpls.
        Parameters:
        - p: determines the 'p'-th cyclotomic polynomial that will be factored
        - facFields: list of primes, p-th cyclotomic polynomial will be factored over
                     fields with these number of elements.
        Return value (three-tuple):
        - qMin: number of elements of the field for which the most sparse irreducible
                factor was found (element of facFields).
        - absMinSup: the support size of this sparsest irreducible factor
        - order: Multiplicative order of qMin in Z/pZ.
        """
        # Only allow primes and intialise
        assert(p in self.primes)
        (qMin, absMinSup, minOrder) = (0, 0, p + 1)
        
        # For every number of elements q in facFields, factor the
        # p-th cyclotomic polynomial over F_q.
        for q in facFields:
            (qFactor, factorsSupp, minSupp) = self.factorCyclo(p, q)
            
            # Calculate multiplicative order of q in Z/pZ
            order = self.multOrder(p, q)
            
            # Check whether case is interesting (that is, when minSupp <= ord_p(q)) and
            # add to list of counterexamples
            if minSupp <= order:
                self.allCounterExmpls.append((p, qFactor, minSupp, order))

                # Update the current 'minimum-information' if we have found
                # a smaller support size or if this is the first interesting case
                if minSupp < absMinSup or absMinSup == 0:
                    (qMin, absMinSup, minOrder) = (qFactor, minSupp, order)

        return (qMin, absMinSup, minOrder)
            
            
    def calcPrimeSparsePairs(self):
        """This function loops over each prime p in self.primeList and factorises
        the p-th cyclotomic polynomial over the fields F_q for q in self.fieldList.
        For each p, this adds (p, qMin, minSup, order) to self.bestCounterExmpls,
        where F_{qMin} is the field that gave the sparsest irreducible factor
        of the p-th cyclotomic polynomial, {minSup} is the support size of that
        factor and order is ord_{p}({qMin}).
        """
        if self.showTrace:
            print(f"Factoring the p-th cyclotomic polynomials over F_q")
            print(f"for p in {self.primeList}")
            print(f"and q in {self.fieldList}\n\n")

        for p in self.primeList:
            assert(p in self.primes)
            
            if self.showTrace:
                print(f"p = {p}", end=('.' if p == self.primeList[-1] else ', '))

            # Only interesting fields F_q are the ones where ord_p(q) \notin {1, 2, p - 1}
            # (for the other cases there is no counterexample, see thesis).
            goodFields = [q for q in self.fieldList
                          if (p != q and (not self.multOrder(p, q) in [1, 2, p - 1]))]
            (qMin, minSup, order) = self.sparsestField(p, goodFields)
            self.bestCounterExmpls.append((p, qMin, minSup, order))


    ## Note: The functions below are meant to be used after initialisation
    def runExperiment(self):
        """This function runs the cyclotomic polynomial factoring experiment."""
        if self.read:
            self.bestCounterExmpls = self.read_csv(self.fileName + 'SPARSEST')
            self.allCounterExmpls  = self.read_csv(self.fileName + 'ALL')
        else:
            # Run experiments
            self.bestCounterExmpls = list()
            self.allCounterExmpls  = list()
            self.calcPrimeSparsePairs()
            
            # Write results to file
            self.write_csv(self.fileName + 'SPARSEST', self.bestCounterExmpls)
            self.write_csv(self.fileName + 'ALL', self.allCounterExmpls)
        return self.bestCounterExmpls, self.allCounterExmpls
    
    
    def searchHighOrder(self, read=None, fileName=None, primeList=None):
        """For each prime p in primeList, this function loops over the first
        10.000 primes for q and factors the p-th cyclotomic polynomial over
        F_q, for the smallest q of order (p - 1)/2."""
        if read is None:
            read = self.read
        if fileName is None:
            fileName = self.fileName
        if primeList is None:
            primeList = self.primeList

        if read:
            return self.read_csv(f'csvFiles/{fileName}OrderTest')

        # Make a list of the first 10000 primes, for each p in primeList,
        # loop over these primes and search for a q of order (p - 1)/2.
        fieldList = [self.primes.unrank(i) for i in range(10000)]
        counterexamples = list()
        print(f'Primes: {primeList}')

        for p in primeList:
            print(f'Testing p = {p}...')
            assert(p in self.primes)
            for indx, q in enumerate(fieldList):
                # Print every 1000-th iteration of the inner-loop.
                if indx % 1000 == 0:
                    print(f'The {indx}-th iteration of q')
                assert(q in self.primes)
                order = self.multOrder(p, q)

                # If q has order (p - 1)/2 in (F_p)^*, see if this entails a counterexample.
                if p != q and order == (p - 1)/2:
                    (qFactor, factors, minSupp) = self.factorCyclo(p, q)
                    if minSupp <= order:
                        print(f'Found counterexample! Prime q = {q} has order = {order} in (F_p)^* for p = {p}.')
                        print('---------------------------------------------------------')
                        counterexamples.append((p, q, minSupp, order))
                    else:
                        print(f'PAY ATTENTION: Prime p = {p} and prime q = {q} gives order = {order}, but')
                        print(f'gives no counterexample when factoring. Sparsity = {minSupp}, factors w/ sparsity = {factors}')
                        print('---------------------------------------------------------')
                    break
        
        # Write results to file
        self.write_csv(f'csvFiles/{fileName}OrderTest', counterexamples)
        return counterexamples
                        
    
    def scatterPlot(self, x, y, colors, xlabel, ylabel, title, ylim=None, savefig=''):
        """This function makes a scatter plot of {x}, {y} and uses colors in {colors}.
        Title of the plot is {title}, axis labels are {xlabel} and {ylabel}.
        - savefig: file path to save the figure. The figure is only saved when savefig is not
                    an empty string."""
        # Assert x, y and colors have the same length.
        assert(len(x) == len(y))
        assert(len(y) == len(colors))
        
        # Plot the data
        fig = plt.figure()
        plt.scatter(x, y, color=colors)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        
        # Set ylim if it is not None
        if ylim:
            plt.ylim(ylim)
        plt.show()
        
        # Save figure is savefig is not empty
        if savefig:
            fig.savefig(savefig, bbox_inches='tight')
    
    
    def plotSparsest(self, savefig=''):
        """This function plots the cardinality of the field F_q that gave the sparsest
        factor against the corresponding primes p (p-th cyclotomic polynomial factored over F_q).
        Parameter:
        - savefig: file path to save the figure. The figure is only saved when savefig is not
                    an empty string."""
        primes = [p for (p, qMin, minSupp, order) in self.bestCounterExmpls]
        fields = [qMin for (p, qMin, minSupp, order) in self.bestCounterExmpls]
        
        # Make point red if invalid counterexample
        colors = ['red' if order >= p or minSupp > order else 'green' for (p, qMin, minSupp, order) in self.bestCounterExmpls]
        
        # Plot the data
        self.scatterPlot(primes, fields, colors, 'Prime', 'Field with lowest sparsity',
                         'Sparsest counterexamples', savefig=savefig)


    def plotAll(self, title='All counterexamples', spec=lambda p, q, sparsity, order:True, savefig=''):
        """This function plots the cardinality of all fields F_q that gave counterexmples
        against the corresponding primes p (p-th cyclotomic polynomial factored over F_q).
        Parameters:
        - title: title of plot.
        - spec: function of p, q, sparsity, and order, returning a boolean indicating whether
                 that counterexample should be plotted. Default: every counterexample is plotted
        - savefig: file path to save the figure. The figure is only saved when savefig is not
                    an empty string.
        """
        primes, fields, colors = ([], [], [])
        for (p, q, spars, order) in self.allCounterExmpls:
            if spec(p, q, spars, order):
                primes.append(p)
                fields.append(q)
                colors.append('red') if order >= p or spars > order else colors.append('green')
        self.scatterPlot(primes, fields, colors, 'Prime', 'Field with low sparsity', title, savefig=savefig)
    

    def plotStrength(self, ylim=None, strengthFunc=(lambda order, spars, p: order - spars), title=None,
                     ylabel='$ord_p(q)$ - sparsity', spec=lambda p, qMin, sparsity, order:True,
                     savefig=''):
        """This function plots the strength of all counterexamples against the primes p and returns
        a list of all counterexamples, sorted by their strength (descending).
        Parameters:
         - ylim: set the yrange
         - strengthFunc: function of the order, sparsity and prime p defining strength
         - title: title of plot
         - ylabel: should be same formula as strengthFunc in words (as a string)
         - spec: function of p, q, sparsity, and order, returning a boolean indicating whether
                 that counterexample should be plotted. Default: every counterexample is plotted.
         - savefig: file path to save the figure. The figure is only saved when savefig is not
                    an empty string.
        """
        prim, strengths, colors = [], [], []
        strengthList = list()
        
        # Calculate strength of all counterexamples.
        for (p, q, sparsity, order) in self.allCounterExmpls:
            if spec(p, q, sparsity, order):
                prim.append(p)
                strength = 0 if sparsity == 0 else strengthFunc(order, sparsity, p)
                strengthList.append((p, q, sparsity, order, strength))
                strengths.append(strength)

                # Red color if invalid, blue if strength = 0 and green otherwise.
                if strength < 0 or p==2 or p==5:
                    colors.append('red')
                elif strength == 0:
                    colors.append('blue')
                else:
                    colors.append('green')
    
        # Make plot
        if not title:
            title = 'Strength of all counterexamples plotted against the prime $p$'
        self.scatterPlot(prim, strengths, colors, 'Prime', f'Strength of counterexample\n{ylabel}', title, ylim, savefig=savefig)
        
        # Sort counterexamples on strength and return list.
        strengthList.sort(key=lambda x:x[-1], reverse=True)
        return [(p, q, sparsity, order, strength, p - 1 - order) for (p, q, sparsity, order, strength) in strengthList]


# -------------------------- STRENGTH ANALYSIS CLASS --------------------------
class StrengthAnalysis:
    """This class implements functions to analyse the strength of a given
    list of counterexamples (sorted on strength)."""
    def __init__(self, cycloFact, strengthList):
        """Initialise StregnthAnalysis object.
        Parameters:
            - cycloFact: CycloFactExperiment object
            - strengthList: a list of tuples (p, q, sparsity, order, strength, p - 1 - order)
                            sorted from highest strength to lowest strength
                            (e.g., output of cycloFact.plotStrength()).
        Attributes:
            - self.cycloFact, self.strengthList (same as parameters above)
            - self.counterFields: dictionary with primes p as keys and list
                                  of tuples (rank, (q, (p - 1)/ord_p(q), sparsity)) that give
                                  counterexamples for this p as value.
            - self.counterExamples: dictionary with primes p as keys and list
                                    of all counterexamples for that p as value.
        """
        assert(len(cycloFact.allCounterExmpls) == len(strengthList))
        self.cycloFact = cycloFact
        self.strengthList = strengthList
        
        # We create a dictionary with as key a prime number p and as
        # value a list of primes q that give counterexamples for that p.
        # We also create a similar dictionary with a list of each
        # counterexample for the key p instead of just the fields q.
        self.counterFields = dict()
        self.counterExamples = dict()
        
        # Initialise dictionaries
        for p in self.cycloFact.primeList:
            # We know there are no counterexamples for p = 2 and p = 5.
            if p == 2 or p == 5:
                continue
            self.counterFields[p] = list()
            self.counterExamples[p] = list()
        
        # Fill in the dictionaries
        for i, (p, q, sparsity, order, strength, diff) in enumerate(self.strengthList):
            self.counterFields[p].append((i, (q, (p - 1)/self.cycloFact.multOrder(p, q), sparsity)))
            self.counterExamples[p].append((i, (p, q, sparsity, order, strength, diff)))
    
    
    def __str__(self):
        """Printing this object results in printing the list of strengths
        with their rankings."""
        string = ""
        string += f"(rank, (p, q, sparsity, order, strength, p - 1 - order))\n"
        string += ',\n'.join([str(elem) for elem in list(enumerate(self.strengthList))])
        return string
    
    
    def __testStrengths(self, question, conditionFunc, firstEncounter=False):
        """This function tests whether a certain condition specified by conditionFunc
        holds for each first encounter of the primes p (thus, it tests whether the condition
        holds for the strongest counterexample of each prime p).
        Parameters:
            - question: a string that is printed at the start of the test.
            - conditionFunc: a function of (p, q, spars, order, strength, diff)
                that returns a boolean. The test 'succeeds' if for each
                counterexample the function returns True.
            - firstEncounter: boolean indicating whether conditionFunc should return
                True for all counterexamples or just for the strongest counterexample
                per prime p.
        Return:
            - list of cases for which the test fails.
        """
        print(f"TEST: {question}")
        print("---------------------------------------------------------------------------")
        # Create a set of the primes that have already been encountered
        # (only used if firstEncounter = True).
        seenPrimes = set()
        conjecture = True
        fails = list()
        
        # Loop over all counterexamples (sorted on strength) and check
        # whether the condition holds.
        for i, (p, q, spars, order, strength, diff) in enumerate(self.strengthList):
            if firstEncounter and p in seenPrimes:
                continue
            seenPrimes.add(p)
            if not conditionFunc(p, q, spars, order, strength, diff):
                conjecture = False
                fails.append((i, (p, q, spars, order, strength, diff)))
                print(f"FAIL: (p, q, spars, ord, strength, p - 1 - ord) = {(p, q, spars, order, strength, diff)}")
                print(f"\nFor the prime {p} (the {i}-th strongest example) the test fails.")
                print(f"{p} has counterexamples for (rank, (q, (p - 1)/ord_p(q), sparsity)) in")
                print(self.counterFields[p])
                print(f"\nAll counterexamples for {p} (sorted on strength) for (p, q, spars, ord, strength, p - 1 - ord) in:")
                print(self.counterExamples[p])
                print('----------------------------------------------------------\n\n')
        res = 'YES!' if conjecture else 'NO...'
        print(f"\nRESULT: {question} {res}")
        return fails
    
    
    def maximalOrderTest(self, firstEncounter=True):
        """This function tests whether counterexamples satisfy ord_p(q) = (p - 1)/2
        (maximal order). The parameter {firstEncounter} specifies whether to check this
        condition for the first encounter of each p or for all counterexamples."""
        word = 'Strongest' if firstEncounter else 'All'
        question = f"{word} counterexamples satisfy ord_p(q) = (p - 1)/2?"
        conditionFunc = lambda p, q, spars, order, strength, diff: order == diff
        return self.__testStrengths(question, conditionFunc, firstEncounter)
    
    
    def smallFieldTest(self):
        """This function tests whether each first encounter of a prime p
        (strongest counterexample for p) corresponds to the smallest
        prime q of all primes q that give a counterexample for p."""
        question = "Strongest counterexample for each p achieved for smallest possible field for p?"
        
        # Sort counterFields on the size of the field F_q.
        sortDict = dict()
        for p in self.cycloFact.primeList:
            # Skip primes without counterexample.
            if p == 2 or p == 5:
                continue
            # self.counterFields[p] is a list of 2-tuples [(rank, (q, (p-1)/ord, spars)), ...]
            sortDict[p] = sorted(self.counterFields[p], key=lambda x:x[1][0])
        conditionFunc = lambda p, q, spars, order, strength, diff, sortDict=sortDict: q == sortDict[p][0][1][0]
        return self.__testStrengths(question, conditionFunc, True)
        
    
    def minimalSparsityTest(self):
        """This function tests whether counterexamples satisfy the relation that
        a lower sparsity corresponds to a stronger counterexample."""
        question = "Strongest counterexample for each p achieved for field with sparsest factor?"

        # Sort counterFields on sparsity
        sortDict = dict()
        for p in self.cycloFact.primeList:
            # Skip primes without counterexample.
            if p == 2 or p == 5:
                continue
            # self.counterFields[p] is a list of 2-tuples [(rank, (q, (p-1)/ord, spars)), ...]
            sortDict[p] = sorted(self.counterFields[p], key=lambda x:x[1][2])
        conditionFunc = lambda p, q, spars, order, strength, diff, sortDict=sortDict: spars == sortDict[p][0][1][2]
        return self.__testStrengths(question, conditionFunc, True)
    
    
    def failToPrimes(self, fails):
        """This function translates the list of failures [(rank, (p, q, spars,..)), ..]
        to a set of prime numbers p for which the test failed."""
        return set([elem[1][0] for elem in fails])
    
    
    def testOccurence(self, question, conditionFunc):
        """This function tests whether a certain condition specified by conditionFunc
        holds for any of the counterexamples for each p.
        Parameters:
            - question: a string that is printed at the start of the test.
            - conditionFunc: a function of (p, q, spars, order, strength, diff)
                that returns a boolean. The test 'succeeds' if for each prime p there is a
                counterexample for which the function returns True.
        """
        print(f"TEST: {question}")
        print("---------------------------------------------------------------------------")
        fails = list()
        
        for p in self.cycloFact.primeList:
            # Skip primes without counterexample.
            if p == 2 or p == 5:
                continue
            # Loop over all counterexamples of p and check whether condition holds.
            for (rank, (p, q, spars, order, strength, diff)) in self.counterExamples[p]:
                if conditionFunc(p, q, spars, order, strength, diff):
                    break
            else:
                print(f"FAIL: prime {p} has no counterexamples that satisfy the condition")
                print(f"Counterexamples for {p} are (rank, (p, q, spars, order, strength, diff)):")
                print(self.counterExamples[p])
                print('----------------------------------------------------------')
                fails.append(p)
        print(f"Found {len(fails)} primes p that fail this test.")
        return fails
    
    
    def compareToKnown(self):
        """This function compares the strength of the counterexamples to the
        counterexamples found by Briet and Labib (ord_p(2) < p - 1).
        The strongest ones for each prime are classified into three groups.
          1. Primes p with ord_p(2) < p - 1 and q = 2 gives strongest counterexample
          2. Primes p with ord_p(2) < p - 1 and q = 2 is not strongest counterexamle
          3. Primes p with ord_p(2) = p - 1.
        The function returns three lists of strongest counterexamples of the prime p,
        one list for each category.
        """
        blStrongest = list()
        blNotStrongest = list()
        noBL = list()

        # Classify each prime to one of the three categories.
        for p in self.cycloFact.primeList:
            # Skip the prime numbers without counterexample.
            if p == 2 or p == 5:
                continue
            if self.cycloFact.multOrder(p, 2) < p - 1:
                # self.counterFields[p] is a list of 2-tuples [(rank, (q, (p-1)/ord, spars)), ...]
                if self.counterFields[p][0][1][0] == 2:
                    blStrongest.append(self.counterExamples[p][0])
                else:
                    blNotStrongest.append(self.counterExamples[p][0])
            else:
                noBL.append(self.counterExamples[p][0])
        summ = len(blStrongest) + len(blNotStrongest) + len(noBL)
        print(f"The {summ} primes analysed can be classified into three groups:")
        print(f"1. For {len(blStrongest)} primes p we have ord_p(2) < p - 1 and q = 2 gives the strongest counterexample")
        print(f"2. For {len(blNotStrongest)} primes p we have ord_p(2) < p - 1 and q = 2 is not the strongest")
        print(f"3. For {len(noBL)} primes p we have ord_p(2) = p - 1.")
        print("----------------------------------------------------------------")
        if len(blStrongest) > 0:
            blStrongest.sort(key=lambda x:x[0])
            print(f"(rank, (p, q, sparsity, order, strength, p - 1 - order)) = {blStrongest[0]} is the strongest one of group 1.")
        if len(blNotStrongest) > 0:
            blNotStrongest.sort(key=lambda x:x[0])
            print(f"(rank, (p, q, sparsity, order, strength, p - 1 - order)) = {blNotStrongest[0]} is the strongest one of group 2.")
        if len(noBL) > 0:
            noBL.sort(key=lambda x:x[0])
            print(f"(rank, (p, q, sparsity, order, strength, p - 1 - order)) = {noBL[0]} is the strongest one of group 3.")
        
        return blStrongest, blNotStrongest, noBL
