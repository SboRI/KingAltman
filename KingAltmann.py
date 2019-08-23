from typing import Optional, List, Tuple, Any
from sympy import symbols
import sympy
from pandas import DataFrame
from random import randint
from functools import reduce



class Wang_algebra():

    @staticmethod
    def wang_product(ele: List[List[Any]]):
        """Alphanumeric multiplication according to wang rules: x * x = 0, y*x + y*x = 0"""

        def wang_multiply(a: List[Any], b: List[Any]):
            res = []
            for el_b in b:
                _res = a + [el_b]
                res.append(_res)
            return res

        elements = list(ele)

        # nothing to multiply
        if len(elements) < 2:
            return elements

        # First round of the algorithm has a flat list [a,b,c] as element 0, helper algorithm needs [a] as input, so put the elements in a list
        res: List[Any] = [[e] for e in elements.pop(0)]
        while len(elements) >= 1:
            term1 = elements.pop(0)
            #print(f'current: {res}, multiplyinng with {term1}')

            _res = []
            for el in res:
                wang_prods = wang_multiply(el, term1)
                for result in wang_prods:
                    #
                    _res.append(result)

                for n_el, el in enumerate(_res):
                    # wang rule 1: x*x = 0
                    if len(set(el)) < len(el):
                        #print(f'delete x*x = 0, el: {el}')
                        _res[n_el] = []

                    # wang rule 2: xy + yx = 0, iterate over all elements and compare if 2 sets of elements are equal.
                    # Skip the comparison of the element to itself
                    for n_el2, el2 in enumerate(_res):
                        if n_el2 == n_el:
                            break
                        if set(el) == set(el2):
                            #print(f'delete xy + xy = 0, els = {el}, {el2}')
                            _res[n_el] = []
                            _res[n_el2] = []

                _res = list(filter(lambda el: len(el) > 0, _res))
            res = _res

            #print(f'result: {res}')

        return res


class Enzymestate():
    def __init__(self, name: str):
        self.enzyme = symbols(name)

    def __str__(self):
        return self.enzyme.__str__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.enzyme == other.enzyme
        return False
    def to_value(self):
        return self.enzyme


class ReactionRate():

    def __init__(self, name: str, pseudoFirstOrderReactant: Optional[str] = None):
        self.name = symbols(name)
        if pseudoFirstOrderReactant:
            self.reactant = symbols(pseudoFirstOrderReactant)
        else:
            self.reactant = None

    def __str__(self):
        if self.reactant:
            return self.name.__str__() + "[" + self.reactant.__str__() + "]"
        else:
            return self.name.__str__()

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.name == other.name
        return False
    
    def to_value(self) -> symbols:
        if self.reactant:
            return self.name * self.reactant
        return self.name


class UnitReaction():
    def __init__(self, from_State: Enzymestate, rate: ReactionRate, to_State: Enzymestate):
        self.from_state = from_State
        self.to_state = to_State
        self.rate = rate

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.from_state == other.from_state and self.to_state == other.to_state and self.rate == other.rate
        return False

    def __str__(self):
        return f"{self.from_state} \t --{self.rate}--> \t {self.to_state}"


class BiDirReaction():
    """Constructs a Non/Bi directional rate from 2 unit reactions"""

    def __init__(self, fwd_reaction: UnitReaction, rev_reaction: UnitReaction):
        self._fwd_reaction = fwd_reaction
        self._rev_reaction = rev_reaction

    def contains_Rate(self, rate: Optional[ReactionRate]):
        if not rate:
            return False
        return (rate == self._fwd_reaction.rate or rate == self._rev_reaction.rate)

    def contains_reaction(self, reaction: Optional[UnitReaction]) -> bool:
        if not reaction:
            return False
        return (reaction == self._fwd_reaction or reaction == self._rev_reaction)

    def produces(self, enzymestate: Enzymestate) -> Optional[UnitReaction]:
        if self._fwd_reaction.to_state == enzymestate:
            return self._fwd_reaction

        if self._rev_reaction.to_state == enzymestate:
            return self._rev_reaction

        return None

    def __str__(self):
        return f"{self._fwd_reaction.rate}//{self._rev_reaction.rate}"


class Reactions():

    def __init__(self):
        self._reactions: List[UnitReaction] = []

        self._enzymeStates: List[Enzymestate] = []
        """list of all UNIQUE enzyme states"""

        self._bidirectionalRates: List[BiDirReaction] = []
        """List of all bidirectional couples, (ki[a], k-1) from E -> k1[a] -> EA, EA -> k-1 -> E)"""

        self._product_forming_complex = []
        """List of all [E]*rate which are producing the product your are looking for, e.g. dp/dt = + [ES/EP]*k6"""

        self._product_consuming_complex = []
        """List of all [E]*rate which are consuming the product your are looking for, e.g. dp/dt = - [E]*k-6 *[P]"""

        self._null_rates = []
        """List of rates and concentrations that are considered to be 0, e.g. for irreversible reactions k-1 = 0, or for v_o conditions [P]=0""" 

    def addReaction(self, reaction: UnitReaction):
        if len(self.reaction_from_Rate(reaction.rate)) > 0:
            raise AttributeError("rate constant already defined")
        self._reactions.append(reaction)
        self._add_enzymeStates(reaction)
        self._add_bidirectionalRates(reaction)

    def _add_enzymeStates(self, u: UnitReaction):
        """adds uniquely an Enzymestate to internal list of all _enzymeStates"""
        for e in [u.from_state, u.to_state]:
            if not e in self._enzymeStates:
                self._enzymeStates.append(e)

    def _add_bidirectionalRates(self, reaction: UnitReaction):
        """ adds a bidirectional Rate (Tuple) if it doesn't exist yet and both directions exist in the network"""

        # if the rate has already been used, it is added a second time, which shouldn't happen
        if any([bi.contains_reaction(reaction) for bi in self._bidirectionalRates]):
            raise AttributeError(
                "Added a unitreaction with the same rate constant a second time")

        # find the reverse reaction
        reverseReaction = self.reverse_Reaction(reaction)
        if not reverseReaction:
            return
        else:
            self._bidirectionalRates.append(
                BiDirReaction(reaction, reverseReaction)
            )
    def add_product_forming_complex(self, e: Enzymestate, r: ReactionRate):
        self._product_forming_complex.append([e, r])

    def add_product_consuming_complex(self, e: Enzymestate, r: ReactionRate):
        self._product_consuming_complex.append([e, r])

    def as_text(self) -> str:
        res = []
        for reac in self._reactions:
            res.append(f'{reac}')

        return '\n'.join(res)

    def reverse_Reaction(self, reaction: UnitReaction) -> Optional[UnitReaction]:
        fwdReaction = reaction
        fwdSub, fwdProd = fwdReaction.from_state, fwdReaction.to_state

        reactions = self.consumed_by(fwdProd)
        reverseReact = list(
            filter(lambda react: react.to_state == fwdSub, reactions))

        if len(reverseReact) == 0:
            return None
        if len(reverseReact) == 1:
            return reverseReact[0]
        else:
            raise AttributeError("more than 1 reverse Reaction defined")

    def reaction_from_Rate(self, rate: ReactionRate) -> List[UnitReaction]:
        """Returns the unitreaction which contains the ReactionRate rate"""
        return list(filter(lambda reac: reac.rate == rate, self._reactions))

    def reaction_from_Reactants(self, from_state: Enzymestate, to_state: Enzymestate) -> Optional[UnitReaction]:
        res: List[UnitReaction] = list(filter(
            lambda reaction: reaction.from_state == from_state and reaction.to_state == to_state, self._reactions))

        if len(res) == 0:
            return None

        if len(res) > 1:
            raise AttributeError(
                "more than 1 reaction with identical reactants defined")

        return res[0]

    def produced_by(self, e: Enzymestate) -> List[UnitReaction]:
        """Returns all UnitReactions that produce input enzyme state {e}"""
        return list(filter(lambda react: react.to_state == e, self._reactions))

    def consumed_by(self, e: Enzymestate) -> List[UnitReaction]:
        """Returns all UnitReactions that consume input enzyme state e"""
        return list(filter(lambda react: react.from_state == e, self._reactions))

    def linear_graph_matrix(self) -> List[List[Optional[BiDirReaction]]]:
        """ Returns the linear graph structure (Qi...Beard 2008, BMC Bioinfo Eq. 4) of the reaction network."""
        kin_matrix = self.kinetic_matrix()

        # res[i][j] BIDIRECTIONAL ReactionRate for ENzymstate i -> Enzymstate j
        res: List[List[Optional[BiDirReaction]]] = [
            [None for _ in self._enzymeStates] for _ in self._enzymeStates]
        for n_es1, _ in enumerate(self._enzymeStates):
            for n_es2, _ in enumerate(self._enzymeStates):
                _res = None if kin_matrix[n_es1][n_es2] == None else list(
                    filter(lambda biDirect: biDirect.contains_Rate(kin_matrix[n_es1][n_es2]), self._bidirectionalRates))[0]
                res[n_es1][n_es2] = _res
        return res

    def kinetic_matrix(self) -> List[List[Optional[ReactionRate]]]:
        """ Returns the kinetic matrix (Qi...Beard 2008, BMC Bioinfo Eq. 5) of the reaction network."""

        # res[i][j] ReactionRate for ENzymstate i -> Enzymstate j
        res: List[List[Optional[ReactionRate]]] = [
            [None for _ in self._enzymeStates] for _ in self._enzymeStates]

        for n_es1, es1 in enumerate(self._enzymeStates):
            for n_es2, es2 in enumerate(self._enzymeStates):
                reac = self.reaction_from_Reactants(es1, es2)
                _res = reac.rate if reac else None
                res[n_es1][n_es2] = _res

        return res

    def kaPatterns(self):
        # make a dictionary to replace bidirectional rates with numbers
        # create a matrix equal to the linear_graph_matrix but with numbers
        ratesDict = dict()
        lin_graph = self.linear_graph_matrix()
        lin_graph_numbers = [[None for _ in self._enzymeStates]
                             for _ in self._enzymeStates]

        for x, rate in enumerate(self._bidirectionalRates):
            rateNumber = x+1
            ratesDict[rateNumber] = rate
            for n, el in enumerate(lin_graph):
                for m, birate in enumerate(el):
                    if birate == rate:
                        lin_graph_numbers[n][m] = rateNumber

        deleteEnzymestate = randint(0, len(lin_graph_numbers) - 1)
        lin_graph_numbers.pop(deleteEnzymestate)

        wang_patterns = [list(filter(lambda x:x != None, el))
                         for el in lin_graph_numbers]
        wang_products = Wang_algebra.wang_product(wang_patterns)

        res = []
        for rates in wang_products:
            _res = []
            for n in rates:
                _res.append(ratesDict[n])
            res.append(_res)
        return res

    def directedPatterns(self, state: Enzymestate) -> List[List[ReactionRate]]:
        # find rate that produces said enzyme form
        kaPatterns = self.kaPatterns()
        def directedPattern(pattern: List[BiDirReaction]):
            pattern = list(pattern)
            
            res: List[UnitReaction] = [] #store the unitreactions that produce target Enzymestate
            target = state
            while len(pattern) > 0:
                reactions_to_delete = []
                for n_, biDir in enumerate(pattern):
                    target_reaction = biDir.produces(target)
                    if target_reaction:
                        res.append(target_reaction)
                        
                for reaction in res:
                    pattern = list(filter(lambda biDirRate: not biDirRate.contains_Rate(reaction.rate), pattern))
                for reaction in res:
                    target = reaction.from_state
                    for n_, biDir in enumerate(pattern):
                        target_reaction = biDir.produces(target)
                        if target_reaction:
                            res.append(target_reaction)
                    for reaction in res:
                        pattern = list(filter(lambda biDirRate: not biDirRate.contains_Rate(reaction.rate), pattern))

            return [el.rate for el in res]

        res = []
        for pattern in kaPatterns:
            res.append(directedPattern(pattern))
            #print(DataFrame(res))

        return res

    def _pattern_to_equation(self, directedPattern: List[List[ReactionRate]]):
        def my_multiply(singlePattern: List[ReactionRate]):
            return reduce(lambda x,y: x*y.to_value(),singlePattern, 1)
        

        return reduce(lambda x, y: x + my_multiply(y), directedPattern, 0)
    
    def numerator(self, enzymstate: Enzymestate):
        pat = self.directedPatterns(enzymstate)
        return self._pattern_to_equation(pat)

    def denominator(self):
        nominators = [self._pattern_to_equation(self.directedPatterns(state)) for state in self._enzymeStates]
        return reduce(lambda x,y: x + y, nominators)
    
    def product_term_summation(self, product_terms: List[List[Any]], with_full_numerator=True):
        def term_helper(enz_rate: List[Any]):
            state=enz_rate[0]
            rate=enz_rate[1]
            numerator = self.numerator(state)
            return rate.to_value() * numerator 
        
        res = None

        if with_full_numerator:
            _res = map(lambda e_r: term_helper(e_r), product_terms)
            res = reduce(lambda x, y: x + y, _res)
        
        else:
            _res = map(lambda e_r: e_r[0].to_value() * e_r[1].to_value(), product_terms)
            res = reduce(lambda x,y: x+y, _res)

        return res.simplify()*symbols("E_{0}")

    def solve_for_product(self, with_full_numerator=True):
        product_forming = self.product_term_summation(self._product_forming_complex, with_full_numerator)
        product_consuming = self.product_term_summation(self._product_consuming_complex, with_full_numerator)
        denominator = self.denominator()
        return ((product_forming - product_consuming)/denominator).simplify()

    def simplify_null_pathways(self):
        eq = self.solve_for_product()
        for el in self._null_rates:
            eq = eq.subs(el, 0)
        return eq




        
# read Reaction mechanism
mechanism = Reactions()

with open("upo.txt", "r") as infile:
    def check_input_items(els):
        if len(els) != 3:
            raise TypeError(
                f"Infile line {count+1}:expected 3 comma separated values per line, instead got:\n{line}")
    
    def check_rate(rate):
        if len(rate) < 1 or len(rate) > 2:
            raise TypeError(
                "Rate costant must be given as 1 parameter or 2 parameters seperated by \";\"")



    
    def sanitize_inputs(line):
        els: List[str] = line.split(',')
        check_input_items(els)
        e1, rate, e2 = els[0].strip(), els[1].split(';'), els[2]
        rate = [x.strip() for x in rate]
        check_rate(rate)
        return e1, rate, e2

 

    for count, line in enumerate(infile):
        is_product_pos = False
        is_product_neg = False
        is_null_pathway = False

        def is_reaction_line():
            res= (is_product_pos or is_product_neg or is_null_pathway)
            return not res
        # skip comments and empty lines
        if line.startswith("#") or len(line.strip()) == 0:
            continue
        if line.startswith('=+:'):
            is_product_pos = True
            line=line.replace('=+:', '')
        if line.startswith('=-:'):
            is_product_neg = True
            line=line.replace('=-:', '')
        if line.startswith('=0:'):
            is_null_pathway = True
            line=line.replace('=0:', '')
    
        if is_product_pos or is_product_neg:
            e1, rate, _ = sanitize_inputs(line)
            es1 = Enzymestate(e1)
            rrate = ReactionRate(*rate)
            if is_product_pos:
                mechanism.add_product_forming_complex(es1, rrate)
            if is_product_neg:
                mechanism.add_product_consuming_complex(es1, rrate)
        
        if is_reaction_line():
            e1, rate, e2 = sanitize_inputs(line)
            es1 = Enzymestate(e1)
            es2 = Enzymestate(e2)
            rrate = ReactionRate(*rate)
            unitReaction = UnitReaction(es1, rrate, es2)
            mechanism.addReaction(unitReaction)

        
        if is_null_pathway:
            for nulls in line.split(','):
                if not len(nulls.strip())==0:
                    mechanism._null_rates.append(symbols(nulls.strip()))
            

        


# reactionR = ReactionRate("k1", "A")
# print(reactionR.to_value().subs(symbols("A"), 1))
# for state in mechanism._enzymeStates:
#     print(f'Enzymestate: {state}')
#     print(mechanism.numerator(state))

# print(f'Denominator:\n{mechanism.denominator()}')

# print("simple")
# sympy.pprint(mechanism.solve_for_product(False))

# print("complex")
# sympy.pprint(mechanism.solve_for_product())

# print(DataFrame(mechanism._null_rates))

print("with zero")
with_zero = mechanism.simplify_null_pathways()
sympy.pprint(with_zero)


k_ib = symbols("K_iB")
k_ma = symbols("K_mA")
k_ma2 = symbols("K_mA2")
k_mb = symbols("K_mB")
k_1, k_2, k_3, k_4, k_5, k_6, k_7, k_8 =symbols(["k_"+str((x+1)*-1) for x in range(8)])
k1, k2, k3, k4, k5, k6, k7, k8 = symbols("k1:9")


eq = with_zero.subs(k8, k_8/k_ib)

eq = eq.subs(k1, (k_1 + k2)/k_ma)
eq = eq.subs(k4, (k_4+k7)/k_ma2)
eq = eq.subs(k3, (k_3+k6)/k_mb)

#vp = symbols("v_p")
#eq = eq.subs(k6, (vp*(k2+k6)/k2))
#vp2 = symbols("v_p2")
#eq = eq.subs(k7, (vp2*(k2+k7)/k2))
eq = eq.factor()

sympy.pprint(eq)


with open("out2.tex", "w") as outfile:
     outfile.write(sympy.latex(eq))




