#%%
from typing import Optional, List, Tuple, Any
from sympy import symbols
from pandas import DataFrame
from random import randint



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

    def directedPatterns(self, state: Enzymestate):
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

        
# read Reaction mechanism
mechanism = Reactions()

with open("fumarase.txt", "r") as infile:
    for count, line in enumerate(infile):
        # skip comments and empty lines
        if line.startswith("#") or len(line.strip()) == 0:
            break
        els: List[str] = line.split(',')
        if len(els) != 3:
            raise TypeError(
                f"Infile line {count+1}:expected 3 comma separated values per line, instead got:\n{line}")

        e1, rate, e2 = els[0].strip(), els[1].split(';'), els[2]

        if len(rate) < 1 or len(rate) > 2:
            raise TypeError(
                "Rate costant must be given as 1 parameter or 2 parameters seperated by \";\"")

        es1 = Enzymestate(e1)
        es2 = Enzymestate(e2)
        rrate = ReactionRate(*[x.strip() for x in rate])

        reaction = UnitReaction(es1, rrate, es2)
        mechanism.addReaction(reaction)

for state in mechanism._enzymeStates:
    print(f'Enzymestate: {state}')
    print(DataFrame(mechanism.directedPatterns(state)))

# l1, l2, l3, l4 = [1, 2, 6], [1, 4], [3, 4, 5], [5, 6]

# reacs = [l1, l2, l3, l4]
# res = Wang_algebra.wang_product(reacs)
# res = [sorted(el) for el in res]
# print(f'{len(res)} products: {res}')


