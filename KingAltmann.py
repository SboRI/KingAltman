from typing import Optional, List, Tuple
from sympy import symbols


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

    def __str__(self):
        return f"{self.from_state} \t --{self.rate}--> \t {self.to_state}"


class Reactions():

    def __init__(self):
        self._reactions: List[UnitReaction] = []

        self._enzymeStates: List[Enzymestate] = []
        """list of all UNIQUE enzyme states"""

        self._bidirectionalRates: List[Tuple[ReactionRate, ReactionRate]] = []
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
        """ adds a bidirectional rate Tuple if it doesn't exist yet and both directions exist in the network"""

        # if the rate has already been used, it is added a second time, which shouldn't happen
        if any([reaction.rate in t for t in self._bidirectionalRates]):
            raise AttributeError(
                "Added a unitreaction with the same rate constant a second time")

        # find the reverse reaction
        reverseReaction = self.reverse_Reaction(reaction)
        if not reverseReaction:
            return
        else:
            self._bidirectionalRates.append(
                (reaction.rate, reverseReaction.rate)
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

    def produced_by(self, e: Enzymestate) -> List[UnitReaction]:
        """Returns all UnitReactions that produce input enzyme state {e}"""
        return list(filter(lambda react: react.to_state == e, self._reactions))

    def consumed_by(self, e: Enzymestate) -> List[UnitReaction]:
        """Returns all UnitReactions that consume input enzyme state e"""
        return list(filter(lambda react: react.from_state == e, self._reactions))

    def linear_graph_matrix(self):
        """ Returns the linear graph structure (Qi...Beard 2008, BMC Bioinfo Eq. 4) of the reaction network."""
        pass

    def kinetic_matrix(self):
        """ Returns the kinetic matrix (Qi...Beard 2008, BMC Bioinfo Eq. 5) of the reaction network."""
        pass


class uniReaction():

    def __init__(self, reaction: Reactions):
        self._bidirectional = reaction


# read Reaction mechanism
mechanism = Reactions()

with open("2substr.txt", "r") as infile:
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

[print(f'{a}, {b}') for a, b in mechanism._bidirectionalRates]

# k_0 = ReactionRate("k1")
# k_1 = ReactionRate("k-1", "A")

# E1, E0 = [Enzymestate(name) for name in ["EA", "E"]]

# r1 = UnitReaction(E0, k_0, E1)
# r2 = UnitReaction(E1, k_1, E0)

# scheme = Reactions()
# scheme.addReaction(r1)
# scheme.addReaction(r2)

# print(scheme.as_text())

# print(k_0.reactant is None)
