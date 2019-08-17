from typing import Optional, List
from sympy import symbols


class Enzymestate():
    def __init__(self, name: str):
        self.enzyme = symbols(name)

    def __str__(self):
        return self.enzyme.__str__()


class ReactionRate():

    def __init__(self, name: str, pseudoFirstOrderReactant: Optional[str] = None):
        self.name = symbols(name)
        if pseudoFirstOrderReactant:
            self.reactant = symbols(pseudoFirstOrderReactant)
        else:
            self.reactant = None

    def __str__(self):
        if self.reactant:
            return self.name.__str__() + self.reactant.__str__()
        else:
            return self.name.__str__()


class UnitReaction():
    def __init__(self, from_State: Enzymestate, rate: ReactionRate, to_State: Enzymestate):
        self.from_state = from_State
        self.to_state = to_State
        self.rate = rate


class Reactions():

    def __init__(self):
        self._reactions: List[UnitReaction] = []

    def addReaction(self, reaction: UnitReaction):
        self._reactions.append(reaction)

    def as_text(self) -> str:
        res = []
        for reac in self._reactions:
            start, rate, end = reac.from_state, reac.rate, reac.to_state
            res.append(f'{start} \t --{rate}--> \t {end} ')

        return '\n'.join(res)


k_0 = ReactionRate("k1")
k_1 = ReactionRate("k-1", "A")

E1, E0 = [Enzymestate(name) for name in ["EA", "E"]]

r1 = UnitReaction(E0, k_0, E1)
r2 = UnitReaction(E1, k_1, E0)

scheme = Reactions()
scheme.addReaction(r1)
scheme.addReaction(r2)

print(scheme.as_text())

print(k_0.reactant is None)
