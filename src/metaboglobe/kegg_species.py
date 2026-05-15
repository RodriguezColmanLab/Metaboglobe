class KeggReactionId:
    """Holds a KEGG reaction ID, for example "rn:R02235". Immutable object."""

    _reaction_id: str  # Of the form "rn:R02235"

    @staticmethod
    def create_from_id(reaction_id: str) -> "KeggReactionId":
        """Creates a KeggReactionId object from a string, for example "R02235". Raises ValueError if reaction_id does
        not start with "R".

        If your string is of the form "rn:R02235", you can just use the constructor KeggReactionId("rn:R02235") instead
        of this method.
        """
        if not reaction_id.startswith("R"):
            raise ValueError("Reaction ID must start with R")
        return KeggReactionId("rn:" + reaction_id)

    def __init__(self, reaction_id: str):
        """Creates a new KeggReactionId object. Raises ValueError if reaction_id does not start with "rn:"."""
        if not reaction_id.startswith("rn:"):
            raise ValueError("Reaction ID must start with rn:")
        self._reaction_id = reaction_id

    @property
    def reaction_id(self) -> str:
        """Just returns the reaction ID, for example "rn:R02235"."""
        return self._reaction_id

    def __str__(self) -> str:
        """Just returns the reaction ID, for example "rn:R02235"."""
        return self._reaction_id

    def __repr__(self) -> str:
        return f"<{self._reaction_id}>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, KeggReactionId) and self._reaction_id == other.reaction_id

    def __hash__(self) -> int:
        return hash(self._reaction_id)


class KeggCompoundId:
    """Holds a KEGG compound ID, for example "cpd:C00504". Immutable object."""

    _compound_id: str  # Of the form "cpd:C00504"

    @staticmethod
    def create_from_id(compound_id: str) -> "KeggCompoundId":
        """Creates a KeggCompoundId object from a string, for example "C00504". Raises ValueError if reaction_id does
        not start with "C".

        If your string is of the form "cpd:C00504", you can just use the constructor KeggCompoundId("cpd:C00504")
        instead of this method.
        """
        if not compound_id.startswith("C"):
            raise ValueError("Compound ID must start with C")
        return KeggCompoundId("cpd:" + compound_id)

    def __init__(self, compound_id: str):
        """Creates a new KeggCompoundId object. Raises ValueError if compound_id does not start with "rn:"."""
        if not compound_id.startswith("cpd:"):
            raise ValueError("Compound ID must start with cpd:")
        self._compound_id = compound_id

    @property
    def compound_id(self) -> str:
        """Just returns the compound ID, for example "cpd:C00504"."""
        return self._compound_id

    def __str__(self) -> str:
        """Just returns the compound ID, for example "cpd:C00504"."""
        return self._compound_id

    def __repr__(self) -> str:
        return f"<{self._compound_id}>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, KeggCompoundId) and self._compound_id == other.compound_id

    def __hash__(self) -> int:
        return hash(self._compound_id)


class KeggPathwayId:
    """Holds a KEGG pathway ID, for example "path:mmu00250". Immutable object."""

    _pathway_id: str  # Of the form "path:mmu00250"

    @staticmethod
    def create_from_id(pathway_id: str) -> "KeggPathwayId":
        """Creates a KeggPathwayId object from a string, for example "mmu00250". Raises ValueError if reaction_id does
        not start with "C".

        If your string is of the form "path:mmu00250", you can just use the constructor KeggPathwayId("path:mmu00250")
        instead of this method.
        """
        if pathway_id.startswith("path:"):
            raise ValueError("Pathway ID must not start with path:")
        return KeggPathwayId("path:" + pathway_id)

    def __init__(self, pathway_id: str):
        """Creates a new KeggPathwayId object. Raises ValueError if pathway_id does not start with "rn:"."""
        if not pathway_id.startswith("path:"):
            raise ValueError("Pathway ID must start with path:")
        self._pathway_id = pathway_id

    @property
    def pathway_id(self) -> str:
        """Just returns the pathway ID, for example "path:mmu00250"."""
        return self._pathway_id

    def __str__(self) -> str:
        """Just returns the pathway ID, for example "path:mmu00250"."""
        return self._pathway_id

    def __repr__(self) -> str:
        return f"<{self._pathway_id}>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, KeggPathwayId) and self._pathway_id == other.pathway_id

    def __hash__(self) -> int:
        return hash(self._pathway_id)
