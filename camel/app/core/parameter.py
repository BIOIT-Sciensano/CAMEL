from typing import Optional

from pydantic import BaseModel, ConfigDict


class Parameter(BaseModel):
    """
    Represents a tool parameter.
    """
    model_config = ConfigDict(extra='forbid')
    name: str
    option: Optional[str] = None
    flag: bool = False
    default: bool = False
    mandatory: bool = False
    value: Optional[str | int | float] = None
    p_index: Optional[int] = 0

    def __repr__(self) -> str:
        """
        Returns the internal representation.
        :return: Internal representation.
        """
        return f"Parameter(name='{self.name}', val='{self.value}')"

    def __str__(self) -> str:
        """
        Returns the parameter as a string (for the command line).
        """
        if self.flag:
            return self.option
        if self.option is not None and self.value is not None:
            return f'{self.option} {self.value}'
        if self.option is None:
            return str(self.value)
        return self.option
