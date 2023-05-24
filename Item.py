from typing import Union


class Item:
    def __init__(self, width: int, height: int, mass: Union[int, float, None] = None):
        self.width = width
        self.height = height
        self.mass = mass

    def __repr__(self) -> str:
        return f'Item({repr(self.width)}, {repr(self.height)}, {repr(self.mass)})'

    def __str__(self) -> str:
        return f'Item({str(self.width)}, {str(self.height)}, {str(self.mass)})'


class Position:
    def __init__(self, x: Union[int, float], y: Union[int, float], rotated: bool):
        self.x = x
        self.y = y
        self.rotated = rotated

    def __repr__(self) -> str:
        return f'Position({repr(self.x)}, {repr(self.y)}, {repr(self.rotated)})'

    def __str__(self) -> str:
        return f'Position({str(self.x)}, {str(self.y)}, {str(self.rotated)})'
