from typing import Iterable, Optional, Tuple, Union
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.axes import Subplot
from matplotlib.figure import Figure

from Item import Item, Position


def draw_packing(pallet_width: int, pallet_height: int, items: Iterable[Item],
                positions: Iterable[Optional[Position]]) -> Tuple[Figure, Subplot]:
    fig, ax = plt.subplots()
    pallet = patches.Rectangle((0, 0), pallet_width, pallet_height, linewidth=5, facecolor='none', edgecolor='blue')
    ax.add_patch(pallet)
    ax.set_xlim(-1, pallet_width + 1)
    ax.set_ylim(-1, pallet_height + 1)
    for (item, position) in zip(items, positions):
        if position is not None:
            w = item.width if not position.rotated else item.height
            h = item.height if not position.rotated else item.width
            item = patches.Rectangle((position.x, position.y), w, h, linewidth=2, facecolor='gray', edgecolor='black')
            ax.add_patch(item)
    return fig, ax


def draw_packing_mass(items: Iterable[Item], positions: Iterable[Optional[Position]],
                      pallet_width: int, pallet_height: int,  x_tol: Union[int, float], y_tol: Union[int, float]) \
        -> Tuple[Figure, Subplot]:
    fig, ax = draw_packing(pallet_width, pallet_height, items, positions)
    cent_x = pallet_width / 2
    cent_y = pallet_height / 2
    ax.add_patch(
            patches.Rectangle((cent_x - x_tol, cent_y - y_tol), 2 * x_tol, 2 * y_tol, facecolor='none', edgecolor='blue')
            )

    x_mass = sum(
            (position.x + (item.width if not position.rotated else item.height) / 2) * item.mass for item, position in
            zip(items, positions) if position is not None) / sum(
            item.mass for item, position in zip(items, positions) if position is not None)

    y_mass = sum(
            (position.y + (item.height if not position.rotated else item.width) / 2) * item.mass for item, position in
            zip(items, positions) if position is not None) / sum(
            item.mass for item, position in zip(items, positions) if position is not None)
    ax.scatter(x_mass, y_mass, color='blue', s=10, zorder=2)

    fig.savefig("resultPacking.png")
    return fig, ax

