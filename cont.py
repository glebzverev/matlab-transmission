import matplotlib.pyplot as plt
import numpy as np

def contour_legend(hCont, txtCont=None, locLeg='best', hOther=None, txtOther=None, posOther='after'):
    """
    Adds a curve description legend to a contour plot.
    """
    if txtCont is None:
        txtCont = []
    if hOther is None:
        hOther = []
    if txtOther is None:
        txtOther = []

    hCont = np.array(hCont)
    hOther = np.array(hOther)

    # Check if all contour objects are in the same axes
    hAx = check_same_axes(hCont)
    if hAx is None:
        raise ValueError("All contour objects should be in the same axes")

    # Get contour curve specification
    curveSpec = contour_curves(hCont)
    curveSpec = flat_to_color(hAx, curveSpec)

    # Get handles to curve specification invisible dummy lines or patches
    hCurve, txtCurve, levMat = get_curve_handles(hAx, curveSpec, txtCont)

    # Draw legend
    hLeg = None
    if locLeg != 'notdrawn':
        hLeg = draw_leg(hAx, locLeg, hCurve, txtCurve, hOther, txtOther, posOther)

    return hLeg, levMat


def draw_leg(hAx, locLeg, hCurve, txtCurve, hOther, txtOther, posOther):
    """
    Draw the legend
    """
    hItem = list(hCurve)
    txtItem = list(txtCurve)
    if hOther:
        if posOther == 'before':
            hItem = list(hOther) + hItem
            txtItem = list(txtOther) + txtItem
        else:
            hItem = hItem + list(hOther)
            txtItem = txtItem + list(txtOther)

    hLeg = hAx.legend(hItem, txtItem, loc=locLeg)
    return hLeg


def get_curve_handles(hAx, curveSpec, txtCont):
    """
    Get handles to curve specification invisible dummy lines or patches
    """
    nCurves = len(curveSpec)
    unLev = np.unique([c[0] for c in curveSpec])
    levMat = []
    for curve in curveSpec:
        currLev = curve[0]
        currType = curve[4]
        if currType == 'curves':
            levMat.append([currLev, currLev])
        elif currType == 'filled':
            currLevInd = np.where(unLev == currLev)[0][0]
            if currLevInd == len(unLev) - 1:
                levMat.append([currLev, np.inf])
            else:
                levMat.append([currLev, unLev[currLevInd + 1]])

    levMat = np.array(sorted(set(tuple(row) for row in levMat)))

    if txtCont:
        if len(txtCont) > len(levMat):
            txtCont = txtCont[:len(levMat)]
        elif len(txtCont) < len(levMat):
            txtCont += [str(lev) for lev in levMat[len(txtCont):]]

    else:
        txtCont = []
        for lev in levMat:
            if lev[0] == lev[1]:
                txtCont.append(str(lev[0]))
            else:
                if not np.isfinite(lev[1]):
                    txtCont.append(f'> {lev[0]}')
                else:
                    txtCont.append(f'{lev[0]} to {lev[1]}')

    hCurve = []
    txtCurve = []
    for curve in curveSpec:
        currLev = curve[0]
        currType = curve[4]
        if currType == 'curves':
            levInd = np.where((levMat[:, 0] == currLev) & (levMat[:, 1] == currLev))[0][0]
            h = hAx.plot(np.nan, np.nan, linestyle=curve[1], linewidth=curve[2], color=curve[3])
        elif currType == 'filled':
            levInd = np.where((levMat[:, 0] == currLev) & (levMat[:, 1] > currLev))[0][0]
            h = hAx.fill_between([], [], color=curve[3])
        hCurve.append(h[0])
        txtCurve.append(txtCont[levInd])

    return hCurve, txtCurve, levMat


def flat_to_color(hAx, curveSpec):
    """
    Set 'flat' to a color specification
    """
    unLev = np.unique([c[0] for c in curveSpec])
    colMap = plt.get_cmap(hAx.get_cmap().name)
    nCol = colMap.N
    colLim = hAx.get_clim()

    if len(unLev) == 1:
        colNo = [nCol // 2] * len(unLev)
    else:
        colNo = np.round(np.interp(unLev, colLim, [1, nCol])).astype(int)

    colNo[unLev < colLim[0]] = 1
    colNo[unLev > colLim[1]] = nCol
    colNo = np.clip(colNo, 1, nCol)

    for i, curve in enumerate(curveSpec):
        if isinstance(curve[3], str) and curve[3] == 'flat':
            levNo = np.where(unLev == curve[0])[0][0]
            curveSpec[i] = list(curve)
            curveSpec[i][3] = colMap(colNo[levNo] - 1)

    curveSpec = [list(x) for x in set(tuple(row) for row in curveSpec)]
    return curveSpec


def contour_curves(hCont):
    """
    Get contour curve specification
    """
    plt.draw()
    curveSpec = []
    for cont in hCont:
        lineStyle = cont.get_linestyle()
        lineWidth = cont.get_linewidth()
        if hasattr(cont, 'collections') and cont.collections:
            type_ = 'filled'
            lineColor = 'flat'
        else:
            type_ = 'curves'
            lineColor = cont.get_edgecolor()

        levVec = cont.levels
        for lev in levVec:
            curveSpec.append([lev, lineStyle, lineWidth, lineColor, type_])

    return curveSpec


def check_same_axes(hVec):
    """
    Check if all handles belong to the same axes
    """
    hAxes = [h.axes for h in hVec]
    if all(ax == hAxes[0] for ax in hAxes):
        return hAxes[0]
    else:
        return None


# Example usage:
X, Y, Z = np.meshgrid(np.linspace(-3, 3, 100), np.linspace(-3, 3, 100), np.linspace(-3, 3, 100))
fig, ax = plt.subplots()
CS = ax.contour(X[:, :, 0], Y[:, :, 0], Z[:, :, 0])
contour_legend(CS.collections)

plt.show()
