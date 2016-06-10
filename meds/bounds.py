class Bounds(object):
    """
    A bounding box

    Based on Bounds from deswl_shapelets
    """
    def __init__(self, rowmin, rowmax, colmin, colmax):
        self.rowmin=rowmin
        self.rowmax=rowmax
        self.colmin=colmin
        self.colmax=colmax

    def contains_point(self, row, col):
        """
        is the point contained within the Bounds?

        parameters
        ----------
        row,col: floats
            The position to check
        """
        return (row <= self.rowmax
                and row >= self.rowmin
                and col <= self.colmax
                and col >= self.colmin)

    def contains_points(self, rows, cols):
        """
        is the point contained within the Bounds?

        parameters
        ----------
        rows,cols: arrays
            The positions to check
        """
        return (  (rows <= self.rowmax)
                & (rows >= self.rowmin)
                & (cols <= self.colmax)
                & (cols >= self.colmin) )


    def contains_bounds(self, bounds):
        """
        Is the input bounds entirely contained?

        parameters
        ----------
        bounds: Bounds
            The bounding box to check
        """
        assert isinstance(bounds,Bounds)
        return (bounds.rowmin >= self.rowmin
                and bounds.rowmax <= self.rowmax
                and bounds.colmin >= self.colmin
                and bounds.colmax <= self.colmax)

    def intersects_bounds(self, bounds):
        """
        Is there non-zero intersection between the two Bounds

        parameters
        ----------
        bounds: Bounds
            The bounding box to check
        """
        assert isinstance(bounds,Bounds)

        return ( not (bounds.rowmin >= self.rowmax)
                and not (bounds.rowmax <= self.rowmin)
                and not (bounds.colmin >= self.colmax)
                and not (bounds.colmax <= self.colmin) )

    def expand_point(self, row, col):
        """
        expand the bounds to include the input point

        parameters
        ----------
        row,col: floats
            The position to check
        """
        if row < self.rowmin:
            self.rowmin = row
        elif row > self.rowmax:
            self.rowmax = row

        if col < self.colmin:
            self.colmin = col
        elif col > self.colmax:
            self.colmax = col

    def expand_points(self, rows, cols):
        """
        expand the bounds to include the input points

        parameters
        ----------
        rows,cols: arrays
            The positions to check
        """

        rowmin = rows.min()
        rowmax = rows.max()
        colmin = cols.min()
        colmax = cols.max()

        if rowmin < self.rowmin:
            self.rowmin = rowmin
        if rowmax > self.rowmax:
            self.rowmax = rowmax

        if colmin < self.colmin:
            self.colmin = colmin
        if colmax > self.colmax:
            self.colmax = colmax



    def expand_bounds(self, bounds):
        """
        expand the bounds to include the input Bounds

        parameters
        ----------
        bounds: Bounds
            The bounding box to check
        """
        assert isinstance(bounds,Bounds)

        if bounds.rowmin < self.rowmin:
            self.rowmin = bounds.rowmin
        elif bounds.rowmax > self.rowmax:
            self.rowmax = bounds.rowmax

        if bounds.colmin < self.colmin:
            self.colmin = bounds.colmin
        elif bounds.colmax > self.colmax:
            self.colmax = bounds.colmax

    def __repr__(self):
        mess="Bounds(rowmin=%g, rowmax=%g, colmin=%g, colmax=%g)"
        mess = mess % (self.rowmin, self.rowmax, self.colmin, self.colmax)
        return mess


