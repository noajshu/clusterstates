class paulis_2d:
    """
    A bit string that represents a Pauli operator supported on a 2D grid.

    Attributes:
        lx (int): Length in the x direction
        ly (int): Length in the y direction
        operator (list of bool): A bit string representing the Pauli operator
        flying (bool): Two bits representing the Pauli operator of the flying
                       qubit
    """

    def __init__(self, lx, ly):
        self.lx = lx
        self.ly = ly
        self.operator = [False] * 2 * lx * ly
        self.flying = [False, False]

    def __repr__(self):
        title = "Pauli operator on a {} x {} grid:\n".format(self.lx, self.ly)
        if self.flying[0]:
            if self.flying[1]:
                title += "Flying qubit: Y.\n"
            else:
                title += "Flying qubit: X.\n"
        else:
            if self.flying[1]:
                title += "Flying qubit: Z.\n"
            else:
                title += "Flying qubit: I.\n"
        content = ""
        for i in range(self.ly):
            for j in range(self.lx):
                if self.operator[j + self.lx * i]:
                    if self.operator[j + self.lx * i + self.lx * self.ly]:
                        content += "Y "
                    else:
                        content += "X "
                else:
                    if self.operator[j + self.lx * i + self.lx * self.ly]:
                        content += "Z "
                    else:
                        content += "I "
            content += "\n"

        return title+content

    def co2idx_x(self, x, y):
        return x + y * self.lx
    
    def co2idx_z(self, x, y):
        return x + y * self.lx + self.lx * self.ly - 1;
    
    def swap(self, x1, y1, x2, y2):
        self.operator[self.co2idx_x(x1,y1)], self.operator[self.co2idx_x(x2,y2)] = self.operator[self.co2idx_x(x2,y2)], self.operator[self.co2idx_x(x1,y1)]
        self.operator[self.co2idx_z(x1,y1)], self.operator[self.co2idx_z(x2,y2)] = self.operator[self.co2idx_z(x2,y2)], self.operator[self.co2idx_z(x1,y1)]

    def cz(self, x1, y1, x2, y2):
        if self.operator[self.co2idx_x(x1,y1)]:
            self.operator[self.co2idx_z(x2,y2)] = not self.operator[self.co2idx_z(x2,y2)]
        if self.operator[self.co2idx_x(x2,y2)]:
            self.operator[self.co2idx_z(x1,y1)] = not self.operator[self.co2idx_z(x1,y1)]

    def swap_flying(self, x, y):
        self.operator[self.co2idx_x(x,y)], self.flying[0] = self.flying[0], self.operator[self.co2idx_x(x,y)]
        self.operator[self.co2idx_z(x,y)], self.flying[1] = self.flying[1], self.operator[self.co2idx_z(x,y)]

    def cz_flying(self, x, y):
        if self.operator[self.co2idx_x(x,y)]:
            self.flying[1] = not self.flying[1]
        if self.flying[0]:
            self.operator[self.co2idx_z(x,y)] = not self.operator[self.co2idx_z(x,y)]
        
    def knit(self, row):

        if row is 0:
            for x in range(self.lx):
                self.cz_flying(x, 0)
                self.swap_flying(x, 0)
            self.cz_flying(self.lx-1, 0)    
        else:
            for x in range(self.lx):
                self.cz_flying(x, row-1)
                self.cz_flying(x, row)
                self.swap_flying(x, row)
            self.cz_flying(self.lx-1, row)

    def knit_all(self):
        for y in range(self.ly):
            self.knit(y)
