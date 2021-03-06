#!/usr/bin/env python3.6
import numpy as np
from numpy import r_
from ase import Atoms
from ase.io import read
import ase.build
from useful_functions import logger_setup
import logging, sys

class lmp_data_hex():
    """Class to contain the LAMMPS data of hexagonal system in studying MoS2/Gr heterostructures"""

    def __init__(self, filename=None):
        """Deinfe the variables of the object"""

        self.comment = "#Empty lammps data"

        self.n_atm = 0
        self.n_type = 0

        self.xlo, self.xhi = 0, 0
        self.ylo, self.yhi = 0, 0
        self.zlo, self.zhi = 0, 0
        self.xy, self.xz, self.yz = 0, 0, 0

        self.mass = {} # Type-mass

        self.style = ""

        self.pos = []
        self.types = []

        if filename:
            with open(filename, 'r') as in_str:
                self.load_stream(in_str)

    def set_cell(self, cell):
        """Set LAMMPS box from cell. VERY VERY specific for our hexagonal geometry"""
        # See get cell definition, must be compliant to that.
        self.xlo, self.ylo, self.zlo = cell[2][0], cell[2][1], cell[0][2]
        self.xhi, self.yhi, self.zhi = cell[0][0], cell[1][1], cell[2][2]
        self.xy, self.xz, self.yz = cell[1][0], 0, 0

    def get_ase(self, types_dic={}):
        """Return a ASE Atoms object

        LAMMPS types to elements via types_dict={type: chem_symbol}"""
        if types_dic:
            symbols = [types_dic[x] for x in self.types]
        else:
            symbols = [x for x in self.types]
        c = Atoms(symbols = symbols,
                  positions=self.pos,
                  cell=self.get_cell(),
                  masses=[self.mass[i] for i in self.types])
        return(c)

    def load_from_ase(self, ase_obj, style="atomic", types_dic={}):
        """Load from a ASE Atoms object. Very specific:
        Cell must be hexagonal
        If not type-number dict {'Symbol': number} is given, atom types are arge given as number following the checmial symbol appearence
        Style can be specified, default atomic.
        Comment line will be Symbol-type mapping and string representation of input Atoms
        Masses are taken from ASE Atoms
        """
        self.set_cell(ase_obj.cell)
        self.pos = ase_obj.positions
        # Get list of atom types.
        types = [atm.symbol for atm in ase_obj]
        if not len(types_dic):
            # Associate a number to each type.
            types_dic = {symbol: i+1
                         for i, symbol in enumerate(set(types))
            }

        for a in ase_obj:
            self.mass[types_dic[a.symbol]] = a.mass
        # Use the association rule to give numberical types to lammps atoms
        self.types = [types_dic[s] for s in types]
        self.n_atm = len(ase_obj)


        self.n_type = len(set(self.types))
        # !! HARDCODED !!
        self.style = style
        self.comment = "#Type-symbol_map: "+" ".join(["%s-%i"%(s,i)
                                                      for s,i in types_dic.items()])+" "+ase_obj.__str__()

    def copy(self):
        """Return a copy of the object"""

        cp = lmp_data()
        cp.comment = self.comment
        cp.n_atm = self.n_atm
        cp.n_type = self.n_type
        cp.xlo, cp.xhi = self.xlo, self.xhi
        cp.ylo, cp.yhi = self.ylo, self.yhi
        cp.zlo, cp.zhi = self.zlo, self.zhi
        cp.xy, cp.xz, cp.yz = self.xy, self.xz, self.yz
        cp.mass =  self.mass
        cp.style = self.style
        cp.pos = self.pos
        cp.types = self.types

        return cp

    def load_stream(self, in_stream):
        """Load from an opened stream"""
        c_log = logger_setup(__name__)
        #c_log.setLevel(logging.DEBUG)

        c = -1
        for l in in_stream.readlines():
            # Increment counter
            c += 1

            #print(c,")", l.strip())
            c_log.debug("On line %i %s" , c, l.strip())
            if c in [1, 4, 9, 10, 11, 11+self.n_type, 11+self.n_type+2, 11+self.n_type+3]:
                c_log.debug("skip")
                #continue
            if c == 0:
                self.comment = l.strip()
                c_log.debug("Set comment")
            if c == 2:
                self.n_atm = int(l.split()[0])
                c_log.debug("Set atm num to %i", self.n_atm)
            if c == 3:
                self.n_type = int(l.split()[0])
                c_log.debug("Set atm types to %i", self.n_type)
            if c == 5:
                self.xlo, self.xhi = [float(x) for x in l.split()[0:2]]
                c_log.debug("Set xlo/hi %f %f", self.xlo, self.xhi)
            if c == 6:
                self.ylo, self.yhi = [float(x) for x in l.split()[0:2]]
                c_log.debug("Set ylo/hi %f %f", self.ylo, self.yhi)
            if c == 7:
                self.zlo, self.zhi = [float(x) for x in l.split()[0:2]]
                c_log.debug("Set zlo/hi %f %f ", self.zlo, self.zhi)
            if c == 8:
                self.xy, self.xz, self.yz = [float(x) for x in l.split()[0:3]]
                c_log.debug("Set xy/xz/yz %f %f %f", self.xy, self.xz, self.yz)
            if c > 11 and c <= 11+self.n_type:
                l = l.split()
                c_type, c_mass = int(l[0]), float(l[1])
                self.mass[c_type] = c_mass
                c_log.debug("Mass type %i = %f", c_type, c_mass)
            if c == 11+self.n_type+2:
                self.style = l.split()[2]
                c_log.debug("Atoms style is %s", self.style)
            if c > 11+self.n_type+3 and c <= 11+self.n_type+3+self.n_atm:
                c_id, c_type = [int(x) for x in l.split()[0:2]]
                c_pos = [float(x) for x in l.split()[2:5]]
                self.types.append(c_type)
                self.pos.append(c_pos)
                c_log.debug("Atm type %i, position %f %f %f",  c_id, *c_pos)

        # Check we are dealing with our standard geometry, rise error otherwise
        if self.xlo != 0 or self.ylo != 0 or self.zlo != 0:
            c_log.error("Input geometry not compatible with expected hexagonal cell: 0 xhi 0 yhi 0 zhi, xy 0 0")
            raise ValueError

    def print_to_stream(self, out_stream=sys.stdout):
        """"""

        # Print comment line and space
        print(self.comment, file=out_stream)
        print("", file=out_stream)

        # Number of atoms and types. Mandatory space
        print(self.n_atm, "atoms", file=out_stream)
        print(self.n_type, "atom types", file=out_stream)
        print("", file=out_stream)

        # Box dimension
        print("%-20.15g %-20.15g xlo xhi" % (self.xlo, self.xhi), file=out_stream)
        print("%-20.15g %-20.15g ylo yhi" % (self.ylo, self.yhi), file=out_stream)
        print("%-20.15g %-20.15g zlo zhi" % (self.zlo, self.zhi), file=out_stream)
        print("%-20.15g %-20.15g %-20.15g xy xz yz" % (self.xy, self.xz, self.yz), file=out_stream)
        print("", file=out_stream)

        # Masses of the system
        print("Masses", file=out_stream)
        print("", file=out_stream)
        for k, v in self.mass.items():
            print("%-10i %20.15f" % (k, v), file=out_stream)
        print("", file=out_stream)

        # Style of atoms
        print("Atoms # %s" % self.style, file=out_stream)
        print("", file=out_stream)

        # Atoms indeces, types and relative positions
        for i, c_pos in enumerate(self.pos):
            #print("%6i %6i %20.15 %20.15f %20.15f" % (i, types[i], *c_pos))
            #print(i) #, self.types[i], *c_pos, file=sys.stderr)
            print("%-6i %-6s %-25.15g  %-25.15g %-25.15g" % (i+1, str(self.types[i]), *c_pos), file=out_stream)

    def print_to_file(self, filename):
        """"""
        with open(filename, 'w') as out_file:
            self.print_to_stream(out_file)

    def get_cell(self):
        """Conversion from LAMMPS (weird) lo-hi system to hexagonal cell (a, b, c)"""
        cell = np.array([[self.xhi, self.ylo, self.zlo],
                         [self.xy, self.yhi, self.zlo],
                         [self.xlo, self.ylo, self.zhi]])
        return cell

#================================================================================================================================

class lmp_data_rect():
    """Class to contain the LAMMPS data of rectangular cell"""

    def __init__(self, filename=None):
        """Deinfe the variables of the object"""

        self.log = logger_setup(__name__)
        self.log.setLevel(logging.INFO)

        self.comment = "#Empty lammps data"

        self.n_atm = 0
        self.n_type = 0

        self.xlo, self.xhi = 0, 0
        self.ylo, self.yhi = 0, 0
        self.zlo, self.zhi = 0, 0
        self.xy, self.xz, self.yz = 0, 0, 0

        self.mass = {} # Type-mass

        self.style = ""

        self.pos = []
        self.types = []

        if filename:
            with open(filename, 'r') as in_str:
                self.load_stream(in_str)

    def set_cell2box(self, cell):
        """Set LAMMPS box from cell. VERY VERY specific for our ortorombic geometry"""
        # See get cell definition, must be compliant to that.
        self.xlo, self.ylo, self.zlo = 0, 0, 0
        self.xhi, self.yhi, self.zhi = cell[0][0], cell[1][1], cell[2][2]
        self.xy, self.xz, self.yz = 0,0,0

    def get_ase(self, types_dic={}, translate_pos=False):
        """Return a ASE Atoms object

        LAMMPS types to elements via types_dict={type: chem_symbol}"""
        if types_dic:
            symbols = [types_dic[x] for x in self.types]
        else:
            symbols = [x for x in self.types]
        # Translate positions in case cell was weird
        v = np.zeros(3)
        if translate_pos:
            v = np.array([self.xlo, self.ylo, self.zlo])
            self.log.info("Translate position of %.2f %.2f %.2f" % tuple(v))

        c = Atoms(symbols = symbols,
                  positions=self.pos-v,
                  cell=self.get_cell(),
                  masses=[self.mass[i] for i in self.types])
        return(c)

    def load_from_ase(self, ase_obj, style="atomic", types_dic={}):
        """Load from a ASE Atoms object. Very specific:
        Cell must be hexagonal
        If not type-number dict {'Symbol': number} is given, atom types are arge given as number following the checmial symbol appearence
        Style can be specified, default atomic.
        Comment line will be Symbol-type mapping and string representation of input Atoms
        Masses are taken from ASE Atoms
        """
        self.set_cell2box(ase_obj.cell)
        self.pos = ase_obj.positions
        # Get list of atom types.
        types = [atm.symbol for atm in ase_obj]
        if not len(types_dic):
            # Associate a number to each type.
            types_dic = {symbol: i+1
                         for i, symbol in enumerate(set(types))
            }

        for a in ase_obj:
            self.mass[types_dic[a.symbol]] = a.mass
        # Use the association rule to give numberical types to lammps atoms
        self.types = [types_dic[s] for s in types]
        self.n_atm = len(ase_obj)


        self.n_type = len(set(self.types))
        # !! HARDCODED !!
        self.style = style
        self.comment = "#Type-symbol_map: "+" ".join(["%s-%i"%(s,i)
                                                      for s,i in types_dic.items()])+" "+ase_obj.__str__()

    def copy(self):
        """Return a copy of the object"""

        cp = lmp_data()
        cp.comment = self.comment
        cp.n_atm = self.n_atm
        cp.n_type = self.n_type
        cp.xlo, cp.xhi = self.xlo, self.xhi
        cp.ylo, cp.yhi = self.ylo, self.yhi
        cp.zlo, cp.zhi = self.zlo, self.zhi
        cp.xy, cp.xz, cp.yz = self.xy, self.xz, self.yz
        cp.mass =  self.mass
        cp.style = self.style
        cp.pos = self.pos.copy()
        cp.types = self.types

        return cp

    def load_stream(self, in_stream):
        """Load from an opened stream"""
        #self.log.setLevel(logging.DEBUG)

        c = -1
        for l in in_stream.readlines():
            # Increment counter
            c += 1

            #print(c,")", l.strip())
            self.log.debug("On line %i %s" , c, l.strip())
            if c in [1, 4, 9, 10, 11, 11+self.n_type, 11+self.n_type+2, 11+self.n_type+3]:
                self.log.debug("skip")
                #continue
            if c == 0:
                self.comment = l.strip()
                self.log.debug("Set comment")
            if c == 2:
                self.n_atm = int(l.split()[0])
                self.log.debug("Set atm num to %i", self.n_atm)
            if c == 3:
                self.n_type = int(l.split()[0])
                self.log.debug("Set atm types to %i", self.n_type)
            if c == 5:
                self.xlo, self.xhi = [float(x) for x in l.split()[0:2]]
                self.log.debug("Set xlo/hi %f %f", self.xlo, self.xhi)
            if c == 6:
                self.ylo, self.yhi = [float(x) for x in l.split()[0:2]]
                self.log.debug("Set ylo/hi %f %f", self.ylo, self.yhi)
            if c == 7:
                self.zlo, self.zhi = [float(x) for x in l.split()[0:2]]
                self.log.debug("Set zlo/hi %f %f ", self.zlo, self.zhi)
            if c == 8:
                self.xy, self.xz, self.yz = [float(x) for x in l.split()[0:3]]
                self.log.debug("Set xy/xz/yz %f %f %f", self.xy, self.xz, self.yz)
            if c > 11 and c <= 11+self.n_type:
                l = l.split()
                c_type, c_mass = int(l[0]), float(l[1])
                self.mass[c_type] = c_mass
                self.log.debug("Mass type %i = %f", c_type, c_mass)
            if c == 11+self.n_type+2:
                self.style = l.split()[2]
                self.log.debug("Atoms style is %s", self.style)
            if c > 11+self.n_type+3 and c <= 11+self.n_type+3+self.n_atm:
                c_id, c_type = [int(x) for x in l.split()[0:2]]
                c_pos = [float(x) for x in l.split()[2:5]]
                self.types.append(c_type)
                self.pos.append(c_pos)
                self.log.debug("Atm type %i, position %f %f %f",  c_id, *c_pos)

        # Check we are dealing with our standard geometry, rise error otherwise
        if self.xy != 0 or self.xz != 0 or self.yz != 0:
            raise ValueError("Input geometry not compatible with expected rectangular cell: xy=0 xz=0 yz=0")

                # Check we are dealing with our standard geometry, rise error otherwise
        if self.xlo != 0 or self.ylo != 0 or self.zlo != 0:
            self.log.warning("Input geometry incompatible with ASE cell (xlo!=0 ylo!=0 zlo!=0). Shift on export.")

    def print_to_stream(self, out_stream=sys.stdout):
        """"""

        # Print comment line and space
        print(self.comment, file=out_stream)
        print("", file=out_stream)

        # Number of atoms and types. Mandatory space
        print(self.n_atm, "atoms", file=out_stream)
        print(self.n_type, "atom types", file=out_stream)
        print("", file=out_stream)

        # Box dimension
        print("%-20.15g %-20.15g xlo xhi" % (self.xlo, self.xhi), file=out_stream)
        print("%-20.15g %-20.15g ylo yhi" % (self.ylo, self.yhi), file=out_stream)
        print("%-20.15g %-20.15g zlo zhi" % (self.zlo, self.zhi), file=out_stream)
        print("%-20.15g %-20.15g %-20.15g xy xz yz" % (self.xy, self.xz, self.yz), file=out_stream)
        print("", file=out_stream)

        # Masses of the system
        print("Masses", file=out_stream)
        print("", file=out_stream)
        for k, v in self.mass.items():
            print("%-10i %20.15f" % (k, v), file=out_stream)
        print("", file=out_stream)

        # Style of atoms
        print("Atoms # %s" % self.style, file=out_stream)
        print("", file=out_stream)

        # Atoms indeces, types and relative positions
        for i, c_pos in enumerate(self.pos):
            #print("%6i %6i %20.15 %20.15f %20.15f" % (i, types[i], *c_pos))
            #print(i) #, self.types[i], *c_pos, file=sys.stderr)
            print("%-6i %-6s %-25.15g  %-25.15g %-25.15g" % (i+1, str(self.types[i]), *c_pos), file=out_stream)

    def print_to_file(self, filename):
        """"""
        with open(filename, 'w') as out_file:
            self.print_to_stream(out_file)

    def get_cell(self):
        """Conversion from LAMMPS (weird) lo-hi system to hexagonal cell (a, b, c)"""
        if self.xlo != 0 or self.ylo != 0 or self.zlo != 0:
            self.log.warning("xlo!=0 ylo!=0 zlo!=0. Cell is shifted to origin. You might want to shift positions of [xlo, ylo, zlo]")
        cell = np.array([[self.xhi-self.xlo, 0,                 0],
                         [0,                 self.xhi-self.ylo, 0],
                         [0,                 0,                 self.zhi-self.zlo]])
        return cell
