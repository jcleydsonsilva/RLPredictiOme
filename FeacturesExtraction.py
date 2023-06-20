class FeacturesExtraction:

    mapVolume = {}
    mapMass = {}
    mapHydro = {}

    seq = ''
    length = 0.0
    Nonpolar_Aliphatic = 0.0
    Aromatic = 0.0
    Polar_Uncharged = 0.0
    Positively_Charged = 0.0
    Negatively_Charged  = 0.0
    mass = 0.0
    volume = 0.0
    hydro = 0.0

    def __init__(self, seq):
        self.seq = seq    # instance variable unique to each instance

    def calculate_proportion(self):
        self.mapVolume['A'] = 88.6
        self.mapVolume['R'] = 173.4
        self.mapVolume['N'] = 114.1
        self.mapVolume['D'] = 111.1
        self.mapVolume['C'] = 108.5
        self.mapVolume['Q'] = 143.8
        self.mapVolume['E'] = 138.4
        self.mapVolume['G'] = 60.1
        self.mapVolume['H'] = 153.2
        self.mapVolume['I'] = 166.7
        self.mapVolume['L'] = 166.7
        self.mapVolume['K'] = 168.6
        self.mapVolume['M'] = 162.9
        self.mapVolume['F'] = 189.9
        self.mapVolume['P'] = 112.7
        self.mapVolume['S'] = 89.0
        self.mapVolume['T'] = 116.1
        self.mapVolume['W'] = 227.8
        self.mapVolume['Y'] = 193.6
        self.mapVolume['V'] = 140.0

        self.mapMass['A'] = 89.0
        self.mapMass['R'] = 174.0
        self.mapMass['N'] = 132.0
        self.mapMass['D'] = 133.0
        self.mapMass['C'] = 121.0
        self.mapMass['Q'] = 146.0
        self.mapMass['E'] = 147.0
        self.mapMass['G'] = 75.0
        self.mapMass['H'] = 155.0
        self.mapMass['I'] = 131.0
        self.mapMass['L'] = 131.0
        self.mapMass['K'] = 146.0
        self.mapMass['M'] = 149.0
        self.mapMass['F'] = 165.0
        self.mapMass['P'] = 115.0
        self.mapMass['S'] = 105.0
        self.mapMass['T'] = 119.0
        self.mapMass['W'] = 204.0
        self.mapMass['Y'] = 181.0
        self.mapMass['V'] = 117.0

        self.mapHydro['A'] = 1.8
        self.mapHydro['R'] = -4.5
        self.mapHydro['N'] = -3.5
        self.mapHydro['D'] = -3.5
        self.mapHydro['C'] = 2.5
        self.mapHydro['Q'] = -3.5
        self.mapHydro['E'] = -3.5
        self.mapHydro['G'] = -0.4
        self.mapHydro['H'] = -3.2
        self.mapHydro['I'] = 4.5
        self.mapHydro['L'] = 3.8
        self.mapHydro['K'] = -3.9
        self.mapHydro['M'] = 1.9
        self.mapHydro['F'] = 2.8
        self.mapHydro['P'] = -1.6
        self.mapHydro['S'] = -0.8
        self.mapHydro['T'] = -0.7
        self.mapHydro['W'] = -0.9
        self.mapHydro['Y'] = -1.3
        self.mapHydro['V'] = 4.2 

        self.length = len(self.seq)

        for aa in self.seq:
            if ( 'G' in aa or 'A' in aa  or 'P'  in aa  or 'V'  in aa  or 'L'  in aa  or 'I'  in aa or 'M'  in aa ):
                self.Nonpolar_Aliphatic += 1
            elif ( 'F'  in aa or 'Y' in aa or 'W'  in aa ):
                self.Aromatic += 1
            elif ( 'S'  in aa or 'T'  in aa  or 'C'  in aa  or 'N' in aa  or 'Q'  in aa  ):
                self.Polar_Uncharged +=1
            elif ( 'K' in aa or 'H'  in aa or 'R'  in aa ):
                self.Positively_Charged +=1
            elif ( 'D' in aa or  'E' in aa):
                self.Negatively_Charged +=1
            else:
                print ('')

            if (aa in self.mapHydro):
                self.hydro += self.mapHydro[aa]
                self.mass += self.mapMass[aa]
                self.volume += self.mapVolume[aa]

        self.Nonpolar_Aliphatic  = 0 if self.Nonpolar_Aliphatic == 0 else self.Nonpolar_Aliphatic / len(self.seq)
        self.Aromatic            = 0 if self.Aromatic  == 0 else self.Aromatic / len(self.seq)
        self.Polar_Uncharged     = 0 if self.Polar_Uncharged == 0 else self.Polar_Uncharged / len(self.seq)
        self.Positively_Charged  = 0 if self.Positively_Charged == 0 else self.Positively_Charged / len(self.seq)
        self.Negatively_Charged  = 0 if self.Negatively_Charged == 0  else self.Negatively_Charged / len(self.seq)
        self.mass                = 0 if self.mass == 0 else self.mass / len(self.seq)
        self.volume              = 0 if self.volume == 0 else self.volume / len(self.seq)

        return 'Ok'
        
    def get_sequence(self):
        return self.seq
   
    def get_length(self):
        return self.length
   
    def get_nonpolar_aliphatic(self):
        return self.Nonpolar_Aliphatic
   
    def get_aromatic(self):
        return self.Aromatic
   
    def get_polar_uncharged(self):
        return self.Polar_Uncharged
   
    def get_positively_charged(self):
        return self.Positively_Charged

    def get_negatively_charged(self):
        return self.Negatively_Charged
   
    def get_mass(self):
        return self.mass
   
    def get_volume(self):
        return self.volume
   
    def get_hydro(self):
        return self.hydro
    
    def get_all_features(self):
        all_features = []
        all_features.append(self.length)
        all_features.append(self.Nonpolar_Aliphatic)
        all_features.append(self.Aromatic)
        all_features.append(self.Polar_Uncharged)
        all_features.append(self.Positively_Charged)
        all_features.append(self.Negatively_Charged)
        all_features.append(self.mass)
        all_features.append(self.volume)
        all_features.append(self.hydro)
        #print (all_features)
        return all_features

   
