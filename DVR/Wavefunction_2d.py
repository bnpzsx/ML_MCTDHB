# coding=utf-8
import utils

#class Wavefunction(object):
class Wavefunction:
    """
    Wavefunction object implementation in the python world

    """

    def __init__(self,**opts):#通过赋值的方式以字典的形式传入参数
        """
        Initialize the wavefunction
        """
        self._tape = None
        self.n = None
        self.m = None
        self.mx = None
        self.my = None
        self.psilen = 0
        self.PSI = []
        
        # start from tape
        tape = opts.get("tape",None)
        if tape is not None:
            self._init_from_tape(tape)

    def _init_from_tape(self,tape):
        self._tape = tape
      
        i = 1
        pnode = []
        subnode = []
        subnode_s = []
        while tape[i] != -2:
            if tape[i] > 1:
                if tape[i+3] == -1:
                    pnode.append(tape[i])
                    pnode.append(tape[i+2])
                    i = i +3
                else:
                    a = tape[i] + i + 3
                    b = tape[i+2:a]
                    subnode.extend(b)
                    i = i + tape[i] + 2
            if tape[i] == 1:
                subnode_s.append(tape[i+2])
                i += 3
            if tape[i] <= 0:
                if tape[i] == -1:  # go down
                    i += 2
                if tape[i] == 0:  # jump up
                    i += 3
                if tape[i] == -2:  # exit
                    break
                                          
        self.n = pnode[0]
        self.m = pnode[1]
        self.mx = subnode[0]
        self.my = subnode[1]
        
    def _PSI_ini(self,ns):
        """
        initializes coefficients of initial wave function of a bosonic layer.
        """
        ns_index_ini = utils.ns_index(ns)
        len_ns = utils.num_combination(self.n + self.m-1,self.n)
        PSI_1 = [0]*len_ns
        PSI_1[ns_index_ini] = 1
        return PSI_1
        
    def _PSI_ml_dict(self):
        k_dict={}
        k = 0
        for y in range(1,self.my+1):
            for x in range(1,self.mx+1):
                k=k+1
                k_dict[x,y] = k
        return k_dict
        
    def _chang_to_dict(self,listspf):
        dictspf = {}
        for ever in enumerate(listspf):
            dictspf[str(ever[1])] = ever[0]+1
        return dictspf
        
    def _PSI_ml(self,startSPF,spfx,spfy):
        """
        initializes coefficients of initial wave function of sublayer.
        """
        len_subPSI_2 = self.mx*self.my
        PSI_2 = []
        
        spfxdict = self._chang_to_dict(spfx)
        spfydict = self._chang_to_dict(spfy)
        k_dict = self._PSI_ml_dict()
        
        for ever in startSPF:
            PSI_2 = [0]*len_subPSI_2
            x_ind = spfxdict[str(ever[0])]
            y_ind = spfydict[str(ever[1])]
            spf_ind = k_dict[x_ind,y_ind]
            PSI_2[spf_ind] = 1
            PSI_2.extend(PSI_2)
        
        return PSI_2
        
    def _chang_to_list(self,list_arrayspf):
        listspf = []
        for ever in list_arrayspf:
            arr = ever.tolist()
            listspf.append(arr)
        return listspf
        
    def _chang_to_list_2(self,list_list_arrayspf):
        listspf_1 = []
        listspf = []
        for list_arrayspf in list_list_arrayspf:
            for ever in list_arrayspf:
                arr = ever.tolist()
                listspf_1.append(arr)
            listspf.append(listspf_1)
        return listspf
        
    def init_coef_sing_spec_B(self,ns,startSPF,spfx,spfy):
        PSI_1 = self._PSI_ini(ns)
        self.PSI.extend(PSI_1)
        
        startSPF = self._chang_to_list_2(startSPF)
        spfx = self._chang_to_list(spfx)
        spfy = self._chang_to_list(spfy)
        
        PSI_2 = self._PSI_ml(startSPF,spfx,spfy)
        self.PSI.extend(PSI_2)
        
        PSI_mx = []
        for spf_x in spfx:
            PSI_mx.extend(spf_x)
        self.PSI.extend(PSI_mx)
        
        PSI_my = []
        for spf_y in spfy:
            PSI_my.extend(spf_y)
        self.PSI.extend(PSI_my)
        
        self.psilen = len(self.PSI)
        
        
    def createWfnFile(self,fname):
        """
        Create a wfn file.
        """
        fout = open(fname,'w')

        fout.write('$tape\n')
        for entry in self._tape:
            fout.write(' %5i\n' % (entry,))
        fout.write('\n')

        fout.write('$psi\n')
        for i in range(self.psilen):
            fout.write(' (%15.12f,%15.12f)\n' % (self.PSI[i].real,self.PSI[i].imag))
        
        fout.close()

#test
if __name__=='__main__':
    tape = (-10,
         2, 1, 3, 
        -1, 1,
         2, 0, 4, 5,
        -1, 1,
         1, 0, 10,
         0,-1, 2,   
         1, 0, 10,  
        -2)
    wfn = Wavefunction(tape=tape)
    print(wfn.n)
    print(wfn.m)
    print(wfn.mx)
    print(wfn.my)

