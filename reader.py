'''
文件读取类
'''
import re
import operator

class base_reader:
    '文件读取器'
    def __init__(self,path):
        self._path=path
        f=self.file=open(path)
        self.eof=False
        self.main()

    def readline(self):
        '读取一行,并且判断是否读到文件尾'
        r=self.file.readline()
        if r=='':
            self.eof=True           
        return r

    def readuntil(self,end=''):
        '读取多行,直到碰到`end`'
        r=[]
        s=self.readline().strip()
        while s!=end and not self.eof:
            r.append(s)
            s=self.readline().strip()
        return r

    def main(self):
        '初始化时调用'
        pass
    
    def __del__(self):
        self.file.close()
    
class psi_reader(base_reader):
    '''用于读取QDTK工具包产生的restart/psi等文件
      - tape 保存在`self.tape`中
      - psi 以{time:psi}的形式在字典`self.psi`中
    restart文件的psi为`self.psi[0]`
    '''
    def main(self):
        self.tape=[]
        self.psi={}
        self.t=0
        while not self.eof:
            s=self.readline().strip()
            if s=='$tape':
                t=self.readuntil('')
                self.tape=[int(i) for i in t]
            elif s=='$time':
                s=self.readline().strip()
                self.t=eval(re.match('[0-9.]+',s).group())
            elif s=='$psi':
                t=self.readuntil('')
                for i in t:
                    if re.match('[a-zA-Z]+',i):
                        raise OSError('Letters found in the psi.')
                c=[complex(*eval(i)) for i in t]
                self.psi[self.t]=c
