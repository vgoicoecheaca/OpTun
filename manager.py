from plotter import Plotter 
from helper import Helper
from config import Config
from branches import Branches

class Manager():
    def __init__(self,**kwargs):
        self.config   = Config()
        self.helper   = Helper(manager=self,**kwargs)
        self.plotter  = Plotter(manager=self,**kwargs)
        self.branches = Branches(manager=self,**kwargs)
