from ovito.modifiers import ExpressionSelectionModifier, ClearSelectionModifier
from ovito.pipeline import ModifierInterface
from ovito.data import DataCollection
from ovito.traits import OvitoObject
from ovito.vis import VectorVis
from traits.api import Float
from ovito.gui import UtilityInterface

import numpy as np

class S2TXA(ModifierInterface):

    alpha = Float(0.1, label="Shear strain threshold")

    slip_direction_vis = OvitoObject(VectorVis,
                                     alignment=VectorVis.Alignment.Base,
                                     flat_shading=False,
                                     color=(1,0,0),
                                     title='Slip direction')
    
    slip_plane_vis = OvitoObject(VectorVis,
                                 alignment=VectorVis.Alignment.Base,
                                 flat_shading=False,
                                 color=(0,1,0),
                                 title='Slip plane')

    def modify(self, data: DataCollection, **kwargs):
        
        if self.alpha <= 0.0:
            raise ValueError("Shear strain threshold must be positive.")
        
        sheared = data.apply(ExpressionSelectionModifier(expression='ShearStrain>%5.8f' %self.alpha))
        
        slip_direction = np.zeros((data.particles.count,3))
        slip_plane = np.zeros((data.particles.count,3))
        tol=1e-2
        
        for pi in range(data.particles.count):
            is_selected = data.particles['Selection'][pi]
            if is_selected:
                Floc = np.transpose(data.particles['Deformation Gradient'][pi].reshape(3,3))
                Sloc = Floc-np.eye(3)
                normS=np.linalg.norm(Sloc)
                Sloc/=normS
                SlocSD = np.dot(Sloc,np.transpose(Sloc))
                SlocSP = np.dot(np.transpose(Sloc),Sloc)
                SlocSP[np.abs(SlocSP[:,:])<tol]=0.
                SlocSP[np.abs(SlocSP[:,:])<tol]=0.
                
                # Eigenvectors computation
                wm,vm = np.linalg.eigh(SlocSD)
                wn,vn = np.linalg.eigh(SlocSP)
                nloc=vn[:,np.argmax(wn)]
                if nloc[2]<0.: nloc*=-1
                mloc=vm[:,np.argmax(wm)]
                slip_direction[pi]=mloc
                slip_plane[pi]=nloc
        data.apply(ClearSelectionModifier())
        data.particles_.create_property("Slip direction", data=slip_direction).vis = self.slip_direction_vis
        data.particles_.create_property("Slip plane", data=slip_plane).vis = self.slip_plane_vis
