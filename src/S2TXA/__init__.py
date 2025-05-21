from ovito.modifiers import ExpressionSelectionModifier, ClearSelectionModifier
from ovito.pipeline import ModifierInterface
from ovito.data import *
from ovito.traits import *
from ovito.vis import *
from traits.api import *
from ovito.gui import UtilityInterface

import numpy as np
from itertools import permutations, product

class S2TXAPoleFigureOverlay(ViewportOverlayInterface):

    group1 = "Positioning"
    alignment =  Map({"Top left": (0.,1.,"north west"), "Top":(0.5, 1., "north"), "Top right":(1.,1., "north east"), "Right":(1., 0.5, "east"), "Bottom right": (1.,0., "south east"), "Bottom": (0.5,0., "south"), "Bottom left":(0.,0., "south west"), "Left":(0.,0.5, "west")}, label="Alignment", ovito_group=group1)
    px = Range(low=-1., high=1., value=0.0, label="X-offset", ovito_unit="percent", ovito_group=group1)
    py = Range(low=-1., high=1., value=0.0, label="Y-offset", ovito_unit="percent", ovito_group=group1)
    w = Range(low=0.05, high=1, value=0.5, label="Width", ovito_unit="percent", ovito_group=group1)
    h = Range(low=0.05, high=1, value=0.5, label="Height", ovito_unit="percent", ovito_group=group1)
    alpha = Range(low=0., high=1., value = 0.5, label="Opacity", ovito_unit="percent", ovito_group=group1)

    group2 = "Figure style"
    title = Str(label="Title", ovito_placeholder="‹auto›", ovito_group=group2)
    color_choice = ColorTrait(default=(0.401, 0.435, 1.0), ovito_group=group2, label = "Unicolor")
    font_size = Range(value = 1., low=0.01, label="Text scaling", ovito_unit="percent", ovito_group=group2)
    
    @staticmethod
    def normalize_vectors(vectors):
        return vectors / np.linalg.norm(vectors, axis=1, keepdims=True)

    @staticmethod
    def stereographic_projection(vectors):
        """
        Perform stereographic projection of 3D unit vectors onto a 2D plane.
        Vectors must be normalized.
        """
        # Stereographic projection from the south pole (0, 0, -1)
        x, y, z = vectors[:, 0], vectors[:, 1], vectors[:, 2]
        denominator = 1 + z  # projection from the south pole to xy-plane
        proj_x = x / denominator
        proj_y = y / denominator
        return proj_x, proj_y
        
    @staticmethod
    def setup_pole_figure(fig, ax, vectors):
        from matplotlib import patches

        n = np.array([0, 0, 1])
        for i in [0.2,0.4,0.6,0.8]:
            tt=np.array([i,0.,0.])
            ttp = tt - np.dot(tt, n) * n
            ctt = patches.Circle((0., 0.), radius=ttp[0], linewidth=2.0, edgecolor='black', fill=False, alpha=0.25)
            ax.add_patch(ctt)
        circle2 = patches.Circle((0., 0.), radius=1.0, linewidth=3.0, color='black', edgecolor='black', fill=False)
        ax.add_patch(circle2)

        ax.set_aspect('equal')

        planes=np.array([[ 0., 0., 1.],
                        [ 1., 1., 1.],
                        [-1., 1., 1.],
                        [ 1.,-1., 1.],
                        [-1.,-1., 1.],
                        [ 1., 0., 0.],
                        [ 0., 1., 0.]])
        planes_x, planes_y = S2TXAPoleFigureOverlay.stereographic_projection(S2TXAPoleFigureOverlay.normalize_vectors(planes))
        ax.scatter(planes_x, planes_y, s=40, color='darkred',zorder=np.inf)

        listplanes=[r'$[001]$',r'$[111]$',r'$[\overline{1}11]$',r'$[1\overline{1}1]$',r'$[\overline{1}\overline{1}1]$',r'$[100]$',r'$[010]$']
        offx=[0.0,0.1,-0.1,0.1,-0.1,0.1,0.]
        offy=[0.1,0.1,0.1,-0.15,-0.15,0.,0.1]
        horal=['center','center','center','center','center','left','center']
        for cnt in range(len(planes)):
            vx=planes_x[cnt]
            vy=planes_y[cnt]
            ax.text(vx+offx[cnt],vy+offy[cnt],listplanes[cnt],fontfamily='serif',weight='extra bold',color='black',horizontalalignment=horal[cnt],transform=ax.transData)

    def render(self, canvas: ViewportOverlayInterface.Canvas, data: DataCollection, frame: int, **kwargs):

        if 'Slip plane' not in data.particles or 'Slip direction' not in data.particles:
            raise RuntimeError(
                "Missing required particle properties: 'Slip plane' and 'Slip direction'.\n"
                "Please make sure to insert the 'S2TXA Modifier' earlier in the pipeline."
            )
        
        slip_planes = data.particles['Slip plane'][...]
        slip_dirs = data.particles['Slip direction'][...]

        with canvas.mpl_figure(pos=(self.alignment_[0] + self.px, self.alignment_[1] + self.py), size=(self.w, self.h), font_scale = self.font_size, anchor=self.alignment_[2], alpha=self.alpha, tight_layout=True) as fig:
            ax = fig.subplots(1,1)
            
            if self.title != "":
                ax.set_title(self.title)
            else:
                ax.set_title('S2TXA Pole Projection (Stereographic)')
            S2TXAPoleFigureOverlay.setup_pole_figure(fig, ax, slip_planes)

            # Normalize vectors to lie on unit sphere
            unit_vectors = S2TXAPoleFigureOverlay.normalize_vectors(slip_planes)
            # Stereographic projection
            proj_x, proj_y = S2TXAPoleFigureOverlay.stereographic_projection(unit_vectors)
            kwargs = {}
            kwargs['color'] = self.color_choice
            kwargs['alpha'] = self.alpha
            ax.scatter(proj_x, proj_y, s=10, **kwargs)

            ax.set_xlim([-1.2, 1.2])
            ax.set_ylim([-1.2, 1.2])
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.axis('off')
            ax.grid(False)
            
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
